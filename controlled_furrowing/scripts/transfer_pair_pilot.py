"""Freeze and run a four-field conditional transfer diagnostic without target outcomes."""

import argparse
import csv
from dataclasses import asdict, replace
from datetime import datetime, timezone
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
from statistics import median
import subprocess
import sys
import time

import numpy as np
import scipy
import torch

from annulus_observer import serial
from quadratic_endpoint_observer import Curve, _state_controls, positive_thickness, subtract
from quadratic_quadrature_study import AxisymmetricModel as CrossModel, evaluate, exact_geometry
from robust_quadratic import Parameters, SolvePolicy
from sparse_quadratic import AxisymmetricModel, IntegrationPolicy, SparsePolicy


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'runs/transfer_pair_pilot_v1'
PREP = OUT / 'preparation'
MEMO = 'runs/exact_integration_cohort_v1/scientific_next_action/receipt.json'
BASE_PROTOCOL = 'runs/joint_gain_extension_v2/protocol.json'
BASE_CALIBRATION = 'runs/joint_gain_extension_v2/direct_results.json'
JOINT_PROTOCOL = 'runs/material_cohort_direct_v1/protocol.json'
JOINT_CALIBRATION = 'runs/material_cohort_direct_v1/results.json'
SMALL = 'data/derived/calibration_small_18um.csv'
LARGE = 'data/derived/spatial_inputs_only.json'
ACTUAL_M = 'evidence/validated_identified_set_actual_M.json'
ACTUAL_M_COVER = 'evidence/validated_identified_set_actual_M_cover.json'
INPUTS = (MEMO, BASE_PROTOCOL, BASE_CALIBRATION, JOINT_PROTOCOL, JOINT_CALIBRATION,
          SMALL, LARGE, ACTUAL_M, ACTUAL_M_COVER, 'PLAN.md')
SOURCE_NAMES = ('transfer_pair_pilot.py', 'sparse_quadratic.py', 'quadratic_quadrature_study.py',
                'robust_quadratic.py', 'axisymmetric_quadratic.py', 'axisymmetric_mixed.py',
                'axisymmetric_prestress.py', 'axisymmetric_mechanics.py',
                'mixed_geometry_certificate.py', 'mixed_geometry_certificate_v2.py',
                'quadratic_endpoint_observer.py', 'annulus_observer.py')
PRIMARY = 'current_disk__spatial_max'
MATERIAL_KEYS = ('apical_area_modulus', 'basal_area_modulus',
                 'apical_surface_tension', 'basal_surface_tension')
LENGTH = Q('33.5')
ROOT_STEPS = 80
SAMPLE_CELLS = 100
SOLVER = SolvePolicy()
SPARSE = SparsePolicy()
MAIN_RULE = IntegrationPolicy(8, 12)
CROSS_RULE = IntegrationPolicy(12, 16)
REVIEW_KEYS = ('driver_sha256', 'protocol_sha256', 'preflight_sha256')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path, **kwargs):
    return json.loads(path.read_text(encoding='utf-8'), **kwargs)


def save(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(payload, allow_nan=False, indent=2) + '\n', encoding='utf-8')
    tmp.replace(path)


def now():
    return datetime.now(timezone.utc).isoformat()


def runtime():
    return {'python': sys.version.split()[0], 'numpy': np.__version__, 'scipy': scipy.__version__,
            'torch': str(torch.__version__), 'torch_device': 'cpu',
            'torch_threads': torch.get_num_threads()}


def hashes(paths):
    return {str(path): sha(ROOT / path) for path in paths}


def select_input(small, large):
    columns = {'record_id', 'reported_diameter_um', 'activation_magnitude'}
    require(len(small) == 17 and len(large) == 31, 'expected 17 small and 31 larger input rows')
    require(all(set(row) == columns for row in large), 'larger input allowlist schema violated')
    require(all(columns <= set(row) for row in small), 'small input columns missing')
    for rows in (small, large):
        require(len({row['record_id'] for row in rows}) == len(rows), 'duplicate source row')
        require(all(Q(row['activation_magnitude']) > 0 and Q(row['reported_diameter_um']) > 0
                    for row in rows), 'nonpositive dose or diameter')
    small_diameters = {Q(row['reported_diameter_um']) for row in small}
    require(len(small_diameters) == 1, 'one reported small diameter required')
    small_diameter = small_diameters.pop()
    small_median = median([Q(row['activation_magnitude']) for row in small])
    largest = max(Q(row['reported_diameter_um']) for row in large)
    require(largest > small_diameter, 'invalid transfer diameter ordering')
    selected = min((row for row in large if Q(row['reported_diameter_um']) == largest),
                   key=lambda row: (abs(Q(row['activation_magnitude']) - small_median), row['record_id']))
    return {'large_record_id': selected['record_id'], 'large_reported_diameter_um_exact': str(largest),
            'matched_small_reported_diameter_um_exact': str(small_diameter),
            'matched_dose_exact': str(Q(selected['activation_magnitude'])),
            'small_training_median_dose_exact': str(small_median)}


def calibration_sources(memo, input_hashes):
    base_protocol = read(ROOT / BASE_PROTOCOL)
    base = read(ROOT / BASE_CALIBRATION)
    joint_protocol = read(ROOT / JOINT_PROTOCOL)
    joint = read(ROOT / JOINT_CALIBRATION)
    require(base['protocol_sha256'] == input_hashes[BASE_PROTOCOL]
            and joint['protocol_sha256'] == input_hashes[JOINT_PROTOCOL], 'calibration replay/protocol mismatch')
    base_params = base_protocol['physical_parameters']
    joint_params = joint_protocol['new_parameters']
    require(base_params == joint_protocol['base_parameters'], 'baseline calibration material mismatch')
    require(base_protocol['integration_policy'] == joint_protocol['integration_policy'] == asdict(MAIN_RULE),
            'historical integration rules differ from declared rule')
    base_variant = next(v for v in base['variants'] if v['variant'] == PRIMARY)
    joint_variant = next(v for v in joint['variants'] if v['variant'] == PRIMARY)
    candidate = memo['screened_candidate_pair']
    base_gain = Q(base_variant['selected_exact_SURROGATE_gain'])
    joint_gain = Q(joint_protocol['new_gain_exact'])
    require(base_gain == Q(candidate['base_gain_exact']) == Q(joint_protocol['base_gain_exact'])
            and joint_gain == Q(candidate['joint_gain_exact']), 'gain provenance mismatch')
    require(joint_protocol['primary_variant'] == candidate['common_conditional_observer'] == PRIMARY,
            'observer provenance mismatch')
    require({k: base_params[k] for k in MATERIAL_KEYS} == candidate['baseline_material']
            and {k: joint_params[k] for k in MATERIAL_KEYS} == candidate['joint_material'],
            'screened material provenance mismatch')
    require(all(base_params[key] == value for key, value in joint_params.items() if key not in MATERIAL_KEYS),
            'non-material calibration geometry or coefficients changed')
    definitions = {}
    for name, parameters, gain, variant, results, protocol in (
        ('baseline', base_params, base_gain, base_variant, BASE_CALIBRATION, BASE_PROTOCOL),
        ('joint', joint_params, joint_gain, joint_variant, JOINT_CALIBRATION, JOINT_PROTOCOL),
    ):
        residuals = variant['all68_residuals']
        require(residuals['records'] == 17 and residuals['constraints'] == len(residuals['constraint_rows']) == 68,
                'calibration replay must retain all 68 small-record summaries')
        definitions[name] = {'parameters': parameters, 'gain_exact': str(gain), 'gain_display': float(gain),
                             'protocol': protocol, 'saved_replay': results,
                             'comparison_scope': '68 small calibration rows, conditional static stored curves'}
    decision = read(ROOT / ACTUAL_M)['primary_decision']
    require(decision['record_count'] == 17 and decision['constraint_count'] == 68
            and len(decision['constraints']) == 68 and len(decision['uniform_constraints_covered']) == 1
            and sum(row['status'] == 'UNRESOLVED' for row in decision['constraints']) == 67,
            'original 67/68 selected-family linkage changed')
    earlier = {'selected_family_record': ACTUAL_M, 'finite_cover_record': ACTUAL_M_COVER,
               'original_constraint_count': 68, 'original_unresolved_constraints': 67,
               'provided_constraint': decision['uniform_constraints_covered'][0],
               'parameter_cell_sha256': decision['parameter_cell_sha256'],
               'inference_context_sha256': decision['inference_context_sha256'],
               'all_branches_status': decision['all_branches_existential_parameter_status'],
               'scope': 'Archived one-channel selected-family application, distinct from the archived pair of stored-curve replays.'}
    return definitions, earlier


def dose_schedule(final_dose, gain):
    schedule = []
    for step in range(4):
        dose = step * final_dose / 3
        requested = dose * gain
        rounded = Q.from_float(float(requested))
        schedule.append({'dose_exact': str(dose), 'requested_strength_exact': str(requested),
                         'solver_strength_binary64_exact': str(rounded),
                         'binary64_rounding_difference_exact': str(rounded - requested)})
    return schedule


def create_protocol():
    input_hashes = hashes(INPUTS)
    memo = read(ROOT / MEMO)
    for path in INPUTS:
        if path in memo['source_and_input_hashes']:
            require(input_hashes[path] == memo['source_and_input_hashes'][path], f'memo binding changed: {path}')
    candidates, earlier = calibration_sources(memo, input_hashes)
    with (ROOT / SMALL).open(encoding='utf-8', newline='') as handle:
        selected = select_input(list(csv.DictReader(handle)), read(ROOT / LARGE, parse_float=str))
    require(all(v == memo['proposed_input_only_pilot'][k] for k, v in selected.items()),
            'outcome-free selection does not reproduce archived receipt')
    shoulder = Q(selected['large_reported_diameter_um_exact']) / 2
    cases = []
    for material, candidate in candidates.items():
        gain = Q(candidate['gain_exact'])
        for name, key in (('small', 'matched_small_reported_diameter_um_exact'),
                          ('large', 'large_reported_diameter_um_exact')):
            diameter = Q(selected[key])
            patch = diameter / (2 * LENGTH)
            case_params = replace(Parameters(**candidate['parameters']), patch_radius=float(patch))
            cases.append({'case_id': f'{material}_{name}', 'candidate': material, 'diameter_level': name,
                          'diameter_um_exact': str(diameter), 'patch_radius_model_exact': str(patch),
                          'patch_radius_model_binary64_exact': str(Q.from_float(case_params.patch_radius)),
                          'parameters': asdict(case_params), 'gain_exact': str(gain),
                          'schedule': dose_schedule(Q(selected['matched_dose_exact']), gain)})
    return {
        'version': 1,
        'scope': 'Two frozen calibration points crossed with two reported optical patch diameters in one conditional P2 finite model.',
        'source_hashes': hashes(f'scripts/{name}' for name in SOURCE_NAMES),
        'input_hashes': input_hashes,
        'runtime': runtime(),
        'input_selection': {**selected, 'rule': 'Largest diameter, then closest to median small dose, then record ID.',
                            'small_intervention': 'Dose-matched counterfactual', 'heldout_response_tables_read': False},
        'cases': cases,
        'calibration_candidates': candidates,
        'original_67_of_68_context': earlier,
        'loading': {'rule': 'Reported dose fractions [0, 1/3, 2/3, 1] shared by all four cases.',
                    'strength': 'Exact archived gain times exact optical dose, rounded once to binary64.',
                    'zero': 'Independent flat-reference zero vector for each case.',
                    'adaptive_inserts': 'Inherited robust bisection on failed strength intervals with a common numerical policy; insertions can differ and are recorded as implied optical doses.',
                    'interpretation': 'Numerical continuation with no stability or elapsed experimental time interpretation.'},
        'mesh': {'rule': 'P2 patch-fitted radial grid and activation-aligned depth layers inherited unchanged. Changing the patch radius regenerates the radial coordinates with the same connectivity and cell counts.',
                 'boundary': 'Axis radial component fixed at R=0. Both displacement components fixed on the outer rim. Apical/basal boundary conditions inherited from variational energy.'},
        'numerical_policy': {'solve': asdict(SOLVER), 'sparse': asdict(SPARSE),
                             'solve_quadrature': asdict(MAIN_RULE), 'cross_force_quadrature': asdict(CROSS_RULE),
                             'checkpoint_gate': 'Reached strength, all floating point residual/force gates, exact stored-node positive J and embedding, exact monotone apical/basal graphs and positive current thickness.'},
        'observable': {'name': 'Q', 'formula': 'z_apical(current R=0) - z_apical(current R=R_shoulder)',
                       'units': 'um', 'positive_direction': 'upward toward apical; Q<0 indicates center below shoulder',
                       'length_unit_um_exact': str(LENGTH), 'shoulder_radius_um_exact': str(shoulder),
                       'root_bisections': ROOT_STEPS, 'full_profile_sample_cells': SAMPLE_CELLS,
                       'thickness_status': 'Positive exact stored-curve thickness; inherited sufficient proof over full graph domain.'},
        'gap_semantics': 'Disjoint exact Q intervals show separation of these serialized P2 curves. Numerical residual/geometry diagnostics do not supply equilibrium, integration, mesh, source or physical-error radii.',
        'claims': {'root_validated': False, 'common_physical_branch_validated': False,
                   'source_compatible_witness': False, 'source_outcome_prediction_scored': False, 'novelty_established': False},
        'launch_gate': {'status': 'PASS', 'requires_all_true_nonempty_checks': True,
                        'binds_sha256': list(REVIEW_KEYS)},
    }


def verify_hashes(frozen, observed):
    require(frozen == observed, 'frozen path inventory or hash changed')


def assert_frozen(protocol):
    verify_hashes(protocol['source_hashes'], hashes(f'scripts/{name}' for name in SOURCE_NAMES))
    verify_hashes(protocol['input_hashes'], hashes(INPUTS))
    require(protocol['runtime'] == runtime(), 'frozen numerical runtime changed')


def verify_review(review, binding):
    require(review.get('status') == 'PASS' and isinstance(review.get('checks'), dict)
            and review['checks'] and all(value is True for value in review['checks'].values()),
            'independent prelaunch review with passing checks required')
    require(all(review.get(key) == binding[key] for key in REVIEW_KEYS),
            'independent prelaunch review is stale or not bound to this freeze')


def load_launch(review_path):
    protocol = read(PREP / 'protocol.json')
    preflight = read(PREP / 'preflight.json')
    assert_frozen(protocol)
    require(protocol == create_protocol(), 'frozen protocol differs from reproducible protocol')
    binding = {'driver_sha256': sha(Path(__file__)), 'protocol_sha256': sha(PREP / 'protocol.json'),
               'preflight_sha256': sha(PREP / 'preflight.json')}
    require(preflight['status'] == 'PASS' and preflight['checks']
            and all(value is True for value in preflight['checks'].values())
            and all(preflight[key] == binding[key] for key in ('driver_sha256', 'protocol_sha256')),
            'matching passed preflight required')
    verify_review(read(review_path), binding)
    return protocol, {**binding, 'prelaunch_review_sha256': sha(review_path)}


def surface_receipt(state, shoulder, full_profile=False):
    receipt = {'checks': {'controls_reconstructed': False, 'both_graphs_strictly_radial': False,
                          'full_common_domain': False, 'positive_current_graph_thickness': False,
                          'signed_Q_enclosed': False}, 'height_units': 'um', 'radius_units': 'um'}
    checks = receipt['checks']
    try:
        controls, height = _state_controls(state, LENGTH)
        receipt.update({'declared_flat_reference_height_um_exact': str(height),
                        'controls_um': {name: [[[str(r), str(z)] for r, z in edge] for edge in controls[key]]
                                        for name, key in (('apical', 'top'), ('basal', 'bottom'))}})
        checks['controls_reconstructed'] = True
        graphs = {label: Curve(controls[key], ROOT_STEPS)
                  for label, key in (('apical', 'top'), ('basal', 'bottom'))}
        derivatives = {name: min(min(c1[0], c1[0] + 2 * c2[0]) for _, _, c1, c2 in curve.curves)
                       for name, curve in graphs.items()}
        receipt['minimum_radius_derivative_um_per_local_parameter_exact'] = {
            name: str(value) for name, value in derivatives.items()}
        require(all(value > 0 for value in derivatives.values()), 'nonpositive minimum radial derivative')
        checks['both_graphs_strictly_radial'] = True
        outer = Q(state['parameters']['outer_radius']) * LENGTH
        require(all(curve.domain == (0, outer) for curve in graphs.values()), 'unexpected axis or rim extent')
        require(0 < shoulder < outer, 'shoulder radius outside open physical domain')
        checks['full_common_domain'] = True
        # Preserve full profiles and Q even if the separate thickness proof fails.
        apical = graphs['apical']
        center = apical.range(Q(0), Q(0))
        shoulder_height = apical.range(shoulder, shoulder)
        q_value = subtract(center, shoulder_height)
        receipt.update({'apical_center_um': serial(center), 'apical_shoulder_um': serial(shoulder_height),
                        'shoulder_radius_um_exact': str(shoulder), 'Q_um': serial(q_value),
                        'Q_um_exact': [str(q_value.lo), str(q_value.hi)]})
        checks['signed_Q_enclosed'] = True
        if full_profile:
            radii = {outer * Q(i, SAMPLE_CELLS) for i in range(SAMPLE_CELLS + 1)} | {shoulder}
            # Include edge endpoints and all exact extrema of each surface height.
            radii.update(r for curve in graphs.values() for r in curve.knots())
            receipt['sampled_full_current_profiles'] = [
                {'radius_um': float(r), 'radius_um_exact': str(r),
                 **{f'{name}_height_um': serial(curve.range(r, r)) for name, curve in graphs.items()}}
                for r in sorted(radii)]
        thickness = positive_thickness(graphs['apical'], graphs['basal'], (0, outer), Q('0.001'), 8192)
        receipt['minimum_current_thickness_um_lower_exact'] = str(thickness)
        checks['positive_current_graph_thickness'] = thickness > 0
    except (ValueError, ArithmeticError, KeyError) as error:
        receipt['observation_error'] = {'type': type(error).__name__, 'message': str(error)}
    receipt['usable_stored_graph'] = all(checks.values())
    return receipt


def build_models(parameters):
    return (AxisymmetricModel(parameters, SOLVER, MAIN_RULE, SPARSE),
            CrossModel(parameters, SOLVER, CROSS_RULE))


def mesh_receipt(model):
    params = model.parameters
    ref = model.reference
    return {'reference_node_count': len(ref), 'triangle_count': len(model.triangles),
            'free_dof_count': len(model.free), 'activation_element_count': int(sum(model.activation_elements)),
            'patch_radius_is_vertex_aligned': bool(params.patch_radius in model.radial),
            'axial_activation_faces_aligned': bool(params.activation_z_min in model.layers
                                                    and params.activation_z_max in model.layers),
            'has_axis_and_outer_rim_fixed_components': bool(np.array_equal(model.fixed[:, 0],
                (ref[:, 0] == 0) | (ref[:, 0] == params.outer_radius)) and
                np.array_equal(model.fixed[:, 1], ref[:, 0] == params.outer_radius)),
            'positive_reference_element_area': bool(np.min(model.area) > 0),
            'reference_coordinate_sha256': hashlib.sha256(ref.astype('<f8').tobytes()).hexdigest(),
            'connectivity_sha256': hashlib.sha256(model.triangles.astype('<i8').tobytes()).hexdigest(),
            'fixed_mask_sha256': hashlib.sha256(model.fixed.tobytes()).hexdigest(),
            'activation_mask_sha256': hashlib.sha256(model.activation_elements.tobytes()).hexdigest()}


def verify_mesh_pair(primary, cross):
    require(mesh_receipt(primary) == mesh_receipt(cross), 'cross quadrature changed the finite field space')


def identity_state(model):
    return {'parameters': asdict(model.parameters), 'activation_strength': 0.,
            'reference_nodes': model.reference.tolist(), 'current_nodes': model.reference.tolist(),
            'triangles': model.triangles.tolist(), 'free_displacement': np.zeros(len(model.free)).tolist(),
            'boundary_top_quadratic_controls': model.reference[model.top_edges].tolist(),
            'boundary_bottom_quadratic_controls': model.reference[model.bottom_edges].tolist()}


def rejects(action):
    try:
        action()
    except ValueError:
        return True
    return False


def synthetic_checks(protocol):
    nominal = protocol['input_selection']
    cases = protocol['cases']
    expected_dose = Q('0.89767374')
    shoulder = Q(protocol['observable']['shoulder_radius_um_exact'])
    checks = {
        'four_unique_named_material_diameter_cells': {case['case_id'] for case in cases}
        == {'baseline_small', 'baseline_large', 'joint_small', 'joint_large'} and len(cases) == 4,
        'actual_17_to_31_input_rule': nominal['large_record_id'] == 'fig4f_col_I'
        and Q(nominal['small_training_median_dose_exact']) == Q('0.78208913')
        and Q(nominal['matched_small_reported_diameter_um_exact']) == Q('18.71024074')
        and Q(nominal['large_reported_diameter_um_exact']) == Q('49.8939753')
        and Q(nominal['matched_dose_exact']) == expected_dose,
        'same_exact_dose_schedule_all_four': all([Q(step['dose_exact']) for step in case['schedule']]
            == [Q(0), expected_dose / 3, 2 * expected_dose / 3, expected_dose] for case in cases),
        'all_rounded_strengths_reproduce': all(Q(step['requested_strength_exact'])
            == Q(step['dose_exact']) * Q(case['gain_exact']) and Q(step['solver_strength_binary64_exact'])
            == Q.from_float(float(Q(step['requested_strength_exact']))) and Q(step['binary64_rounding_difference_exact'])
            == Q(step['solver_strength_binary64_exact']) - Q(step['requested_strength_exact'])
            for case in cases for step in case['schedule']),
        'fixed_material_and_gain_within_each_candidate': all(case['gain_exact']
            == protocol['calibration_candidates'][case['candidate']]['gain_exact'] and
            {k: v for k, v in case['parameters'].items() if k != 'patch_radius'} ==
            {k: v for k, v in protocol['calibration_candidates'][case['candidate']]['parameters'].items()
             if k != 'patch_radius'} for case in cases),
        'gains_differ_without_recalibration': len({case['gain_exact'] for case in cases}) == 2,
        'exact_diameter_to_radius_and_units': shoulder == Q('24.94698765') and all(
            2 * LENGTH * Q(case['patch_radius_model_exact']) == Q(case['diameter_um_exact'])
            and Q(case['patch_radius_model_binary64_exact']) == Q.from_float(case['parameters']['patch_radius'])
            for case in cases),
    }
    fixture = lambda label, diameter, dose: {'record_id': label, 'reported_diameter_um': diameter,
                                           'activation_magnitude': dose}
    small = [fixture(f's{i:02}', '2', '1') for i in range(17)]
    large = [fixture(f'l{i:02}', '3', '1') for i in range(28)] + [
        fixture('lower_dose_lexically_later', '10', '0.5'), fixture('higher_dose', '10', '1.5'),
        fixture('distance_runner_up', '10', '2')]
    checks['selection_prioritizes_size_then_exact_distance_then_id'] = (
        select_input(small, large)['large_record_id'] == 'higher_dose')
    checks['changed_input_schema_blocks_outcomes'] = rejects(
        lambda: select_input(small, [*large[:-1], {**large[-1], 'heldout_response': 'blocked'}]))
    checks['changed_input_duplicate_record_blocks'] = rejects(lambda: select_input(small, [*large[:-1], large[-2]]))
    with_file_change = {**protocol['input_hashes'], LARGE: 'changed'}
    with_source_change = {**protocol['source_hashes'], 'scripts/robust_quadratic.py': 'changed'}
    checks['changed_input_hash_blocks'] = rejects(lambda: verify_hashes(protocol['input_hashes'], with_file_change))
    checks['changed_dependency_hash_blocks'] = rejects(lambda: verify_hashes(protocol['source_hashes'], with_source_change))
    binding = {key: key for key in REVIEW_KEYS}
    review = {'status': 'PASS', 'checks': {'fixture': True}, **binding}
    verify_review(review, binding)
    for key in REVIEW_KEYS:
        checks[f'changed_review_{key}_blocks'] = rejects(lambda: verify_review({**review, key: 'changed'}, binding))
    checks['missing_or_failed_review_checks_block'] = all(rejects(lambda: verify_review({**review, 'checks': v}, binding))
                                                          for v in ({}, {'fixture': False}, {'fixture': 1}))

    scale = LENGTH
    # In model units r=t+t² and z=10−2t+t². At R=(3/4)L this has t=1/2.
    first = [(0, 10), (Q(3, 4), Q(37, 4)), (2, 9)]
    curve = Curve([[(scale * r, scale * z) for r, z in first]], ROOT_STEPS)
    contrast = subtract(curve.range(0, 0), curve.range(scale * Q(3, 4), scale * Q(3, 4)))
    checks['signed_center_minus_shoulder_with_exact_scale'] = contrast.lo == contrast.hi == scale * Q(3, 4)
    shifted = Curve([[(scale * r, scale * (z - 10)) for r, z in first]], ROOT_STEPS)
    checks['signed_contrast_translation_invariant'] = (
        subtract(shifted.range(0, 0), shifted.range(scale * Q(3, 4), scale * Q(3, 4))) == contrast)
    checks['current_radius_rejects_outside_domain'] = rejects(lambda: curve.range(3 * scale, 3 * scale))
    checks['folded_p2_surface_rejected'] = rejects(lambda: Curve([[(0, 0), (2, 1), (1, 2)]]))
    irrational = Curve([[(0, 0), (3, 2), (8, 4)]], ROOT_STEPS)
    bound = irrational.range(Q(2), Q(2))
    radius = lambda z: z + z * z / 4  # r=4t+4t² and z=4t
    checks['irrational_radius_inverse_is_enclosed'] = (bound.lo < bound.hi
        and radius(bound.lo) <= 2 <= radius(bound.hi) and bound.hi - bound.lo <= Q(4, 2**ROOT_STEPS))
    checks['signed_interval_subtraction'] = difference(['-4', '-3'], ['1', '2']) == ['-6', '-4']
    return checks


def preflight(protocol):
    checks = synthetic_checks(protocol)
    meshes = {}
    shoulder = Q(protocol['observable']['shoulder_radius_um_exact'])
    for case in protocol['cases']:
        case_id = case['case_id']
        model, cross = build_models(Parameters(**case['parameters']))
        verify_mesh_pair(model, cross)
        mesh = mesh_receipt(model)
        state = identity_state(model)
        geometry = exact_geometry(state)
        surfaces = surface_receipt(state, shoulder, full_profile=case_id == 'baseline_large')
        main = evaluate(model, state)
        alt = evaluate(cross, state)
        checks[f'{case_id}_active_aligned_mesh_and_fixed_bc'] = (
            mesh['patch_radius_is_vertex_aligned'] and mesh['axial_activation_faces_aligned'] and
            mesh['has_axis_and_outer_rim_fixed_components'] and mesh['positive_reference_element_area']
            and 0 < mesh['activation_element_count'] < mesh['triangle_count'])
        checks[f'{case_id}_positive_J_embedding_and_zero_shape'] = (
            geometry['all_exact_J_positive'] and geometry['global_embedding_certified']
            and surfaces['usable_stored_graph'] and surfaces['Q_um_exact'] == ['0', '0'])
        checks[f'{case_id}_initial_zero_residual_both_rules'] = (
            main['scaled_force_linf'] < SOLVER.force_tolerance and alt['scaled_force_linf'] < SOLVER.force_tolerance)
        meshes[case_id] = {**mesh, 'Q8_residual': main['scaled_force_linf'], 'Q12_residual': alt['scaled_force_linf'],
                           'exact_J_lower': geometry['J_lower_exact'], 'exact_J_upper': geometry['J_upper_exact'],
                           'Q_um_exact': surfaces.get('Q_um_exact'), 'cross_rule_same_mesh': True}
        if case_id == 'baseline_large':
            profile = surfaces.get('sampled_full_current_profiles', [])
            radii = [Q(row['radius_um_exact']) for row in profile]
            checks['actual_identity_profile_covers_both_full_surfaces'] = (
                len(radii) >= SAMPLE_CELLS + 1 and radii[0] == 0 and radii[-1] == 3 * LENGTH and shoulder in radii
                and all(row['apical_height_um']['lower']['display'] == float(LENGTH)
                        and row['basal_height_um']['upper']['display'] == 0 for row in profile))
    for key in ('connectivity_sha256', 'fixed_mask_sha256', 'activation_mask_sha256', 'free_dof_count'):
        checks[f'common_mesh_topology_and_support_{key}'] = len({mesh[key] for mesh in meshes.values()}) == 1
    checks['only_diameter_regenerates_reference_grid'] = all(
        meshes[f'baseline_{level}']['reference_coordinate_sha256'] == meshes[f'joint_{level}']['reference_coordinate_sha256']
        for level in ('small', 'large')) and (
        meshes['baseline_small']['reference_coordinate_sha256'] != meshes['baseline_large']['reference_coordinate_sha256'])
    return {'status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks, 'actual_initial_meshes': meshes,
            'driver_sha256': sha(Path(__file__)), 'protocol_sha256': sha(PREP / 'protocol.json'),
            'heldout_response_tables_read': False, 'nonlinear_jobs_launched': 0,
            'scope': 'Synthetic checks, original calibration linkage, actual zero-node field and residual evaluation. No nonlinear optimizer or continuation invoked.'}


def prepare():
    require(not (PREP / 'protocol.json').exists(), 'preserve the existing frozen protocol')
    protocol = create_protocol()
    save(PREP / 'protocol.json', protocol)
    result = preflight(protocol)
    save(PREP / 'preflight.json', result)
    print(json.dumps({'status': result['status'], 'checks': len(result['checks']),
                      'driver_sha256': result['driver_sha256'], 'protocol_sha256': result['protocol_sha256'],
                      'preflight_sha256': sha(PREP / 'preflight.json')}), flush=True)
    require(result['status'] == 'PASS', 'preflight contains failed checks')


def checkpoint(model, cross, state, step, shoulder, last):
    strength = float(Q(step['solver_strength_binary64_exact']))
    main = evaluate(model, state)
    alternate = evaluate(cross, state)
    geometry = exact_geometry(state)
    surfaces = surface_receipt(state, shoulder, full_profile=last)
    checks = {
        'reached_requested_strength': state['reached_target'] is True and state['activation_strength'] == strength,
        'solver_converged': state['converged'] is True,
        'Q8_replayed_residual_below_limit': main['scaled_force_linf'] < SOLVER.force_tolerance,
        'Q12_cross_residual_below_limit': alternate['scaled_force_linf'] < SOLVER.force_tolerance,
        'exact_stored_J_positive': geometry['all_exact_J_positive'] is True and Q(geometry['J_lower_exact']) > 0,
        'exact_stored_global_embedding': geometry['global_embedding_certified'] is True,
        'profile_checks': surfaces['usable_stored_graph'] is True,
    }
    return {'status': 'PASS' if all(checks.values()) else 'INCOMPLETE', 'checks': checks,
            'Q8_mechanical_evaluation': main, 'Q12_cross_evaluation': alternate,
            'geometry': geometry, 'surface': surfaces, 'root_validated': False}


def relative(path):
    return path.relative_to(ROOT).as_posix()


def implied_dose(accepted, gain):
    return [{'from_dose_exact': str(Q.from_float(a['from_strength']) / gain),
             'to_dose_exact': str(Q.from_float(a['to_strength']) / gain),
             'from_dose_display': a['from_strength'] / float(gain),
             'to_dose_display': a['to_strength'] / float(gain)} for a in accepted]


def run_worker(review_path, case_id):
    protocol, bindings = load_launch(review_path)
    planned = {item['case_id']: item for item in protocol['cases']}
    case = planned[case_id]
    folder = OUT / 'jobs' / case_id
    # The supervisor owns subprocess launch and OS-exit receipts.
    require(folder.is_dir() and not (folder / 'worker.json').exists(), 'expected fresh supervisor-owned case output')
    start = time.monotonic()
    report = {'status': 'RUNNING', 'worker_pid': os.getpid(), 'started_at_utc': now(),
              'case': case, 'bindings': bindings, 'source_hashes': protocol['source_hashes'],
              'input_hashes': protocol['input_hashes'], 'checkpoints': [], 'root_validated': False}
    save(folder / 'worker.json', report)
    try:
        model, cross = build_models(Parameters(**case['parameters']))
        verify_mesh_pair(model, cross)
        report['mesh'] = mesh_receipt(model)
        shoulder = Q(protocol['observable']['shoulder_radius_um_exact'])
        gain = Q(case['gain_exact'])
        previous = None
        for index, step in enumerate(case['schedule']):
            assert_frozen(protocol)
            stage = folder / f'checkpoint_{index:02d}'
            last = index == len(case['schedule']) - 1
            state = model.advance(float(Q(step['solver_strength_binary64_exact'])), previous,
                                  output_directory=stage / 'solver')
            field_path = stage / 'field.json'
            save(field_path, state)
            receipt = checkpoint(model, cross, state, step, shoulder, last or not state['reached_target'])
            verify_path = stage / 'verification.json'
            save(verify_path, receipt)
            continuation = state['robust_continuation']
            row = {
                **step, 'index': index, 'field': relative(field_path), 'field_sha256': sha(field_path),
                'verification': relative(verify_path), 'verification_sha256': sha(verify_path),
                'actual_strength_binary64_exact': str(Q.from_float(state['activation_strength'])),
                'actual_dose_implied_by_binary_strength_exact': str(Q.from_float(state['activation_strength']) / gain),
                'converged': state['converged'], 'optimizer_success': state['continuation'][-1]['optimizer_success'],
                'numerical_diagnostic_status': receipt['status'], 'checks': receipt['checks'],
                'Q8_residual': receipt['Q8_mechanical_evaluation']['scaled_force_linf'],
                'Q12_residual': receipt['Q12_cross_evaluation']['scaled_force_linf'],
                'exact_J_lower': receipt['geometry']['J_lower_exact'], 'exact_J_upper': receipt['geometry']['J_upper_exact'],
                'attempts': len(continuation['attempts']),
                'rejected_attempts': sum(not attempt['accepted'] for attempt in continuation['attempts']),
                'accepted_path_implied_optical_doses': implied_dose(continuation['accepted_increments'], gain),
                'continuation_stop_reason': continuation['stop_reason'],
                'Q_um_exact': receipt['surface'].get('Q_um_exact'),
            }
            report['checkpoints'].append(row)
            report['last_retained_field'] = row
            reached = receipt['status'] == 'PASS'
            if last or not reached:
                # Full serialized curves are always in verification.json; preserve a
                # dense profile even when a different gate stops this checkpoint.
                if not last and state['reached_target']:
                    receipt['surface'] = surface_receipt(state, shoulder, full_profile=True)
                profile_path = folder / 'profiles.json'
                save(profile_path, {'case_id': case_id, 'planned_final_checkpoint': last,
                                    'diagnostic_gate_passed': reached, 'field_sha256': sha(field_path),
                                    **receipt['surface']})
                report.update({'profiles': relative(profile_path), 'profiles_sha256': sha(profile_path)})
            if not reached:
                report['status'] = 'CHECKPOINT_INCOMPLETE'
            save(folder / 'worker.json', report)
            print(json.dumps({'case_id': case_id, 'checkpoint': index, 'status': receipt['status'],
                              'Q8_residual': row['Q8_residual'], 'Q12_residual': row['Q12_residual']}), flush=True)
            if not reached:
                break
            previous = state
        assert_frozen(protocol)
        if report['status'] == 'RUNNING':
            report['status'] = 'COMPLETE_CONDITIONAL_STORED_FIELD'
    except Exception as error:
        report['status'] = 'WORKER_ERROR'
        report['error'] = {'type': type(error).__name__, 'message': str(error)}
        raise
    finally:
        report['worker_finished_at_utc'] = now()
        report['elapsed_seconds'] = time.monotonic() - start
        save(folder / 'worker.json', report)
    return 0 if report['status'] == 'COMPLETE_CONDITIONAL_STORED_FIELD' else 2


def launch_case(review_path, case, bindings):
    case_id = case['case_id']
    folder = OUT / 'jobs' / case_id
    folder.mkdir(parents=True, exist_ok=False)
    sidecar = folder / 'job.json'
    start = time.monotonic()
    command = [sys.executable, '-B', str(Path(__file__).resolve()), 'worker', '--review', str(review_path),
               '--case', case_id]
    receipt = {'case_id': case_id, 'bindings': bindings, 'command': command,
               'started': False, 'terminal': False, 'started_at_utc': None, 'pid': None, 'exit_code': None}
    try:
        with (folder / 'stdout.log').open('wb') as stdout, (folder / 'stderr.log').open('wb') as stderr:
            process = subprocess.Popen(command, cwd=ROOT, stdout=stdout, stderr=stderr)
            receipt.update({'status': 'PROCESS_STARTED', 'started': True, 'pid': process.pid, 'started_at_utc': now()})
            save(sidecar, receipt)
            print(json.dumps({'case_id': case_id, 'started': True, 'pid': process.pid}), flush=True)
            code = process.wait()
    except OSError as error:
        if receipt['started']:
            raise
        receipt.update({'status': 'LAUNCH_FAILED', 'terminal': True, 'ended_at_utc': now(),
                        'error': {'type': type(error).__name__, 'message': str(error)}})
        save(sidecar, receipt)
        return receipt
    receipt.update({'status': 'PROCESS_EXITED', 'terminal': True, 'exit_code': code, 'ended_at_utc': now(),
                    'elapsed_seconds': time.monotonic() - start})
    result = folder / 'worker.json'
    if result.is_file():
        receipt.update({'worker_receipt': relative(result), 'worker_receipt_sha256': sha(result),
                        'worker_status': read(result)['status']})
    save(sidecar, receipt)
    print(json.dumps({'case_id': case_id, 'terminal': True, 'exit_code': code,
                      'worker_status': receipt.get('worker_status')}), flush=True)
    return receipt


def difference(after, before):
    a, b = map(Q, after)
    c, d = map(Q, before)
    return [str(a - d), str(b - c)]


def interval_receipt(bounds):
    return {'exact': bounds, 'display_um': [float(Q(value)) for value in bounds]}


def summarize(jobs):
    completed = {}
    for job in jobs:
        if job['exit_code'] == 0 and job.get('worker_status') == 'COMPLETE_CONDITIONAL_STORED_FIELD':
            report = read(ROOT / job['worker_receipt'])
            completed[job['case_id']] = report['checkpoints'][-1]
    q = {case_id: interval_receipt(row['Q_um_exact']) for case_id, row in completed.items()}
    diameter_induced = {name: interval_receipt(difference(q[f'{name}_large']['exact'], q[f'{name}_small']['exact']))
                        for name in ('baseline', 'joint') if f'{name}_small' in q and f'{name}_large' in q}
    pair = {'status': 'INCOMPLETE', 'joint_minus_baseline_large_Q_um': None}
    if 'baseline_large' in q and 'joint_large' in q:
        large_gap = difference(q['joint_large']['exact'], q['baseline_large']['exact'])
        pair.update({'status': 'DISJOINT_STORED_CURVE_INTERVALS' if Q(large_gap[0]) > 0 or Q(large_gap[1]) < 0
                     else 'OVERLAP_STORED_CURVE_INTERVALS', 'joint_minus_baseline_large_Q_um': interval_receipt(large_gap)})
    return {'status': 'COMPLETE_FOUR_CONDITIONAL_FIELDS' if len(completed) == 4 else 'INCOMPLETE_PILOT',
            'terminal_case_receipts': jobs, 'diagnostic_fields_accepted': completed,
            'signed_Q_um': q, 'signed_diameter_change_large_minus_small_Q_um': diameter_induced,
            'large_diameter_gap_test': pair, 'source_outcomes_read_or_scored': False,
            'root_validated': False, 'source_compatible_witness': False, 'novelty_established': False,
            'scope': 'Intervals enclose observable evaluation for stored-node curves only. There are no propagated root, integration, continuum, or source uncertainty bounds.'}


def run(review_path):
    review_path = review_path.resolve()
    protocol, bindings = load_launch(review_path)
    (OUT / 'jobs').mkdir(exist_ok=False)
    job_path = OUT / 'job.json'
    started = time.monotonic()
    jobs = []
    receipt = {'status': 'RUNNING', 'pid': os.getpid(), 'started_at_utc': now(), 'bindings': bindings}
    save(job_path, receipt)
    for case in protocol['cases']:
        assert_frozen(protocol)
        jobs.append(launch_case(review_path, case, bindings))
        save(job_path, {**receipt, 'terminal_cases': jobs})
    assert_frozen(protocol)
    result = summarize(jobs)
    output = OUT / 'results.json'
    save(output, {**result, 'bindings': bindings, 'source_hashes': protocol['source_hashes'],
                  'input_hashes': protocol['input_hashes']})
    save(job_path, {**receipt, 'status': result['status'], 'ended_at_utc': now(), 'terminal_cases': jobs,
                    'elapsed_seconds': time.monotonic() - started, 'results_sha256': sha(output)})
    print(json.dumps({'status': result['status'], 'results_sha256': sha(output)}), flush=True)
    return 0 if result['status'] == 'COMPLETE_FOUR_CONDITIONAL_FIELDS' else 2


def main():
    torch.set_num_threads(1)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('prepare', 'run', 'worker'))
    parser.add_argument('--review', type=Path)
    parser.add_argument('--case', choices=('baseline_small', 'baseline_large', 'joint_small', 'joint_large'))
    args = parser.parse_args()
    if args.stage == 'prepare':
        prepare()
        return 0
    if args.review is None:
        parser.error('scheduled execution requires an independent prelaunch review')
    if args.stage == 'run':
        return run(args.review)
    if args.case is None:
        parser.error('internal worker requires one case')
    return run_worker(args.review, args.case)


if __name__ == '__main__':
    sys.exit(main())
