"""Four exact-spatial-FE singleton targets from the reviewed nonlinear point caches."""
import argparse
import ast
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
import json
import os
from pathlib import Path
from time import monotonic
import traceback
from unittest.mock import patch

from flint import ctx as arb_context
import numpy as np

import exact_cohort_integrator as E
import exact_fe_target_safe as T
import root_body_geometry as G
import transfer_ambiguity_cached_pilots as CP
import transfer_ambiguity_observer_sweep as sweep


SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[1]
BASE = f'{sweep.PRIMARY}/certification_design'
CACHE = f'{BASE}/cached_pilots_v1'
META = f'{CACHE}/independent_actual_preparation'
REPLAY = f'{META}/actual_replay'
OUT = ROOT / BASE / 'exact_targets_v1'
PREP = OUT / 'preparation'
GEOMETRY = OUT / 'geometry'
CACHE_ADMISSION = f'{CACHE}/terminal_admission.json'
CACHE_RESULT = f'{CACHE}/actual/results.json'
CACHE_TERMINAL = f'{CACHE}/actual/terminal.json'
CACHE_REVIEW = f'{REPLAY}/receipt.json'
CACHE_REVIEW_TERMINAL = f'{REPLAY}/terminal.json'
CACHE_REVIEW_OS = f'{REPLAY}/os_exit_observation.json'
BODY_REVIEW = 'evidence/root_body_geometry_independent_review.json'
AXIS_REVIEW = 'evidence/exact_axis_edge_physical_independent_review.json'
RADII = (Q(1, 2**36), Q(1, 2**34))
CHARTS = {'off_axis': 132, 'axis_vertex': 6, 'axis_edge_physical': 6}
RANK = 611
PILOTS = tuple(f'{endpoint}_{stage}' for endpoint in sweep.ENDPOINTS for stage in ('small_17', 'large_03'))
PINS = {
    'scripts/transfer_ambiguity_cached_pilots.py': '1b22ebc307c0cee51bf4193ff3999242b1cfe278d07c10cc63619f2795489c58',
    'scripts/exact_fe_target_safe.py': '889b61cfa794d432021fbd4670b781169033d4b8451554b18e7f7a569d9d3954',
    'scripts/exact_cohort_integrator.py': '092dc51161e708cdba3791ec47a7f8e843dd0e1afc866b6f65f286b8657d1745',
    'scripts/exact_axis_edge_physical.py': 'ebe0cbe021e04bf51aaf8985df6bfbe6f6fd20e1fc60194acbd1fd968654199f',
    'scripts/exact_uniform_surface_bridge_B_safe.py': '61cf47d2b1a9e7ee46f9fa680c24b4fdf7b4e31e4889379c56951ab6a7455ee9',
    'scripts/root_body_geometry.py': '18013c6c50c5f679eb6c4abc04e3cf9ce447ee18d0bf1efbf6dcfbd6806f49c2',
    'scripts/mixed_geometry_certificate_v2.py': '165d3b914b1ab9bad2e6591a7df63c09e004f33a742cfa5f4430584746edfbd2',
    'scripts/mixed_geometry_certificate.py': 'fb0d49d680286ff49d0dc6713baca41fbebbf724442bff9731292e1441b03fd8',
    BODY_REVIEW: 'f7463b263d39294d576417372469f5580c651e5070ca2ea627aaef86cdcbd6f2',
    AXIS_REVIEW: '8eaa76bd6ee3a6980815e303aba619c3f21243191c44035e0b1c3dabb26edffd',
    CACHE_ADMISSION: '406abdfc4dcc2f2a09e7081362d4ec77bd727d78f1c003224fa2cecfc61265d7',
    CACHE_RESULT: '3138dd8c6c72dd086a9864d6ac3e662d8e1a70ecbd17d401c688687e0beb0757',
    CACHE_TERMINAL: '8c67e43fa12024ef5e8eb15a3b101cd8891d0d77d58b47bb0de9ad8b797887f4',
    CACHE_REVIEW: '4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17',
    CACHE_REVIEW_TERMINAL: '08da2c2f5754d3f8bf8b392157228f52ed316ac2e9d6af53fc8745181b658a8a',
    CACHE_REVIEW_OS: '65858f55c1281d3953b775f147c7475a8604891f37cf56ed0fdd47697e6df6e3',
    f'{META}/protocol.json': 'c476eb5c5ae406f2515e79475705efcb62ddf820fad3dc13679d08838b9f5382',
    f'{META}/preflight.json': 'd3407cddd0ae54ae9fa4004d7d088019b62222b970c29648a1e0a28f0fb302b1',
    f'{META}/replay_frozen.py': '3f71e5cae5ee314e260e3bef6c478a923fb7c28f446c184c53d912d6ae296a94',
}
SCOPE = ('Four independently specified fixed-parameter axisymmetric P2 exact-spatial-FE potential targets. '
         'The singleton is in material-parameter space. Each explicit free-coordinate ball is a proposed '
         'local geometry domain around the exact binary-coordinate sum X+inj(u). '
         'The cached proof concerns quadrature only. Whole-box spatial derivative, exact-integral root, '
         'local minimum, continuum and observational claims require separate work.')
require, save, read, pinned, sha = sweep.require, sweep.save, sweep.read, sweep.pinned, sweep.sha


@dataclass(frozen=True)
class SingletonContext:
    pilot_id: str
    record_id: str
    row: dict
    target: T.Target
    S: np.ndarray
    cell: dict
    base: dict
    dose: Q
    drive: Q

    @property
    def R(self):
        return self.target.a['preconditioner']

    @property
    def u(self):
        return self.target.a['u']

    @property
    def free(self):
        return self.target.a['free']


def relative(path):
    return path.resolve().relative_to(ROOT).as_posix()


def exact(value):
    require(type(value) in (str, int, Q), 'non-Boolean exact rational spelling required')
    return Q(value)


def doc(value):
    return json.loads(json.dumps(value, default=str))


def array_sha(value):
    return sha256(value.tobytes(order='C')).hexdigest()


def link(links, name, digest):
    path = (ROOT / name).resolve()
    require(path.is_relative_to(ROOT) and path.relative_to(ROOT).as_posix() == name
            and links.get(name, digest) == digest and sha(path) == digest,
            f'changed or mixed authority bytes: {name}')
    links[name] = digest


def bootstrap():
    arb_context.prec = 640
    T.core.A.cache_clear()
    require(E.core is T.core and E.T is T and E.SUR.T is T and arb_context.prec == 640,
            'corrected Target, shared mathematical aliases and Arb640 required')


def local_source_closure():
    seeds = [SOURCE, Path(E.__file__), Path(T.__file__), Path(G.__file__), Path(CP.__file__), Path(sweep.__file__)]
    output = {}
    while seeds:
        path = seeds.pop().resolve()
        name = relative(path)
        if name in output:
            continue
        output[name] = sha(path)
        tree = ast.parse(path.read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            modules = ([x.name for x in node.names] if isinstance(node, ast.Import) else
                       [node.module] if isinstance(node, ast.ImportFrom) else [])
            for module in modules:
                filename = ROOT / 'scripts' / f"{(module or '').split('.')[0]}.py"
                if filename.is_file():
                    seeds.append(filename)
    return dict(sorted(output.items()))


def authority():
    bootstrap()
    links = {}
    for name, h in PINS.items():
        link(links, name, h)
    for name in (BODY_REVIEW, AXIS_REVIEW):
        require(read(ROOT / name)['status'].startswith('PASS'), 'reviewed lower mathematical kernel required')
    _, inherited, four = CP.admitted_targets()
    for name, h in {**CP.PINS, **inherited}.items():
        link(links, name, h)
    admission, result, terminal, review, replay_end, replay_os = (
        read(ROOT / name) for name in (CACHE_ADMISSION, CACHE_RESULT, CACHE_TERMINAL, CACHE_REVIEW,
                                       CACHE_REVIEW_TERMINAL, CACHE_REVIEW_OS))
    require(admission['status'] == 'PASS_ROOT_OBSERVED_CACHED_PILOT_TERMINALS'
            and type(admission['observed_exit_code']) is int and admission['observed_exit_code'] == 0
            and len(admission['sha256']) == 23 and admission['cache_payloads_opened'] is False
            and admission['supervisor_pid'] == terminal['pid']
            and terminal['status'] == 'FINISHED' and terminal['results_sha256'] == PINS[CACHE_RESULT]
            and result['status'] == 'COMPLETE_FOUR_CACHED_QUADRATURE_PILOT_ATTEMPTS'
            and terminal['bindings'] == result['bindings'], 'original successful cached supervisor required')
    require(review['status'] == 'PASS_INDEPENDENT_FOUR_CACHED_LOCAL_ROOTS'
            and sweep.all_true(review['checks'])
            and review['root_certificates_cached_quadrature_only'] == 4
            and review['root_certificates_exact_spatial_integral'] == 0
            and review['parents']['parent_admission_sha256'] == PINS[CACHE_ADMISSION]
            and replay_end['status'] == 'FINISHED' and type(replay_end['exit_code']) is int
            and replay_end['exit_code'] == 0 and replay_end['receipt_sha256'] == PINS[CACHE_REVIEW]
            and replay_os['status'] == 'PASS_OBSERVED_SERIAL_INDEPENDENT_REPLAY_TRUE_OS_EXIT'
            and type(replay_os['observed_exit_code']) is int and replay_os['observed_exit_code'] == 0
            and replay_end['pid'] == replay_os['pid'] == review['parents']['pid']
            and replay_os['terminal_sha256'] == PINS[CACHE_REVIEW_TERMINAL]
            and replay_os['receipt_sha256'] == PINS[CACHE_REVIEW],
            'independent cached proof/exit does not match the admitted parents')
    for name, digest in review['producer_and_payload_sha256'].items():
        link(links, name, digest)
    require(all(review['producer_and_payload_sha256'][p] == h for p,h in admission['sha256'].items()),
            'independent cached proof does not cover all 23 terminal pins')
    require(tuple(p['pilot_id'] for p in four) == tuple(p['pilot_id'] for p in result['pilots'])
            == tuple(p['pilot_id'] for p in review['field_results']) == PILOTS,
            'expected four distinct small/large singleton point identities')
    data = {}
    for fresh, produced, independent in zip(four, result['pilots'], review['field_results'], strict=True):
        pid = fresh['pilot_id']
        require(independent['status'] == 'PASS_INDEPENDENT_CACHED_LOCAL_ROOT'
                and produced['status'] == 'CANDIDATE_CERTIFIED_CACHED_QUADRATURE_LOCAL_ROOT'
                and all(fresh[k] == independent[k] for k in ('field', 'field_sha256', 'record_id'))
                and links[produced['cache_packet']] == produced['cache_packet_sha256']
                and links[f'{CACHE}/actual/{pid}/results.json'] == produced['results_sha256'],
                'mixed reviewed pilot field, fresh cache or worker result')
        per_path = f'{REPLAY}/{pid}.json'
        link(links, per_path, independent['receipt_sha256'])
        link(links, fresh['field'], fresh['field_sha256'])
        link(links, fresh['verification'], fresh['verification_sha256'])
        data[pid] = (fresh, produced, independent, per_path)
    for name, h in local_source_closure().items():
        link(links, name, h)
    return dict(sorted(links.items())), data


SHAPES = {'reference': (325, 2), 'free': (611,), 'u': (611,), 'triangles': (144, 6),
          'shape': (64, 6), 'shape_grad': (144, 64, 6, 2), 'qr': (144, 64),
          'volume_weights': (144, 64), 'active': (144, 1), 'energy_scale': (1,),
          'preconditioner': (611, 611)}
for _face in (0, 1):
    for _name, _shape in {'edges': (12, 3), 'N': (12, 3), 'dN': (12, 3),
                          'r0': (12, 12), 'l0': (12,), 'weights': (12, 12)}.items():
        SHAPES[f'surface_{_face}_{_name}'] = _shape


def validate_cache_shapes(a):
    require(set(a) == set(SHAPES) and all(a[k].shape == s for k, s in SHAPES.items()),
            'exact reviewed 325/144/611 Q8/S12 cache shapes required')
    indices = ('free', 'triangles', 'surface_0_edges', 'surface_1_edges')
    require(all(a[k].dtype.kind in ('i', 'u') for k in indices)
            and all(a[k].dtype.kind == 'f' and a[k].dtype.itemsize == 8 and np.isfinite(a[k]).all()
                    for k in a if k not in indices),
            'integer topology and finite binary64 cache tensors required')


def exact_parameters(packet, pilot):
    parameters = packet['parameters']
    drive = CP.exact_drive(pilot['step'], pilot['material'])
    require(parameters == pilot['material']['parameters']
            and packet['material_parameter_storage_exact'] == pilot['material']['parameter_storage_exact']
            and all(packet[k] == v for k, v in drive.items())
            and T.rational(parameters['mu']) == 1 and T.rational(parameters['bulk']) == 5000,
            'fresh candidate material and exact nonrounded physical drive required')
    base = {k: T.rational(parameters[k]) for k in T.ORDER[:-1]}
    base['common_gain'] = exact(drive['common_gain_exact'])
    require(all(v >= 0 for v in base.values()) and base['common_gain'] > 0
            and base['apical_area_modulus'] > 0 and base['basal_area_modulus'] > 0,
            'nonnegative tensions and positive material moduli and gain required')
    dose = exact(drive['optical_dose_exact'])
    strength = exact(drive['activation_strength_exact'])
    require(dose * base['common_gain'] == strength != exact(drive['activation_strength_binary_exact']),
            'exact candidate dyadic-coefficient/nondyadic-drive distinction required')
    return base, dose, strength


def verify_state(target, state, packet, pilot, base, drive):
    a = target.a
    x, tri = a['reference'], a['triangles']
    params = packet['parameters']
    require(state['parameters'] == params and state['integration_policy'] == packet['integration_policy']
            and state['reached_target'] is True and state['converged'] is True
            and T.rational(state['activation_strength']) == exact(packet['activation_strength_binary_exact'])
            and np.array_equal(state['reference_nodes'], x) and np.array_equal(state['triangles'], tri)
            and np.array_equal(state['free_displacement'], a['u']),
            'saved primary field is not the precise frozen center')
    rim = T.rational(params['outer_radius'])
    height = T.rational(params['height'])
    patch_r, zlo, zhi = [T.rational(params[k]) for k in ('patch_radius', 'activation_z_min', 'activation_z_max')]
    require(T.rational(params['reference_curvature']) == 0 and height == 1 and rim == 3
            and (params['nr'], params['nz']) == (12, 6) and 0 < patch_r < rim and 0 <= zlo < zhi <= height
            and T.rational(state['reference_volume']) == target.energy_scale
            == exact(packet['normalization_binary_exact']) > 0
            and float(a['volume_weights'].sum()) == float(target.energy_scale),
            'flat physical geometry or positive binary energy normalization drift')
    exact_x = [[T.rational(v) for v in row] for row in x]
    desired_free = [2*n+k for n,(r,_) in enumerate(exact_x) for k in (0,1)
                    if r != rim and (r != 0 or k == 1)]
    require(desired_free == a['free'].tolist() == pilot['free_map']
            and len(desired_free) == RANK and all(0 <= r <= rim and 0 <= z <= height for r,z in exact_x),
            'actual ordered fixed-axis/rim mesh projection required')
    mask = []
    for corners in tri[:, :3]:
        plane = [exact_x[i] for i in corners]
        rc, zc = [sum(p[k] for p in plane)/3 for k in (0, 1)]
        active = int(rc < patch_r and zlo < zc < zhi)
        if active:
            require(all(r <= patch_r and zlo <= z <= zhi for r,z in plane),
                    'an active material triangle crosses the frozen patch or depth boundary')
        mask.append(active)
    require(mask == a['active'][:, 0].tolist() == state['activation_elements']
            and [i for i, hit in enumerate(mask) if hit] == pilot['active_elements'] and sum(mask) == 24
            and base['common_gain'] * exact(packet['optical_dose_exact']) == drive,
            'the geometric 24-element mask or requested candidate activation drifted')
    centers = target.nodes()
    require(len(state['current_nodes']) == 325 and all(len(row) == 2 for row in state['current_nodes'])
            and all(T.rational(float(c)) == T.rational(v) for r,t in zip(centers,state['current_nodes'],strict=True)
                    for c,v in zip(r,t,strict=True)),
            'serialized current nodes must be binary64 rounding of X + injection(u)')
    gap = max(abs(c - T.rational(v)) for r,t in zip(centers,state['current_nodes'],strict=True)
              for c,v in zip(r,t,strict=True))
    require(gap == exact(pilot['xi_dimensionless_exact']), 'reviewed binary addition gap drift')
    for face, z in enumerate((Q(0), height)):
        edges = a[f'surface_{face}_edges'].tolist()
        require(all(len(set(e)) == 3 and min(e) >= 0 and max(e) < 325 for e in edges)
                and all(a[-1] == b[0] for a,b in zip(edges,edges[1:]))
                and exact_x[edges[0][0]] == [0,z] and exact_x[edges[-1][-1]] == [rim,z]
                and all(exact_x[j][1] == z for e in edges for j in e)
                and all(exact_x[e[0]][0] < exact_x[e[1]][0] < exact_x[e[2]][0] for e in edges),
                'both ordered conforming horizontal P2 material trace maps required')


def objects(pilot_id, data):
    require(pilot_id in PILOTS, 'singleton pilot ID required')
    pilot, result, independent, proof_path = data[pilot_id]
    name = f'transfer_exact_{pilot_id}'
    packet = pinned(ROOT / result['cache_packet'], result['cache_packet_sha256'])
    proof = pinned(ROOT / proof_path, independent['receipt_sha256'])
    state = pinned(ROOT / pilot['field'], pilot['field_sha256'])
    require(packet['status'] == 'FROZEN_FRESH_EXACT_BINARY_QUADRATURE_CACHE'
            and packet['pilot_id'] == proof['pilot_id'] == pilot_id
            and packet['record_id'] == proof['record_id'] == pilot['record_id']
            and packet['input_state'] == proof['field'] == pilot['field']
            and packet['input_state_sha256'] == proof['field_sha256'] == pilot['field_sha256']
            and packet['field_verification_sha256'] == pilot['verification_sha256']
            and packet['free_dof'] == RANK and packet['integration_policy']['volume_duffy_gauss_order'] == 8
            and packet['integration_policy']['surface_gauss_order'] == 12
            and packet['frozen_binary_cache']['sha256'] == proof['cache_array_contents_sha256']
            and proof['status'] == 'PASS_INDEPENDENT_CACHED_LOCAL_ROOT'
            and proof['rounded_strength_used_in_proof'] is False
            and proof['cache_packet_sha256'] == result['cache_packet_sha256']
            and proof['site_counts'] == CP.TERMS and exact(proof['radius_exact']) == RADII[0]
            and 0 <= exact(proof['independent_posterior_error_upper_exact']) < RADII[0] < RADII[1],
            'same independent current cache, source field, exact drive and narrow cached proof required')
    base, dose, drive = exact_parameters(packet, pilot)
    require(drive == exact(proof['exact_strength']), 'native independently verified exact point strength differs')
    target = T.Target(packet)
    validate_cache_shapes(target.a)
    require(packet['frozen_binary_cache']['array_shapes'] == {k: list(v) for k,v in SHAPES.items()},
            'embedded cache shape metadata drift')
    verify_state(target, state, packet, pilot, base, drive)
    cell = {k: [str(base[k]), str(base[k])] for k in T.ORDER}
    cell_hash = E._digest(cell)
    contract = {'status': 'FROZEN_FRESH_SINGLETON_EXACT_SPATIAL_FE_TARGET',
                'record_id': name, 'pilot_id': pilot_id, 'physical_source_record_id': pilot['record_id'],
                'source_field': pilot['field'], 'source_field_sha256': pilot['field_sha256'],
                'source_geometry_verification': pilot['verification'],
                'source_geometry_verification_sha256': pilot['verification_sha256'],
                'point_cache': result['cache_packet'], 'point_cache_sha256': result['cache_packet_sha256'],
                'frozen_binary_cache_sha256': packet['frozen_binary_cache']['sha256'],
                'cached_root_independent_receipt': proof_path, 'cached_root_independent_receipt_sha256': independent['receipt_sha256'],
                'binary_preconditioner_sha256': array_sha(target.a['preconditioner']),
                'cached_root_radius_exact': proof['radius_exact'],
                'cached_root_error_exact': proof['independent_posterior_error_upper_exact'],
                'observer_design_goal_exact': pilot['e_goal_dimensionless_exact'],
                'parameter_cell': cell, 'parameter_cell_sha256': cell_hash,
                'exact_requested_drive': CP.exact_drive(pilot['step'], pilot['material']),
                'material_parameters_binary_exact': packet['material_parameter_storage_exact'],
                'energy_scale_binary_exact': str(target.energy_scale),
                'energy_scale_convention': 'Exact 2*pi divided by this positive frozen binary normalization in spatial FE.',
                'root_center': 'Exact fraction X+inj(u) in the 611 selected coordinates.',
                'geometry_radii_exact': list(map(str,RADII)), 'scope': SCOPE,
                'exact_spatial_root_asserted': False, 'larger_radius_cached_derivative_asserted': False}
    S = np.zeros((RANK,5), dtype='<f8')
    pred = {'status': 'ZERO_AFFINE_PREDICTOR_FOR_FIXED_SINGLETON',
            'record_id': name, 'parameter_cell_sha256': cell_hash, 'shape': [RANK,5],
            'little_endian_zero_binary64_sha256': array_sha(S), 'identically_zero': True}
    return target, S, base, dose, drive, contract, pred


def make_context(parts, paths):
    target, S, base, dose, drive, contract, prediction = parts
    S.flags.writeable = False
    for a in target.a.values():
        a.flags.writeable = False
    cell = {k: (base[k],base[k]) for k in T.ORDER}
    row = {'record_id': contract['record_id'], 'source_record_id': contract['physical_source_record_id'],
           'dose_exact': str(dose), 'parameter_order': list(T.ORDER),
           'parameter_cell': contract['parameter_cell'],
           'parameter_cell_sha256': contract['parameter_cell_sha256'],
           'base_parameters_exact': {k: str(Q(v)) for k,v in base.items()},
           'point_cache_file': contract['point_cache'], 'point_cache_sha256': contract['point_cache_sha256'],
           'frozen_binary_cache_sha256': contract['frozen_binary_cache_sha256'],
           'model_observation_contract_sha256': paths['model_contract_sha256'],
           'predictor_sha256': paths['predictor_sha256']}
    ctx = SingletonContext(contract['pilot_id'], contract['record_id'], row, target, S, cell, base, dose, drive)
    check_context(ctx)
    return ctx


def check_context(ctx):
    require(type(ctx) is SingletonContext and ctx.pilot_id in PILOTS
            and ctx.record_id == f'transfer_exact_{ctx.pilot_id}'
            and set(ctx.cell) == set(ctx.base) == set(T.ORDER)
            and all(ctx.cell[k] == (ctx.base[k], ctx.base[k]) for k in T.ORDER)
            and ctx.S.shape == (RANK, 5) and np.all(ctx.S == 0)
            and ctx.dose * ctx.base['common_gain'] == ctx.drive
            and exact(ctx.target.packet['activation_strength_exact']) == ctx.drive,
            'authority must be a fresh reviewed singleton with exact requested drive')
    return E._context_identity(ctx)


def exact_node_box(ctx, radius):
    identity = check_context(ctx)
    r = exact(radius)
    require(r in RADII, 'only the two frozen positive dyadic geometry radii are declared')
    X = [list(map(T.rational,row)) for row in ctx.target.a['reference']]
    Y = [row.copy() for row in X]
    rad = [[Q(0),Q(0)] for _ in X]
    for k,u in zip(ctx.free,ctx.u,strict=True):
        n, comp = divmod(int(k),2)
        Y[n][comp] += T.rational(u)
        rad[n][comp] = r
    require(Y == ctx.target.nodes() and len(Y) == 325 and sum(v == r for row in rad for v in row) == RANK,
            'all 650 explicit coordinate hulls must match the frozen exact injection and free ball')
    return doc({**identity, 'reference_nodes_exact': X, 'current_centers_exact': Y,
                'coordinate_radii_exact': rad, 'coordinate_hulls_exact': [
                    [[v-e,v+e] for v,e in zip(row,errs,strict=True)] for row,errs in zip(Y,rad,strict=True)],
                'predictor_hulls_exact': [[Q(0),Q(0)] for _ in range(RANK)],
                'triangles': ctx.target.a['triangles'].tolist(), 'free': ctx.free.tolist(),
                'corrector_radius_exact': str(r), 'root_membership_assumed': False,
                'any_spatial_derivative_bound_assumed': False})


def routes(ctx, rule):
    check_context(ctx)
    current = ctx.target.nodes()
    bulk, surface = [], []
    for i in range(144):
        # Test both the canonical metadata and the corrected physical-axis routing.
        canonical = ctx.target.element(i, nodes=current)
        fields, ids = E._bulk_fields(ctx.target, current, i)
        chart = fields['chart']['kind']
        require(canonical['global_indices'] == fields['global_indices']
                and canonical['free_local'] == fields['free_local']
                and canonical['active'] == fields['active'], 'physical chart changed local projection')
        bulk.append({'element': i, 'chart': chart, 'nodes': ids, 'free_local': fields['free_local'],
                     'global_indices': fields['global_indices'], 'active': str(fields['active']),
                     'auxiliary_element_rule_sha256': E._bulk_rule(ctx.target,rule,i)})
    for f in (0,1):
        for e in range(12):
            fields = ctx.target.surface(f,e,nodes=current)
            E._surface_rule(ctx.target,rule,f,e)
            surface.append({'face': f, 'edge': e, 'nodes': ctx.target.a[f'surface_{f}_edges'][e].tolist(),
                            'free_local': fields['free_local'], 'global_indices': fields['global_indices'],
                            'Rlo_exact': str(fields['Rlo']), 'Rhi_exact': str(fields['Rhi']),
                            'trace_rule_sha256': E._digest(doc([rule['surface_s'][f],
                                                              rule['surface_weights'][f][e]]))})
    counts = {k: sum(r['chart'] == k for r in bulk) for k in CHARTS}
    require(counts == CHARTS and sum(row['active'] == '1' for row in bulk) == 24
            and {i for row in bulk for i in row['global_indices']} == set(range(RANK)),
            'complete physical three-chart partition and 611 free projection required')
    return {'status': 'FROZEN_LOCAL_CHART_RULE_AND_PROJECTION_IDENTITIES_NO_INTEGRALS',
            'identity': check_context(ctx), 'bulk': bulk, 'surface': surface, 'bulk_chart_counts': counts}


@contextmanager
def no_integrals():
    with ExitStack() as stack:
        methods = ((E,'bulk_component'),(E,'surface_component'),(E,'prepare_bulk_selector'),
                   (E.U,'uniform_cache'),(E.SUR,'run_edge'),(E.SUR,'cache_edge_correction'),
                   (T.core,'assemble_point'),(T.core,'local_terms'))
        for module, method in methods:
            stack.enter_context(patch.object(module,method,side_effect=AssertionError(
                f'force, moment or cached-AD work forbidden during target/geometry phases: {method}')))
        yield


def receipts(identifier, linked, fresh):
    folder = PREP / identifier
    paths = {}
    for label, filename, content in fresh:
        path = folder / filename
        h = save(path, content)
        paths[label] = relative(path)
        paths[f'{label}_sha256'] = h
    return {**linked, **paths}


def reject(action):
    try:
        action()
    except ValueError:
        return True
    return False


def check_fixtures(ctx, entries, data):
    from dataclasses import replace
    nonsingle = dict(ctx.cell, common_gain=(ctx.base['common_gain'],ctx.base['common_gain']+1))
    pilot = data[ctx.pilot_id][0]
    packet = ctx.target.packet
    bad = dict(ctx.target.a)
    bad['preconditioner'] = np.zeros((610,611))
    bad_drive = dict(packet,activation_strength_exact=packet['activation_strength_binary_exact'])
    other = next(i for i in PILOTS if i != ctx.pilot_id)
    return {
        'all_four_full_650_coordinate_hulls_and_rules_created': len(entries) == 4 and len(PILOTS) == 4,
        'reject_bool_nonnumeric_and_unprepared_radius': all(
            reject(lambda: exact_node_box(ctx,r)) for r in (True,0,'NaN',Q(1,2**35))),
        'reject_old_cohort_or_nonsingleton_and_wrong_point': (
            reject(lambda: check_context(object())) and reject(lambda: check_context(replace(ctx,cell=nonsingle)))
            and reject(lambda: check_context(replace(ctx,pilot_id=other)))),
        'reject_nonzero_S_mixed_shape_or_binary_strength_as_drive': (
            reject(lambda: check_context(replace(ctx,S=np.ones((RANK,5)))))
            and reject(lambda: validate_cache_shapes(bad)) and reject(lambda: exact_parameters(bad_drive,pilot))
            and reject(lambda: exact_parameters(packet,data[other][0]))),
        'radii_are_distinct_and_only_smaller_has_cached_proof': RADII[0] * 4 == RADII[1],
    }


def prepare():
    PREP.mkdir(parents=True, exist_ok=False)
    with (PREP/'transfer_ambiguity_exact_targets_frozen.py').open('xb') as handle:
        handle.write(SOURCE.read_bytes())
    try:
        dependencies, data = authority()
        entries = []
        with (no_integrals(),
              patch.object(G,'certify_box',side_effect=AssertionError('actual geometry forbidden before separate gate')),
              patch.object(E.SUR,'uniform_G_bounds',side_effect=AssertionError('actual trace geometry forbidden before gate'))):
            for pid in PILOTS:
                parts = objects(pid,data)
                bound = receipts(pid, {'pilot_id': pid}, (
                    ('model_contract','model_observation_contract.json',parts[5]), ('predictor','predictor.json',parts[6])))
                ctx = make_context(parts,bound)
                rule = ctx.target.auxiliary_rule()
                route = routes(ctx,rule)
                bound.update(receipts(pid, {}, (('row','row.json',ctx.row),('rule','auxiliary_rule.json',doc(rule)),
                    ('routes','routes.json',route),
                    *[(f'box_{i}',f'coordinate_box_{i}.json',exact_node_box(ctx,r)) for i,r in enumerate(RADII)])))
                entries.append(bound)
            checks = check_fixtures(ctx,entries,data)
        checks['no_actual_geometry_or_force_output'] = not GEOMETRY.exists()
        require(sweep.all_true(checks), 'singleton construction checks failed')
        manifest_name = PREP/'target_manifest.json'
        manifest_h = save(manifest_name, {'status': 'FOUR_SINGLETON_TARGETS_PREPARED',
                       'pilot_order': list(PILOTS), 'radii_exact': list(map(str,RADII)), 'targets': entries})
        protocol_h = save(PREP/'protocol.json', {
            'status': 'FROZEN_FRESH_SINGLETON_TARGETS_GEOMETRY_REVIEW_PENDING', 'scope': SCOPE,
            'source_sha256': sha(SOURCE), 'input_hashes': dependencies,
            'target_manifest_sha256': manifest_h, 'target_manifest': relative(manifest_name),
            'precision_bits': 640, 'radii_exact': list(map(str,RADII)),
            'geometry_phase': 'New whole-body and trace-domain certifications on all eight independently declared boxes.',
            'large_radius_policy': 'The larger radius has no new cached or exact spatial derivative premise.',
            'further_component_launch_gate': 'Require independent actual geometry admission and a separately sealed integration protocol.',
            'new_spatial_integrals': 0, 'new_cached_AD_assemblies': 0})
        preflight_h = save(PREP/'preflight.json', {'status': 'PASS', 'checks': checks, 'source_sha256': sha(SOURCE),
                       'protocol_sha256': protocol_h, 'target_manifest_sha256': manifest_h,
                       'actual_geometry_jobs': 0, 'new_spatial_integrals': 0, 'new_cached_AD_assemblies': 0})
        print(json.dumps({'source_sha256': sha(SOURCE),'protocol_sha256': protocol_h,
                          'target_manifest_sha256':manifest_h,'preflight_sha256':preflight_h}),flush=True)
    except BaseException as error:
        save(PREP/'failure.json',{'error':str(error),'traceback':traceback.format_exc()})
        raise


def prepared():
    dependencies, data = authority()
    protocol = read(PREP/'protocol.json')
    require(protocol['source_sha256'] == sha(SOURCE)
            == sha(PREP/'transfer_ambiguity_exact_targets_frozen.py')
            and protocol['input_hashes'] == dependencies, 'current frozen singleton dependencies differ')
    manifest = pinned(ROOT/protocol['target_manifest'],protocol['target_manifest_sha256'])
    require(manifest['status'] == 'FOUR_SINGLETON_TARGETS_PREPARED'
            and tuple(manifest['pilot_order']) == tuple(row['pilot_id'] for row in manifest['targets']) == PILOTS
            and list(map(exact,manifest['radii_exact'])) == list(RADII),
            'four unchanged explicit singleton identities and two exact radii required')
    return protocol, manifest, data


def open_context(pilot_id, manifest, data):
    require(pilot_id in PILOTS,'unknown singleton pilot')
    entry, = [r for r in manifest['targets'] if r['pilot_id'] == pilot_id]
    parts = objects(pilot_id,data)
    for name, value in (('model_contract',parts[5]),('predictor',parts[6])):
        require(pinned(ROOT/entry[name],entry[name+'_sha256']) == value,'generated point authority drift')
    ctx = make_context(parts,entry)
    rule = ctx.target.auxiliary_rule()
    require(pinned(ROOT/entry['row'],entry['row_sha256']) == ctx.row
            and pinned(ROOT/entry['rule'],entry['rule_sha256']) == doc(rule),
            'frozen explicit point metadata or auxiliary rule drift')
    return ctx, rule, entry


def load_context(pilot_id):
    """Open a source-bound singleton context. This does not admit components."""
    _, manifest, data = prepared()
    return open_context(pilot_id,manifest,data)[0]


def geometry_box(ctx, rule, entry, i, folder):
    r = RADII[i]
    key = f'box_{i}'
    box = pinned(ROOT/entry[key],entry[key+'_sha256'])
    require(box == exact_node_box(ctx,r),'geometry must use this exact current singleton box')
    route = pinned(ROOT/entry['routes'],entry['routes_sha256'])
    require(route == routes(ctx,rule),'record-specific frozen chart/trace proof route changed')
    body = G.certify_box(box['reference_nodes_exact'],box['current_centers_exact'],
                         box['triangles'],box['coordinate_radii_exact'])
    body_h = save(folder/'body_geometry.json',body)
    status = 'UNRESOLVED_FRESH_BODY_GEOMETRY'
    traces = []
    good = body['status'] == 'CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING'
    if good:
        nodes = E.SUR.common_nodes(ctx.target,box)
        for row in route['bulk']:
            fields, ids = E._bulk_fields(ctx.target,nodes,row['element'])
            require(ids == row['nodes'] and fields['chart']['kind'] == row['chart']
                    and fields['global_indices'] == row['global_indices'],'whole-box chart map drift')
            E._body_rho(body['triangles'][row['element']])
        for row in route['surface']:
            f,e = row['face'],row['edge']
            try:
                fields = ctx.target.surface(f,e,nodes=nodes)
                domain = E.SUR.uniform_G_bounds(fields,max_depth=8)
                site = [E.SUR.local_site_domain(ctx.target,box,f,e,k) for k in range(12)]
                require(domain['status'] == 'CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN'
                        and fields['global_indices'] == row['global_indices'],'fresh surface domain unresolved')
                traces.append({**row,'domain':domain,'cache_site_geometry':site})
            except ValueError as error:
                traces.append({**row,'status':'UNRESOLVED_FRESH_SURFACE_DOMAIN','error':str(error)})
        status = 'PASS_FRESH_GEOMETRY_ONLY_UNREVIEWED' if all('domain' in t for t in traces) else 'UNRESOLVED_FRESH_SURFACE_DOMAIN'
    return {'status': status,'pilot_id':ctx.pilot_id,'record_id':ctx.record_id,'radius_exact':str(r),
            'box_sha256':entry[key+'_sha256'],'routes_sha256':entry['routes_sha256'],'rule_sha256':entry['rule_sha256'],
            'model_contract_sha256':entry['model_contract_sha256'],'body_geometry_sha256':body_h,
            'trace_geometry_only':traces,'local_component_selector_frozen':False,
            'new_force_or_hessian_work':False,'new_cached_AD_assembly':False,'root_asserted':False}


def geometry(review_sha):
    protocol, manifest, data = prepared()
    expected = {'source_sha256':sha(SOURCE),'protocol_sha256':sha(PREP/'protocol.json'),
                'preflight_sha256':sha(PREP/'preflight.json'),
                'target_manifest_sha256':protocol['target_manifest_sha256'],
                'cached_root_independent_review_sha256':PINS[CACHE_REVIEW]}
    accepted = pinned(PREP/'prelaunch_independent_review.json',review_sha)
    require(accepted.get('status') == 'PASS' and sweep.all_true(accepted.get('checks'))
            and all(accepted.get(k) == v for k,v in expected.items()),'fresh independent geometry prelaunch gate required')
    pf = pinned(PREP/'preflight.json',expected['preflight_sha256'])
    require(pf['status']=='PASS' and sweep.all_true(pf['checks'])
            and all(pf[k]==expected[k] for k in ('source_sha256','protocol_sha256','target_manifest_sha256'))
            and pf['actual_geometry_jobs']==pf['new_spatial_integrals']==pf['new_cached_AD_assemblies']==0,
            'preflight differs from the frozen cold geometry plan')
    GEOMETRY.mkdir(exist_ok=False)
    start = monotonic()
    binding = {**expected,'prelaunch_review_sha256':review_sha}
    save(GEOMETRY/'started.json',{'pid':os.getpid(),'bindings':binding})
    try:
        rows=[]
        with no_integrals():
            for pid in PILOTS:
                ctx,rule,entry = open_context(pid,manifest,data)
                for i in range(len(RADII)):
                    folder = GEOMETRY/pid/f'box_{i}'
                    result = geometry_box(ctx,rule,entry,i,folder)
                    h = save(folder/'geometry.json',result)
                    rows.append({'pilot_id':pid,'status':result['status'],'file':relative(folder/'geometry.json'),'sha256':h})
        require(len(rows)==8,'all eight geometry cases required')
        results_sha = save(GEOMETRY/'results.json',{'status':'COMPLETE_FRESH_GEOMETRY_ATTEMPTS_UNREVIEWED',
                           'bindings':binding,'geometry_cases':rows,'seconds':monotonic()-start,
                           'independent_actual_geometry_review_required':True,'force_or_H_jobs_authorized':False})
        save(GEOMETRY/'terminal.json',{'status':'FINISHED','pid':os.getpid(),'exit_code':0,
                                       'results_sha256':results_sha,'bindings':binding})
        print(json.dumps({'geometry_results_sha256':results_sha,'geometry_cases':len(rows)}),flush=True)
    except BaseException as e:
        save(GEOMETRY/'failure.json',{'error':str(e),'traceback':traceback.format_exc(),'bindings':binding})
        save(GEOMETRY/'terminal.json',{'status':'FAILED','pid':os.getpid(),'exit_code':1,'bindings':binding})
        raise


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    stages = cli.add_subparsers(dest='stage',required=True)
    stages.add_parser('prepare')
    gate = stages.add_parser('geometry')
    gate.add_argument('--review-sha256',required=True)
    args = cli.parse_args()
    if args.stage == 'prepare':
        prepare()
    else:
        geometry(args.review_sha256)

