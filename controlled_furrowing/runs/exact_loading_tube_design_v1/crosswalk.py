"""Read-only archive and fixed-model anchor identity check, no kernel imports."""

from base64 import b64decode
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZipFile


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
B = Path('runs/transfer_ambiguity_nonlinear_probe_v1')
EX = B / 'certification_design/exact_targets_v1'
W = EX / 'remaining32_bridge_v1'
END = 'max_negative_contrast_change'
NAMES = tuple(END + '_' + code for code in ('small_16', 'small_17'))
PINNED = {
    'PLAN.md': '993240af5d30c8e5b72ab31a1278a524314d38ede48486b677b021d84bcde500',
    B / 'field_manifest.json': '308ca4c05e476ba6a14cbc7b8d34018983b8f11263e2c224488672b35eb38ac0',
    B / 'independent_actual_preparation/actual_replay/receipt.json':
        'c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd',
    W / 'first_wave_reduction_v1/independent_actual_preparation/actual_replay/receipt.json':
        'e73ea05a4703256235a46ce7b914298d3f8258e9fef029f4dbe952b04bb9798a',
    W / 'first_wave_reduction_v1/independent_actual_preparation/actual_replay/os_exit_observation.json':
        'aa3c90101c64ccefe1403af2c296bcb22cafad3b2e35d73539aa7e3c5227f281',
    EX / 'bridge_reduction_v1/independent_actual_preparation/all_four_readback.json':
        '859f507ad3762cfef44d53672bbe15a43592ce54d6e08ba8828bf50a925ec340',
}


def require(ok, why):
    if not ok:
        raise ValueError(why)


def hash_bytes(data):
    return sha256(data).hexdigest()


def text(path, hashes, pinned=None):
    p = (ROOT / path).resolve()
    name = p.relative_to(ROOT).as_posix()
    raw = p.read_bytes()
    got = hash_bytes(raw)
    require(got == (pinned or got) == hashes.setdefault(name, got), f'authority mismatch {name}')
    return raw


def pick(rows, **match):
    [row] = [r for r in rows if all(r[k] == v for k, v in match.items())]
    return row


def archive(protocol, manifest, audited, hashes):
    fields = manifest['fields']
    require(len(fields) == audited['checkpoints'] == 44
            and len({(r['endpoint_id'], r['checkpoint']) for r in fields}) == 44,
            'unique complete archive and independent receipt')
    schedules = {}
    for candidate in protocol['endpoints']:
        end = candidate['endpoint_id']
        schedules[end] = {}
        for geom in ('small', 'large'):
            rows = [r for r in fields if r['endpoint_id'] == end and r['checkpoint'].startswith(geom+'_')]
            plan = ([candidate['small_zero']] + candidate['small_records'] if geom == 'small'
                    else candidate['large_checkpoints'])
            require(len(rows) == len(plan) == (18 if geom == 'small' else 4), 'whole exact geometry schedule')
            gain = Q(candidate['materials'][geom]['gain_exact'])
            doses = []
            for j, (row, step) in enumerate(zip(rows, plan, strict=True)):
                label = 'small_zero' if (geom, j) == ('small', 0) else f'{geom}_{j:02d}'
                dose = Q(step['dose_exact'])
                require(row['checkpoint'] == label and Q(row['step']['dose_exact']) == dose
                        and Q(row['step']['requested_strength_exact']) == Q(step['requested_strength_exact']) == gain*dose
                        and Q(row['actual_dose_implied_by_binary_strength_exact'])*gain
                        == Q(row['actual_strength_binary64_exact']), 'nominal exact mapping and stored binary solver drive')
                for kind in ('field', 'verification'):
                    path = row[kind]
                    require(audited['all_hashes'][path] == row[kind+'_sha256'], 'independently reviewed archive digest')
                    text(path, hashes, row[kind+'_sha256'])
                doses.append(dose)
            require(doses[0] == 0 and all(a < b for a,b in zip(doses,doses[1:])), 'separate ascending dose schedule')
            schedules[end][geom] = {'count': len(rows), 'interval_count': len(rows)-1,
                                    'roles': dict(Counter(r['role'] for r in rows)), 'doses_exact': list(map(str,doses))}
    require(sum(g['interval_count'] for e in schedules.values() for g in e.values()) == 40, 'four disjoint schedules')
    return schedules


def cache_array_hashes(packet, contract, archive_row):
    drive = contract['exact_requested_drive']
    require(packet['input_state_sha256'] == archive_row['field_sha256']
            and packet['normalization_binary_exact'] == contract['energy_scale_binary_exact']
            and all(Q(packet[key]) == Q(drive[key]) for key in
                    ('optical_dose_exact', 'common_gain_exact', 'activation_strength_exact')),
            'cache packet binding to exact nominal load')
    cache = packet['frozen_binary_cache']
    raw = b64decode(cache['data'], validate=True)
    require(hash_bytes(raw) == cache['sha256'] == contract['frozen_binary_cache_sha256']
            and len(raw) == cache['bytes'], 'exact embedded cache digest')
    with ZipFile(BytesIO(raw)) as zipfile:
        return {name: hash_bytes(zipfile.read(name)) for name in zipfile.namelist()}


def check(hashes):
    def obj(path, pin=None):
        return json.loads(text(path, hashes, pin))

    for name, pin in PINNED.items():
        text(name, hashes, pin)
    protocol = obj(B / 'preparation/protocol.json')
    manifest = obj(B / 'field_manifest.json')
    audited = obj(B / 'independent_actual_preparation/actual_replay/receipt.json')
    require(audited['status'] == 'PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY'
            and manifest['independent_primary_review_sha256'] == hashes[(B/'independent_actual_preparation/actual_replay/receipt.json').as_posix()],
            'independent archived results')
    schedules = archive(protocol, manifest, audited, hashes)

    k,l = NAMES
    new = W / 'first_wave_reduction_v1/independent_actual_preparation/actual_replay'
    wave, os_doc = (obj(new / name) for name in ('receipt.json', 'os_exit_observation.json'))
    require(os_doc['status'] == 'PASS_ROOT_OBSERVED_FIRST_WAVE_INDEPENDENT_GLOBAL_EXIT'
            and type(os_doc['observed_exit_code']) is int and os_doc['observed_exit_code'] == 0
            and os_doc['receipt_sha256'] == hashes[(new/'receipt.json').as_posix()], 'genuine independent K completion')
    kr = pick(wave['case_results'], case_id=k)
    ka = obj(kr['native_file'], kr['native_sha256'])
    require(os_doc['native_cases_sha256'][kr['native_file']] == kr['native_sha256']
            and ka['status'] == 'PASS_INDEPENDENT_FIRST_WAVE_EXACT_SPATIAL_ROOT'
            and ka['arithmetic']['strict_selfmap'] is True
            and ka['observer_and_optional_minimum']['strict_minimum_proved'] is True, 'K exact spatial minimum')
    old = obj(EX / 'bridge_reduction_v1/independent_actual_preparation/all_four_readback.json')
    lr = pick(old['cases'], pilot_id=l)
    la = obj(lr['independent_proof_file'], lr['independent_proof_sha256'])
    lo = obj(lr['independent_os_observation_file'], lr['independent_os_observation_sha256'])
    require(old['status'] == 'PASS_FOUR_INDEPENDENT_FIXED_INPUT_EXACT_SPATIAL_ROOTS_AND_MINIMA'
            and type(lo['observed_exit_code']) is int and lo['observed_exit_code'] == 0
            and la['target_identity'] == lr['target_identity'] and lr['selfmap']['strict_selfmap'] is True
            and Q(lr['whole_box_normalized_H_min_exact']) > 0, 'L exact spatial minimum')

    kt, lt = ka['target'], lr['target_identity']
    kpath = EX / 'remaining32_targets_v1/preparation' / k / 'model_contract.json'
    lpath = EX / 'preparation' / l / 'model_observation_contract.json'
    kc, lc = obj(kpath, kt['model_contract_sha256']), obj(lpath, lt['model_contract_sha256'])
    arches = [pick(manifest['fields'], endpoint_id=END, checkpoint=c) for c in ('small_16', 'small_17')]
    ks, ls = (obj(row['field'], row['field_sha256']) for row in arches)
    mat = pick(protocol['endpoints'], endpoint_id=END)['materials']['small']
    require(kc['archived_field'] == arches[0]['field'] and lc['source_field'] == arches[1]['field']
            and kt['original_observer_field_sha256'] == arches[0]['field_sha256']
            and lt['field_sha256'] == arches[1]['field_sha256']
            and all(ks[field] == ls[field] for field in
                    ('parameters','reference_nodes','triangles','activation_elements','integration_policy'))
            and ks['parameters'] == mat['parameters'] and sum(ks['activation_elements']) == 24
            and len(ks['reference_nodes']) == 325 and len(ks['triangles']) == 144,
            'same frozen nominal geometry and elementconstant activation')
    norm = Q(kc['energy_scale_binary_exact'])
    require(norm == Q(lc['energy_scale_binary_exact']) > 0
            and kc['material_parameters_binary_exact'] == lc['material_parameters_binary_exact'] == mat['parameter_storage_exact']
            and kc['parameter_cell'] == lc['parameter_cell'] and all(lo==hi for lo,hi in kc['parameter_cell'].values()),
            'identical material singleton and positive exact energy normalization')
    gain = Q(mat['gain_exact'])
    for row, contract, target, load_key, cache_key, Rkey in zip(
            arches, (kc,lc), (kt,lt), ('exact_drive','drive_exact'), ('point_cache_sha256','cache_packet_sha256'),
            ('binary_R_sha256','binary_preconditioner_sha256'), strict=True):
        drive = contract['exact_requested_drive']
        require(Q(drive['optical_dose_exact']) == Q(row['step']['dose_exact'])
                and Q(drive['common_gain_exact']) == gain == Q(contract['parameter_cell']['common_gain'][0])
                and Q(target[load_key]) == Q(drive['activation_strength_exact']) == gain*Q(row['step']['dose_exact'])
                and Q(target['energy_scale_binary_exact']) == norm
                and target[cache_key] == contract['point_cache_sha256'] and target['R_binary_sha256'] == contract[Rkey],
                'each independent exact root uses the same force law at its exact dose')
    packets = [obj(c['point_cache'], c['point_cache_sha256']) for c in (kc,lc)]
    ac, bc = (cache_array_hashes(p,c,r) for p,c,r in zip(packets,(kc,lc),arches,strict=True))
    changed = sorted(n for n in set(ac) | set(bc) if ac.get(n) != bc.get(n))
    require(changed == ['preconditioner.npy','u.npy'] and {'reference.npy','free.npy','active.npy'} <= set(ac),
            'all frozen numerical maps, mask, weights and reference agree, except R and center')
    kroute, lroute = (obj(f.with_name('routes.json')) for f in (kpath,lpath))
    require(kroute['bulk'] == lroute['bulk'] and kroute['surface'] == lroute['surface']
            and [bool(int(r['active'])) for r in kroute['bulk']] == ks['activation_elements'],
            'same mathematical element and edge maps and masks')
    geom = obj(EX / 'remaining32_geometry_v1/preparation/protocol.json')
    trial = pick(geom['trials'], case_id=k)
    kbox = obj(trial['trial_box'], kt['trial_box_sha256'])
    lbpath = lpath.with_name('coordinate_box_1.json')
    lbox = obj(lbpath, old['current_file_sha256'][lbpath.as_posix()])
    require(la['source_and_actual_sha256'][lbpath.as_posix()] == hashes[lbpath.as_posix()]
            and Q(kbox['corrector_radius_exact']) == Q(lbox['corrector_radius_exact']) == Q(1,2**34)
            and kbox['free'] == lbox['free'] and len(kbox['free']) == 611
            and kbox['triangles'] == lbox['triangles'] == ks['triangles']
            and kbox['reference_nodes_exact'] == lbox['reference_nodes_exact'] == [[str(Q(v)) for v in r]
                                                                                 for r in ks['reference_nodes']],
            'same exact proof reference/free injection and closed radius')
    exact_source = 'scripts/exact_fe_target_safe.py'
    model = text(exact_source, hashes).decode()
    older = old['current_file_sha256'][exact_source]
    newer = obj(EX / 'remaining32_targets_v1/preparation/protocol.json')['local_source_closure_sha256'][exact_source]
    kernel = text('scripts/axisymmetric_mixed.py', hashes).decode()
    quadratic = text('scripts/axisymmetric_quadratic.py', hashes).decode()
    require(hashes[exact_source] == older == newer and 'chi*g*dose/2*(q_b^2+q_d^2)' in model
            and '(2*pi()/S_binary)' in model and '(r_R,r_Z,z_R,z_Z,r/R)' in model
            and 'fm[:,:,:,1]' in kernel and 'self.active*activation_strength/2' in kernel
            and "self.parameters.bulk/2*result['logJ']**2" in quadratic,
            'same exact target and vertical activation source')
    doses = tuple(Q(r['step']['dose_exact']) for r in arches)
    require(doses == (Q('1.28898317'),Q('1.51057344')) and
            schedules[END]['small']['doses_exact'][-2:] == list(map(str,doses)), 'exact adjacent nominal dose interval')
    return schedules, {
        'case_order': list(NAMES), 'exact_dose_interval': list(map(str,doses)), 'exact_dose_gap': str(doses[1]-doses[0]),
        'same_gain_exact': str(gain), 'exact_requested_drives': [str(gain*i) for i in doses],
        'common_binary_energy_normalization_exact': str(norm), 'exact_material_parameters': kc['material_parameters_binary_exact'],
        'exact_active_direction': 'e_Z', 'active_slots': [i for i,v in enumerate(ks['activation_elements']) if v],
        'same_binary_arrays_sha256': {n:h for n,h in ac.items() if n not in changed}, 'different_binary_arrays': changed,
        'K_spatial_proof_sha256': kr['native_sha256'], 'L_spatial_proof_sha256': lr['independent_proof_sha256'],
        'K_box_sha256': kt['trial_box_sha256'], 'L_box_sha256': hashes[lbpath.as_posix()],
        'endpoint_posterior_error_upper_exact': [ka['arithmetic']['posterior_error_upper_exact'],
                                                 lr['selfmap']['posterior_error_upper_exact']],
        'K_original_observer_shift_preserved_exact': kt['archived_center_shift_exact'],
        'joint_dose_interval_certified': False,
    }


def main():
    hashes = {}
    schedules,pilot = check(hashes)
    for doc in ('ENERGY-GAIN-CELL','ENERGY-GAIN-CELL-REVIEW','ENERGY-GAIN-GEOMETRY','VALIDATED-IDENTIFIED-SET',
                'COHORT-ENERGY-STABILITY','ROOT-BODY-GEOMETRY','ROOT-BODY-GEOMETRY-REVIEW','JOINT-OBSERVATION-CONTRACT',
                'PHYSICAL-INFERENCE-READINESS','NODAL-OBSERVER-ERROR-TRANSFER','L-CENTERED-ROOT'):
        text(f'docs/{doc}.md', hashes)
    for path in ('scripts/root_body_geometry.py','scripts/algebraic_root_feasibility.py','scripts/exact_centered_bulk_L.py',
                 'scripts/exact_centered_cohort_bulk.py','scripts/transfer_ambiguity_exact_targets.py',
                 'scripts/transfer_ambiguity_remaining32_exact_targets.py','scripts/transfer_ambiguity_remaining32_exact_bridge.py',
                 B / 'report/REPORT.md'):
        text(path, hashes)
    result = {'status': 'PASS_FROZEN_ARCHIVE_AND_ANCHOR_IDENTITY_FOR_LOADING_TUBE_DESIGN',
              'method_sha256': hash_bytes((OUT/'METHOD.md').read_bytes()),
              'script_sha256': hash_bytes(Path(__file__).read_bytes()), 'schedules': schedules,
              'archived_checkpoint_count': 44, 'archived_intercheckpoint_intervals':40, 'pilot':pilot,
              'input_sha256_by_study_relative_path':dict(sorted(hashes.items())),
              'numerical_kernel_invocations':0, 'new_geometry_or_force_AD_integrations_or_solves':0,
              'joint_path_certified':False, 'full_five_gate_complete':False}
    with (OUT/'receipt.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(json.dumps({'receipt_sha256':hash_bytes((OUT/'receipt.json').read_bytes()),
                      'input_count':len(hashes), 'joint_path_certified':False}))


if __name__ == '__main__':
    main()
