"""Read-only validation and visualization of the unresolved L enclosure.

This composes existing reviewed artifacts; it does not evaluate a force or
assert a root. Exact rational arithmetic decides the displayed predicates.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'runs/exact_integration_cohort_v1'
PINS = {
    'records/optoGEF2_small_col_L/original_cube_v1/assembled.json':
        '13e06c3f23c0e840a9b3043f8b09d40757cff5b84831f728c4b413ae9a8bc438',
    'reduction_preparation/optoGEF2_small_col_L_reduction_review.json':
        '5a55e0d139edb39206ab3999681d9858d6b0162f73c7e5103d7a5bed3ebc4514',
    'reduction_preparation/optoGEF2_small_col_L_signed_endpoint_diagnostics.json':
        'd251d4b797877e2a5d588b33f33483a7fe6bc8f3af2c1b27d4a7358590ebc619',
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def write_once(path, data):
    with path.open('x', encoding='utf8', newline='\n') as stream:
        json.dump(data, stream, indent=2, default=str, allow_nan=False)
        stream.write('\n')


def main():
    bindings = {}
    docs = []
    for relative, expected in PINS.items():
        path = BASE / relative
        assert digest(path) == expected, relative
        bindings[path.relative_to(ROOT).as_posix()] = expected
        docs.append(json.loads(path.read_bytes()))
    producer, review, diagnostic = docs
    # Existing independent endpoint checks remain shared premises. Revalidate
    # every artifact identity instead of relying on a stale PASS string.
    for relative, expected in review['source_input_and_component_hashes'].items():
        assert digest(ROOT / relative) == expected, relative
        bindings[relative] = expected
    assert review['status'] == 'PASS_REDUCTION_UNRESOLVED_SUFFICIENT_ROOT_TEST'
    assert review['record_id'] == producer['record_id'] == diagnostic['record_id']
    assert review['component_count'] == 168 and review['force_entries'] == 611
    assert review['upper_Hessian_entries'] == 6676
    eta0, z0, rf, b, r = [F(review[key]) for key in (
        'cached_eta_upper_exact', 'cached_Z_upper_exact', 'RF_uniform_upper_exact',
        'b_upper_exact', 'normal_radius_exact')]
    eta, z = eta0 + rf, z0 + b
    lhs = eta + z * r
    hypothetical = eta0 + z * r
    required_b_if_rf_zero = 1 - z0 - eta0 / r
    assert eta == F(review['eta_upper_exact'])
    assert z == F(review['z_upper_exact'])
    assert lhs == F(review['selfmap_lhs_exact'])
    assert 0 <= z < 1 and lhs > r and hypothetical > r
    assert required_b_if_rf_zero < b
    assert review['sufficient_conditions_pass'] is False
    assert review['root_error_upper_exact'] is None and review['energy'] is None
    assert diagnostic['RF_endpoint_sha256'] == review['RF_endpoint_sha256']
    assert F(diagnostic['certified_RF_infinity_norm_lower_exact']) == 0
    assert diagnostic['coordinates_strictly_excluding_zero'] == 0
    assert len(producer['preconditioned_force_defect']) == 611
    assert all(F(x['lower_exact']) <= 0 <= F(x['upper_exact'])
               for x in producer['preconditioned_force_defect'])
    checks = {
        'every_independent_input_and_component_digest_is_current': True,
        'complete_168_packet_611_force_6676_upper_Hessian_scope': True,
        'independent_exact_eta_z_and_selfmap_arithmetic': True,
        'contraction_below_one_but_selfmap_fails': True,
        'zero_force_hypothesis_still_fails_at_current_b': True,
        'all_611_producer_signed_defect_intervals_contain_zero': True,
        'independent_signed_endpoint_digest_and_zero_lower_bound': True,
        'no_root_error_or_energy_conclusion_returned': True,
    }
    out = BASE / 'checkpoint_L_v1'
    out.mkdir(exist_ok=True)
    image = out / 'root_test_limits.png'
    assert not image.exists(), 'write-once figure already exists'
    with plt.style.context('dark_background'):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
        fig.patch.set_facecolor('#0c1420')
        for ax in axes:
            ax.set_facecolor('#111e2d')
            ax.spines[['top', 'right']].set_visible(False)
            ax.grid(axis='x', alpha=.16)
        axes[0].barh(['Contraction z'], [float(z)], color='#63c8ba')
        axes[0].axvline(1, color='white', ls='--', lw=1)
        axes[0].set_xlim(0, 1.15)
        axes[0].set_title('Contraction condition passes')
        axes[0].set_xlabel('Required: z < 1')
        axes[0].text(float(z) + .03, 0, f'{float(z):.6f}', va='center')
        vals = [float(lhs/r), float(hypothetical/r)]
        labels = ['Current bound / radius', 'If force defect were zero']
        axes[1].barh(labels, vals, color=['#e3a06b', '#9faec4'])
        axes[1].set_xscale('log')
        axes[1].axvline(1, color='white', ls='--', lw=1)
        axes[1].set_xlim(.6, 50000)
        axes[1].set_title('Self-map condition remains unresolved')
        axes[1].set_xlabel('Required: bound / radius < 1 (log scale)')
        for i, value in enumerate(vals):
            axes[1].text(value * 1.15, i, f'{value:.5g}', va='center')
        fig.suptitle('Record L: limits of this enclosure, not failure of equilibrium')
        fig.savefig(image, dpi=160, facecolor=fig.get_facecolor())
        plt.close(fig)
    bindings[Path(__file__).resolve().relative_to(ROOT).as_posix()] = digest(Path(__file__))
    receipt = {
        'status': 'PASS_COMPOSITION_OF_UNRESOLVED_L_ENCLOSURE',
        'checks': checks, 'source_and_input_hashes': bindings,
        'values_exact': {'eta0': eta0, 'z0': z0, 'rf': rf, 'b': b, 'r': r,
                         'eta': eta, 'z': z, 'lhs': lhs,
                         'lhs_over_r': lhs/r, 'zero_rf_lhs': hypothetical,
                         'zero_rf_lhs_over_r': hypothetical/r,
                         'necessary_b_upper_for_strict_test_at_zero_rf': required_b_if_rf_zero},
        'figure': image.relative_to(ROOT).as_posix(), 'figure_sha256': digest(image),
        'arithmetic': 'Exact rational predicates; floating values used only to plot.',
        'scope': 'Existing reviewed local enclosure engines and endpoint reduction are premises. '
                 'No independent force integration, nonexistence theorem, physical identified set, '
                 'PDE accuracy, biological validation, or novelty claim.',
    }
    target = ROOT / 'evidence/exact_cohort_L_checkpoint.json'
    write_once(target, receipt)
    print(json.dumps({'status': receipt['status'], 'bindings': len(bindings),
                      'receipt_sha256': digest(target),
                      'figure_sha256': digest(image),
                      'lhs_over_r': float(lhs/r)}))


if __name__ == '__main__':
    main()
