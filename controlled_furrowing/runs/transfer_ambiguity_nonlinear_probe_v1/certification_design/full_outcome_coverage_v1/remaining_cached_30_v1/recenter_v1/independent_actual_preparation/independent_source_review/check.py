"""Cold peer checks. No actual cache arrays, step vectors, or mechanics jobs."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import sys
from unittest.mock import patch

HERE = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("reviewed_replay", HERE / "replay.py")
m = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = m
spec.loader.exec_module(m)
p_sha = "3361790dace50108802693832741382d6d6b2fc2da0b74a6de9ec66e2615d450"
f_sha = "e0b7c3426bc995e9808284e6c501723c29f210cdc29cc9b981dec4418bdfcdc5"
s_sha = "15632d9f96a94f8454abbd59774448ce8a8d0da890af28d51b4b71d868173642"
a_sha = "d5a40643b35b9e38c1b94939199eb26d2b88bddc3bb07a54ce06f0263273d0e0"
a_path = m.PRODUCED / "root_terminal_admission.json"
b = dict(source_sha256=s_sha, protocol_sha256=p_sha, preflight_sha256=f_sha,
         producer_root_admission_sha256=a_sha)
assert m.sha(m.SOURCE) == m.sha(HERE / "replay_frozen.py") == s_sha
p = m.pinned(HERE / "protocol.json", p_sha)
pf = m.pinned(HERE / "preflight.json", f_sha)
original_pinned = m.pinned


def metadata_only(path, *args, **kwargs):
    assert path.name not in {"exact_binary_step.json", "recentered_cache_packet.json",
                             "cache_packet.json", "field.json", "verification.json"}
    return original_pinned(path, *args, **kwargs)


with patch.object(m, "pinned", metadata_only):
    assert m.static_protocol() == p
    a, records, pids = m.admitted(a_path, a_sha, p)
    assert tuple(record[1]["case_id"] for record in records) == m.CASES
    assert len(a["sha256"]) == 17 and len(set(pids)) == 2
    assert m.rejects(lambda: m.admitted(a_path, "1" * 64, p))
assert m.yes(m.fixtures(p))
import algebraic_root_feasibility as core
import review_algebraic_root as jet
A = core.A


def contains(x, y):
    l, h = m.signed_pair(core.pack(x))
    a, b = m.signed_pair(core.pack(y))
    return l <= a <= b <= h


def overlaps(x, y):
    l, h = m.signed_pair(core.pack(x))
    a, b = m.signed_pair(core.pack(y))
    return l <= b and a <= h


def audited_chart(kind, center, direction, material):
    center, direction, material = [tuple(map(A, a)) for a in (center, direction, material)]
    closed, automatic = ((core.bulk_local, jet.bulk_AD) if kind == "bulk"
                         else (core.surface_local, jet.surface_AD))
    weight = A(Q(7, 5))
    gc, hc, jc = closed(center, *material)
    g0, h0, j0, _ = automatic(center, *material)
    assert overlaps(jc, j0)
    assert all(overlaps(x, y) for x, y in zip(gc, g0, strict=True))
    for left, right in zip(hc, h0, strict=True):
        assert all(overlaps(x, y) for x, y in zip(left, right, strict=True))
    term = dict(kind=kind, q=center, B=[[x] for x in direction], indices=[0],
                weight=weight, parameters=material)
    _, point, sites, _ = core.assemble_point(iter([term]), 1)
    radius = A(Q(1, 64))
    variation, _, _ = core.variation_rows(sites, 1, radius)
    v = m.signed_pair(core.pack(variation[0, 0]))[1]
    for u in (Q(-1,64), Q(-1,128), Q(0), Q(1,128), Q(1,64)):
        q = [x + a * A(u) for x, a in zip(center, direction, strict=True)]
        _, hessian, _, _ = automatic(q, *material)
        sample = weight * sum(a * hessian[i][j] * b for i, a in enumerate(direction)
                              for j, b in enumerate(direction))
        assert contains(core.arb(0, A(v).upper()), sample - point[0, 0])


audited_chart("bulk", (Q(5,4), Q(1,8), Q(-1,4), Q(9,8), Q(7,8)),
              (Q(1,4), Q(-1,2), Q(3,4), Q(1,8), Q(-1,4)), (Q(3,2), Q(7,3), Q(1,5)))
audited_chart("surface", (Q(5,4), Q(1,8), Q(-3,4)), (Q(-1,2), Q(1,8), Q(3,4)),
              (Q(3,2), Q(2,3), Q(5,7)))
for target in (Q(1) + Q(1,2**53), -Q(1) - Q(3,2**53), Q(2)**-1075,
               Q(3)*Q(2)**-1075, Q(-3)*Q(2)**-1075):
    assert Q.from_float(m.independent_round(target)) == Q.from_float(float(target))
d = m.classify(Q(1,16), Q(1,4), Q(1,8), Q(1,4), Q(1,128), Q(1,8))
assert Q(d["physical_nodal_coordinate_error_um_exact"]) == Q(67,2) * (Q(1,8)+Q(1,12)+Q(1,128))
assert not m.classify(Q(1,16),Q(1,4),Q(1,8),Q(5,24),Q(0),Q(1,8))["strict_old_observer_goal_pass"]
checks = dict(
    current_replay_snapshot_protocol_and_preflight_match=True,
    real_external_terminal_metadata_only_two_cases_17_artifacts=True,
    archived_operator_tensor_and_exact_singleton_drive_identity_reviewed=True,
    nearest_even_signed_original_center_and_subnormal_ties=True,
    independent_bulk_surface_gradient_hessian_and_variation_math_reviewed=True,
    synthetic_nondiagonal_bulk_and_surface_jet_fixtures=True,
    strict_shift_and_archived_addition_gap_verified=True,
    replay_does_not_use_producer_derivative_enclosures_for_proof=True,
)
r = dict(status=m.REVIEW, **b, checks=checks, reviewed_child_pids=pids,
         scope="Cold source and synthetic mathematics. Shared frozen-array site map generator and existing closed-form mechanics kernels are pinned. Actual recentered cache and numerical proof are still subject to replay.",
         actual_producer_packets_opened=0, actual_new_center_arrays_read_or_written=0,
         actual_site_derivative_jobs=0)
m.reviewed_phase(p, pf, r, b)
print(m.save(Path(__file__).resolve().parent / "receipt.json", r))
