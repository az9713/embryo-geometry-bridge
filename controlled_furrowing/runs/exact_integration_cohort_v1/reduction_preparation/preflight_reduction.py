"""Frozen preparation for the cohort assembly endpoint review."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import review_cohort_exact_record as E


def tests():
    checks = {}

    def check(name, value):
        if not value:
            raise AssertionError(name)
        checks[name] = True

    def reject(name, action, expected=None):
        try:
            action()
        except ValueError as error:
            check(name, expected is None or expected in str(error))
        else:
            raise AssertionError(f"{name}: accepted invalid case")

    def pack(lo, hi=None):
        return {"lower_exact": str(lo), "upper_exact": str(lo if hi is None else hi)}

    def entry(i, j, lo, hi=None):
        return {"i": i, "j": j, "enclosure": pack(lo, hi)}

    g, H = [E.ZERO] * 3, {}
    a = ([2, 0], [pack(-3, -2), pack(-1, 2)], [entry(0, 0, -4, 1), entry(0, 1, 2, 3), entry(1, 1, -1, 1)])
    b = ([1, 2], [pack(-3, 2), pack(2, 3)], [entry(0, 0, 1, 2), entry(0, 1, -1, 2), entry(1, 1, 1, 3)])
    E._scatter(g, H, *a)
    E._scatter(g, H, *b)
    check("negative and zero-crossing force endpoints scatter signed",
          g == [(Q(-1), Q(2)), (Q(-3), Q(2)), (Q(-1), Q(1))])
    check("unordered global columns store one symmetric upper entry",
          H[(0, 2)] == (Q(2), Q(3)) and H[(2, 2)] == (Q(-3), Q(4)))
    rows = E._h_rows(3, H)
    check("absolute H row count reflects off-diagonal twice and diagonal once", rows == [4, 4, 9])
    reject("missing local upper Hessian rejected",
           lambda: E._scatter([E.ZERO] * 3, {}, a[0], a[1], a[2][:-1]), "coverage")
    reject("duplicate local Hessian coordinate rejected",
           lambda: E._scatter([E.ZERO] * 3, {}, a[0], a[1], a[2] + a[2][:1]), "upper-triangular")
    reject("negative local Hessian index rejected",
           lambda: E._scatter([E.ZERO] * 3, {}, a[0], a[1], [entry(-1, 0, 0)]), "upper-triangular")
    reject("Boolean local Hessian index rejected",
           lambda: E._scatter([E.ZERO] * 3, {}, a[0], a[1], [entry(False, 0, 0)]), "upper-triangular")
    reject("duplicate global position rejected",
           lambda: E._scatter([E.ZERO] * 3, {}, [0, 0], a[1], a[2]), "projection")
    reject("unordered packed force endpoints rejected",
           lambda: E._scatter([E.ZERO] * 3, {}, a[0], [pack(1, 0), a[1][1]], a[2]), "ordered")
    reject("nondyadic Arb endpoint rejected", lambda: E._packed(pack("1/3")), "dyadic")
    reject("float endpoint rejected", lambda: E._packed({"lower_exact": 0.1, "upper_exact": "1"}), "exact rational")

    R = [[0.5, -0.25, 0.0], [1.5, 0.0, -0.5], [0.0, 0.125, 2.0]]
    rf, rh, norm = E._dyadic_reduce(R, g, rows)
    check("signed integer-dyadic force products match hand calculation",
          rf == [(Q(-1), Q(7, 4)), (Q(-2), Q(7, 2)), (Q(-19, 8), Q(9, 4))])
    check("positive integer-dyadic row products and exact matrix norm",
          rh == [Q(3), Q(21, 2), Q(37, 2)] and norm == Q(17, 8))
    reject("negative H absolute row rejected", lambda: E._dyadic_reduce(R, g, [Q(0), Q(-1), Q(0)]), "direction")
    reject("wrong matrix dimension rejected", lambda: E._dyadic_reduce(R[:-1], g, rows), "dimension")

    ok = E._strict_test(Q(1, 8), Q(1, 4), Q(1, 8), Q(1, 8), Q(1))
    check("strict predicate and error computed from exact eta and z", ok == (Q(1, 4), Q(3, 8), Q(5, 8), True, Q(2, 5)))
    boundary = E._strict_test(Q(1, 4), Q(1, 4), Q(0), Q(1, 4), Q(1, 2))
    check("closed equality does not satisfy strict self-map test", boundary == (Q(1, 4), Q(1, 2), Q(1, 2), False, None))
    check("z equal one is unresolved with no root error",
          E._strict_test(Q(0), Q(1, 4), Q(0), Q(3, 4), Q(1))[-2:] == (False, None))
    reject("zero radius rejected", lambda: E._strict_test(Q(0), Q(0), Q(0), Q(0), Q(0)), "radius")

    attempts = [{"degree": 2, "force_tail_upper_exact": "2", "H_entry_tail_upper_exact": "1/2"},
                {"degree": 3, "force_tail_upper_exact": "1/2", "H_entry_tail_upper_exact": "1/2"}]
    E._order_history(3, attempts, 4, Q(1), Q(1))
    check("earliest two-budget selection accepted", True)
    reject("partial order coverage rejected", lambda: E._order_history(3, attempts[1:], 4, Q(1), Q(1)), "history")
    bad = deepcopy(attempts)
    bad[0]["force_tail_upper_exact"] = "1"
    reject("passing earlier order rejected", lambda: E._order_history(3, bad, 4, Q(1), Q(1)), "first order")
    bad = deepcopy(attempts)
    bad[-1]["H_entry_tail_upper_exact"] = "2"
    reject("final order must meet both budgets", lambda: E._order_history(3, bad, 4, Q(1), Q(1)), "first order")
    reject("negative tail endpoint rejected", lambda: E._tail(pack(-1, 0), Q(1)), "nonnegative")
    reject("tail over budget rejected", lambda: E._tail(pack(0, 2), Q(1)), "budget")
    reject("inward rho rejected", lambda: E._rho(Q(1, 8), Q(1, 8) + Q(1, 2**49)), "outward")

    manifest = [{"file": "parts/a", "sha256": "a" * 64}, {"file": "parts/b", "sha256": "b" * 64}]
    check("complete exact component manifest accepted", E._manifest("parts", manifest, ["a", "b"]) is manifest)
    reject("partial component manifest rejected", lambda: E._manifest("parts", manifest[:1], ["a", "b"]), "complete")
    reject("duplicate component across manifest rejected", lambda: E._manifest("parts", manifest[:1] * 2, ["a", "b"]), "distinct")
    reject("component substitution path rejected", lambda: E._manifest("other", manifest, ["a", "b"]), "canonical")
    reject("nonterminal process rejected",
           lambda: E._terminal({"status": "RUNNING", "pid": 10, "exit_code": 0, "results_sha256": "0" * 64}), "terminal")
    reject("Boolean pseudo-success exit rejected",
           lambda: E._terminal({"status": "FINISHED", "pid": 10, "exit_code": False, "results_sha256": "0" * 64}), "terminal")
    reject("changed terminal result hash rejected",
           lambda: E._terminal({"status": "FINISHED", "pid": 10, "exit_code": 0, "results_sha256": "0" * 64}, "1" * 64), "binding")
    reject("suspect short digest rejected", lambda: E._sha("0" * 63), "SHA-256")
    reject("noncanonical record rejected", lambda: E._record_dir("../optoGEF2_small_col_L"), "record ID")
    reject("inward producer endpoint rejected", lambda: E._subset((Q(-1), Q(0)), (Q(-2), Q(0)), "toy"), "containment")
    check("negative and positive endpoint containment accepted", E._contains((Q(-3), Q(2)), (Q(-2), Q(1))))

    # End-to-end arithmetic fixture uses zero H with an explicitly complete 1x1 upper triangle.
    pretend = {
        "force_integral_minus_cache": [pack(-Q(1, 8), Q(1, 8))],
        "preconditioned_force_defect": [pack(-Q(1, 8), Q(1, 8))],
        "H_absolute_row_bounds": [pack(0)],
        "H_integral_minus_cache_upper": [entry(0, 0, 0)],
        "RF_defect_uniform_upper_exact": "1/8", "b_upper_exact": "0",
        "eta_upper_exact": "1/4", "z_upper_exact": "1/4", "selfmap_lhs_exact": "1/2", "r_exact": "1",
        "sufficient_conditions_pass": True, "root_assertion_independently_reviewed": False,
        "status": "CANDIDATE_UNIFORM_EXACT_FE_ROOT_UNREVIEWED", "root_error_upper_exact": "1/3",
    }
    force = [(Q(-1, 8), Q(1, 8))]
    verify = lambda d: E._verify_producer(d, force, {(0, 0): E.ZERO}, force, [Q(0)], Q(1, 8),
                                         Q(0), Q(1, 8), Q(1, 4), Q(1))
    check("end-to-end synthetic producer enclosure and scalar agreement", verify(pretend))
    widened = {**pretend, "H_absolute_row_bounds": [pack(Q(1, 16))],
               "b_upper_exact": "1/16", "z_upper_exact": "5/16", "selfmap_lhs_exact": "9/16",
               "root_error_upper_exact": "4/11"}
    check("producer's widened scalar Hessian row upper bound dominates the exact row", verify(widened))
    reject("negative producer Hessian row rejected",
           lambda: verify({**pretend, "H_absolute_row_bounds": [pack(-Q(1, 16), 0)]}), "H row")
    reject("incorrect reported epsilon rejected", lambda: verify({**pretend, "root_error_upper_exact": "1/4"}), "arithmetic")
    reject("incorrect reported b rejected", lambda: verify({**pretend, "b_upper_exact": "-1"}), "domination")
    return checks


def run(freeze):
    checks = tests()
    source = Path(E.__file__).resolve()
    script = Path(__file__).resolve()
    bindings = {}
    for name, digest in E.PREPARATION_PINS.items():
        E._bound(name, digest, bindings, parse=False)
    context = E._bound(E.BASE + "/context_preparation/protocol.json",
                       bindings[E.BASE + "/context_preparation/protocol.json"], bindings)
    E._bound(context["source_file"], context["source_sha256"], bindings, parse=False)
    for p in (source, script):
        bindings[p.relative_to(ROOT).as_posix()] = sha256(p.read_bytes()).hexdigest()
    checks["historical cached-energy and current frozen context/integrator SHA bindings"] = True
    protocol = {
        "version": 1, "status": "FROZEN_ENDPOINT_REVIEW_PREPARATION_ONLY",
        "utc": datetime.now(timezone.utc).isoformat(), "input_hashes": bindings,
        "actual_admission_requires_external_primary_pins": E.PRIMARY_PINS,
        "review_policy": "Freeze 168 expected packet paths and digests from two successful terminal part manifests. "
                         "Run the reducer only with the resulting external admission digest.",
        "arithmetic": "Exact rational signed scatter of every force and upper Hessian endpoint, exact absolute row counts, "
                      "integer-dyadic signed Rg and |R|H-row contractions, then exact scalar strict predicate.",
        "producer_comparison": "Containment for signed F, H and RF intervals. Upper-bound domination for H absolute "
                               "rows and norm scalars, which are bounds already computed on widened Arb intervals.",
        "energy_policy": "Reuse the independently checked same-record cached SPD proof and exact binary R binding. "
                         "H_Q+t(H_I-H_Q) stays nonsingular by the finite infinity-norm bound, preserves inertia, "
                         "and yields a positive normalized curvature bound through an inverse norm.",
        "shared_premises": "Explicit frozen context loader, independently reviewed local integral/cache engines "
                           "and independently reviewed original-cube body and surface domains.",
        "real_component_payloads_opened": False, "real_record_review_run": False,
        "superseded_preparations": {
            path.relative_to(OUT).as_posix(): sha256(path.read_bytes()).hexdigest()
            for path in sorted(OUT.glob("development_revision_*/*")) if path.is_file()
        },
    }
    report = {
        "status": "PASS_REDUCTION_PREPARATION_ONLY", "checks": checks, "check_count": len(checks),
        "source_sha256": bindings[source.relative_to(ROOT).as_posix()],
        "real_record_review_run": False, "component_integration_run": False,
    }
    if freeze:
        snapshot = OUT / "review_cohort_exact_record_frozen.py"
        if any(path.exists() for path in (snapshot, OUT / "protocol.json", OUT / "preflight.json")):
            raise ValueError("write-once reduction outputs already exist")
        with snapshot.open("xb") as f:
            f.write(source.read_bytes())
        report["protocol_sha256"] = E._write(OUT / "protocol.json", protocol)
        report["preflight_sha256"] = E._write(OUT / "preflight.json", report)
    print(json.dumps({key: value for key, value in report.items() if key != "checks"}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    run(parser.parse_args().freeze)
