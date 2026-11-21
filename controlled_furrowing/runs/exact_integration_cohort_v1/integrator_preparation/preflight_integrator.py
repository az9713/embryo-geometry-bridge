"""Prepare the explicit-record integrator with bounded synthetic checks only."""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import monotonic
from types import SimpleNamespace
from unittest.mock import patch

from flint import arb, ctx
import numpy as np


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import exact_cohort_context as C
import exact_cohort_integrator as I


def toy(corners, active):
    corners = [list(map(Q, row)) for row in corners]
    nodes = corners + [[(corners[i][c] + corners[j][c]) / 2 for c in (0, 1)]
                       for i, j in ((0, 1), (1, 2), (2, 0))]
    free = [2 * i + axis for i, row in enumerate(nodes) for axis in (0, 1) if axis or row[0]]
    target = I.T.Target.__new__(I.T.Target)
    target.a = {
        "triangles": np.array([[0, 1, 2, 3, 4, 5]], dtype=np.int64),
        "reference": np.array(nodes, dtype=np.float64),
        "active": np.array([[active]], dtype=np.float64),
        "shape": np.array([[0, 0, 0, 0.5, 0.5, 0]] * 64),
        "volume_weights": np.array([[0.125] * 64]),
        "surface_0_N": np.array([[0.0, 1.0, 0.0]] * 12),
        "surface_0_weights": np.array([[0.125] * 12]),
    }
    target.index = {value: i for i, value in enumerate(free)}
    target.energy_scale = Q(7)
    target.packet = {"parameters": {"mu": 1, "bulk": 5000}}
    cell = {key: (Q(1), Q(2)) for key in I.T.ORDER}
    cell["common_gain"] = (Q(3, 4), Q(5, 4))
    cell_strings = {key: list(map(str, bounds)) for key, bounds in cell.items()}
    context = SimpleNamespace(
        target=target, record_id="synthetic", dose=Q(2, 3), cell=cell,
        row={
            "record_id": "synthetic", "dose_exact": "2/3", "parameter_cell": cell_strings,
            "parameter_cell_sha256": I._digest(cell_strings),
            "model_observation_contract_sha256": "synthetic",
            "point_cache_sha256": "synthetic", "predictor_sha256": "synthetic",
        },
    )
    rule = {
        "volume_barycentric": [[Q(1, 4), Q(1, 2), Q(1, 4)]] * 64,
        "volume_weights": [[Q(1, 56)] * 64],
        "surface_s": [[Q(1, 2)] * 12],
        "surface_weights": [[[Q(1, 56)] * 12]],
    }
    body = {
        "positive_J_for_entire_box": True,
        "J_bound": {"lower_exact": "9/10", "upper_exact": "11/10"},
    }
    return context, nodes, rule, body


def tests():
    checks = {}

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks[name] = True

    def reject(name, action, text=None):
        try:
            action()
        except ValueError as error:
            check(name, text is None or text in str(error))
        else:
            raise AssertionError(f"{name}: accepted malformed input")

    forbidden = AssertionError("preflight entered an integration or all-site cache routine")
    with patch.object(I.P.Functional, "__init__", side_effect=forbidden), \
            patch.object(I.U, "series", side_effect=forbidden), \
            patch.object(I.U, "uniform_cache", side_effect=forbidden):
        for chart, corners, active in (
            ("off_axis", [(1, 0), (3, 0), (1, 2)], 0),
            ("axis_vertex", [(0, 0), (2, 0), (2, 3)], 1),
            ("axis_edge_physical", [(0, 0), (2, 0), (0, 3)], 1),
        ):
            c, nodes, rule, body = toy(corners, active)
            target = c.target
            f, _ = I._bulk_fields(target, nodes, 0)
            expected_type = I.AXbridge.PhysicalFunctional if chart == "axis_edge_physical" else I.P.Functional
            check(chart + ": uses intended reviewed functional type", I._functional_type(f) is expected_type)
            duffy = target.element(0, nodes)
            check(chart + ": same free projection as reviewed original chart",
                  f["global_indices"] == duffy["global_indices"] and f["free_local"] == duffy["free_local"])
            check(chart + ": canonical identity J exactly one",
                  I.JET.bulk_jacobian_jet(f["q"], f["B"])["value"] == I.T.ONE)
            drive = I._bulk_materials(c, f)
            check(chart + ": exact active mask and gain times dose",
                  drive.contains(I.U.A(Q(1, 2))) and drive.contains(I.U.A(Q(5, 6))) if active else drive.is_zero())

            force_budget, H_budget = Q(1, 10**9), Q(1, 10**7)
            s = I.prepare_bulk_selector(c, nodes, rule, body, 0, force_budget, H_budget)
            check(chart + ": selector records explicit element and chart",
                  s["chart"] == chart and s["global_indices"] == f["global_indices"]
                  and s["precision_bits"] == 640 and s["active"] == str(active))
            scaled = Q(1, 10) * 2**48
            check(chart + ": rho ceiling computed exactly", Q(s["rho_exact"]) == Q((scaled.numerator + 4) // 5, 2**48)
                  and Q(s["rho_exact"]) >= Q(s["raw_rho_exact"]) == Q(1, 10))
            check(chart + ": selected first order passing both budgets",
                  s["degree"] > 2 and [row["degree"] for row in s["order_attempts"]] == list(range(2, s["degree"] + 1))
                  and all(Q(row["force_tail_upper_exact"]) > force_budget
                          or Q(row["H_entry_tail_upper_exact"]) > H_budget for row in s["order_attempts"][:-1])
                  and Q(s["force_tail_upper_exact"]) <= force_budget and Q(s["H_entry_tail_upper_exact"]) <= H_budget)
            mass = Q(s["geometry_two_mass_upper_exact"])
            area = Q(f["det"], 2)
            R_avg = sum(Q(x[0]) for x in corners) / 3
            check(chart + ": exact two-mass bound agrees with independent triangle formula",
                  2 * Q(3) * area * R_avg / 7 + Q(64, 56) < mass
                  < 2 * Q(22, 7) * area * R_avg / 7 + Q(64, 56))
            first, second = (Q(s[key]) for key in ("phi_prime_tail_upper_exact", "phi_second_tail_upper_exact"))
            g, h = Q(s["max_DJ_norm_upper_exact"]), Q(s["max_D2J_norm_upper_exact"])
            check(chart + ": exact force and Hessian scalar majorants",
                  Q(s["force_tail_upper_exact"]) == mass * first * g
                  and Q(s["H_entry_tail_upper_exact"]) == mass * (second * g * g + first * h))
            uncertain = arb(0, I.U.A(Q(1, 2)))
            check(chart + ": uncertain zero-crossing coefficient contributes",
                  I._norm({(0, 0): uncertain}, f) >= Q(1, 2))
            if chart == "off_axis":
                check(chart + ": inverse radius uses positive lower endpoint",
                      I._norm({(-1, 0): Q(2)}, f) == 2)
            else:
                reject(chart + ": inverse radius rejected",
                       lambda: I._norm({(-1, 0): Q(1)}, f), "exponents")
            reject(chart + ": negative axial power rejected",
                   lambda: I._norm({(1, -1): Q(1)}, f), "exponents")
            reject(chart + ": Boolean power rejected",
                   lambda: I._norm({(True, 0): Q(1)}, f), "exponents")

            build = lambda **kw: I.prepare_bulk_selector(
                c, kw.get("nodes", nodes), kw.get("rule", rule), kw.get("body", body),
                kw.get("index", 0), kw.get("force", force_budget), kw.get("hessian", H_budget),
                kw.get("maximum", 32))
            altered = deepcopy(s)
            altered["global_indices"] = altered["global_indices"][::-1]
            reject(chart + ": altered free projection rejected before integration",
                   lambda: I.bulk_component(c, nodes, rule, altered), "projection")
            altered = deepcopy(s)
            altered["rho_exact"] = "0"
            reject(chart + ": altered outward rho rejected before integration",
                   lambda: I.bulk_component(c, nodes, rule, altered), "metadata changed")
            shifted = deepcopy(nodes)
            shifted[-1][1] += Q(1, 10**5)
            reject(chart + ": new nodal field rejected under old selector",
                   lambda: I.bulk_component(c, shifted, rule, s), "nodes or rule changed")

        reject("zero force budget rejected", lambda: build(force=0), "positive")
        reject("negative Hessian budget rejected", lambda: build(hessian=-1), "positive")
        reject("infinite budget rejected", lambda: build(force=float("inf")), "finite")
        reject("Boolean budget rejected", lambda: build(hessian=True), "non-Boolean")
        reject("Boolean maximum order rejected", lambda: build(maximum=True), "integer")
        reject("maximum order above 32 rejected", lambda: build(maximum=33), "2..32")
        reject("Boolean element rejected", lambda: build(index=False), "element")
        reject("nonpositive body rejected",
               lambda: build(body={**body, "positive_J_for_entire_box": False}), "positive-J")
        reject("analytic cutoff rejected",
               lambda: build(body={**body, "J_bound": {"lower_exact": "1/1000000", "upper_exact": "2"}}), "cutoff")
        reject("divergent rho rejected",
               lambda: build(body={**body, "J_bound": {"lower_exact": "1/2", "upper_exact": "2"}}), "converge")
        try:
            build(maximum=2)
        except I.TailBudgetUnresolved as error:
            check("bounded unresolved selector preserves attempted order and tails",
                  len(error.attempts) == 1 and error.attempts[0]["degree"] == 2
                  and Q(error.attempts[0]["force_tail_upper_exact"]) > force_budget)
        else:
            raise AssertionError("bounded unresolved order accepted")

        altered_rule = deepcopy(rule)
        altered_rule["volume_weights"][0][0] *= 2
        reject("changed bulk weight rejected",
               lambda: build(rule=altered_rule), "bulk weight")
        altered_rule = deepcopy(rule)
        altered_rule["volume_barycentric"][0] = [Q(1, 3)] * 3
        reject("changed bulk site rejected",
               lambda: build(rule=altered_rule), "bulk site")
        c.row["dose_exact"] = "1"
        reject("changed dose context rejected", build, "exact dose")
        c.row["dose_exact"] = "2/3"
        c.cell["common_gain"] = (Q(1), Q(2))
        reject("changed active parameter cell rejected", build, "unchanged cell")
        c.cell["common_gain"] = (Q(3, 4), Q(5, 4))
        ctx.prec = 512
        try:
            reject("changed caller precision rejected", build, "640-bit")
        finally:
            ctx.prec = I.PRECISION_BITS
            I.core.A.cache_clear()

    # A single matched toy map exercises the safe cache-consistency primitive only.
    eye = [[Q(int(i == j)) for j in range(5)] for i in range(5)]
    strain = list(map(I.U.A, [1, Q(1, 20), 0, 1, 1]))
    same = I.U.C.consistency_bound(eye, eye, strain, list(range(5)), 1, 5000, I.U.A(1))
    check("matched synthetic bulk maps yield exact zero delta q and Hessian difference",
          all(value.is_zero() for value in same["delta_q"] + list(same["difference"].entries())))

    case = {"face": 0, "edge": 0, "order": 53, "record_id": c.record_id}
    with patch.object(I.SUR, "run_edge", return_value={"scope": "synthetic scope"}) as endpoint:
        result = I.surface_component(c, nodes, rule, case)
        check("surface adapter forwards explicit target, nodes, case and larger order",
              endpoint.call_count == 1 and all(a is b for a, b in zip(endpoint.call_args.args[:5],
              (c.target, nodes, c.row, case, rule))) and endpoint.call_args.args[5] == 53)
        check("surface result binds cohort record, dose and cell",
              result["record_id"] == c.record_id and result["dose_exact"] == "2/3"
              and result["parameter_cell_sha256"] == c.row["parameter_cell_sha256"]
              and result["status"] == "UNIFORM_COHORT_SURFACE_COMPONENT_UNREVIEWED")
        reject("cross-record surface case rejected",
               lambda: I.surface_component(c, nodes, rule, {**case, "record_id": "different"}), "different record")
        reject("Boolean surface order rejected",
               lambda: I.surface_component(c, nodes, rule, {**case, "order": True}), "integer")
        altered = deepcopy(rule)
        altered["surface_weights"][0][0] = [Q(1)] * 12
        reject("changed surface weights rejected",
               lambda: I.surface_component(c, nodes, altered, case), "weights required")
        check("rejected surface cases never reach kernel", endpoint.call_count == 1)
    return checks


def _read_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def run(freeze):
    start = monotonic()
    ctx.prec = I.PRECISION_BITS
    I.core.A.cache_clear()
    checks = tests()
    source = Path(I.__file__).resolve()
    context_protocol = ROOT / "runs/exact_integration_cohort_v1/context_preparation/protocol.json"
    context_preflight = context_protocol.with_name("preflight.json")
    context_binding = json.loads(context_protocol.read_bytes())
    context_receipt = json.loads(context_preflight.read_bytes())
    check_paths = [Path(module.__file__).resolve() for module in (
        C, I, I.AX, I.AXbridge, I.U, I.SUR, I.P, I.T, I.M, I.JET, I.core, I.U.C, I.U.C.target_core,
        I.SUR.S, I.SUR.L, I.SUR.L.target)]
    check_paths += [context_protocol, context_preflight, Path(__file__).resolve()]
    paths = {path.relative_to(ROOT).as_posix(): _read_hash(path) for path in sorted(set(check_paths))}
    if context_binding["source_sha256"] != _read_hash(Path(C.__file__)) or \
            context_receipt["protocol_sha256"] != _read_hash(context_protocol) or \
            context_receipt["status"] != "PASS_CONTEXT_PREPARATION_ONLY":
        raise AssertionError("changed context preparation")
    checks["frozen context protocol and loader receipt bindings agree"] = True
    protocol = {
        "version": 1,
        "status": "FROZEN_INTEGRATOR_PREPARATION_REVIEW_PENDING",
        "utc": datetime.now(timezone.utc).isoformat(),
        "source_hashes_and_inputs": paths,
        "precision_bits": I.PRECISION_BITS,
        "caller_precision_policy": "Set flint.ctx.prec=640 and clear core.A's cached conversions before constructing interval nodes.",
        "bulk_policy": "Use the same exact affine/P2 target, closed C and record dose. Select on certified positive element J "
                       "with exact dyadic outward rho; order 2 through the caller's maximum (at most 32). "
                       "Both caller-provided positive tail budgets must pass.",
        "geometry_mass_policy": "Compute only the zeroth material mass using det*sum(R_corner)/6 and the exact normalized old "
                                "quadrature weights. Scale only the material term by 2*pi()/energy_scale.",
        "chart_policy": "Count actual exact radial zeros in each corner triple. Two axis corners take reviewed physical "
                        "ordinary polynomial q/B moments. One axis corner uses its regularized Duffy chart. "
                        "The off-axis chart keeps its positive-radius Laurent moments.",
        "cache_policy": "Bind explicit rule to exact cached sites and normalized weights. The same Y nodal box, B columns, "
                        "material parameters and drive feed both finite and integral evaluations.",
        "execution_precondition": "Caller independently reviews source, selected record, domain proof, selector/case inputs "
                                  "and proof margin before any actual component run.",
        "process_memory_policy": "Use one record per process and clear functional instance caches after each component.",
        "surface_policy": "Forward explicit selected-record target, nodes, rule, row and declared order to shared safe uniform "
                          "trace kernel. Whole-edge G and local cache pressure-domain bounds remain case inputs.",
        "integral_force_or_Hessian_expansion_executed": False,
        "full_bulk_or_surface_cache_sum_executed": False,
        "root_or_assembly_result": False,
    }
    report = {
        "status": "PASS_INTEGRATOR_PREPARATION_ONLY",
        "source_sha256": _read_hash(source),
        "checks": checks,
        "check_count": len(checks),
        "elapsed_seconds": monotonic() - start,
        "precision_bits": I.PRECISION_BITS,
        "actual_record_force_H_moment_expansions": False,
        "actual_bulk_or_surface_cache_sums": False,
        "assembly_or_root_assertion": False,
    }
    if freeze:
        files = [OUT / name for name in ("exact_cohort_integrator_frozen.py", "protocol.json", "preflight.json")]
        if any(path.exists() for path in files):
            raise ValueError("write-once preparation output exists")
        with files[0].open("xb") as out:
            out.write(source.read_bytes())
        with files[1].open("x", encoding="utf8", newline="\n") as out:
            out.write(json.dumps(protocol, indent=2) + "\n")
        report["protocol_sha256"] = _read_hash(files[1])
        with files[2].open("x", encoding="utf8", newline="\n") as out:
            out.write(json.dumps(report, indent=2) + "\n")
        report["preflight_sha256"] = _read_hash(files[2])
    print(json.dumps({key: value for key, value in report.items() if key != "checks"}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    run(parser.parse_args().freeze)
