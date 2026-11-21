"""Freeze and test the record-neutral bulk adapter without reading active domains."""
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from flint import ctx as precision


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import exact_scaled_cohort_bulk as A
import exact_cohort_context as C


E, D, S = A.E, A.D, A.S
BASE = D.BASE
REMAINING = BASE + "/remaining15_preparation"
PILOT = str(Path(A.SCALED_KERNEL).parent.as_posix())
SYNTHETIC = BASE + "/integrator_preparation/preflight_integrator.py"
PIN = {
    A.SCALED_KERNEL: "a173b3759fb66fe66e5f3e442c1bcce9a5a185e67515649a94499f82d1d67f60",
    PILOT + "/actual_independent_review.json": "f87cc324c6469f0717cf399ebd4f695302ca8495698296ea886987b7cf45e9d4",
    BASE + "/integrator_preparation/protocol.json": "cced8eb0a3b7b15d96fe0c4f43da042c807ecd109bcaf1dd37e0a0de1e979096",
    REMAINING + "/protocol.json": "36f58fe4a7b4c100add489b1608b33a264fd35b3d2f30d23d9e03fe2f1a640c0",
    REMAINING + "/preflight.json": "2c7a39c9d2e225afe4cd0f25e8172d5f2c8d12c8820e2a10d4f527bb19143f3d",
    REMAINING + "/manifest.json": "93a2f5858a450ca9a8db9f9bc944b918e2a7446abba614f97de5e47fe434495d",
    "scripts/exact_cohort_context.py": "c22612a1de1ccfd34e8b0f0859f6413c4da09c765ba4945d060340e93e928f9f",
    "scripts/exact_cohort_integrator.py": "092dc51161e708cdba3791ec47a7f8e843dd0e1afc866b6f65f286b8657d1745",
    "scripts/review_cohort_exact_record.py": "67db748a5d40f31471ddd3d6c84fae280ff7ead650e6eb4e950a4766e54e0285",
    "evidence/exact_cohort_integrator_independent_review.json": "1d233c61e5b7283ddfeb856e206a064d28cb2731246f13bff716e9a4063dbfa8",
}
EXCLUDED = {"optoGEF2_small_col_B", "optoGEF2_small_col_L"}


def _relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _inputs():
    hashes = {}
    d = {p: D._bound(p, hid, hashes, parse=p.endswith(".json")) for p, hid in PIN.items()}
    review = d[PILOT + "/actual_independent_review.json"]
    D.require(review["status"] == "PASS_ACTUAL_SCALED_E21_ENDPOINT_REVIEW" and D._yes(review["checks"])
              and _sha(S.__file__) == PIN[A.SCALED_KERNEL], "reviewed scaled primitives must match actual L21 review")
    D._bound(review["support_file"], review["support_sha256"], hashes, parse=False)
    integrator = d[BASE + "/integrator_preparation/protocol.json"]
    for path, digest in integrator["source_hashes_and_inputs"].items():
        D._bound(path, digest, hashes, parse=False)
    pre, manifest = d[REMAINING + "/preflight.json"], d[REMAINING + "/manifest.json"]
    D.require(pre["status"] == "PASS_REMAINING15_PREPARATION_ONLY" and D._yes(pre["checks"])
              and pre["manifest_sha256"] == PIN[REMAINING + "/manifest.json"]
              and manifest["status"] == "READY_FOR_INDEPENDENT_RECORD_DOMAIN_PRELAUNCH_REVIEWS",
              "previous frozen remaining15 preparation required")
    rows = manifest["records"]
    ids = [v["record_id"] for v in rows]
    D.require(len(ids) == len(set(ids)) == 15 and ids == manifest["expected_record_ids"]
              and set(ids).isdisjoint(EXCLUDED), "precisely the remaining 15 unique selected records required")
    return hashes, rows, integrator["source_hashes_and_inputs"][SYNTHETIC]


def _array_key(array):
    return {"dtype": array.dtype.str, "shape": array.shape,
            "binary_sha256": sha256(array.tobytes(order="C")).hexdigest()}


def _route_profile(context):
    a = context.target.a
    route = [sum(a["reference"][int(node), 0] == 0 for node in tri[:3]) for tri in a["triangles"]]
    mapped = {name: [i for i, r in enumerate(route) if r == count]
              for count, name in enumerate(("off_axis", "axis_vertex", "axis_edge_physical"))}
    D.require(len(route) == 144 and list(map(len, mapped.values())) == [132, 6, 6]
              and mapped["axis_vertex"] == list(range(6)) and mapped["axis_edge_physical"] == list(range(72, 78)),
              "unchanged 132 and 6+6 element chart routes required")
    return mapped


def _record_profiles(rows, hashes):
    # Only the frozen L cache's geometry and rule are used for equality. No L
    # field or any element force output supplies a different record's data.
    baseline = C.load_context("optoGEF2_small_col_L")
    for p, h in baseline.input_hashes.items():
        D._bound(p, h, hashes, parse=False)
    common = {key: _array_key(a) for key, a in baseline.target.a.items() if key not in ("u", "preconditioner")}
    output = []
    for row in rows:
        record = row["record_id"]
        ctx = C.load_context(record)
        for p, h in ctx.input_hashes.items():
            D._bound(p, h, hashes, parse=False)
        for p in row["files"].values():
            D._bound(p["file"], p["sha256"], hashes, parse=False)
        own = {key: _array_key(a) for key, a in ctx.target.a.items() if key not in ("u", "preconditioner")}
        D.require(own == common and ctx.row == row["selected_context"]
                  and _array_key(ctx.u)["binary_sha256"] == row["u_binary_sha256"]
                  and _array_key(ctx.R)["binary_sha256"] == row["R_binary_sha256"]
                  and _array_key(ctx.S)["binary_sha256"] == row["S_binary_sha256"]
                  and ctx.old_radius == Q(row["selected_corrector_radius_exact"])
                  and ctx.dose == Q(row["exact_dose"]) and ctx.drive == Q(row["exact_drive"]),
                  "same static geometry/rule with this selected record's own point, R and S required")
        output.append({
            "record_id": record, "selected_context": ctx.row, "choice_of_chart": _route_profile(ctx),
            "static_geometry_and_auxiliary_cache_bitidentity_sha256": D._digest(own),
            "own_u_R_S_binary_sha256": {n: _array_key(v)["binary_sha256"] for n, v in (
                ("u", ctx.u), ("R", ctx.R), ("S", ctx.S))},
            "old_radius_exact": ctx.old_radius, "exact_dose": ctx.dose, "exact_drive": ctx.drive,
            "existing_original_force_tail_budget_exact": row["force_tail_budget_exact"],
            "existing_original_H_entry_tail_budget_exact": row["H_entry_tail_budget_exact"],
        })
    return output


def tests(toy):
    checks = {}

    def ck(name, value):
        D.require(value, name)
        checks[name] = True

    def reject(name, action):
        try:
            action()
        except ValueError:
            checks[name] = True
        else:
            raise AssertionError(name)

    for route, vertices in (
        ("off_axis", [(1, 0), (3, 0), (1, 2)]),
        ("axis_vertex", [(0, 0), (2, 0), (2, 3)]),
        ("axis_edge_physical", [(0, 0), (2, 0), (0, 3)]),
    ):
        context, nodes, rule, body = toy(vertices, active=1)
        free = sorted(context.target.index)
        context.target.a["free"] = np.asarray(free, dtype=np.int64)
        context.target.a["u"] = np.zeros(len(free), dtype=np.float64)
        selector = E.prepare_bulk_selector(context, nodes, rule, body, 0, Q(1, 10**9), Q(1, 10**7))
        proof = A.prepare_component(context, nodes, rule, body, selector, 0)
        expected_route = "scaled_positive_R_affine_material" if route == "off_axis" else route
        ck(route + ": exact chart routing", proof["chart"] == route and proof["integration_chart"] == expected_route)
        if route == "off_axis":
            ck("off-axis pure symbolic strain/derivative identity checks", proof["scaled_polynomial_identity_count"] > 100
               and proof["scaled_exact_Qstar_node_count"] == 64)
        explicit = {"context": context, "nodes": nodes, "rule": rule, "selector": selector}
        identity = {k: v for k, v in E._bulk_setup(context, nodes, rule, 0)[0].items()}
        with (
            patch.object(E, "bulk_component", return_value=identity) as old,
            patch.object(A, "_scaled_component", return_value=identity) as scaled,
        ):
            packet = A.bulk_component(context, nodes, rule, body, selector, 0, proof)
            ck(route + ": intended matching kernel alone called",
               (old.call_count, scaled.call_count) == ((0, 1) if route == "off_axis" else (1, 0))
               and packet["integration_chart"] == proof["integration_chart"])
            endpoint = scaled if route == "off_axis" else old
            ck(route + ": explicit original context and rule passed",
               all(a is b for a, b in zip(endpoint.call_args.args[:3], (explicit["context"], explicit["nodes"], explicit["rule"]))))
            reject(route + ": changed body J coverage rejected", lambda: A.bulk_component(context, nodes, rule,
                {**body, "J_bound": {"lower_exact": "4/5", "upper_exact": "11/10"}}, selector, 0, proof))
            reject(route + ": changed symbolic proof rejected", lambda: A.bulk_component(context, nodes, rule,
                body, selector, 0, {**proof, "original_selector_sha256": "f" * 64}))
            other = SimpleNamespace(target=context.target, cell=context.cell, dose=context.dose,
                                    record_id="different_record", row={**context.row, "record_id": "different_record"})
            reject(route + ": other record cannot reuse selector", lambda: A.bulk_component(other, nodes, rule,
                body, selector, 0, proof))
    s = {"max_DJ_norm_upper_exact": "2", "max_D2J_norm_upper_exact": "3", "phi_prime_tail_upper_exact": "1/10",
         "phi_second_tail_upper_exact": "1/5", "force_tail_budget_exact": "1", "hessian_tail_budget_exact": "4"}
    f, H = A._tail_bounds(s, Q(2), [Q(9), Q(1)], {(0, 0): Q(8), (0, 1): Q(1), (1, 1): Q(1)})
    ck("per-entry new chart bounds may fall back to original qualified majorants", f == [Q(2, 5), Q(1, 5)]
       and H == {(0, 0): Q(11, 5), (0, 1): Q(1), (1, 1): Q(3, 5)})
    reject("force tail over original budget rejected",
           lambda: A._tail_bounds({**s, "force_tail_budget_exact": "1/20"}, Q(2), [Q(9)], {(0, 0): Q(8)}))
    return checks


def run():
    hashes, rows, synthetic_digest = _inputs()
    D._bound(SYNTHETIC, synthetic_digest, hashes, parse=False)
    sys.path.insert(0, str((ROOT / SYNTHETIC).parent))
    import preflight_integrator as synthetic
    D.require(_sha(synthetic.__file__) == synthetic_digest, "same existing pure chart fixture required")
    precision.prec = 640
    E.core.A.cache_clear()
    forbidden = AssertionError("a real component/moment/cache sum was attempted during adapter preparation")
    with (
        patch.object(S.Functional, "__init__", side_effect=forbidden),
        patch.object(S, "_eval_parts", side_effect=forbidden),
        patch.object(E.P.Functional, "__init__", side_effect=forbidden),
        patch.object(E.U, "uniform_cache", side_effect=forbidden),
        patch.object(E.U.C, "exact_bulk_maps", side_effect=forbidden),
        patch.object(E.M, "triangle_moment", side_effect=forbidden),
    ):
        profiles = _record_profiles(rows, hashes)
        checks = tests(synthetic.toy)
    checks.update({
        "all_15_static_geometry_and_rule_bitidentical_to_L": len(profiles) == 15 and len({
            r["static_geometry_and_auxiliary_cache_bitidentity_sha256"] for r in profiles}) == 1,
        "all_15_distinct_selected_anchor_states": len({p["own_u_R_S_binary_sha256"]["u"] for p in profiles}) == 15,
        "all_original_private_state_hashes_and_doses_loaded": True,
        "all_2160_routes_verified_from_static_geometry": sum(sum(map(len, p["choice_of_chart"].values())) for p in profiles) == 2160,
        "no_actual_body_or_scale_moment_cache_operations_run": True,
    })
    D.require(D._yes(checks), "every remaining15 neutral chart guard required")
    src = Path(A.__file__).resolve()
    script = Path(__file__).resolve()
    names = ("exact_scaled_cohort_bulk_frozen.py", "protocol.json", "preflight.json")
    D.require(not any((OUT / name).exists() for name in names), "preserve write-once generic scaled preparation")
    with (OUT / names[0]).open("xb") as f:
        f.write(src.read_bytes())
    protocol = {
        "status": "FROZEN_REMAINING15_SCALED_BULK_ADAPTER_REVIEW_PENDING", "source_file": _relative(src),
        "source_sha256": _sha(src), "preflight_source_file": _relative(script), "preflight_source_sha256": _sha(script),
        "input_hashes": hashes, "utc": datetime.now(timezone.utc).isoformat(),
        "record_profiles": profiles, "precision_bits": 640,
        "api": "prepare_component(context,nodes,rule,body_element,selector,element) creates a pure bound route proof. "
               "bulk_component(context,nodes,rule,body_element,selector,element,proof) verifies identity and evaluates its route.",
        "scope": "Record-neutral component helper on each record's original independently reviewed box, closed C, "
                 "rule, order, rho and force/Hessian tail budgets. No cross-record field data; off-axis scaled "
                 "material series and signed moments from reviewed pure primitives; axis cases delegate unchanged.",
        "job_gate": "Actual record phases require true terminal successful geometry domains and independent per-record "
                    "domain/operator review, including a freeze of derived route proofs; no old B/L gate is inherited.",
        "actual_component_runs_started": False,
        "failed_attempts": {_relative(p): _sha(p) for p in sorted(OUT.glob("development_revision_*/*")) if p.is_file()},
    }
    psha = D._write(OUT / names[1], protocol)
    report = {
        "status": "PASS_REMAINING15_SCALED_ADAPTER_PREPARATION_ONLY", "source_sha256": _sha(src), "protocol_sha256": psha,
        "checks": checks, "check_count": len(checks), "record_count": len(profiles), "actual_component_runs_started": False,
    }
    presha = D._write(OUT / names[2], report)
    print(json.dumps({"status": report["status"], "source_sha256": _sha(src), "protocol_sha256": psha,
                      "preflight_sha256": presha, "records": len(profiles), "checks": len(checks)}), flush=True)


if __name__ == "__main__":
    run()
