"""Independent local-adapter review. Integration and full cache sums are blocked."""
import ast
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

from flint import arb, ctx
import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
PREP = OUT.parent / "integrator_preparation"
SOURCE = "scripts/exact_cohort_integrator.py"
EVIDENCE = ROOT / "evidence/exact_cohort_integrator_independent_review.json"
FROZEN = {
    SOURCE: "092dc51161e708cdba3791ec47a7f8e843dd0e1afc866b6f65f286b8657d1745",
    "protocol.json": "cced8eb0a3b7b15d96fe0c4f43da042c807ecd109bcaf1dd37e0a0de1e979096",
    "preflight.json": "a6152aa1223d34c0c06a6e6f51285d4056791d2c4554ebd608d05dd2d091d9e8",
    "evidence/exact_cohort_context_independent_review.json": "898fb6828c16fc7f0db410c9d677373cf1696db5d0d418c08ef0004a30e215c5",
}
BODY = {"positive_J_for_entire_box": True, "J_bound": {"lower_exact": "19/20", "upper_exact": "21/20"}}


def digest(data):
    return sha256(data).hexdigest()


def tree_function(source, name):
    fn, = [node for node in ast.parse(source).body
           if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name == name]
    return fn


def call_names(node):
    return {ast.unparse(child.func) for child in ast.walk(node) if isinstance(child, ast.Call)}


def source_checks(I, sources):
    wrapper = sources[SOURCE]
    forbidden = {"read", "write", "open", "read_bytes", "read_text", "run", "prepare", "RECORD",
                 "OUT", "ROOT", "BOXFILE", "GEOM", "SELECTED", "point_packet_file"}
    assert not {node.id for node in ast.walk(ast.parse(wrapper)) if isinstance(node, ast.Name)} & forbidden
    assert not {node.attr for node in ast.walk(ast.parse(wrapper)) if isinstance(node, ast.Attribute)} & forbidden
    bulk = tree_function(wrapper, "bulk_component")
    assert {"U.series", "JET.contract_volumetric_defect", "U.uniform_cache",
            "U.RuleTarget", "functional.apply", "method.cache_clear"} <= call_names(bulk)
    surf = tree_function(wrapper, "surface_component")
    calls = call_names(surf)
    assert {"SUR.run_edge", "SUR.S.SurfaceFunctional.integral.cache_clear",
            "SUR.S.SurfaceFunctional.defect.cache_clear"} <= calls
    assert len([n for n in ast.walk(surf) if isinstance(n, ast.Call) and ast.unparse(n.func) == "SUR.run_edge"]) == 1
    assert "P.Functional" in ast.unparse(tree_function(wrapper, "_functional_type"))
    assert "AXbridge.PhysicalFunctional" in ast.unparse(tree_function(wrapper, "_functional_type"))
    run_edge = tree_function(sources["scripts/exact_uniform_surface_bridge_B_safe.py"], "run_edge")
    assert not call_names(run_edge) & {"read", "write", "sha", "Path"}
    assert {"uniform_series", "S.SurfaceFunctional", "cache_edge_correction",
            "JET.contract_volumetric_defect", "verify_frozen_domain",
            "S.binomial_tail", "material_boxes"} <= call_names(run_edge)
    for file, name in [
        ("scripts/exact_uniform_bulk_defect_B_safe.py", "uniform_cache"),
        ("scripts/exact_uniform_bulk_defect_B_safe.py", "series"),
        ("scripts/exact_axis_edge_physical.py", "physical_element"),
        ("scripts/exact_axis_edge_physical.py", "physical_fields"),
        ("scripts/exact_cache_consistency.py", "exact_bulk_maps"),
        ("scripts/exact_cache_consistency.py", "consistency_bound"),
        ("scripts/exact_surface_cache_consistency.py", "actual_maps"),
        ("scripts/exact_surface_cache_consistency.py", "consistency_bound"),
        ("scripts/exact_uniform_surface_bridge_B_safe.py", "cache_edge_correction"),
    ]:
        assert not {name.split(".")[-1] for name in call_names(tree_function(sources[file], name))} & forbidden
    assert I.JET.T is I.T and I.SUR.T is I.T and I.AX.T is I.T
    assert I.SUR.core is I.core and I.AX.core is I.core


def synthetic(I, corners, active):
    vertices = [tuple(map(Q, row)) for row in corners]
    nodes = vertices + [tuple((a + b) / 2 for a, b in zip(vertices[j], vertices[k]))
                        for j, k in ((0, 1), (1, 2), (2, 0))]
    target = I.T.Target.__new__(I.T.Target)
    target.a = {
        "triangles": np.array([range(6)], dtype=np.int64),
        "reference": np.array(nodes, dtype=np.float64),
        "shape": np.tile([0., 0., 0., 1/4, 1/2, 1/4], (64, 1)),
        "volume_weights": np.tile([1/8, 1/16], (1, 32)),
        "active": np.array([[active]], dtype=np.float64),
        "surface_0_N": np.tile([1/4, 1/2, 1/4], (12, 1)),
        "surface_0_weights": np.tile([1/32, 1/16, 1/8], (1, 4)),
    }
    target.index = {v: i for i, v in enumerate(2 * j + c for j, row in enumerate(nodes)
                                                for c in (0, 1) if c == 1 or row[0])}
    target.energy_scale = Q(11, 4)
    target.packet = {"parameters": {"mu": 1, "bulk": 5000}}
    cell = {k: (Q(1), Q(3)) for k in I.T.ORDER}
    cell["common_gain"] = (Q(7, 9), Q(11, 9))
    cell_doc = {k: list(map(str, v)) for k, v in cell.items()}
    cell_hash = digest(json.dumps(cell_doc, sort_keys=True, separators=(",", ":")).encode())
    context = SimpleNamespace(record_id="independent_fixture", target=target, dose=Q(9, 7), cell=cell,
                              row={"record_id": "independent_fixture", "dose_exact": "9/7",
                                   "parameter_cell": cell_doc, "parameter_cell_sha256": cell_hash,
                                   "model_observation_contract_sha256": "contract", "point_cache_sha256": "cache",
                                   "predictor_sha256": "predictor"})
    rule = {"volume_barycentric": [(Q(1, 4), Q(3, 8), Q(3, 8))] * 64,
            "volume_weights": [[Q(1, 22), Q(1, 44)] * 32],
            "surface_s": [[Q(1, 2)] * 12],
            "surface_weights": [[[Q(1, 88), Q(1, 44), Q(1, 22)] * 4]]}
    return context, [list(row) for row in nodes], rule


def assert_rejected(run, message=None):
    try:
        run()
    except ValueError as error:
        assert message is None or message in str(error), str(error)
    else:
        raise AssertionError("malformed fixture reached the numerical kernel")


def mock_bulk(I, chart, corners, active):
    C, Y, rule = synthetic(I, corners, active)
    target = C.target
    selector = I.prepare_bulk_selector(C, Y, rule, BODY, 0, Q(1, 10**12), Q(1, 10**10))
    assert selector["chart"] == chart
    f, ids = I._bulk_fields(target, Y, 0)
    assert ids == list(range(6))
    assert f["global_indices"] == target.element(0, Y)["global_indices"]
    n = len(f["global_indices"])
    ij = tuple(combinations_with_replacement(range(n), 2))
    rho = Q(selector["rho_exact"])
    assert rho * 2**48 == 14073748835533 and Q(1, 20) <= rho < Q(1, 20) + Q(1, 2**48)
    mass = Q(selector["geometry_two_mass_upper_exact"])
    raw_twice_area = abs((corners[1][0] - corners[0][0]) * (corners[2][1] - corners[0][1])
                         - (corners[2][0] - corners[0][0]) * (corners[1][1] - corners[0][1]))
    volume_geometry = Q(raw_twice_area, 6) * sum(row[0] for row in corners)
    aux = Q(24, 11)
    assert mass > aux + 2 * Q(314159, 100000) * volume_geometry / target.energy_scale
    assert mass < aux + 2 * Q(31416, 10000) * volume_geometry / target.energy_scale
    for trial in selector["order_attempts"]:
        f_tail, H_tail = (Q(trial[k]) for k in ("force_tail_upper_exact", "H_entry_tail_upper_exact"))
        fits = f_tail <= Q(1, 10**12) and H_tail <= Q(1, 10**10)
        assert fits == (trial["degree"] == selector["degree"])
    tail1, tail2 = (Q(selector[k]) for k in ("phi_prime_tail_upper_exact", "phi_second_tail_upper_exact"))
    assert tail1 > 0 and tail2 > 0
    forward = {"gradient": [arb(1 + i) for i in range(n)],
               "hessian_upper": {(i, j): arb(17 + i + 2 * j) for i, j in ij}}
    cache_g = [arb(-3 - 2 * i) for i in range(n)]
    cache_H = {(i, j): arb(9 - i - 3 * j) for i, j in ij}
    sample = (Q(3, 7), Q(2, 9)) if chart != "off_axis" else (Q(4, 3), Q(2, 9))
    as_ball = I.U.A
    checked = []

    def sample_poly(poly):
        return sum((as_ball(c) * as_ball(sample[0]**i * sample[1]**j) for (i, j), c in poly.items()), arb(0))

    def stub_functional(passed_target, fields, index, passed_rule):
        assert passed_target is target and passed_rule is rule and index == 0
        assert fields["global_indices"] == f["global_indices"]
        checked.append("functional")
        return SimpleNamespace(apply=sample_poly, defect=lambda *_: None, requested_defect_keys=set())

    def stub_cache(wrapped_target, fields, index, nodes, drive):
        assert wrapped_target.target is target and wrapped_target.auxiliary_rule() is rule
        assert nodes is Y and index == 0 and fields["free_local"] == f["free_local"]
        assert (drive.contains(as_ball(1)) and drive.contains(as_ball(Q(11, 7)))) if active else drive.is_zero()
        checked.append("cache")
        return cache_g, cache_H, {"all64_convex_strain_lines_positive": True}

    with patch.object(I, "_functional_type", return_value=stub_functional), \
            patch.object(I.U, "series", return_value=("mocked series A", "mocked series B")), \
            patch.object(I.JET, "contract_volumetric_defect", return_value=forward), \
            patch.object(I.U, "uniform_cache", side_effect=stub_cache):
        result = I.bulk_component(C, Y, rule, selector)
        assert checked == ["functional", "cache"]
        altered = {**selector, "rho_exact": "0"}
        assert_rejected(lambda: I.bulk_component(C, Y, rule, altered), "metadata changed")
        altered = deepcopy(rule)
        altered["volume_weights"][0][7] = Q(0)
        assert_rejected(lambda: I.bulk_component(C, Y, altered, selector), "weight")
        altered_nodes = deepcopy(Y)
        altered_nodes[1][1] += Q(1, 5)
        assert_rejected(lambda: I.bulk_component(C, altered_nodes, rule, selector), "nodes or rule changed")
        assert checked == ["functional", "cache"]

    def exact_sample(poly):
        return sum((Q(c) * sample[0]**i * sample[1]**j for (i, j), c in poly.items()), Q(0))

    strains = [exact_sample(p) for p in f["q"]]
    B = [[exact_sample(p) for p in row] for row in f["B"]]
    g = [1, Q(1), 1, Q(1), 1]
    assert len(result["force_integral_minus_cache"]) == n
    assert [(h["i"], h["j"]) for h in result["H_integral_minus_cache_upper"]] == list(ij)
    endpoints = (Q(1), Q(11, 7)) if active else (Q(0), Q(0))

    def packed_covers(doc, value):
        lo, hi = Q(doc["lower_exact"]), Q(doc["upper_exact"])
        assert lo <= value <= hi

    for drive in endpoints:
        g[1] = g[3] = Q(1) + drive
        for i, enclosure in enumerate(result["force_integral_minus_cache"]):
            force = sum(g[k] * strains[k] * B[k][i] for k in range(5)) + (1 + i) + (-3 - 2 * i)
            packed_covers(enclosure, force)
        for h in result["H_integral_minus_cache_upper"]:
            i, j = h["i"], h["j"]
            Hessian = sum(g[k] * B[k][i] * B[k][j] for k in range(5)) + (17 + i + 2 * j) + (9 - i - 3 * j)
            packed_covers(h["enclosure"], Hessian)
    assert all(info.currsize == 0 for info in (
        I.P.Functional.integral.cache_info(), I.P.Functional.defect.cache_info(),
        I.P.Functional.node_power.cache_info(), I.AXbridge.PhysicalFunctional.integral.cache_info()))
    return C, Y, rule


def run():
    if not __debug__:
        raise RuntimeError("Review replay requires enabled assertions")
    bindings = {}

    def raw(name, expected=None):
        path = ROOT / name
        assert path.resolve().is_relative_to(ROOT)
        data = path.read_bytes()
        val = digest(data)
        assert expected is None or val == expected, name
        assert bindings.setdefault(name, val) == val
        return data

    prep = PREP.relative_to(ROOT).as_posix()
    protocol = json.loads(raw(f"{prep}/protocol.json", FROZEN["protocol.json"]))
    preflight = json.loads(raw(f"{prep}/preflight.json", FROZEN["preflight.json"]))
    assert preflight["status"] == "PASS_INTEGRATOR_PREPARATION_ONLY"
    assert all(preflight["checks"].values()) and len(preflight["checks"]) == preflight["check_count"]
    assert preflight["protocol_sha256"] == FROZEN["protocol.json"]
    assert preflight["source_sha256"] == FROZEN[SOURCE]
    assert all(preflight[key] is False for key in ("actual_record_force_H_moment_expansions",
                                                   "actual_bulk_or_surface_cache_sums",
                                                   "assembly_or_root_assertion"))
    assert raw(SOURCE, FROZEN[SOURCE]) == raw(f"{prep}/exact_cohort_integrator_frozen.py", FROZEN[SOURCE])
    sources = {name: raw(name, h).decode() for name, h in protocol["source_hashes_and_inputs"].items()}
    context_review_file = "evidence/exact_cohort_context_independent_review.json"
    context_review = json.loads(raw(context_review_file, FROZEN[context_review_file]))
    assert context_review["status"] == "PASS" and context_review["records"] == 17
    assert protocol["precision_bits"] == preflight["precision_bits"] == 640
    sys.path.insert(0, str(ROOT / "scripts"))
    import exact_cohort_context as C
    import exact_cohort_integrator as I
    ctx.prec = 640
    I.core.A.cache_clear()
    source_checks(I, sources)
    checks = {"frozen_source_protocol_preflight_and_dependencies_match": True,
              "reviewed_numeric_dispatch_B_runner_and_receipt_paths_unreachable": True}
    c = C.load_context("optoGEF2_small_col_L")
    for name, h in c.input_hashes.items():
        raw(name, h)
    charts = {}
    y = c.target.nodes()
    for i, element in enumerate(c.target.a["triangles"]):
        actual, _ = I._bulk_fields(c.target, y, i)
        chart = actual["chart"]["kind"]
        expected = ["off_axis", "axis_vertex", "axis_edge_physical"][
            sum(c.target.a["reference"][n, 0] == 0 for n in element[:3])]
        assert chart == expected
        duffy = c.target.element(i, y)
        assert actual["global_indices"] == duffy["global_indices"] and actual["free_local"] == duffy["free_local"]
        charts[chart] = charts.get(chart, 0) + 1
    assert charts == {"off_axis": 132, "axis_vertex": 6, "axis_edge_physical": 6}
    checks["actual_L_connectivity_selects_all_144_expected_charts_and_free_maps"] = True

    fixtures = [("off_axis", [(1, 0), (3, 0), (1, 2)], 0),
                ("axis_vertex", [(0, 0), (2, 0), (2, 3)], 1),
                ("axis_edge_physical", [(0, 0), (2, 0), (0, 3)], 1)]
    for chart, corners, active in fixtures:
        fixture, nodes, rule = mock_bulk(I, chart, corners, active)
        checks[f"{chart}: independently mocked signed wrapper arithmetic and guards"] = True
    fields, _ = I._bulk_fields(fixture.target, nodes, 0)
    uncertain = arb(0, I.U.A(Q(1, 2)))
    assert I._norm({(0, 0): uncertain}, fields) >= Q(1, 2)
    assert (0, 0) in I.T.constant(uncertain) and (0, 0) not in I.T.constant(arb(0))
    bad_axis = deepcopy(nodes)
    bad_axis[0][0] = uncertain
    assert_rejected(lambda: I._bulk_fields(fixture.target, bad_axis, 0), "exactly zero")
    checks["interval_zero_crossing_preserved_and_axis_requires_exact_zero"] = True

    order = 53
    case = {"face": 0, "edge": 0, "order": order, "record_id": fixture.record_id}
    with patch.object(I.SUR, "run_edge", return_value={"tag": "mocked", "order": order, "scope": "old label"}) as edge:
        result = I.surface_component(fixture, nodes, rule, case)
        assert edge.call_count == 1
        assert all(x is y for x, y in zip(edge.call_args.args,
                                         (fixture.target, nodes, fixture.row, case, rule)))
        assert edge.call_args.args[5] == order
        assert result["tag"] == "mocked" and result["record_id"] == fixture.record_id
        assert result["dose_exact"] == "9/7" and result["order"] == order
        assert result["scope"].startswith("One surface edge on the reviewed explicit")
        assert_rejected(lambda: I.surface_component(fixture, nodes, rule, {**case, "record_id": "other"}), "different record")
        altered = deepcopy(rule)
        altered["surface_s"][0][0] = Q(2, 3)
        assert_rejected(lambda: I.surface_component(fixture, nodes, altered, case), "surface sites")
        ctx.prec = 512
        try:
            assert_rejected(lambda: I.surface_component(fixture, nodes, rule, case), "640-bit")
        finally:
            ctx.prec = 640
            I.core.A.cache_clear()
        assert edge.call_count == 1
    checks["surface_forwards_explicit_target_box_cell_rule_order_and_rejects_invalid_wrapper_inputs"] = True
    assert ctx.prec == 640
    assert I.SUR.S.SurfaceFunctional.integral.cache_info().currsize == 0
    assert I.SUR.S.SurfaceFunctional.defect.cache_info().currsize == 0
    checks["caller_precision_unchanged_and_functional_cache_cleanup_present"] = True
    replay_file = Path(__file__).relative_to(ROOT).as_posix()
    replay_hash = digest(raw(replay_file))
    for name, expected in tuple(bindings.items()):
        raw(name, expected)
    evidence = {
        "status": "PASS",
        "scope": "Generic local force and Hessian enclosure engine on caller-supplied reviewed domains. No actual nonzero moment defect, full 64-site correction, assembly or root process was run.",
        "source_sha256": FROZEN[SOURCE], "protocol_sha256": FROZEN["protocol.json"],
        "preflight_sha256": FROZEN["preflight.json"], "review_source_file": replay_file,
        "review_source_sha256": replay_hash, "checks": checks, "actual_L_chart_counts": charts,
        "input_hashes": dict(sorted(bindings.items())),
        "caller_owned_gates": [
            "Use the independently reviewed record/cell/dose, actual common box and certified positive-J body element for the selected index. A selector's internal J consistency is not proof of that domain.",
            "Bind each surface case, order, local site guards and whole-edge G/q0/rho to the exact box; record_id is optional in the numeric wrapper. Fresh proof margins and global assembly remain separate.",
            "Set 640-bit Arb precision and clear core.A before constructing interval nodes. The adapter checks the current precision; the caller owns startup cache hygiene. Use the declared isolated process policy.",
        ],
        "actual_force_H_moment_defects_run": False,
        "actual_full_cache_sums_run": False,
        "matrix_assembly_run": False,
    }
    if "--freeze" in sys.argv:
        receipt = OUT / "replay.json"
        assert not receipt.exists() and not EVIDENCE.exists()
        with receipt.open("x", encoding="utf8", newline="\n") as stream:
            stream.write(json.dumps(evidence, indent=2) + "\n")
        evidence["review_receipt_file"] = receipt.relative_to(ROOT).as_posix()
        evidence["review_receipt_sha256"] = digest(receipt.read_bytes())
        with EVIDENCE.open("x", encoding="utf8", newline="\n") as stream:
            stream.write(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps({"status": evidence["status"], "checks": len(checks), "actual_L_chart_counts": charts}))


if __name__ == "__main__":
    run()
