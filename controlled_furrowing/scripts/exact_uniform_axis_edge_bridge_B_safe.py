"""Six same-target uniform bulk replacements in physical axis-edge coordinates.

prepare freezes selectors/tests only. Actual six integrations/cache sums require
an independent exact-bound prelaunch review. No physical source is modified.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from types import SimpleNamespace
import argparse
import hashlib
import json
import time

from flint import arb, arb_mat, ctx

import exact_uniform_bulk_defect_B_safe as U
import exact_axis_edge_physical as AX

P, T, M, JET, core = U.P, U.T, U.M, U.JET, U.core
ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(__file__)
OUT = ROOT / "runs/exact_integration_bridge_v1/uniform_axis_edge_B_safe_v2"
REP = ROOT / "runs/exact_integration_bridge_v1/axis_edge_physical_preparation"
SELECTED = tuple(range(72, 78))
PRECISION = 512
FORCE_TAIL_TARGET = F(1, 10**22)
H_TAIL_TARGET = F(1, 10**16)
SCOPE = "SIX full-axis-edge contributions on SAME common B nodal domain/C, exact same target in physical polynomial coordinates. Not all144 replacement, assembled611 force/H bound, root certificate, PDE validation or a new model/concept."


def read(p):
    return json.loads(Path(p).read_text(encoding="utf8"))


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(p, value):
    P.write(p, value)


def A(value):
    if isinstance(value, arb):
        if not value.is_finite():
            raise ValueError("finite interval coefficient required")
        return value
    return core.A(T.rational(value))


def abs_coefficient(value):
    if isinstance(value, arb):
        if not value.is_finite():
            raise ValueError("finite polynomial coefficient required")
        return F(str(core.upper(value).fmpq()))
    return abs(T.rational(value))


def pnorm(poly, fields):
    """Ordinary physical polynomial majorant; no inverse-R routine at zero."""
    if fields["chart"]["kind"] != "axis_edge_physical":
        raise ValueError("explicit full-axis physical polynomial chart required")
    if any(isinstance(i, bool) or isinstance(j, bool) or not isinstance(i, int) or not isinstance(j, int) or i < 0 or j < 0 for i, j in poly):
        raise ValueError("ordinary nonnegative powers required, including at R=0")
    Rhi = max(v[0] for v in fields["corners"])
    Zmax = max(abs(v[1] - fields["zcenter"]) for v in fields["corners"])
    if Rhi <= 0:
        raise ValueError("positive material radial extent required")
    return sum((abs_coefficient(v) * Rhi**i * Zmax**j for (i, j), v in poly.items()), F(0))


class PhysicalFunctional(P.Functional):
    """Same signed moments/auxiliary rule, with ordinary physical integrals.

    The base constructor routes this explicit new chart through physical
    (R,zeta) auxiliary points. Dynamic dispatch uses this integral override
    even for the constructor's material mass.
    """
    @lru_cache(maxsize=None)
    def integral(self, a, b):
        if self.f["chart"]["kind"] != "axis_edge_physical" or any(isinstance(k, bool) or not isinstance(k, int) or k < 0 for k in (a, b)):
            raise ValueError("physical axis-edge ordinary moments only")
        return M.normalized_bulk_moment(self.vertices, a, b, self.t.energy_scale, self.f["zcenter"], ctx.prec)


def geometry_mass_upper(target, fields, index, rule):
    """Positive-measure geometry only; no actual field-defect moments in prep."""
    integral_R = fields["det"] * sum(v[0] for v in fields["corners"]) / 6
    auxiliary = sum(rule["volume_weights"][index], F(0))
    if integral_R <= 0 or auxiliary <= 0:
        raise ValueError("positive material and auxiliary measures required")
    return F(str(core.upper(target.normalization() * A(integral_R) + A(auxiliary)).fmpq()))


def choose(target, fields, index, rule, rho):
    rho = F(rho)
    if not 0 <= rho < 1:
        raise ValueError("same certified fresh-body scalar rho required")
    jet = JET.bulk_jacobian_jet(fields["q"], fields["B"])
    mass = geometry_mass_upper(target, fields, index, rule)
    g = max(pnorm(v, fields) for v in jet["gradient"])
    h = max(pnorm(v, fields) for v in jet["hessian_upper"].values())
    attempts = []
    for n in range(2, 19):
        t1, t2 = M.volumetric_derivative_tail_bounds(n, rho, 1, 5000)
        fb, hb = mass * t1 * g, mass * (t2 * g * g + t1 * h)
        attempts.append({"degree": n, "force_tail_upper_exact": str(fb), "H_entry_tail_upper_exact": str(hb)})
        if fb <= FORCE_TAIL_TARGET and hb <= H_TAIL_TARGET:
            return {"element": index, "rho_exact": str(rho), "degree": n, "force_tail_upper_exact": str(fb), "H_entry_tail_upper_exact": str(hb), "active": str(fields["active"]), "geometry_two_mass_upper_exact": str(mass), "order_attempts": attempts}
    raise ValueError("same declared degree2..18 tail policy unresolved")


def packed_overlap(a, b):
    return max(F(a["lower_exact"]), F(b["lower_exact"])) <= min(F(a["upper_exact"]), F(b["upper_exact"]))


def evaluate(target, index, Y, rule, selector, protocol):
    start = time.monotonic()
    fields = AX.physical_element(target, index, nodes=Y)
    if fields["global_indices"] != selector["global_indices"]:
        raise ValueError("same physical free-column map required")
    jet = JET.bulk_jacobian_jet(fields["q"], fields["B"])
    functional = PhysicalFunctional(target, fields, index, rule)
    rho, n = F(selector["rho_exact"]), selector["degree"]
    drive = U.gain_drive(protocol, fields["active"])
    polyA, polyB = U.series(T.sub(jet["value"], T.ONE), n)
    truncation = JET.contract_volumetric_defect(jet, polyA, polyB, functional.defect)
    mass = core.upper(functional.mass_integral + functional.mass_rule)
    t1, t2 = M.volumetric_derivative_tail_bounds(n, rho, 1, 5000)
    norms = [pnorm(p, fields) for p in jet["gradient"]]
    diag = [A(1), A(1) + drive, A(1), A(1) + drive, A(1)]
    d = len(fields["global_indices"])
    fg, HH, ft, ht = [], {}, [], {}
    for i in range(d):
        shear = {}
        for k in range(5):
            shear = T.add(shear, T.scale(T.mul(fields["q"][k], fields["B"][k][i]), diag[k]))
        fg.append(truncation["gradient"][i] + functional.apply(shear))
        ft.append(mass * A(t1 * norms[i]))
    for i, j in jet["hessian_upper"]:
        shear = {}
        for k in range(5):
            shear = T.add(shear, T.scale(T.mul(fields["B"][k][i], fields["B"][k][j]), diag[k]))
        HH[i, j] = truncation["hessian_upper"][i, j] + functional.apply(shear)
        ht[i, j] = mass * A(t2 * norms[i] * norms[j] + t1 * pnorm(jet["hessian_upper"][i, j], fields))
    cf, ch, cache_info = U.uniform_cache(U.RuleTarget(target, rule), fields, index, Y, drive)
    def inflate(value, tail):
        return value + arb(0, core.upper(tail))
    Fpacked = [core.pack(inflate(value, ft[i]) + cf[i]) for i, value in enumerate(fg)]
    Hpacked = [{"i": i, "j": j, "enclosure": core.pack(inflate(value, ht[i, j]) + ch[i, j]), "two_mass_tail_upper": core.pack(ht[i, j])} for (i, j), value in HH.items()]
    pointfile = ROOT / selector["point_packet_file"]
    if sha(pointfile) != selector["point_packet_sha256"]:
        raise ValueError("changed consumed SAFE point packet")
    point = read(pointfile)
    pointH = {(v["i"], v["j"]): v["enclosure"] for v in point["H_integral_minus_cache_upper_triangle"]}
    diagnostics = {
        "base_point_free_map_matches": point["global_indices"] == fields["global_indices"],
        "same_force_entry_count": len(Fpacked) == len(point["force_integral_minus_cache"]),
        "same_Hessian_entry_keys": {(x["i"], x["j"]) for x in Hpacked} == set(pointH),
        "all_base_point_force_entries_overlap": len(Fpacked) == len(point["force_integral_minus_cache"]) and all(packed_overlap(a, b) for a, b in zip(Fpacked, point["force_integral_minus_cache"])),
        "all_base_point_H_entries_overlap": all((v["i"], v["j"]) in pointH and packed_overlap(v["enclosure"], pointH[v["i"], v["j"]]) for v in Hpacked),
    }
    # Necessary consistency diagnostics, not a substitute for the uniform
    # construction or an exact integral-root theorem.
    return {
        "status": "SAME_TARGET_PHYSICAL_AXIS_EDGE_UNIFORM_COMPONENT_UNREVIEWED" if all(diagnostics.values()) else "UNRESOLVED_POINT_OVERLAP_DIAGNOSTIC",
        "element": index, "chart": fields["chart"]["kind"], "active": str(fields["active"]),
        "degree": n, "rho_exact": str(rho), "global_indices": fields["global_indices"],
        "force_integral_minus_cache": Fpacked,
        "H_integral_minus_cache_upper": Hpacked,
        "force_two_mass_tail_upper": [core.pack(v) for v in ft],
        "cache_consistency": cache_info, "point_overlap_diagnostics_only": diagnostics,
        "representation": "Exact physical P2 q/B, affine and ordinary nonnegative powers; physical J degree3. Same current nodal coefficients, free map, material domain, parameters, measure and cache.",
        "unique_moment_defects": len(functional.requested_defect_keys),
        "seconds": time.monotonic() - start,
        "scope": SCOPE,
    }


def validate_review(review, source_hash, protocol_hash, preflight_hash):
    if not isinstance(review, dict) or review.get("status") != "PASS" or not isinstance(review.get("checks"), dict) or not review["checks"] or not all(v is True for v in review["checks"].values()) or review.get("source_sha256") != source_hash or review.get("protocol_sha256") != protocol_hash or review.get("preflight_sha256") != preflight_hash:
        raise ValueError("nonempty literal-True exact-bound independent review required")


def tests():
    checks = {}
    def ck(name, value):
        if value is not True:
            raise AssertionError(name)
        checks[name] = True
    vertices = [[F(0), F(0)], [F(2), F(0)], [F(0), F(3)]]
    current = vertices + [[F(1), F(0)], [F(1), F(3, 2)], [F(0), F(3, 2)]]
    selected = [k for k in range(12) if k % 2 or current[k // 2][0] > 0]
    fields = AX.physical_fields(vertices, current, selected)
    fields["global_indices"] = list(range(len(selected)))
    toy_target = SimpleNamespace(energy_scale=F(7), normalization=lambda: 2 * arb.pi() / A(7))
    rule = {"volume_barycentric": [[F(1, 3)] * 3], "volume_weights": [[F(1, 11)]]}
    physical_fun = PhysicalFunctional(toy_target, fields, 0, rule)
    oldfields = T.element_fields(vertices, current, selected)
    old_fun = P.Functional(toy_target, oldfields, 0, rule)
    ck("identity_J_exact", JET.bulk_jacobian_jet(fields["q"], fields["B"])["value"] == T.ONE)
    ck("identity_material_mass_matches_Duffy", physical_fun.mass_integral.overlaps(old_fun.mass_integral))
    ck("identity_auxiliary_mass_identical", core.pack(physical_fun.mass_rule) == core.pack(old_fun.mass_rule))
    ck("physical_auxiliary_coordinates_exact", physical_fun.points[0][0] == (F(2, 3), F(0)))
    chart = oldfields["chart"]
    R, Z = chart["R"], T.sub(chart["Z"], T.constant(fields["zcenter"]))
    for i in range(4):
        for j in range(4):
            mapped = T.mul(T.power(R, i), T.power(Z, j))
            duffy = sum((v * M.axis_chart_moment("axis_edge", a, b, chart["R1"], fields["det"]) for (a, b), v in mapped.items()), F(0))
            physical = M.triangle_moment(vertices, i + 1, j, fields["zcenter"])
            ck(f"exact_synthetic_measure_identity_{i}_{j}", not physical.logarithms and physical.constant == duffy)
    ck("geometry_mass_formula_encloses_functional_masses", A(geometry_mass_upper(toy_target, fields, 0, rule)) >= core.upper(physical_fun.mass_integral + physical_fun.mass_rule))
    uncertain = A(F(-1, 3)).union(A(F(1, 2)))
    ck("interval_physical_polynomial_norm", pnorm({(2, 1): uncertain}, fields) >= F(4))
    ck("zero_crossing_polynomial_coefficient_retained", (0, 0) in T.constant(arb(0, 1)))
    for name, action in (
        ("negative_R_power_rejected", lambda: pnorm({(-1, 0): F(1)}, fields)),
        ("negative_zeta_power_rejected", lambda: physical_fun.integral(0, -1)),
        ("negative_R_moment_rejected", lambda: physical_fun.integral(-1, 0)),
        ("Boolean_moment_rejected", lambda: physical_fun.integral(True, 0)),
        ("nonfinite_coefficient_rejected", lambda: pnorm({(0, 0): arb("nan")}, fields)),
    ):
        try:
            action()
            ck(name, False)
        except ValueError:
            ck(name, True)
    valid = {"status": "PASS", "checks": {"one": True}, "source_sha256": "s", "protocol_sha256": "p", "preflight_sha256": "f"}
    validate_review(valid, "s", "p", "f")
    ck("valid_prelaunch_guard", True)
    for label, values in (("empty", {}), ("false", {"one": False}), ("integer", {"one": 1}), ("list", [True])):
        try:
            validate_review({**valid, "checks": values}, "s", "p", "f")
            ck("bad_prelaunch_rejected_" + label, False)
        except ValueError:
            ck("bad_prelaunch_rejected_" + label, True)
    try:
        validate_review(valid, "wrong", "p", "f")
        ck("wrong_source_prelaunch_rejected", False)
    except ValueError:
        ck("wrong_source_prelaunch_rejected", True)
    return checks


def prepare():
    if (OUT / "protocol.json").exists():
        raise ValueError("preserve frozen preparation")
    ctx.prec = PRECISION
    core.A.cache_clear()
    checks = tests()
    repr_protocol = read(REP / "protocol.json")
    repr_preflight = read(REP / "preflight.json")
    repr_tests = read(REP / "tests.json")
    if repr_preflight["source_sha256"] != sha(Path(AX.__file__)) or repr_preflight["protocol_sha256"] != sha(REP / "protocol.json") or repr_preflight["tests_sha256"] != sha(REP / "tests.json") or repr_tests["source_sha256"] != sha(Path(AX.__file__)) or repr_tests["protocol_sha256"] != sha(REP / "protocol.json") or not repr_tests["checks"] or not all(v is True for v in repr_tests["checks"].values()) or [x["element"] for x in repr_protocol["cases"]] != list(SELECTED):
        raise ValueError("exact-source six-element physical representation preparation required")
    geom_review = read(U.REVIEW)
    if not geom_review["status"].startswith("PASS") or not geom_review["checks"] or not all(v is True for v in geom_review["checks"].values()):
        raise ValueError("independently verified same fresh common geometry required")
    record = read(ROOT / P.RECORD)
    box = read(U.GEOM / "coordinate_boxes.json")
    body = read(U.GEOM / "certification_run/whole_body_geometry.json")
    if geom_review["body_output_sha256"] != sha(U.GEOM / "certification_run/whole_body_geometry.json") or geom_review["producer_results_sha256"] != sha(U.GEOM / "certification_run/results.json") or geom_review["protocol_sha256"] != sha(U.GEOM / "protocol.json"):
        raise ValueError("geometry review must bind exact current body/results/protocol")
    if box["parameter_cell_sha256"] != record["parameter_cell_sha256"] or repr_protocol["parameter_cell_sha256"] != record["parameter_cell_sha256"] or sha(U.GEOM / "coordinate_boxes.json") != repr_protocol["common_geometry_box_sha256"] or box["root_membership_assumed"] is not False:
        raise ValueError("SAME explicit common domain and C required")
    target = T.Target(read(ROOT / record["point_cache_file"]))
    Y, rule = U.field_nodes(box), target.auxiliary_rule()
    if [i for i, ids in enumerate(target.a["triangles"]) if sum(target.a["reference"][ids[:3], 0] == 0) == 2] != list(SELECTED):
        raise ValueError("same six axis-edge selectors required")
    point_nodes = target.nodes()
    checks["saved_base_point_is_in_same_common_nodal_box"] = all(F(c) - F(r) <= point_nodes[i][k] <= F(c) + F(r) for i, (cr, rr) in enumerate(zip(box["current_centers_exact"], box["coordinate_radii_exact"])) for k, (c, r) in enumerate(zip(cr, rr)))
    checks["base_point_drive_is_in_same_shared_gain_cell"] = F(record["parameter_cell"]["common_gain"][0]) <= F(record["base_strength_exact"]) / F(record["dose_exact"]) <= F(record["parameter_cell"]["common_gain"][1])
    bindings = {}
    def bind(file, expected=None):
        file = Path(file).resolve()
        digest = sha(file)
        if expected is not None and digest != expected:
            raise ValueError("changed input: " + str(file))
        name = file.relative_to(ROOT.resolve()).as_posix()
        if name in bindings and bindings[name] != digest:
            raise ValueError("conflicting input")
        bindings[name] = digest
    for name, digest in read(U.GEOM / "protocol.json")["source_and_input_hashes"].items():
        bind(ROOT / name, digest)
    for file in (SOURCE, Path(AX.__file__), Path(U.__file__), Path(P.__file__), Path(U.C.__file__), Path(T.__file__), Path(M.__file__), Path(JET.__file__), Path(core.__file__), U.REVIEW, U.GEOM / "protocol.json", U.GEOM / "certification_run/results.json", U.GEOM / "certification_run/whole_body_geometry.json", U.GEOM / "coordinate_boxes.json", ROOT / P.RECORD, ROOT / record["point_cache_file"], ROOT / "evidence/exact_cache_consistency_independent_review_v2.json", ROOT / "evidence/exact_safe_polynomial_kernel_independent_review.json", REP / "protocol.json", REP / "preflight.json", REP / "tests.json", REP / "PROOF.md"):
        bind(file)
    selectors = []
    for index in SELECTED:
        fields = AX.physical_element(target, index, nodes=Y)
        bounds = body["triangles"][index]["J_bound"]
        raw_rho = max(abs(F(bounds["lower_exact"]) - 1), abs(F(bounds["upper_exact"]) - 1))
        scaled = raw_rho * 2**48
        rho = F((scaled.numerator + scaled.denominator - 1) // scaled.denominator, 2**48)
        if not rho >= raw_rho or not 0 <= rho < 1 or F(bounds["lower_exact"]) <= F(1, 10**6):
            raise ValueError("same outward fresh-body analytic scalar radius required")
        selector = choose(target, fields, index, rule, rho)
        selector["global_indices"] = fields["global_indices"]
        packet = ROOT / f"runs/exact_integration_bridge_v1/bulk_assembly_B_safe_v2/elements/{index:03d}.json"
        bind(packet)
        point = read(packet)
        if point["global_indices"] != fields["global_indices"] or len(point["force_integral_minus_cache"]) != len(fields["global_indices"]) or len(point["H_integral_minus_cache_upper_triangle"]) != len(fields["global_indices"]) * (len(fields["global_indices"]) + 1) // 2:
            raise ValueError("matching complete SAFE point diagnostic packet required")
        selector.update({"point_packet_file": packet.relative_to(ROOT).as_posix(), "point_packet_sha256": sha(packet), "physical_q_degrees": [AX.degree(p) for p in fields["q"]], "physical_J_degree": AX.degree(JET.bulk_jacobian_jet(fields["q"], fields["B"])["value"])})
        selectors.append(selector)
        print(json.dumps({"prepared_axis_edge": index, "degree": selector["degree"], "rho_display_only": float(rho)}), flush=True)
    if not all(v is True for v in checks.values()):
        raise AssertionError(checks)
    protocol = {
        "source_and_input_hashes": bindings, "precision_bits": PRECISION,
        "parameter_cell": record["parameter_cell"], "parameter_cell_sha256": record["parameter_cell_sha256"], "dose_exact": record["dose_exact"],
        "geometry_domain_id": box["domain_id"], "geometry_boxes_sha256": sha(U.GEOM / "coordinate_boxes.json"),
        "selectors": selectors,
        "rho_tail_policy": "SAME fresh body per-element rho, ceil(rawrho*2**48)/2**48; degrees2..18; same force1e-22/Hentry1e-16 targets; no geometry/parameter change.",
        "integration_policy": "PhysicalFunctional overrides only integral using normalized_bulk_moment on ordinary nonnegative physical powers, including R=0. Inherited auxiliary points/defect/apply unchanged. Radial measure once; no Duffy measure/inverseR pnorm.",
        "cache_policy": "Reuse reviewed U.uniform_cache, SAME fullY exact coefficientwise deltaq, all64 sites, original global map/old-normalized weights.",
        "point_overlap_scope": "Diagnostic only, all6 consumed SAFE point packets hash bound; no root proof inferred from interval overlap.",
        "representation_preparation_protocol_sha256": sha(REP / "protocol.json"),
        "root_certified": False, "scope": SCOPE,
    }
    write(OUT / "protocol.json", protocol)
    preflight = {"status": "PASS_PHYSICAL_AXIS_EDGE_BRIDGE_PREPARATION_ONLY", "checks": checks, "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json"), "selected_elements": list(SELECTED), "selected_degrees": [v["degree"] for v in selectors], "actual_six_force_H_integrals_or_cache_sums_run": False, "preparation_geometry_only_measure_formula": True, "scope": SCOPE}
    write(OUT / "preflight.json", preflight)
    print(json.dumps({"status": preflight["status"], "checks": len(checks), "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json"), "preflight_sha256": sha(OUT / "preflight.json"), "degrees": preflight["selected_degrees"], "bindings": len(bindings)}), flush=True)


def run():
    p = read(OUT / "protocol.json")
    validate_review(read(OUT / "prelaunch_independent_review.json"), sha(SOURCE), sha(OUT / "protocol.json"), sha(OUT / "preflight.json"))
    for name, digest in p["source_and_input_hashes"].items():
        if sha(ROOT / name) != digest:
            raise ValueError("changed frozen input: " + name)
    ctx.prec = p["precision_bits"]
    core.A.cache_clear()
    record = read(ROOT / P.RECORD)
    box = read(U.GEOM / "coordinate_boxes.json")
    if sha(U.GEOM / "coordinate_boxes.json") != p["geometry_boxes_sha256"] or box["parameter_cell_sha256"] != p["parameter_cell_sha256"]:
        raise ValueError("unchanged SAME domain and C required")
    target = T.Target(read(ROOT / record["point_cache_file"]))
    Y, rule = U.field_nodes(box), target.auxiliary_rule()
    write(OUT / "started.json", {"utc": P.now(), "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json")})
    outputs = []
    for selector in p["selectors"]:
        result = evaluate(target, selector["element"], Y, rule, selector, p)
        file = OUT / f"element_{selector['element']:03d}.json"
        write(file, result)
        outputs.append({"element": selector["element"], "file": file.name, "sha256": sha(file), "status": result["status"]})
        print(json.dumps({"axis_edge_uniform_complete": selector["element"], "seconds": result["seconds"], "status": result["status"]}), flush=True)
    for name, digest in p["source_and_input_hashes"].items():
        if sha(ROOT / name) != digest:
            raise ValueError("input changed during six integrations: " + name)
    write(OUT / "results.json", {"status": "SIX_SAME_TARGET_PHYSICAL_AXIS_EDGE_COMPONENTS_UNREVIEWED", "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json"), "elements": outputs, "all_frozen_bindings_current": True, "all_point_diagnostics_pass": all(v["status"] == "SAME_TARGET_PHYSICAL_AXIS_EDGE_UNIFORM_COMPONENT_UNREVIEWED" for v in outputs), "assembled611_bound_or_root_proved": False, "scope": SCOPE})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("prepare", "run", "test"))
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    elif args.command == "run":
        run()
    else:
        ctx.prec = PRECISION
        core.A.cache_clear()
        checks = tests()
        print(json.dumps({"status": "PASS", "checks": len(checks)}))
