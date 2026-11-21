"""Preparation-only exact physical-P2 representation for six full-axis edges.

The physical model, nodal box, free coordinates and parameters are unchanged.
No moment integration, cache sums or root theorem is implemented here.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import datetime
import hashlib
import json

import numpy as np
from flint import arb, ctx

import exact_surface_bridge_safe as SAFE

T, JET, core = SAFE.T, SAFE.JET, SAFE.core
ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(__file__)
OUT = ROOT / "runs/exact_integration_bridge_v1/axis_edge_physical_preparation"
PREPARATION = ROOT / SAFE.PREP
BOXFILE = ROOT / "runs/exact_integration_bridge_v1/retube_B_geometry/coordinate_boxes.json"
PRECISION = 640


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    return json.loads(Path(p).read_text(encoding="utf8"))


def write(p, value):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        raise ValueError("preserve write-once output: " + str(p))
    p.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf8")


def coefficient(value):
    if isinstance(value, arb):
        if not value.is_finite():
            raise ValueError("finite current coordinate required")
        return value
    return T.rational(value)


def divide_R(poly):
    """Exact polynomial division; overlapping-zero coefficients are not zero."""
    poly = T.clean(poly)
    if any(i < 1 or j < 0 for i, j in poly):
        raise ValueError("exact polynomial divisibility by physical R required")
    return {(i - 1, j): value for (i, j), value in poly.items()}


def degree(poly):
    return max((i + j for i, j in poly), default=0)


def physical_fields(corners, current, free_local):
    """Same coherent physical P2 map; only the polynomial coordinates change.

    Every full-axis radial nodal value must be EXACT zero. Divide rational
    off-axis radial basis functions by R BEFORE applying interval coordinates,
    instead of forming and cancelling an interval radial polynomial at R=0.
    """
    g = T.material_basis(corners)
    if sum(p[0] == 0 for p in g["corners"]) != 2:
        raise ValueError("exactly two axis corner vertices required; not axis-vertex cells")
    if len(current) != 6 or any(len(row) != 2 for row in current):
        raise ValueError("six two-coordinate current nodes required")
    current = [[coefficient(v) for v in row] for row in current]
    free_local = list(free_local)
    if not free_local or len(set(free_local)) != len(free_local) or any(isinstance(k, bool) or not isinstance(k, (int, np.integer)) or not 0 <= k < 12 for k in free_local):
        raise ValueError("valid unique nonempty local free columns required")
    phi = g["phi"]
    gradients = [[T.diff(p, axis) for axis in (0, 1)] for p in phi]
    axis_nodes = [j for j, p in enumerate(phi) if any(i == 0 for i, _ in p)]
    if len(axis_nodes) != 3:
        raise ValueError("full axis P2 edge must have three axis material nodes")
    for j in axis_nodes:
        if not T.is_exact_zero(current[j][0]):
            raise ValueError("axis radial current coordinate must be exactly zero, not merely contain zero")
        if 2 * j in free_local:
            raise ValueError("axis radial coordinates must be constrained, not free")
    hoop_bases = [None if j in axis_nodes else divide_R(p) for j, p in enumerate(phi)]
    q = [
        T.linear_combination([p[0] for p in gradients], [p[0] for p in current]),
        T.linear_combination([p[1] for p in gradients], [p[0] for p in current]),
        T.linear_combination([p[0] for p in gradients], [p[1] for p in current]),
        T.linear_combination([p[1] for p in gradients], [p[1] for p in current]),
        T.linear_combination([p for p in hoop_bases if p is not None], [current[j][0] for j in range(6) if hoop_bases[j] is not None]),
    ]
    B = [[] for _ in range(5)]
    for k in free_local:
        j, component = divmod(int(k), 2)
        radial = component == 0
        column = [gradients[j][0] if radial else {}, gradients[j][1] if radial else {}, {} if radial else gradients[j][0], {} if radial else gradients[j][1], hoop_bases[j] if radial else {}]
        if any(p is None for p in column):
            raise ValueError("free radial axis basis reached")
        for row, p in zip(B, column):
            row.append(p)
    polynomials = q + [p for row in B for p in row]
    if any(i < 0 or j < 0 for p in polynomials for i, j in p) or any(degree(p) > 1 for p in polynomials):
        raise ValueError("physical axis-edge q/B must be affine ordinary polynomials")
    return {
        **g, "q": q, "B": B, "free_local": list(map(int, free_local)),
        "axis_local_nodes": axis_nodes,
        "chart": {"kind": "axis_edge_physical", "R": T.X, "jacobian": T.ONE},
        "variables": ("R", "zeta"), "measure_polynomial": T.X,
        "integration_contract": "Physical affine-reference triangle including R=0; ordinary NONNEGATIVE-exponent moments only, radial measure counted once. No inverted-R off-axis norm routine is admissible at this axis.",
    }


def physical_element(target, index, nodes=None):
    if isinstance(index, bool) or not isinstance(index, (int, np.integer)) or not 0 <= index < 144:
        raise ValueError("explicit valid integer element required")
    index = int(index)
    ids = list(map(int, target.a["triangles"][index]))
    reference = target.a["reference"][ids[:3]]
    if sum(p[0] == 0 for p in reference) != 2:
        raise ValueError("only full-axis-edge elements accepted")
    nodes = target.nodes() if nodes is None else nodes
    dofs = [2 * node + component for node in ids for component in (0, 1)]
    free = [j for j, k in enumerate(dofs) if k in target.index]
    fields = physical_fields(reference, [nodes[j] for j in ids], free)
    fields.update({"element": index, "nodes": ids, "global_indices": [target.index[dofs[j]] for j in free], "active": T.rational(target.a["active"][index, 0])})
    return fields


def hull(value):
    value = coefficient(value)
    if isinstance(value, arb):
        return F(str(value.lower().fmpq())), F(str(value.upper().fmpq()))
    return value, value


def pack_poly(p):
    return [{"R_power": i, "zeta_power": j, "coefficient_exact_bounds": list(map(str, hull(v)))} for (i, j), v in sorted(p.items())]


def nodal_box(box):
    centers = [[F(x) for x in row] for row in box["current_centers_exact"]]
    radii = [[F(x) for x in row] for row in box["coordinate_radii_exact"]]
    if len(centers) != 325 or len(radii) != 325 or any(len(row) != 2 for row in centers + radii) or any(x < 0 for row in radii for x in row):
        raise ValueError("same finite ordered 325-node box required")
    return [[core.A(c - r).union(core.A(c + r)) for c, r in zip(cr, rr)] for cr, rr in zip(centers, radii)]


def tests():
    checks = {}
    def ck(name, value):
        if value is not True:
            raise AssertionError(name)
        checks[name] = True
    ref = [[F(0), F(0)], [F(1), F(0)], [F(0), F(1)]]
    current = ref + [[F(1, 2), F(0)], [F(1, 2), F(1, 2)], [F(0), F(1, 2)]]
    selected = [k for k in range(12) if k % 2 or current[k // 2][0] > 0]
    field = physical_fields(ref, current, selected)
    ck("identity_q_exact", field["q"] == [T.ONE, {}, {}, T.ONE, T.ONE])
    ck("identity_J_exact", JET.bulk_jacobian_jet(field["q"], field["B"])["value"] == T.ONE)
    ck("explicit_new_chart_not_offaxis_inverse_chart", field["chart"]["kind"] == "axis_edge_physical")
    tiny = arb(0, core.A(F(1, 10**80)))
    ck("nonexact_zero_crossing_retained", (0, 0) in T.clean({(0, 0): tiny}))
    ck("R_multiple_zero_crossing_division_retained", (0, 0) in divide_R({(1, 0): tiny}))
    near = T.sub(T.constant(arb(1, core.A(F(1, 10**80)))), T.ONE)
    ck("near_cancelled_coefficient_retained", bool(near) and not near[(0, 0)].is_zero())
    for name, action in (
        ("zero_crossing_constant_not_divisible", lambda: divide_R({(0, 0): tiny})),
        ("exact_nonzero_constant_not_divisible", lambda: divide_R(T.ONE)),
        ("nonfinite_coefficient", lambda: divide_R({(1, 0): arb("nan")})),
        ("axis_vertex_rejected", lambda: physical_fields([[F(0), F(0)], [F(1), F(0)], [F(1), F(1)]], current, selected)),
        ("axis_radial_free_rejected", lambda: physical_fields(ref, current, [0] + selected)),
        ("Boolean_free_index_rejected", lambda: physical_fields(ref, current, [True])),
    ):
        try:
            action()
            ck(name, False)
        except ValueError:
            ck(name, True)
    perturbed = [row.copy() for row in current]
    perturbed[0][0] = tiny
    try:
        physical_fields(ref, perturbed, selected)
        ck("overlapping_zero_axis_coordinate_rejected", False)
    except ValueError:
        ck("overlapping_zero_axis_coordinate_rejected", True)
    perturbed = [row.copy() for row in current]
    perturbed[1][0] = core.A(1).union(core.A(F(101, 100)))
    uniform = physical_fields(ref, perturbed, selected)
    ck("offaxis_interval_coordinate_propagates_to_hoop", any(not v.is_exact() for v in uniform["q"][4].values()))
    return checks


def prepare():
    if (OUT / "protocol.json").exists():
        raise ValueError("preserve frozen preparation")
    ctx.prec = PRECISION
    core.A.cache_clear()
    checks = tests()
    prep, box = read(PREPARATION), read(BOXFILE)
    row = next(x for x in prep["records"] if x["record_id"] == "optoGEF2_small_col_B")
    if box["parameter_cell_sha256"] != row["parameter_cell_sha256"] or box["root_membership_assumed"] is not False:
        raise ValueError("same C, explicit common set without root membership required")
    if sha(ROOT / row["point_cache_file"]) != row["point_cache_sha256"]:
        raise ValueError("unchanged B cache required")
    target = T.Target(read(ROOT / row["point_cache_file"]))
    if box["reference_nodes_exact"] != [[str(T.rational(x)) for x in row] for row in target.a["reference"]] or box["triangles"] != target.a["triangles"].tolist():
        raise ValueError("same reference and connectivity required")
    point = target.nodes()
    uniform = nodal_box(box)
    selected = [i for i, ids in enumerate(target.a["triangles"]) if sum(target.a["reference"][ids[:3], 0] == 0) == 2]
    if len(selected) != 6:
        raise ValueError("exactly six full-axis-edge elements required")
    bindings = {}
    def bind(path, expected=None):
        path = Path(path).resolve()
        digest = sha(path)
        if expected is not None and digest != expected:
            raise ValueError("changed source/input: " + str(path))
        key = path.relative_to(ROOT.resolve()).as_posix()
        if key in bindings and bindings[key] != digest:
            raise ValueError("conflicting dependency")
        bindings[key] = digest
    for group in (prep["source_hashes"], prep["input_hashes"]):
        for name, digest in group.items():
            bind(ROOT / name, digest)
    for path in (SOURCE, Path(T.__file__), Path(JET.__file__), Path(core.__file__), PREPARATION, BOXFILE, ROOT / "evidence/exact_safe_polynomial_kernel_independent_review.json", ROOT / "runs/exact_integration_bridge_v1/uniform_bulk_pilot_B_safe_v2_independent_results/axis_edge_basis_probe.py", ROOT / "runs/exact_integration_bridge_v1/uniform_bulk_pilot_B_safe_v2_independent_results/axis_edge_basis_identity.json"):
        bind(path)
    rng = np.random.default_rng(719243)
    cases = []
    for index in selected:
        old = target.element(index, nodes=point)
        physical = physical_element(target, index, nodes=point)
        interval_field = physical_element(target, index, nodes=uniform)
        chart = old["chart"]
        Rchart, Zchart = chart["R"], T.sub(chart["Z"], T.constant(old["zcenter"]))
        for k in range(5):
            checks[f"same_point_q_chart_identity_{index}_{k}"] = T.compose(physical["q"][k], Rchart, Zchart) == old["q"][k]
            for j in range(len(physical["free_local"])):
                checks[f"same_free_basis_chart_identity_{index}_{k}_{j}"] = T.compose(physical["B"][k][j], Rchart, Zchart) == old["B"][k][j]
        jp = JET.bulk_jacobian_jet(physical["q"], physical["B"])
        jo = JET.bulk_jacobian_jet(old["q"], old["B"])
        checks[f"same_point_J_chart_identity_{index}"] = T.compose(jp["value"], Rchart, Zchart) == jo["value"]
        for j, p in enumerate(jp["gradient"]):
            checks[f"same_point_Jgradient_identity_{index}_{j}"] = T.compose(p, Rchart, Zchart) == jo["gradient"][j]
        for (i, j), p in jp["hessian_upper"].items():
            checks[f"same_point_JHessian_identity_{index}_{i}_{j}"] = T.compose(p, Rchart, Zchart) == jo["hessian_upper"][i, j]
        allpolys = interval_field["q"] + [p for row in interval_field["B"] for p in row]
        checks[f"all_uniform_q_B_affine_nonnegative_powers_{index}"] = all(degree(p) <= 1 and all(i >= 0 and j >= 0 for i, j in p) for p in allpolys)
        ju = JET.bulk_jacobian_jet(interval_field["q"], interval_field["B"])
        checks[f"uniform_J_degree_at_most_three_{index}"] = degree(ju["value"]) <= 3
        checks[f"all_uniform_J_derivative_powers_nonnegative_{index}"] = all(i >= 0 and j >= 0 for p in [ju["value"]] + ju["gradient"] + list(ju["hessian_upper"].values()) for i, j in p)
        ids = physical["nodes"]
        centers = [[F(v) for v in box["current_centers_exact"][node]] for node in ids]
        radii = [[F(v) for v in box["coordinate_radii_exact"][node]] for node in ids]
        # Every q coefficient is linear in nodal coordinates. Verify its two
        # exact extremizers over the original full local box, beyond samples.
        full_basis = physical_fields(physical["corners"], centers, physical["free_local"])
        for component in range(5):
            support = set(interval_field["q"][component])
            for monomial in support:
                lo = hi = F(0)
                for local_node in range(6):
                    phi = physical["phi"][local_node]
                    if component == 4:
                        basis = {} if local_node in physical["axis_local_nodes"] else divide_R(phi)
                        coord = 0
                    else:
                        basis = T.diff(phi, component % 2)
                        coord = component // 2
                    value = basis.get(monomial, F(0))
                    c, r = centers[local_node][coord], radii[local_node][coord]
                    endpoints = value * (c - r), value * (c + r)
                    lo += min(endpoints)
                    hi += max(endpoints)
                observed = interval_field["q"][component][monomial]
                hlo, hhi = hull(observed)
                checks[f"linear_q_coefficient_box_extrema_{index}_{component}_{monomial}"] = hlo <= lo <= hi <= hhi
        sample_count = 0
        for sample in range(33):
            sampled = [row.copy() for row in point]
            for j, node in enumerate(ids):
                sampled[node] = [centers[j][c] + (F(0) if sample == 0 else F(int(rng.choice((-1, 1))))) * radii[j][c] for c in range(2)]
            candidate = physical_element(target, index, nodes=sampled)
            for component in range(5):
                outer = interval_field["q"][component]
                for monomial, value in candidate["q"][component].items():
                    lo, hi = hull(outer.get(monomial, F(0)))
                    checks[f"sampled_point_coefficient_contains_{index}_{sample}_{component}_{monomial}"] = lo <= value <= hi
            sample_count += 1
        if not all(v is True for k, v in checks.items() if f"_{index}_" in k or k.endswith("_" + str(index))):
            raise AssertionError("axis representation identity/enclosure check failed: " + str(index))
        cases.append({
            "element": index, "nodes": ids, "global_indices": physical["global_indices"],
            "free_local": physical["free_local"], "axis_local_nodes": physical["axis_local_nodes"],
            "physical_q_interval_coefficients": [pack_poly(p) for p in interval_field["q"]],
            "physical_B_exact_coefficients": [[pack_poly(p) for p in row] for row in interval_field["B"]],
            "point_q_degrees": [degree(p) for p in physical["q"]],
            "physical_J_degree": degree(ju["value"]), "Duffy_J_degree": degree(jo["value"]),
            "physical_J_support": len(ju["value"]), "Duffy_J_support": len(jo["value"]),
            "sampled_common_box_fields": sample_count,
            "integration_chart_kind": "axis_edge_physical", "moments_evaluated": False,
        })
        print(json.dumps({"prepared_axis_edge": index, "physical_J_degree": degree(ju["value"]), "Duffy_J_degree": degree(jo["value"]), "sample_fields": sample_count}), flush=True)
    if not all(v is True for v in checks.values()):
        raise AssertionError("representation preparation check failed")
    protocol = {
        "version": 1, "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "source_and_input_hashes": bindings, "selected_record": row,
        "common_geometry_box_file": BOXFILE.relative_to(ROOT).as_posix(), "common_geometry_box_sha256": sha(BOXFILE),
        "parameter_cell_sha256": row["parameter_cell_sha256"], "common_geometry_domain_id": box["domain_id"],
        "precision_bits": PRECISION, "cases": cases,
        "fixed_seed": 719243, "samples_per_element": 33,
        "integration_contract": "Use exact physical triangle ordinary nonnegative-exponent moments with R measure exactly once. Explicit axis_edge_physical chart; do not call inverted-R off-axis norms or generic axis-vertex formulas.",
        "same_target": "Frozen binary current nodal coefficients, exact P2 basis on affine binary-corner reference, same free columns, SAME common box/C. Symbolic exact-R basis division before intervalization only.",
        "root_membership_assumed": False, "actual_moment_integrals_run": False, "actual_cache_jobs_run": False,
        "scope": "Representation and enclosure preparation for six full-axis-edge elements only; not axis-vertex cells, an improved error bound, a root, a physical model change or a completed five-gate result.",
    }
    write(OUT / "protocol.json", protocol)
    write(OUT / "tests.json", {"status": "PASS_AXIS_EDGE_PHYSICAL_PREPARATION_ONLY", "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json"), "checks": checks})
    preflight = {"status": "PASS_AXIS_EDGE_PHYSICAL_PREPARATION_ONLY", "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json"), "tests_sha256": sha(OUT / "tests.json"), "checks": {"all_literal_true": all(v is True for v in checks.values()), "all_six_axis_edges": len(cases) == 6, "no_actual_moments": True, "no_actual_cache_jobs": True, "same_common_box_and_C": True}, "test_count": len(checks), "scope": protocol["scope"]}
    write(OUT / "preflight.json", preflight)
    print(json.dumps({"status": preflight["status"], "checks": len(checks), "source_sha256": sha(SOURCE), "protocol_sha256": sha(OUT / "protocol.json"), "preflight_sha256": sha(OUT / "preflight.json"), "bindings": len(bindings)}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("test", "prepare"))
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    else:
        ctx.prec = PRECISION
        core.A.cache_clear()
        checks = tests()
        print(json.dumps({"status": "PASS", "checks": len(checks)}))
