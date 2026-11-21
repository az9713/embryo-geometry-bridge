"""Independent preparation replay of all 144 centered L bulk force components."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement, product
import json
from math import ceil
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

from flint import arb

OUT = Path(__file__).resolve().parent
STAGE = OUT.parent
ROOT = STAGE.parents[3]
PREP = STAGE / "preparation"
SOURCE = "scripts/exact_centered_bulk_L.py"
EXPECTED = {
    "source_sha256":"a83c89ffc302be341df220fcb97c5be82f191a70e317fea9af839d2917a87492",
    "protocol_sha256":"466048ec456499aa904dd50daa744245aff2eecb6fc2a30fb15bec3d97d23f9b",
    "preflight_sha256":"2ed81adb1f8dd7e52d7c49685be49e508af4ee395911c0e292a74c67148632e0",
    "selectors_sha256":"8fbb0e8844d7f47f119b795e453032ae46d24a53ea489d2ce9fe2119576a9bf6",
}


def canon(value):
    return sha256(json.dumps(value,default=str,separators=(",",":"),sort_keys=True).encode()).hexdigest()


def bounds(row):
    a,b = (Q(row[k]) for k in ("lower_exact","upper_exact"))
    assert a <= b
    return a,b


def overlap(a,b):
    x,y = max(a[0],b[0]),min(a[1],b[1])
    assert x <= y
    return x,y


def matrix(packet,n):
    rows = {(entry["i"],entry["j"]):bounds(entry["enclosure"]) for entry in packet["H_integral_minus_cache_upper"]}
    keys = set(combinations_with_replacement(range(n),2))
    assert len(rows) == len(packet["H_integral_minus_cache_upper"]) == len(keys) and rows.keys() == keys
    return rows


def form_derivatives(corners,l):
    a,b = (corners[j][0]-corners[0][0] for j in (1,2))
    c,d = (corners[j][1]-corners[0][1] for j in (1,2))
    det = a*d-b*c
    Lz = ((b-a)/det,-b/det,a/det)
    return [(4*l[i]-1)*Lz[i] for i in range(3)] + [
        4*(Lz[i]*l[j]+l[i]*Lz[j]) for i,j in ((0,1),(1,2),(2,0))]


def material_coordinates(fields,l):
    vertices = fields["corners"]
    if fields["chart"]["kind"] == "axis_vertex":
        i,j = [i for i,p in enumerate(vertices) if p[0]>0]
        s = l[i]+l[j]
        return s,l[j]/s
    return sum(l[i]*v[0] for i,v in enumerate(vertices)),sum(l[i]*v[1] for i,v in enumerate(vertices))-fields["zcenter"]


def covers(C,value,expected):
    lo,hi = C.E.SUR.rational_hull(value)
    assert lo <= expected <= hi


def synthetic_stiffness(C):
    base = [(Q(0),Q(-1)),(Q(2),Q(1)),(Q(2),Q(3))]
    kinds = {"axis_vertex":base,"off_axis":[(r+1,z) for r,z in base],
             "axis_edge_physical":[(Q(0),Q(-1)),(Q(2),Q(1)),(Q(0),Q(3))]}
    dose = Q(5,7)
    for chart,verts in kinds.items():
        ref = verts + [tuple((a+b)/2 for a,b in zip(verts[i],verts[j]))
                       for i,j in ((0,1),(1,2),(2,0))]
        now = [(R*Q(11,9)+R*Z/30, Z+Q(3,7)+R/17) for R,Z in ref]
        free = [k for k in range(12) if (k%2 == 1 or ref[k//2][0] > 0) and k not in (3,6)]
        fields = C.E.AX.physical_fields(verts,now,free) if chart == "axis_edge_physical" else C.T.element_fields(verts,now,free)
        Y = [v for p in now for v in p]
        cells = [(Q(2,7),Q(3,7),Q(2,7)),(Q(1,4),Q(1,4),Q(1,2))]
        weights = [Q(1,7),Q(-2,9)]
        site_points = [material_coordinates(fields,l) for l in cells]
        gradients = [form_derivatives(fields["corners"],l) for l in cells]
        def exact_signed_functional(poly):
            return C.A(sum(w*C.T.evaluate(poly,*site) for w,site in zip(weights,site_points)))
        fun = SimpleNamespace(apply=exact_signed_functional)
        maps = []
        for qidx,d in enumerate(gradients):
            star = [[Q(0)]*12 for _ in range(5)]
            cached = [[Q(0)]*12 for _ in range(5)]
            for i,v in enumerate(d):
                for parity,m in enumerate((1,3)):
                    star[m][2*i+parity] = v
                    cached[m][2*i+parity] = v*(1+Q((-1)**(i+qidx+parity),18+i))
            maps.append({"star_full":star,"cache_full":cached,"auxiliary_normalized_weight_exact":Q(qidx+1,17)})
        calculated = C.affine_gain_derivative(dose,1,fields,fun,maps,list(map(C.A,Y)))
        dropped = [Q(0) if i not in free else y for i,y in enumerate(Y)]

        def exact_derivative(Y0):
            output = []
            for i in free:
                parity = i%2
                j = 1 if parity == 0 else 3
                def energy(Y1):
                    dot = lambda first,second:sum(v*w for v,w in zip(first,second))
                    square = lambda row:dot(row,Y1)**2
                    material = sum(w*dot(d,Y1[parity::2])**2 for w,d in zip(weights,gradients))
                    finite = sum(m["auxiliary_normalized_weight_exact"] *
                                 (square(m["star_full"][j])-square(m["cache_full"][j])) for m in maps)
                    return dose*(material+finite)/2
                plus,minus = Y0.copy(),Y0.copy()
                plus[i] += 1
                minus[i] -= 1
                output.append((energy(plus)-energy(minus))/2)
            return output

        exact = exact_derivative(Y)
        assert any(a != b for a,b in zip(exact,exact_derivative(dropped)))
        assert any(v > 0 or v < 0 for v in exact)
        for v,expect in zip(calculated,exact):
            covers(C,v,expect)
        fail = SimpleNamespace(apply=lambda p:(_ for _ in ()).throw(AssertionError("inactive element used moments")))
        assert all(v.is_zero() for v in C.affine_gain_derivative(dose,0,fields,fail,[],list(map(C.A,Y))))


def synthetic_two_segments(C):
    force = lambda x,y,g:((2+g)*x-3*y+x*x, -3*x+(5-g)*y)
    x0,y0,g0 = Q(1,2),Q(-1,4),Q(7,8)
    X,Y,G = (Q(1,4),Q(3,4)),(Q(-3,4),Q(1,4)),(Q(1,2),Q(1))
    H = {(0,0):(2+G[0]+2*X[0],2+G[1]+2*X[1]),(0,1):(Q(-3),Q(-3)),(1,1):(5-G[1],5-G[0])}
    ball = lambda a,b:C.A(a).union(C.A(b))
    center,_,_ = C.centered_force(list(map(C.A,force(x0,y0,g0))),
                                  [ball(X[0]-x0,X[1]-x0),ball(Y[0]-y0,Y[1]-y0)],H,
                                  [ball(*X),ball(-Y[1],-Y[0])],ball(G[0]-g0,G[1]-g0))
    for x,y,g in product(X,Y,G):
        for have,want in zip(center,force(x,y,g)):
            covers(C,have,want)
    crossing = arb(0,C.A(Q(1,16)))
    assert not C.T.is_exact_zero(crossing)
    assert C.F._series_first({(1,0):crossing},6).keys() == {(i,0) for i in range(7)}


def all_selectors(C,data,stored,read):
    seen = Counter()
    active = Counter()
    checked_F,checked_H = 0,0
    for i,wanted in enumerate(stored):
        proof,fields,jet,nodes,case,H,F,radii = C._element(data,i)
        assert canon(proof) == canon(wanted) and proof["element"] == i
        names = [int(v) for v in data.context.target.a["triangles"][i]]
        assert nodes == names == wanted["local_nodes"]
        chart = fields["chart"]["kind"]
        axis_vertices = sum(Q(fields["corners"][j][0]) == 0 for j in range(3))
        assert chart == ("off_axis","axis_vertex","axis_edge_physical")[axis_vertices]
        seen[chart] += 1
        active[chart] += int(wanted["active"])
        global_index = {int(dof):i for i,dof in enumerate(data.context.free)}
        dofs = [2*node+axis for node in nodes for axis in (0,1)]
        local = [j for j,x in enumerate(dofs) if x in global_index]
        idx = [global_index[dofs[j]] for j in local]
        n = len(idx)
        assert idx == proof["global_indices"] == fields["global_indices"]
        assert local == proof["free_local"] == fields["free_local"]
        assert tuple(radii) == tuple(Q(data.box["coordinate_radii_exact"][node][axis])
                                     for j,(node,axis) in enumerate(product(nodes,range(2))) if j in local)
        lo,hi = bounds(data.body["triangles"][i]["J_bound"])
        assert lo > Q(1,10**6) and lo <= hi
        r_raw = max(abs(1-lo),abs(1-hi))
        assert Q(wanted["rho_exact"]) == Q(ceil(r_raw*2**48),2**48) < 1
        assert wanted["order"] == case["degree"]
        assert Q(wanted["point_force_tail_max_exact"]) <= Q(wanted["original_force_tail_budget_exact"])
        assert all(C.T.evaluate(p,*material_coordinates(fields,data.rule["volume_barycentric"][i%64])) == v
                   for p,v in zip(C.axial_gradient(fields),form_derivatives(fields["corners"],data.rule["volume_barycentric"][i%64])))
        filename = f"{C.ORIGINAL}/bulk/element_{i:03d}.json"
        old = read(filename,data.old_files[filename])
        f0 = list(map(bounds,old["force_integral_minus_cache"]))
        h0 = matrix(old,n)
        if i in data.scaled_files:
            linked = data.scaled_files[i]
            scaled = read(linked["file"],linked["sha256"])
            assert scaled["same_Y_cache_consistency"] == old["cache_consistency"]
            h1 = matrix(scaled,n)
            f1 = list(map(bounds,scaled["force_integral_minus_cache"]))
            f0 = [overlap(v,w) for v,w in zip(f0,f1)]
            h0 = {key:overlap(v,h1[key]) for key,v in h0.items()}
        assert len(f0) == len(F) == n and f0 == list(map(bounds,F))
        assert h0 == H and len(H) == n*(n+1)//2
        checked_F += n
        checked_H += len(H)
    assert dict(seen) == {"off_axis":132,"axis_vertex":6,"axis_edge_physical":6}
    assert sum(active.values()) == 24
    return seen,active,checked_F,checked_H


def run():
    if not __debug__:
        raise RuntimeError("Review requires enabled assertions")
    bindings = {}

    def read(name,digest=None,parse=True):
        file = (ROOT/name).resolve()
        assert file.is_relative_to(ROOT)
        path = file.relative_to(ROOT).as_posix()
        value = file.read_bytes()
        found = sha256(value).hexdigest()
        assert digest is None or digest == found, path
        assert bindings.setdefault(path,found) == found, path
        return json.loads(value) if parse else value

    prefix = PREP.relative_to(ROOT).as_posix()
    source = read(SOURCE,EXPECTED["source_sha256"],False)
    assert source == read(f"{prefix}/exact_centered_bulk_L_frozen.py",EXPECTED["source_sha256"],False)
    p,pre,selectors = [read(f"{prefix}/{key}.json",EXPECTED[f"{key}_sha256"]) for key in ("protocol","preflight","selectors")]
    assert p["source_sha256"] == pre["source_sha256"] == EXPECTED["source_sha256"]
    assert p["selectors_sha256"] == pre["selectors_sha256"] == EXPECTED["selectors_sha256"]
    assert pre["protocol_sha256"] == EXPECTED["protocol_sha256"]
    assert pre["status"] == "PASS_CENTERED_BULK_L_PREPARATION_ONLY"
    assert len(pre["checks"]) == 20 and all(value is True for value in pre["checks"].values())
    assert pre["actual_record_force_integrals_run"] is pre["actual_cache_site_sums_run"] is False
    assert [s["element"] for s in selectors] == list(range(144))
    for path,digest in p["input_hashes"].items():
        read(path,digest,False)
    assert len(p["failed_attempts"]) == 2
    for path,digest in p["failed_attempts"].items():
        read(path,digest,False)
    history = read(f"{prefix}/development_revision_0/attempt.json")
    assert history["status"] == "REJECTED_INPUT_LINK_CHECK" and history["moment_or_cache_sums_started"] is False
    assert history["source_sha256"] == p["failed_attempts"][f"{prefix}/development_revision_0/exact_centered_bulk_L.py"]

    sys.path.insert(0,str(ROOT/"scripts"))
    import exact_centered_bulk_L as C
    fail = AssertionError("actual integration or cache evaluation entered preparation audit")
    with patch.object(C.E.P.Functional,"__init__",side_effect=fail), \
            patch.object(C.M,"triangle_moment",side_effect=fail), \
            patch.object(C.U.C,"exact_bulk_maps",side_effect=fail), \
            patch.object(C.F,"_force_parts",side_effect=fail):
        data = C._source_inputs()
        assert data.hashes == p["input_hashes"]
        assert canon(data.center) == p["point_center_sha256"]
        assert data.context.base["common_gain"] == Q(p["exact_base_gain"])
        assert data.context.dose == Q(p["exact_source_dose"]) > 0
        assert data.gain_bounds == tuple(map(Q,p["closed_gain_shift_exact"]))
        assert data.gain_bounds[0] <= 0 <= data.gain_bounds[1]
        assert data.context.drive == data.context.dose*Q(p["exact_base_gain"])
        assert data.context.row["parameter_cell_sha256"] == p["same_parameter_cell_sha256"]
        assert data.protocol["coordinate_box_sha256"] == p["original_coordinate_box_sha256"]
        counts,active,count_F,count_H = all_selectors(C,data,selectors,read)
        synthetic_stiffness(C)
        synthetic_two_segments(C)
    assert counts == p["chart_counts"]
    own = Path(__file__).relative_to(ROOT).as_posix()
    read(own,parse=False)
    for path,digest in tuple(bindings.items()):
        read(path,digest,False)
    report = {
        **EXPECTED,"status":"PASS",
        "record_id":data.context.record_id,
        "scope":"Independent prelaunch for 144 L centered bulk-force components on the frozen original nodal box and parameter cell. Shared reviewed moment, AD, cache and Hessian packets are premises. No actual force, H, or cache sums.",
        "checks":{
             "frozen_input_hashes_source_snapshot_rejection_history_and_corrected_replay_support_binding":True,
             "exact_same_immutable_L_original_center_box_C_positive_dose_gain0_and_in_cell_segments":True,
             "all_144_element_chart_geometry_original_rho_order_budget_and_signed_force_H_sources":True,
             "all_144_free_and_full_nodal_maps_authoritative_projection_and_independent_axial_stiffness_sample":True,
             "independent_exact_energy_difference_oracle_all_charts_no_factor_two_and_constrained_coordinates_present":True,
             "same_Y_original_cache_frozen_scaled_witnesses_and_inactive_gain_exact_zero":True,
             "mixed_H_and_affine_gain_follow_two_convex_segments_signed_corner_oracle":True,
             "zero_crossing_series_coefficients_retained_and_actual_cache_integration_entrypoints_blocked":True,
        },
        "chart_counts":counts,"active_by_chart":active,"force_endpoint_count":count_F,"H_endpoint_count":count_H,
        "actual_force_integrals_run":False,"actual_cache_site_sums_run":False,
        "input_hashes":dict(sorted(bindings.items())),"review_source_file":own,"review_source_sha256":bindings[own],
    }
    if "--freeze" in sys.argv:
        receipt,gate = OUT/"replay.json",PREP/"prelaunch_independent_review.json"
        assert not receipt.exists() and not gate.exists()
        with receipt.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2)+"\n")
        report["support_file"] = receipt.relative_to(ROOT).as_posix()
        report["support_sha256"] = sha256(receipt.read_bytes()).hexdigest()
        with gate.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"checks":len(report["checks"]),"chart_counts":counts,
                      "active_by_chart":active,"signed_F":count_F,"signed_H":count_H}),flush=True)


if __name__ == "__main__":
    run()
