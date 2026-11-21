"""Independent packet endpoints and direct reduction readback. No moment or cache re-evaluation."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import ceil
from pathlib import Path, PurePosixPath
import sys

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
ROOT = OUT.parents[2]
BASE = OUT.parent.relative_to(ROOT).as_posix()
PREP = OUT.relative_to(ROOT).as_posix()
SOURCE = "scripts/exact_cohort_record_reduction.py"
GATE = PREP+"/prelaunch_independent_review.json"
PINS = {
    SOURCE:"8e3177434299fbffeee2fa5161dd4025fc0e9cd0b7f46be85fceaeccae204ed2",
    GATE:"66c9a028bc45c204d5c889c4c698fef28511c63ab9df6bddf0cd02e8452d425a",
    PREP+"/protocol.json":"772f47398fd7a376fab274bb1d34c9cf7e3fe0bdd661fdd0b1f8b790eccdcc28",
    PREP+"/preflight.json":"a598a545502f09287f909cf8d2753757da1c02863fe2619d06841b8ce1e79496",
}
IDLIST = ["optoGEF2_small_col_"+k for k in "CDEFGHIJKMNOPQR"]
ZERO = (Q(0),Q(0))


def need(test,label):
    if not test:
        raise ValueError(label)


def sha(b):
    return sha256(b).hexdigest()


def stable(x):
    return sha(json.dumps(x,default=str,sort_keys=True,separators=(",",":")).encode())


def safe(path):
    p = (ROOT/path).resolve()
    need(p.is_relative_to(ROOT) and p.relative_to(ROOT).as_posix() == path,"canonical local path")
    return p


def read(path, expected, payload=False):
    need(payload or "/components/" not in path,"no packet opens before admission and freeze")
    b = safe(path).read_bytes()
    need(sha(b) == expected,f"changed input {path}")
    return json.loads(b)


def write(path,obj):
    data = (json.dumps(obj,default=str,indent=2,allow_nan=False)+"\n").encode()
    with path.open("xb") as dst:
        dst.write(data)
    return sha(data)


def q(x):
    need(type(x) in (str,int,Q),"exact rational required")
    return Q(x)


def pair(obj, exact_binary=True):
    lo, hi = q(obj["lower_exact"]),q(obj["upper_exact"])
    need(lo <= hi and (not exact_binary or all(t.denominator.bit_count() == 1 for t in (lo,hi))),
         "strict exact ordered endpoints")
    return lo,hi


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def magnitude(a):
    return max(map(abs,a))


def yes(x):
    return isinstance(x,dict) and bool(x) and all(v is True for v in x.values())


def cast(I,c):
    return tuple(sorted((I[0]*c,I[1]*c)))


def scatter(g,H,ids,force,hessian):
    n = len(ids)
    need(n == len(set(ids)) == len(force) and n > 0 and all(type(i) is int and 0 <= i < len(g) for i in ids),
         "local free map")
    expected = set(combinations_with_replacement(range(n),2))
    used = set()
    for local, global_id in enumerate(ids):
        g[global_id] = add(g[global_id],pair(force[local]))
    for elem in hessian:
        k = (elem["i"],elem["j"])
        need(all(type(i) is int for i in k) and k in expected and k not in used,"unique local upper entries")
        used.add(k)
        i,j = sorted(ids[c] for c in k)
        H[i,j] = add(H.get((i,j),ZERO),pair(elem["enclosure"]))
    need(used == expected,"complete local upper matrix")


def transposed_rows(R,g,H):
    hrows = [Q(0)]*len(g)
    for (i,j),v in H.items():
        for k in {i,j}:
            hrows[k] += magnitude(v)
    RF, b, rnorm = [], [], []
    for r in R:
        rr = [Q.from_float(float(v)) for v in r]
        values = [cast(v,c) for v,c in zip(g,rr)]
        RF.append((sum(x[0] for x in values),sum(x[1] for x in values)))
        b.append(sum(abs(x)*v for x,v in zip(rr,hrows)))
        rnorm.append(sum(map(abs,rr)))
    return RF,hrows,b,max(rnorm)


def synth():
    g,H = [ZERO,ZERO],{}
    interval = lambda a,b:{"lower_exact":str(a),"upper_exact":str(b)}
    h = [{"i":i,"j":j,"enclosure":interval(-2,-1)} for i,j in ((0,0),(0,1),(1,1))]
    scatter(g,H,[1,0],[interval(-3,-2),interval(1,2)],h)
    scatter(g,H,[1,0],[interval(2,3),interval(-1,0)],h)
    RF,rows,b,norm = transposed_rows([[1.0,-0.5],[0.0,0.5]],g,H)
    need(g == [(Q(0),Q(2)),(Q(-1),Q(1))] and RF == [(Q(-1,2),Q(5,2)),(Q(-1,2),Q(1,2))]
         and rows == [8,8] and b == [12,4] and norm == Q(3,2),"independent signed arithmetic")
    try:
        scatter([ZERO,ZERO],{},[1,0],[interval(0,0)]*2,h[:-1])
    except ValueError:
        return
    raise AssertionError("missing Hessian accepted")


def rho_check(outward,raw):
    need(q(outward) == Q(ceil(raw*2**48),2**48) and 0 <= raw <= q(outward) < 1,"outward exact rho")


def tail(actual,budget,exact=None):
    lo,hi = pair(actual)
    if exact is None:
        need(0 <= lo <= hi <= budget,"original nonnegative bounded tail")
        return hi
    value = q(exact)
    need(0 <= value <= budget and lo <= value <= hi,"exact scaled tail and rounded enclosure")
    return value


def packet_meta(ctx,p,body,case,proof,v,part,index):
    a = ctx.target.a
    node_ids = a["triangles"][index] if part == "bulk" else a[f"surface_{index//12}_edges"][index%12]
    dofs = [2*int(n)+k for n in node_ids for k in (0,1)]
    projected = [i for i,k in enumerate(dofs) if k in ctx.target.index]
    ids = [ctx.target.index[dofs[i]] for i in projected]
    expected = {"record_id":ctx.record_id,"dose_exact":str(ctx.dose),"precision_bits":640,"global_indices":ids,
        **{k:ctx.row[k] for k in ("parameter_cell_sha256","predictor_sha256","point_cache_sha256","model_observation_contract_sha256")}}
    need(all(v.get(k) == val for k,val in expected.items()) and case["record_id"] == ctx.record_id
         and case["global_indices"] == ids,"same record, cell, dose and projected fields")
    fmax,hmax = q(p["force_tail_budget_exact"]),q(p["H_entry_tail_budget_exact"])
    force = v["force_integral_minus_cache"]
    h = v["H_integral_minus_cache_upper"]
    need(len(v["force_two_mass_tail_upper"]) == len(force) == len(ids),"all local forces/tails")
    if part == "bulk":
        geometry = body["triangles"][index]
        lo,hi = pair(geometry["J_bound"],False)
        raw = max(abs(lo-1),abs(hi-1))
        axis = sum(Q.from_float(float(a["reference"][n,0])) == 0 for n in node_ids[:3])
        chart = ("off_axis","axis_vertex","axis_edge_physical")[axis]
        route = ("scaled_positive_R_affine_material",chart,chart)[axis]
        actual = v["cache_consistency"]
        need(geometry["positive_J_for_entire_box"] is True and lo > Q(1,10**6)
             and list(map(q, case["J_bounds_exact"])) == [lo,hi] and raw == q(case["raw_rho_exact"])
             and v["status"] == "UNIFORM_COHORT_BULK_COMPONENT_UNREVIEWED"
             and v["chart"] == case["chart"] == chart and v["integration_chart"] == proof["integration_chart"] == route
             and v["record_component_proof_sha256"] == stable(proof) and proof["original_selector_sha256"] == stable(case)
             and proof["body_element_sha256"] == stable(geometry)
             and v["element"] == case["element"] == index and v["free_local"] == case["free_local"] == projected
             and v["active"] == case["active"] == str(int(a["active"][index,0]))
             and all(v[k] == case[k] == proof[k] for k in ("nodal_hulls_sha256","auxiliary_element_rule_sha256"))
             and v["degree"] == case["degree"] == proof["unchanged_degree"]
             and v["rho_exact"] == case["rho_exact"] == proof["unchanged_rho_exact"]
             and actual["all64_convex_strain_lines_positive"] is True
             and pair(actual["minimum_cache_convex_line_J_lower"])[0] > Q(1,10**6)
             and pair(actual["max_coefficientwise_delta_q_upper"])[0] >= 0
             and fmax == q(case["force_tail_budget_exact"]) and hmax == q(case["hessian_tail_budget_exact"]),
             f"same element proof and cache-line analyticity {index}")
        rho_check(v["rho_exact"],raw)
        if not axis:
            e = v["force_two_mass_tail_upper_exact"]
            need(len(e) == len(ids),"exact force tails")
            fs = [tail(x,fmax,y) for x,y in zip(v["force_two_mass_tail_upper"],e)]
            hs = [tail(z["two_mass_tail_upper"],hmax,z["two_mass_tail_upper_exact"]) for z in h]
            need(max(fs) == q(proof["scaled_force_tail_max_exact"]) and max(hs) == q(proof["scaled_H_entry_tail_max_exact"]),
                 "scaled norm maxima match the frozen same-state proofs")
    else:
        face,edge = divmod(index,12)
        dom = case["domain"]
        lo,hi = map(q,dom["G_bounds_exact"])
        raw = (hi-lo)/(hi+lo)
        need(0 < lo <= hi and q(dom["hoop_lower_exact"]) > 0 and raw == q(dom["raw_rho_exact"])
             and v["status"] == "UNIFORM_COHORT_SURFACE_COMPONENT_UNREVIEWED"
             and v["face"] == case["face"] == face and v["edge"] == case["edge"] == edge
             and v["order"] == case["order"] == proof["unchanged_order"] and proof["case_sha256"] == stable(case)
             and q(v["q0_exact"]) == q(dom["q0_exact"]) == (lo+hi)/2
             and v["rho_exact"] == dom["rho_exact"] == proof["unchanged_rho_exact"],
             f"same certified surface expansion {face} {edge}")
        rho_check(v["rho_exact"],raw)
        site_weights = [Q.from_float(float(w))/ctx.target.energy_scale for w in a[f"surface_{face}_weights"][edge]]
        need(len(site_weights) == len(v["actual_cache_sites"]) == len(case["cache_site_domains"]) == 12
             and len(v["force_auxiliary_minus_cache"]) == len(ids)
             and pair(v["material_mass"])[0] > 0
             and pair(v["auxiliary_mass"])[0] <= sum(site_weights) <= pair(v["auxiliary_mass"])[1],
             "all surface cache entries and same positive two masses")
        for i,(x,y,w) in enumerate(zip(v["actual_cache_sites"],case["cache_site_domains"],site_weights)):
            need(x["quadrature_index"] == y["quadrature_index"] == i and y["global_indices"] == ids
                 and y["selected_local"] == projected and y["nodes"] == list(map(int,node_ids))
                 and q(x["weight_applied_exactly_once"]) == q(y["auxiliary_normalized_weight_exact"]) == w > 0
                 and q(y["convex_line_exact"][0][0]) > 0 and q(y["tangent_norm_squared_lower_exact"]) > 0
                 and all(pair(x[k])[0] >= 0 for k in ("unweighted_F_max","unweighted_H_max")),
                 "same reviewed surface cache maps, weights and strain-line domains")
        for x in v["force_auxiliary_minus_cache"]:
            pair(x)
        for x in h:
            pair(x["auxiliary_minus_cache"])
    if part != "bulk" or axis:
        for x in v["force_two_mass_tail_upper"]:
            tail(x,fmax)
        for x in h:
            tail(x["two_mass_tail_upper"],hmax)
    return ids


def freeze():
    synth()
    for path,digest in PINS.items():
        need(sha(safe(path).read_bytes()) == digest,"frozen replay inputs")
    src = Path(__file__).read_bytes()
    with (HERE/"replay_frozen.py").open("xb") as f:
        f.write(src)
    p = write(HERE/"protocol.json",{"status":"FROZEN_INDEPENDENT_DIRECT_PACKET_AND_ENDPOINT_REPLAY",
        "source_sha256":sha(src),"input_pins":PINS,"records":IDLIST,
        "required_external_pins":["observation file and digest containing exact per-part wait/terminal/manifest",
                                  "producer admission/terminal/results digests observed by retained process controller"],
        "scope":"Read all admitted packets, verify identity and tail/domain/measure metadata, then independent signed Fraction scatter and R endpoint arithmetic. Reviewed component kernels and cached-Q proof/geometry are shared premises. No component and cache arithmetic replay."})
    f = write(HERE/"preflight.json",{"status":"PASS_PREPARATION_ONLY","source_sha256":sha(src),"protocol_sha256":p,
        "checks":{"independent_signed_reordered_scatter_and_R":True,"malformed_missing_H_rejected":True,
                  "frozen_input_pins_current":True,"no_actual_component_reads":True,"no_moments_or_cache_sums":True}})
    print(json.dumps({"source_sha256":sha(src),"protocol_sha256":p,"preflight_sha256":f}))


def observed(record, obs_path, obs_sha):
    obs = read(obs_path,obs_sha)
    item, = [v for v in obs["records"] if v["record_id"] == record]
    need(set(item["parts"]) == {"bulk","surface"},"per-part observation")
    return item["parts"]


def part_files(record, links, pins, part):
    folder = f"{BASE}/records/{record}/scaled_integrals_v1/actual/{part}"
    manifest = read(folder+"/results.json",pins["manifest_sha256"])
    terminal = read(folder+"/terminal.json",pins["terminal_sha256"])
    supervisor = PurePosixPath(pins["supervisor_dir"])
    need(supervisor.parent.as_posix() == f"{BASE}/remaining15_integrals_preparation/launches"
         and len(supervisor.name)==64,"canonical supervisor")
    piece = f"{supervisor}/{record}_{part}"
    os_term = read(piece+".os_terminal.json",pins["os_terminal_sha256"])
    pid = pins["pid"]
    need(type(pid) is int and pid > 0 and os_term["pid"] == terminal["pid"] == pid
         and type(os_term["exit_code"]) is type(terminal["exit_code"]) is int
         and os_term["exit_code"] == terminal["exit_code"] == 0
         and terminal["status"] == "FINISHED" and terminal["results_sha256"] == pins["manifest_sha256"]
         and manifest["status"] == "COMPONENT_PART_COMPLETE_UNREVIEWED" and manifest["part"] == part
         and all(manifest.get(k)==v for k,v in links.items()),"explicit successful local worker and independently observed true OS exit")
    files = manifest["files"]
    if part == "bulk":
        expected = [f"{folder}/components/element_{i:03d}.json" for i in range(144)]
        need([x["element"] for x in files] == list(range(144)),"all bulk indices")
    else:
        expected = [f"{folder}/components/face{f}_edge{e:02d}.json" for f in range(2) for e in range(12)]
        need([(x["face"],x["edge"]) for x in files] == [divmod(i,12) for i in range(24)],"all surface indices")
    need([x["file"] for x in files] == expected and len(files) == pins["components"], "canonical unique packet names")
    return files


def replay(args):
    here = HERE.relative_to(ROOT).as_posix()
    doc = read(here+"/protocol.json",args.protocol_sha256)
    pre = read(here+"/preflight.json",args.preflight_sha256)
    src = Path(__file__).read_bytes()
    need((HERE/"replay_frozen.py").read_bytes() == src and doc["source_sha256"] == pre["source_sha256"] == sha(src)
         and pre["status"] == "PASS_PREPARATION_ONLY" and yes(pre["checks"])
         and pre["protocol_sha256"] == args.protocol_sha256 and doc["input_pins"] == PINS,
         "frozen review and sealed prep")
    for path,digest in PINS.items():
        need(sha(safe(path).read_bytes()) == digest,"changed frozen prerequisite")
    synth()
    record = args.record
    need(record in IDLIST and not (HERE/f"{record}.json").exists(),"unique approved record replay")
    obs = observed(record,args.observation,args.observation_sha256)
    upstream = PREP+"/"+record
    producer = read(upstream+"/results.json",args.results_sha256)
    terminal = read(upstream+"/terminal.json",args.terminal_sha256)
    admission = read(PREP+f"/{record}_input_admission.json",args.admission_sha256)
    need(terminal["status"] == "FINISHED" and type(terminal["exit_code"]) is int and terminal["exit_code"] == 0
         and terminal["results_sha256"] == args.results_sha256 and terminal["pid"] > 0
         and producer["admission_sha256"] == args.admission_sha256 and producer["record_id"] == record
         and admission["record_id"] == record and admission["local_component_payloads_opened"] is False
         and admission["status"] == "ALL_168_EXPECTED_PACKETS_BOUND_BEFORE_LOCAL_READ"
         and admission["review_sha256"] == PINS[GATE], "same completed producer input admission")
    sys.path.insert(0,str(ROOT/"scripts"))
    import exact_cohort_record_reduction as A
    runtime, plan, jobs, expect, hashes = A._boot()
    ctx,p,domains,body,routes,links,stage = A._load_record(runtime,plan,jobs,expect,record,hashes)
    loaded = {}
    for part,pins in obs.items():
        need(all(pins[k] == v for k,v in admission["external_part_pins"][part].items()),"producer input terminal set")
        loaded[part] = part_files(record,links,pins,part)
        need(loaded[part] == admission["expected_packets"][part],"same canonical admitted packets")
    # First read of a component occurs only after both independently observed parts have been pinned.
    f,H,counts,hcounts = [ZERO]*611,{},[0]*611,[0]*611
    tail_receipt = {"bulk":0,"surface":0,"axis_vertex":0,"axis_edge_physical":0,"scaled":0}
    for part in ("bulk","surface"):
        cases = domains["bulk_selectors" if part=="bulk" else "surface_cases"]
        for i, item in enumerate(loaded[part]):
            v = read(item["file"],item["sha256"],payload=True)
            ids = packet_meta(ctx,p,body,cases[i],routes[part][i],v,part,i)
            scatter(f,H,ids,v["force_integral_minus_cache"],v["H_integral_minus_cache_upper"])
            for j in ids:
                counts[j] += 1
                hcounts[j] += len(ids)
            tail_receipt[part] += 1
            if part == "bulk":
                tail_receipt[v["integration_chart"] if v["chart"]!="off_axis" else "scaled"] += 1
    need(counts == p["force_incidence"] and hcounts == p["H_row_incidence"]
         and len(H) == 6676 and tail_receipt == {"bulk":144,"surface":24,"axis_vertex":6,"axis_edge_physical":6,"scaled":132},
         "complete force/H and all routes")
    RF,hrows,br,norm = transposed_rows(ctx.R,f,H)
    d,b = max(map(magnitude,RF)),max(br)
    eta0,z0,r = [q(p[k]) for k in ("cached_eta_upper_exact","cached_Z_upper_exact","normal_radius_exact")]
    data = read(GATE,PINS[GATE])["old_proof_and_spd_bindings"][record]
    need(z0 == q(data["cached_Z_upper_exact"]) and r == q(data["exact_radius"]) == ctx.old_radius
         and p["cached_bound_provider"] == data["cached_provider"]
         and norm == q(p["R_infinity_norm_exact"]) and 0 <= z0 < 1 and eta0+z0*r < r,"same independent old-Q premise and R")
    eta,z = eta0+d,z0+b
    lhs = eta+z*r
    passes = z < 1 and lhs < r
    error = eta/(1-z) if passes else None
    expected = {"RF_defect_uniform_upper_exact":d,"b_upper_exact":b,"old_cached_eta_upper_exact":eta0,
                "old_cached_Z_upper_exact":z0,"old_radius_exact":r,"total_eta_upper_exact":eta,"total_Z_upper_exact":z,
                "strict_selfmap_lhs_exact":lhs,"R_infinity_norm_exact":norm}
    need(all(q(producer[k])==v for k,v in expected.items())
         and producer["sufficiency_passed"] is passes and producer["root_assertion_independently_reviewed"] is False
         and producer["old_cached_bound_provider"] == data["cached_provider"]
         and producer["same_closed_C_sha256"] == ctx.row["parameter_cell_sha256"]
         and (producer["conditional_root_error_exact"] is None if error is None else q(producer["conditional_root_error_exact"]) == error)
         and producer["status"] == ("CONDITIONAL_UNREVIEWED_DIRECT_SUFFICIENT_TEST_PASSES" if passes
                                   else "UNREVIEWED_DIRECT_SUFFICIENT_TEST_UNRESOLVED"),
         "exact scalar direct test readback")
    need([tuple(map(q,row)) for row in producer["global_F_intervals_exact"]] == f
         and [tuple(map(q,row)) for row in producer["preconditioned_RF_intervals_exact"]] == RF
         and producer["global_H_upper_intervals_exact"] == [[i,j,str(lo),str(hi)] for (i,j),(lo,hi) in sorted(H.items())]
         and list(map(q,producer["H_absolute_row_bounds_exact"])) == hrows
         and list(map(q,producer["abs_R_H_row_bounds_exact"])) == br,"all exact 611F, 6676H, signed RF and row endpoints")
    curvature = None
    if passes:
        # Cached endpoint stays SPD through the symmetric homotopy since ||I-RH_t|| < 1.
        cached = read(data["cached_spd_audit_file"],data["cached_spd_audit_sha256"])
        need(cached["record_id"] == record and cached["status"] == "PASS" and yes(cached["checks"])
             and q(cached["uniform_Hessian_lower_bound_exact"]) > 0
             and q(cached["normal_root_radius_exact"]) == r and q(cached["independent_fulltube_Z_upper_exact"]) == z0
             and q(cached["exact_R_infinity_norm"]) == norm,"same-record cached positive endpoint")
        curvature = (1-z)/norm
        pc = producer["conditional_curvature"]
        need(pc["strict_local_minimum"] is True and q(pc["normalized_curvature_lower_exact"]) == curvature
             and pc["cached_independent_review_sha256"] == data["cached_spd_audit_sha256"],"homotopy curvature exact readback")
    else:
        need(producer["conditional_curvature"] is None,"unresolved direct test")
    receipt = {"status":"PASS_DIRECT_ROOT_AND_MINIMUM" if passes else "PASS_DIRECT_REDUCTION_UNRESOLVED",
               "record_id":record,"root_certified":passes,"strict_local_minimum_certified":passes,
               "source_sha256":sha(src),"protocol_sha256":args.protocol_sha256,"preflight_sha256":args.preflight_sha256,
               "producer_pins":{"observation_file":args.observation,"observation_sha256":args.observation_sha256,
                   "admission_sha256":args.admission_sha256,"terminal_sha256":args.terminal_sha256,"results_sha256":args.results_sha256},
               "packet_route_counts":tail_receipt,"F_entries":611,"H_upper_entries":6676,"exact_scalars":expected,
               "root_error_upper_exact":error,"normalized_curvature_lower_exact":curvature,
               "same_closed_C_sha256":ctx.row["parameter_cell_sha256"],
               "RF_entries_containing_zero":sum(a<=0<=b for a,b in RF),
               "scope":"All admitted component endpoints, reviewed chart, domain, tail, and surface cache-weight metadata, and independently recomputed exact signed reducer. Shared independently reviewed component-enclosure engines and cached-Q proofs. No component, moment, or cache-sum rerun."}
    h = write(HERE/f"{record}.json",receipt)
    print(json.dumps({"status":receipt["status"],"record_id":record,"review_sha256":h,
                      "rf":float(d),"b":float(b),"lhs":float(lhs),"r":float(r),"z":float(z)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    cli = parser.add_subparsers(dest="mode",required=True)
    cli.add_parser("prepare")
    p = cli.add_parser("replay")
    p.add_argument("record")
    for flag in ("protocol-sha256","preflight-sha256","observation","observation-sha256",
                 "admission-sha256","terminal-sha256","results-sha256"):
        p.add_argument("--"+flag,required=True)
    args = parser.parse_args()
    freeze() if args.mode == "prepare" else replay(args)

