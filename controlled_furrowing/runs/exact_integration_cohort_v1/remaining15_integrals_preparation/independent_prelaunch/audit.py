"""Review the frozen remaining-cohort launch and per-record route packets."""

import ast
from fractions import Fraction as Q
from hashlib import sha256
import importlib
import json
from pathlib import Path
import sys
from unittest.mock import Mock, patch

OUT = Path(__file__).resolve().parent
BATCH = OUT.parent
ROOT = OUT.parents[3]
BASE = BATCH.parent.relative_to(ROOT).as_posix()
SOURCE = "scripts/remaining_cohort_integrals.py"
P0 = BATCH.relative_to(ROOT).as_posix()
IDS = ["optoGEF2_small_col_" + c for c in "CDEFGHIJKMNOPQR"]
AXIS = [*range(6), *range(72,78)]
PARTS = ("bulk","surface")
PIN = {
    "source_sha256":"c9654ab19be69aa70a0100ede0919bc8a87040728ca38487758f55ad3047dcfc",
    "protocol_sha256":"ce18ae97d8e9a31db1bb12460371700eb0b17a5cb90e05daf0248584e7d8585e",
    "preflight_sha256":"a6b5939d7e7aefecc150f45067764e30d87e9b8551355a2b5e223e4d0e992d50",
    "manifest_sha256":"e4cbc8826b763d235e24ff827d684a00526b79e18bef63b51bc74f727dcdd363",
}
ADAPTER = {"file":f"{BASE}/remaining15_scaled_preparation/independent_review.json",
           "sha256":"821b3730305b1d55da9a1cb119780a1717adf99b953539e2357535d81e597fe2"}
RRESULT = f"{BASE}/remaining15_domains_review/results.json"
RRESULT_SHA = "161e67f2bf1d2d01cafe303d82971edc6a3cc48387d8d664dc6a84d23c5023af"
RTERM_SHA = "7e52139b4425346a69b7ca338969842d6c1aa8b5048017ae846f00c51e4a5170"


def insist(x, message):
    if not x:
        raise ValueError(message)


def sha(raw):
    return sha256(raw).hexdigest()


def stable(x):
    return sha(json.dumps(x, default=str, sort_keys=True, separators=(",",":")).encode())


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


def raw(path):
    f = (ROOT/path).resolve()
    insist(str(path) == f.relative_to(ROOT).as_posix(), f"unexpected path {path}")
    return f.read_bytes()


def read(path, pin=None, parse=True):
    b = raw(path)
    if pin is not None:
        insist(sha(b) == pin, f"changed {path}")
    return json.loads(b) if parse else b


def yes(values):
    return isinstance(values,dict) and bool(values) and all(v is True for v in values.values())


def frac(x):
    return Q.from_float(float(x))


def exact_hist(last, attempts, budgetf, budgeth, cap):
    insist(2 <= last <= cap and [x["degree"] for x in attempts] == list(range(2,last+1)),
           "consecutive attempts and reviewed order")
    for item in attempts:
        f, h = Q(item["force_tail_upper_exact"]), Q(item["H_entry_tail_upper_exact"])
        insist(f >= 0 and h >= 0 and ((f <= budgetf and h <= budgeth) == (item["degree"] == last)),
               "selected minimal order respects both positive budgets")


def launch_shape(L, runtime, protocol, adapter):
    tree = ast.parse(raw(SOURCE))
    methods = {n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
    boot_source = ast.unparse(methods["_boot"])
    insist(boot_source.index("checked(n, h, parse=False)") < boot_source.index("importlib.import_module("),
           "verify every frozen project source before import")
    worker = ast.unparse(methods["_worker"])
    launch = ast.unparse(methods["launch"])
    insist(launch.index("_gate(") < launch.index("subprocess.Popen(") < launch.index("process.wait()")
           < launch.index("os_terminal.json") < launch.index("require(code == 0")
           and launch.count("subprocess.Popen(") == 1 and "shell=True" not in launch
           and "terminal['pid'] == process.pid" in launch and "_terminal(terminal, result_hash)" in launch,
           "explicit pre-gate and serial true OS supervisor")
    insist(worker.index("_gate(") < worker.index("_boot()") < worker.index(".mkdir(parents=True, exist_ok=False)")
           < min(worker.index("runtime.I.bulk_component("),worker.index("runtime.E.surface_component("))
           and "data.hashes == protocol['domain_and_record_input_hashes']" in worker
           and "runtime.I.bulk_component(data.ctx, data.nodes, data.rule" in worker
           and "runtime.E.surface_component(data.ctx, data.nodes, data.rule, case)" in worker
           and "_surface_proof(runtime, data, face, edge, case) == proof" in worker
           and "except BaseException as error" in worker,
           "rechecked worker input bindings and exclusive output before one exact-record callback")
    insist(runtime.pinned == protocol["transitive_sources"]
           and adapter["supplemental_runtime_source_hashes"] == {L.AUX_FORCE:runtime.pinned[L.AUX_FORCE]}
           and runtime.S.F is runtime.F and runtime.S.T is runtime.I.T is runtime.E.T
           and runtime.S.A is runtime.F.A is runtime.I.A
           and runtime.S.core is runtime.I.core is runtime.E.core,
           "bound supplemental source and resolved runtime aliases")


def body_and_basic(data, leaf, row, terminal, record):
    ctx, d, p = data.ctx, data.domains, data.protocol
    box_path = next(k for k in data.hashes if k.endswith("/coordinate_box.json"))
    box = read(box_path,data.hashes[box_path])
    body_path = next(k for k in data.hashes if k.endswith("/body_geometry.json"))
    body = read(body_path,data.hashes[body_path])
    actual = {k:v for k,v in data.hashes.items() if k.startswith(rel(data.original)+"/") or k.endswith(".os_terminal.json")}
    insist(actual == leaf["input_sha256"] and len(actual) == 8 and body == data.body
           and ctx.row == row["selected_context"] == p["row"] and terminal["exit_code"] == 0
           and d["body_sha256"] == data.hashes[body_path] and d["coordinate_box_sha256"] == data.hashes[box_path]
           and leaf["same_C_sha256"] == ctx.row["parameter_cell_sha256"] == box["parameter_cell_sha256"]
           and leaf["proof_center_state"] == ctx.row["proof_center_state_file"]
           and Q(leaf["dose_exact"]) == ctx.dose and ctx.drive == ctx.base["common_gain"]*ctx.dose
           and Q(box["corrector_radius_exact"]) == Q(p["normal_radius_exact"]) == ctx.old_radius
           and p["root_asserted"] is False and box["root_membership_assumed"] is False,
           f"selected context and exact reviewed original-domain leaves {record}")
    assert ctx.record_id == record
    return box


def route_packets(runtime, data, proofs, box, static, report):
    ctx, E, S = data.ctx, runtime.E, runtime.S
    a, d, p = ctx.target.a, data.domains, data.protocol
    triangles = a["triangles"].tolist()
    names = ["off_axis","axis_vertex","axis_edge_physical"]
    route = ["scaled_positive_R_affine_material", *names[1:]]
    free = {int(k):j for j,k in enumerate(ctx.free)}
    dose = str(ctx.dose)
    budgetf, budgeth = Q(p["force_tail_budget_exact"]),Q(p["H_entry_tail_budget_exact"])
    insist(budgetf > 0 and budgeth > 0 and p["precision_bits"] == 640
           and set(proofs) == {"bulk","surface"} and len(proofs["bulk"]) == len(triangles) == 144
           and len(proofs["surface"]) == 24, "complete work units and exact positive budgets")
    # These are actual outward Arb hulls passed by the worker. Check enclosure of the exact selected domain.
    hulls = []
    for j, pair in enumerate(data.nodes):
        hull = [E.SUR.rational_hull(t) for t in pair]
        for axis,(l,h) in enumerate(hull):
            y,r = Q(box["current_centers_exact"][j][axis]), Q(box["coordinate_radii_exact"][j][axis])
            insist(l <= y-r <= y+r <= h and (axis != 0 or not (frac(a["reference"][j,0]) == 0) or l == h == 0),
                   "outward same-Y hull and exact radial axis condition")
        hulls.append([list(map(str,x)) for x in hull])
    assert sum(map(len,hulls)) == 650
    ids0 = {"record_id":ctx.record_id,"dose_exact":dose,"precision_bits":640,
            **{k:ctx.row[k] for k in ("parameter_cell_sha256","model_observation_contract_sha256",
                                     "point_cache_sha256","predictor_sha256")}}
    local_sites = []
    for row in a["shape"]:
        n = list(map(frac,row))
        xi,eta = n[1]+(n[3]+n[4])/2, n[2]+(n[4]+n[5])/2
        local_sites.append([1-xi-eta,xi,eta])
    insist(local_sites == data.rule["volume_barycentric"] and len(local_sites) == 64, "unchanged local auxiliary sites")
    counts = dict.fromkeys(route,0)
    total_symbolic = 0
    for i, (nodes, proof, selector, geometry) in enumerate(zip(triangles,proofs["bulk"],d["bulk_selectors"],data.body["triangles"])):
        local_full = [2*node+j for node in nodes for j in range(2)]
        local_free = [j for j,k in enumerate(local_full) if k in free]
        global_free = [free[local_full[j]] for j in local_free]
        corners = tuple(tuple(map(frac,a["reference"][node])) for node in nodes[:3])
        axis = sum(v[0] == 0 for v in corners)
        expect_route = route[axis]
        counts[expect_route] += 1
        insist((axis != 0) == (i in AXIS), "axis chart and same record projection")
        weights = [frac(w)/ctx.target.energy_scale for w in a["volume_weights"][i]]
        rule_hash = stable([[list(map(str,x)) for x in local_sites],list(map(str,weights))])
        local_hash = stable([hulls[j] for j in nodes])
        expect = {**ids0,"element":i,"chart":names[axis],"free_local":local_free,"global_indices":global_free,
                  "active":str(frac(a["active"][i,0])),"auxiliary_element_rule_sha256":rule_hash,
                  "nodal_hulls_sha256":local_hash}
        insist(all(proof.get(k) == v == selector.get(k) for k,v in expect.items())
              and list(data.rule["volume_weights"][i]) == weights and min(weights) > 0
              and proof["local_nodes"] == nodes and proof["integration_chart"] == expect_route
              and proof["body_element_sha256"] == stable(geometry) and proof["original_selector_sha256"] == stable(selector),
              f"per-element original geometry/rule/projection {ctx.record_id} {i}")
        lo,hi = [Q(geometry["J_bound"][k]) for k in ("lower_exact","upper_exact")]
        rawrho = max(abs(lo-1),abs(hi-1))
        rho = Q(selector["rho_exact"])
        insist(type(selector["degree"]) is int and selector["degree"] == proof["unchanged_degree"]
               and Q(proof["unchanged_rho_exact"]) == rho
               and list(map(Q,selector["J_bounds_exact"])) == [lo,hi]
               and Q(1,10**6) < lo <= hi and geometry["positive_J_for_entire_box"] is True
               and Q(selector["raw_rho_exact"]) == rawrho and rawrho <= rho < 1
               and (rho*2**48).denominator == 1 and rho-Q(1,2**48) < rawrho
               and Q(proof["original_force_tail_budget_exact"]) == Q(selector["force_tail_budget_exact"]) == budgetf
               and Q(proof["original_H_entry_tail_budget_exact"]) == Q(selector["hessian_tail_budget_exact"]) == budgeth,
               f"same positive J domain, outward dyadic rho, degree and budgets {ctx.record_id} {i}")
        exact_hist(selector["degree"],selector["order_attempts"],budgetf,budgeth,p["max_bulk_order"])
        n = len(local_free)
        symbolic = 5*(1+n)+1+n+n*(n+1)//2 if not axis else 0
        insist(proof["scaled_polynomial_identity_count"] == symbolic, "full symbolic q/B/J/DJ/D2J identity count")
        total_symbolic += symbolic
        if axis:
            insist("scaled_chart_geometry" not in proof,"axis dispatch retains the reviewed original kernel")
            continue
        insist(0 <= Q(proof["scaled_force_tail_max_exact"]) <= Q(selector["force_tail_upper_exact"]) <= budgetf
               and 0 <= Q(proof["scaled_H_entry_tail_max_exact"]) <= Q(selector["H_entry_tail_upper_exact"]) <= budgeth,
               "bounded minima of chart and original norm majorants")
        if corners not in static:
            chart = S.chart(corners)
            static[corners] = (stable(S._metadata(chart)),stable([S.chart_point(chart,site) for site in local_sites]))
        svalue, loc = static[corners]
        insist(stable(proof["scaled_chart_geometry"]) == svalue and proof["scaled_exact_Qstar_node_count"] == 64
                and proof["scaled_Qstar_node_sha256"] == loc, "same physical 64-site signed functional coordinates")
    surface_sites = 0
    for i,(case,proof) in enumerate(zip(d["surface_cases"],proofs["surface"])):
        face,edge = divmod(i,12)
        face_nodes = a[f"surface_{face}_edges"][edge]
        mapped = [free[k] for node in face_nodes for k in (2*int(node),2*int(node)+1) if k in free]
        sites = [frac(n[2])+frac(n[1])/2 for n in a[f"surface_{face}_N"]]
        weights = [frac(w)/ctx.target.energy_scale for w in a[f"surface_{face}_weights"][edge]]
        causal = case["domain"]
        rlo,rhi = map(Q,causal["G_bounds_exact"])
        rr = (rhi-rlo)/(rhi+rlo)
        rho = Q(causal["rho_exact"])
        insist(list(data.rule["surface_s"][face]) == sites and all(0<s<1 for s in sites)
               and list(data.rule["surface_weights"][face][edge]) == weights and min(weights) > 0
               and [v["quadrature_index"] for v in case["cache_site_domains"]] == list(range(12))
               and case["record_id"] == ctx.record_id and (case["face"],case["edge"]) == (face,edge)
               and case["global_indices"] == mapped
               and proof == {"face":face,"edge":edge,"case_sha256":stable(case),"global_indices":mapped,
                             "unchanged_order":case["order"],"unchanged_rho_exact":causal["rho_exact"]},
               f"same reviewed surface cases and cache maps {ctx.record_id} {face} {edge}")
        insist(0 < rlo <= rhi and Q(causal["q0_exact"]) == (rlo+rhi)/2
               and Q(causal["raw_rho_exact"]) == rr and 0 <= rr <= rho < 1
               and (rho*2**48).denominator == 1 and rho-Q(1,2**48) < rr,
               "preserved positive surface G, q0 and outward rho")
        exact_hist(case["order"],case["order_attempts"],budgetf,budgeth,p["max_surface_order"])
        surface_sites += len(sites)
    report[ctx.record_id] = {"routes":counts,"symbolic_identities":total_symbolic,"surface_edges":24,"surface_sites":surface_sites}


def mocked_supervisor(L):
    """Exercise launch supervision with fake waitable processes and output fixtures only."""
    rec = IDS[0]
    part = "bulk"
    base = OUT/"supervisor_mocks"
    base.mkdir(exist_ok=False)
    results = {}
    for tag, code, wrong_pid in (("success",0,False),("child_error",11,False),("pid_mismatch",0,True)):
        work = base/tag
        work.mkdir()
        record_dir = work/"records"/rec
        actual = record_dir/"scaled_integrals_v1"/"actual"/part
        manifest, review = "a"*64, "b"*64
        job = [rec+":"+part]
        callback = Mock()
        p = Mock()
        p.pid = 321
        p.wait.return_value = code

        def fake_start(args, **kwargs):
            insist(checker.called and kwargs["cwd"] == ROOT and args[-1] == job[0] and args[4] == "worker",
                   "gate runs before execution of one explicitly selected job")
            callback(args)
            if code == 0:
                actual.mkdir(parents=True)
                (actual/"results.json").write_bytes(b"{}\n")
                (actual/"terminal.json").write_text(json.dumps({"status":"FINISHED","exit_code":0,
                    "pid":p.pid+int(wrong_pid),"results_sha256":sha((actual/"results.json").read_bytes())}),encoding="utf8")
            return p

        with (patch.object(L,"OUT",work),
              patch.object(L,"_record_dir",return_value=record_dir),
              patch.object(L,"_gate",return_value=({}, {(rec,part):None})) as checker,
              patch.object(L.subprocess,"Popen",side_effect=fake_start)):
            failed = False
            try:
                L.launch(manifest,review,job)
            except ValueError:
                failed = True
            insist(callback.call_count == p.wait.call_count == 1, "one wait on exactly the launched worker")
            launch_dir, = list((work/"launches").iterdir())
            found = json.loads((launch_dir/f"{rec}_{part}.os_terminal.json").read_bytes())
            terminal = json.loads((launch_dir/"terminal.json").read_bytes())
            success = code == 0 and not wrong_pid
            insist(found["pid"] == p.pid and found["exit_code"] == code and success == (not failed)
                   and terminal["status"] == ("FINISHED" if success else "FAILED")
                   and (launch_dir/"failure.json").exists() is not success,
                   "OS and application terminal agreement or retained failure")
            if success:
                try:
                    L._unused(rec, part)
                except ValueError:
                    results["already_started_part_rejected"] = True
                else:
                    raise AssertionError("already started case was not rejected")
            results["mocked_os_"+tag] = True
    return results


def run():
    read(SOURCE,PIN["source_sha256"],False)
    read(f"{P0}/remaining_cohort_integrals_frozen.py",PIN["source_sha256"],False)
    protocol, pre, manifest = [read(f"{P0}/{n}.json",PIN[f"{n}_sha256"]) for n in ("protocol","preflight","manifest")]
    prep_started = read(f"{P0}/started.json","b654094604dda75ba3c6bd50c4a86f1c75cde782f7ad8b3219e8d5b80632cfe0")
    prep_term = read(f"{P0}/terminal.json","77311167bba761e6c54e9ba536506c4fc0001a4bd964865b4ee0b39dfc7353cc")
    rresult = read(RRESULT,RRESULT_SHA)
    rterm = read(f"{BASE}/remaining15_domains_review/terminal.json",RTERM_SHA)
    adapter = read(ADAPTER["file"],ADAPTER["sha256"])
    insist(protocol["records"] == [r["record_id"] for r in manifest["records"]]
           == [r["record_id"] for r in rresult["records"]] == IDS
           and protocol["precision_bits"] == 640
           and protocol["status"] == "FROZEN_REMAINING15_INTEGRAL_LAUNCHER_PREPARATION"
           and manifest["status"] == "READY_FOR_INDEPENDENT_INTEGRAL_PRELAUNCH" and not manifest["rejected_attempts"]
           and pre["status"] == "PASS_PREPARATION_ONLY" and yes(pre["checks"])
           and all(v["force_hessian_or_cache_sums_run"] is False for v in (protocol,pre,manifest))
           and rresult["status"] == "PASS_REMAINING15_DOMAINS_REPLAY"
           and prep_started["phase"] == "PURE_ROUTE_PREPARATION" and prep_started["pid"] == prep_term["pid"]
           and prep_started["protocol_sha256"] == PIN["protocol_sha256"]
           and all(all(v[k] == PIN[k] for k in v if k in PIN) for v in (protocol,pre,manifest))
           and adapter["status"] == "PASS" and yes(adapter["checks"]), "frozen pair of preparations and reviewed prerequisites")
    for path,pin in {**protocol["pins"],**protocol["transitive_sources"]}.items():
        read(path,pin,False)
    sys.path.insert(0,str(ROOT/"scripts"))
    L = importlib.import_module("remaining_cohort_integrals")
    insist(rel(L.__file__) == SOURCE and protocol["pins"] == L.PIN and tuple(L.PARTS) == PARTS,"same launcher source")
    fun = L._boot()
    launch_shape(L,fun,protocol,adapter)
    L._terminal(prep_term,PIN["manifest_sha256"])
    L._terminal(rterm,RRESULT_SHA)
    jobs, stats, static = [], {}, {}
    forbidden = AssertionError("actual integral/force/H/cache work during route review")
    packs = L._domains(fun)
    with (patch.object(fun.I,"bulk_component",side_effect=forbidden),
          patch.object(fun.E,"surface_component",side_effect=forbidden),
          patch.object(fun.S,"material_moment",side_effect=forbidden),
          patch.object(fun.S,"_certificate",side_effect=forbidden),
          patch.object(fun.F,"_asset",side_effect=forbidden),
          patch.object(fun.E.core,"local_terms",side_effect=forbidden),
          patch.object(fun.E.U,"uniform_cache",side_effect=forbidden)):
        for record, pack, m, i in zip(IDS,packs,manifest["records"],rresult["records"],strict=True):
            row, terminal, gate = pack
            insist(row["record_id"] == terminal["record_id"] == gate["record_id"] == m["record_id"] == i["record_id"] == record,
                   "record alignment")
            leaf = read(i["file"],i["sha256"])
            insist(leaf["record_id"] == record and leaf["status"] == "PASS_RECORD_DOMAIN_REPLAY"
                   and yes(leaf["checks"])
                   and [leaf[k] for k in ("bulk_elements","surface_edges","surface_cache_sites")] == [144,24,288],
                   "all independent positive domain proofs")
            record_root = f"{BASE}/records/{record}/scaled_integrals_v1/preparation"
            insist(set(m["files"]) == {record_root+"/"+n+".json" for n in ("protocol","preflight","proofs")}
                   and [m[k] for k in ("bulk","surface")] == [144,24], "canonical three-file route set")
            ps, fs, bs = [f"{record_root}/{name}.json" for name in ("protocol","preflight","proofs")]
            pr, fl, proofs = [read(n,m["files"][n]) for n in (ps,fs,bs)]
            insist(pr["source_sha256"] == fl["source_sha256"] == PIN["source_sha256"]
                   and pr["batch_protocol_sha256"] == PIN["protocol_sha256"]
                   and pr["record_id"] == fl["record_id"] == record
                   and pr["independently_reviewed_actual_domains_path"] == i["file"]
                   and pr["original_axis_indices"] == AXIS and pr["required_parts"] == list(PARTS)
                   and pr["proofs_sha256"] == fl["proofs_sha256"] == m["files"][bs]
                   and fl["protocol_sha256"] == m["files"][ps] and fl["status"] == "PASS_ROUTE_PREPARATION_ONLY"
                   and yes(fl["checks"]) and fl["force_hessian_or_cache_sums_run"] is False,
                   f"frozen route set, domain review leaf and no work for {record}")
            data = L._record_data(fun,*pack)
            insist(data.hashes == pr["domain_and_record_input_hashes"], "current context matches frozen record inputs")
            box = body_and_basic(data,leaf,row,terminal,record)
            route_packets(fun,data,proofs,box,static,stats)
            jobs.extend({"record_id":record,"part":part,"domain_review_sha256":i["sha256"]} for part in PARTS)
            print(json.dumps({"record_id":record,"independent_proof_packet_readback":"PASS"}),flush=True)
    insist(len(jobs) == 30 and len(static) == 132 and len(stats) == 15
           and sum(s["symbolic_identities"] for s in stats.values()) == 297360
           and sum(sum(s["routes"].values()) for s in stats.values()) == 2160
           and sum(s["surface_edges"] for s in stats.values()) == 360
           and sum(s["surface_sites"] for s in stats.values()) == 4320,"required actual frozen route counts")
    checks = mocked_supervisor(L)
    checks.update(dict.fromkeys(("frozen_source_preflight_manifest_and_terminal","complete_source_bootstrap_and_aliases",
        "exclusive_worker_output_before_component_callbacks","all_15_independent_domain_leaf_bytes_and_actual_terminal_identity",
        "same_selected_point_C_radius_exact_node_box_and_normalization","all_2160_bulk_and_360_surface_proof_metadata_maps",
        "all_297360_symbolic_identities_reported_with_correct_per_element_counts",
        "original_J_G_rho_minimal_order_and_budgets","all_4320_reviewed_surface_cache_domains_bound",
        "no_L_component_substitution","checked_gate_before_true_child_spawn"), True))
    receipt = {"status":"PASS",**PIN,"adapter_review":ADAPTER,"jobs":jobs,"checks":checks,
               "preparation_terminal":{"pid":prep_term["pid"],"file":f"{P0}/terminal.json"},
               "domain_review_results":{"file":RRESULT,"sha256":RRESULT_SHA},
               "record_proof_counts":stats,"source_review_sha256":sha(Path(__file__).read_bytes()),
               "review_ran_actual_integration":False,
               "concurrency":"Each launch invocation runs one child at a time. The retained controller bounds active supervisors to 2. Concurrent duplicate supervisors can spawn two prechecking workers, but the exclusive output-directory claim prevents duplicate component execution.",
               "terminal_scope":"An approved component requires the supervisor's wait-derived OS terminal and matching successful child PID and result digest. Preserve failed work and logs.",
               "authorization_scope":"Only the 30 explicitly listed same-record bulk/surface parts. Actual evidence remains unreviewed until independently admitted and replayed."}
    gate = BATCH/"prelaunch_independent_review.json"
    bytes_out = (json.dumps(receipt,indent=2,allow_nan=False)+"\n").encode()
    digest = sha(bytes_out)
    check = L.checked

    def staged(path, expected, parse=True):
        if path == rel(gate):
            insist(expected == digest and parse, "readback gate SHA")
            return receipt
        return check(path,expected,parse)

    cmds = [j["record_id"]+":"+j["part"] for j in jobs]
    with patch.object(L,"checked",side_effect=staged):
        _, selected = L._gate(PIN["manifest_sha256"],digest,cmds)
        insist(set(selected) == {(j["record_id"],j["part"]) for j in jobs}, "all explicit jobs admitted at pure gate")
        for invalid in (cmds[:1]+cmds[:1],["optoGEF2_small_col_L:bulk"],["all:bulk"]):
            try:
                L._gate(PIN["manifest_sha256"],digest,invalid)
            except ValueError:
                continue
            raise AssertionError("invalid actual job accepted")
    with gate.open("xb") as handle:
        handle.write(bytes_out)
    _, selected = L._gate(PIN["manifest_sha256"],digest,cmds)
    insist(len(selected) == 30, "gate durability readback")
    print(json.dumps({"status":"PASS","gate_sha256":digest,"jobs":len(selected),"checks":len(checks)}))


if __name__ == "__main__":
    run()

