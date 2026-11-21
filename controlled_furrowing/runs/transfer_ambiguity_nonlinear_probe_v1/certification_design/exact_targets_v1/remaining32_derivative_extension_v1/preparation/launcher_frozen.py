"""Cold preparation and separately gated whole-box geometry/cached variation."""
import argparse
from collections import Counter
from contextlib import contextmanager
from fractions import Fraction
import os
import re
from tempfile import TemporaryDirectory
from time import sleep
from pathlib import Path
import subprocess
import sys
import traceback
from unittest.mock import patch

import algebraic_root_feasibility as core
import review_algebraic_root as maps
import transfer_ambiguity_remaining32_exact_targets as A

SOURCE=Path(__file__).resolve()
ROOT=A.ROOT
GOUT=A.X.OUT/"remaining32_geometry_v1"
DOUT=A.X.OUT/"remaining32_derivative_extension_v1"
OUTS={"geometry":GOUT,"variation":DOUT}
TARGET=A.PREP/"manifest.json"
TARGET_REVIEW=A.OUT/"independent_prelaunch/review.json"
T_SHA="dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242"
T_SOURCE="d1a9c9e4839ff85a9926b6fd5ecc11e4cd8e0831c7e6263ecfbe068ef22ea3f4"
T_PROTOCOL="a12669fd1e50cbd7435d016ba60223791c8a70119f622048c28e24ef5a9f456b"
T_PREFLIGHT="86ddb06c75c2deedea29d62bb7e03da6142d55b41a0291e1567c9e4dba6aab4f"
G=A.X.G
E=A.E
Q=A.Q
RADIUS=A.R34
N=611
COUNT={"bulk":9216,"surface":288}
MARGINS={"bulk_det_min","bulk_hoop_min","surface_radius_min","surface_tangent_squared_min","normalized_weight_min"}
ZERO_AD="EXISTING_MATCHED_WHOLE_CUBE_CACHED_ZQ"
FRESH={"EXISTING_MATCHED_SMALL_CUBE_LARGER_TRIAL_REQUIRES_NEW_AD","FRESH_VARIATION_REQUIRED_BEFORE_SPATIAL_LAUNCH"}
COLD={"new_geometry_certificates":0,"new_AD_derivatives":0,"new_spatial_integrations":0,"new_solves":0}
SCOPE="Whole closed same-center 2^-34 cubes, preserving independently archived observer fields. Separate exact whole-body analytic geometry and finite cached derivative bounds only."
STAGE_REVIEW="PASS_INDEPENDENT_REMAINING32_GEOMETRY_VARIATION_PRELAUNCH"
HOLD=GOUT/"independent_prelaunch/hold_4dd6.json"
HOLD_SHA="8859b9f22e4781f4533e641cbd3270eb5d7723020fd16ad43a86c0869cad9c5e"
KERNEL={
 Path(core.__file__).resolve():"7321d38bf17cf354cc03be2871f3fbaa45bcebbcbfca7be90d0f867b80078acf",
 Path(maps.__file__).resolve():"6af2311d0667a28af263b6fb6237063481aa8e65209ae9722c2b0246a8cff27a",
 Path(G.__file__).resolve():"18013c6c50c5f679eb6c4abc04e3cf9ce447ee18d0bf1efbf6dcfbd6806f49c2",
 Path(E.__file__).resolve():"092dc51161e708cdba3791ec47a7f8e843dd0e1afc866b6f65f286b8657d1745",
 Path(E.SUR.__file__).resolve():"61cf47d2b1a9e7ee46f9fa680c24b4fdf7b4e31e4889379c56951ab6a7455ee9",
 Path(E.AX.__file__).resolve():"ebe0cbe021e04bf51aaf8985df6bfbe6f6fd20e1fc60194acbd1fd968654199f",
}
sha,read,save,require=A.sha,A.read,A.X.save,A.require
rel=A.rel


def positive_int(x):
    return type(x) is int and x>0


def zero(x):
    return type(x) is int and x==0


def yes(d):
    return isinstance(d,dict) and bool(d) and all(v is True for v in d.values())


def input_phase():
    require(sha(A.SOURCE)==T_SOURCE,"exact adapter source")
    p=read(A.PREP/"protocol.json",T_PROTOCOL)
    f=read(A.PREP/"preflight.json",T_PREFLIGHT)
    m=read(TARGET,T_SHA)
    require(p["source_sha256"]==f["source_sha256"]==m["source_sha256"]==T_SOURCE
        and f["status"]=="PASS" and yes(f["checks"]) and m["protocol_sha256"]==f["protocol_sha256"]==T_PROTOCOL
        and f["manifest_sha256"]==T_SHA and m["geometry_or_variation_actual_launch_authorized"] is False
        and tuple(e["case_id"] for e in m["targets"])==A.IDS,"current complete cold target preparation required")
    closure=p["local_source_closure_sha256"]
    for n,h in closure.items():
        require(sha(A.canonical(n))==h,"transitive frozen target source changed")
    require(core.ctx is maps.ctx and all(sha(k)==h for k,h in KERNEL.items()),"fixed geometry and cached kernels")
    return m,{**closure,**{rel(k):h for k,h in KERNEL.items()}}


def target_gate(path,h,m):
    require(path.resolve()==TARGET_REVIEW.resolve(),"canonical target context review required")
    v=read(path,h)
    require(v["status"]=="PASS_INDEPENDENT_REMAINING32_TARGET_CONTEXTS" and yes(v["checks"])
            and v["source_sha256"]==T_SOURCE and v["protocol_sha256"]==T_PROTOCOL
            and v["manifest_sha256"]==T_SHA and v["preflight_sha256"]==T_PREFLIGHT
            and tuple(v["ordered_case_ids"])==A.IDS and v["derivative_policy_counts"]==m["derivative_policy_counts"],
            "matched independent all32 target context review")
    expected=("case_id","record_id","model_contract_sha256","box_sha256","routes_sha256","rule_sha256",
              "planned_geometry_radius_exact","derivative_policy")
    checks=v["per_case"]
    require(isinstance(checks,list) and len(checks)==32 and all(
        {k:a[k] for k in expected}=={k:b[k] for k in expected}
        for a,b in zip(m["targets"],checks,strict=True)), "every reviewed target input must match")
    return v


def box_pair(old,new):
    require(set(old)==set(new),"same exact coordinate box schema")
    different={"coordinate_radii_exact","coordinate_hulls_exact","corrector_radius_exact"}
    require(all(old[k]==new[k] for k in old if k not in different),"box growth cannot change center, R map or target")
    require(Q(old["corrector_radius_exact"])==A.R36 and Q(new["corrector_radius_exact"])==RADIUS
            and old["root_membership_assumed"] is new["root_membership_assumed"] is False,
            "bound enlargement from preserved smaller box to new trial")
    for b,r in ((old,A.R36),(new,RADIUS)):
        free=b["free"]
        require(type(free) is list and len(free)==len(set(free))==N
                and all(type(k) is int and 0<=k<650 for k in free),"typed free-coordinate injection required")
        f=set(free)
        for j,(y,ds,hulls) in enumerate(zip(b["current_centers_exact"],b["coordinate_radii_exact"],
                                          b["coordinate_hulls_exact"],strict=True)):
            for i,(x,d,ends) in enumerate(zip(y,ds,hulls,strict=True)):
                d,x=Q(d),Q(x)
                require(d==(r if 2*j+i in f else 0) and tuple(map(Q,ends))==(x-d,x+d),
                        "literal exact enlarged cube coordinate endpoint mismatch")


def pair(value):
    return Q(value["lower_exact"]),Q(value["upper_exact"])


def row_bounds(m):
    return [pair(core.pack(sum((core.upper(m[i,j]) for j in range(m.ncols())),core.arb(0))))[1]
            for i in range(m.nrows())]


def majorize(abs_r,v):
    require(v.ncols()==1 and abs_r.nrows()==abs_r.ncols()==v.nrows()
            and all(pair(core.pack(v[i,0]))[0]>=0 for i in range(v.nrows())),
            "nonnegative exact R-weighted row majorant required")
    return row_bounds(abs_r*v)


def stage_cases(stage,rows):
    return tuple(x["case_id"] for x in rows if stage=="geometry" or x["derivative_policy"] in FRESH)


def pack_trial(e):
    c=read(A.canonical(e["model_contract"]),e["model_contract_sha256"])
    preserved=e["box"]
    preserved_h=e["box_sha256"]
    old=read(A.canonical(preserved),preserved_h)
    b,bh=preserved,preserved_h
    nread=[]
    if Q(e["planned_geometry_radius_exact"])==A.R36:
        ctx,_,current=A.open_context(e["case_id"],T_SHA)
        require(current==e and c["derivative_policy"] in FRESH,"same original small cube")
        new=A.node_box(ctx,RADIUS)
        box_pair(old,new)
        path=GOUT/"preparation/enlarged_trial_boxes"/(e["case_id"]+".json")
        path.parent.mkdir(exist_ok=True)
        b,bh=rel(path),save(path,new)
        nread=[c["archived_field"],c["archived_verification"],c["original_point_cache"],c["point_cache"]]
        if c["recenter_exact_step_sha256"]:
            nread.append(rel(A.canonical(c["point_cache"]).with_name("exact_binary_step.json")))
    require(Q(read(A.canonical(b),bh)["corrector_radius_exact"])==RADIUS,"matched planned bridge r34")
    return {**e,"trial_box":b,"trial_box_sha256":bh,"trial_radius_exact":str(RADIUS),
            "archived_center_shift_exact":c["archived_shift_exact"],"center_binary_sha256":c["center_binary_sha256"],
            "R_binary_sha256":c["binary_R_sha256"],"exact_drive":c["exact_requested_drive"]["activation_strength_exact"],
            "point_cache_sha256":c["point_cache_sha256"],"energy_scale_binary_exact":c["energy_scale_binary_exact"],
            "native_cached_bound_usable_at_trial_radius":c["derivative_policy"]==ZERO_AD,
            "native_geometry_claimed":False,"opened_preexisting_context_files_at_prep":sorted(set(nread))}


def current_context(trial):
    ctx,rule,e=A.open_context(trial["case_id"],T_SHA)
    box=read(A.canonical(trial["trial_box"]),trial["trial_box_sha256"])
    old=read(A.canonical(e["box"]),e["box_sha256"])
    require(all(trial[k]==v for k,v in e.items())
            and box==A.node_box(ctx,RADIUS) and Q(trial["trial_radius_exact"])==RADIUS
            and Q(trial["exact_drive"])==ctx.drive and A.a_sha(ctx.R)==trial["R_binary_sha256"]
            and A.a_sha(ctx.u)==trial["center_binary_sha256"]
            and Q(trial["energy_scale_binary_exact"])==ctx.target.energy_scale
            and ctx.row["point_cache_sha256"]==trial["point_cache_sha256"],"same exact proof center, drive R and trial box")
    if Q(e["planned_geometry_radius_exact"])==A.R36:
        box_pair(old,box)
    else:
        require(old==box and e["box_sha256"]==trial["trial_box_sha256"],"native matched planned larger cube unchanged")
    route=read(A.canonical(e["routes"]),e["routes_sha256"])
    return ctx,rule,box,route


def calc_geometry(trial,folder):
    with A.X.no_integrals():
        ctx,_,box,route=current_context(trial)
        body=G.certify_box(box["reference_nodes_exact"],box["current_centers_exact"],
                           box["triangles"],box["coordinate_radii_exact"])
        bh=save(folder/"body_geometry.json",body)
        bulk=[]
        surface=[]
        status="UNRESOLVED_WHOLE_BODY_GEOMETRY"
        if body["status"]=="CERTIFIED_ENTIRE_NODAL_BOX_EMBEDDING":
            require(len(body["triangles"])==len(route["bulk"])==144 and len(route["surface"])==24,
                    "all actual bulk and trace routes required")
            nodes=E.SUR.common_nodes(ctx.target,box)
            for row in route["bulk"]:
                i=row["element"]
                f,ids=E._bulk_fields(ctx.target,nodes,i)
                require(ids==row["nodes"] and f["chart"]["kind"]==row["chart"]
                        and f["free_local"]==row["free_local"] and f["global_indices"]==row["global_indices"]
                        and Q(row["active"])==f["active"],"whole-box chart/projection unchanged")
                try:
                    lo,hi,raw,rho=E._body_rho(body["triangles"][i])
                    bulk.append({**row,"J_bounds_exact":list(map(str,(lo,hi))),
                                 "raw_rho_exact":str(raw),"rho_exact":str(rho)})
                except ValueError as err:
                    bulk.append({**row,"status":"UNRESOLVED_BULK_DOMAIN","error":str(err)})
            for row in route["surface"]:
                f,e=row["face"],row["edge"]
                fields=ctx.target.surface(f,e,nodes=nodes)
                require(fields["free_local"]==row["free_local"] and fields["global_indices"]==row["global_indices"],
                        "whole-box surface projection unchanged")
                try:
                    d=E.SUR.uniform_G_bounds(fields,max_depth=8)
                    sites=[E.SUR.local_site_domain(ctx.target,box,f,e,k) for k in range(12)]
                    require(d["status"]=="CERTIFIED_UNIFORM_HOOP_G_AND_BINOMIAL_DOMAIN","whole-box surface trace domain")
                    surface.append({**row,"domain":d,"cache_site_geometry":sites})
                except ValueError as err:
                    surface.append({**row,"status":"UNRESOLVED_SURFACE_DOMAIN","error":str(err)})
            status=("PASS_CANDIDATE_WHOLE_BOX_GEOMETRY_AWAITS_INDEPENDENT_REVIEW"
                    if all("rho_exact" in x for x in bulk) and all("domain" in x for x in surface)
                    else "UNRESOLVED_WHOLE_BOX_ANALYTIC_SUBDOMAIN")
    return {"status":status,"whole_body_sha256":bh,"bulk_geometry_only":bulk,"trace_geometry_only":surface,
            "scope":SCOPE,"new_AD_derivatives":0,"new_spatial_integrations":0,"root_asserted":False}
def calc_variation(trial,folder):
    ctx,_,_,_=current_context(trial)
    require(trial["derivative_policy"] in FRESH and trial["native_cached_bound_usable_at_trial_radius"] is False,
            "only the ten matched fresh cached variation cases")
    core.ctx.prec=256
    core.A.cache_clear()
    maps.A.cache_clear()
    counts=dict.fromkeys(COUNT,0)
    def terms():
        for site in maps.terms_from_arrays(ctx.target.a,ctx.target.packet["parameters"],ctx.drive):
            counts[site["kind"]]+=1
            yield site
    flags={"scope":SCOPE,"new_spatial_integrations":0,"exact_spatial_derivatives_computed":False,
           "root_asserted":False,"cached_root_retested_but_not_admitted":True,
           "independent_closed_cube_derivative_review_required":True,"precision_bits":256}
    try:
        g,H,sites,min_point=core.assemble_point(terms(),N)
        require(counts==COUNT and len(sites)==sum(COUNT.values())
                and pair(core.pack(min_point))[0]>Fraction(1,10**6),"complete positive cached point")
        rvalues=list(map(A.T.rational,ctx.R.ravel()))
        R=core.arb_mat(N,N,list(map(core.A,rvalues)))
        RA=core.arb_mat(N,N,[core.A(abs(x)) for x in rvalues])
        identity=core.arb_mat(N,N)
        for i in range(N):
            identity[i,i]=1
        rf=R*g
        signed=[core.pack(rf[i,0]) for i in range(N)]
        eta=max(max(map(abs,pair(x))) for x in signed)
        center_rows=row_bounds(identity-R*H)
        variation,min_j,domain=core.variation_rows(sites,N,core.A(RADIUS))
        margins=domain["strict_lower_margins"]
        require(domain["term_counts"]==COUNT and set(margins)==MARGINS
                and pair(core.pack(min_j))[0]>Fraction(1,10**6)
                and all(pair(x)[0]>0 for x in margins.values()),"complete positive cached whole cube")
        scale=min(pair(core.pack(pt["parameters"][2]))[0] for pt in sites if pt["kind"]=="surface")
        minimum_surface=core.A(scale)*core.A(pair(margins["surface_radius_min"])[0])*core.A(
            pair(margins["surface_tangent_squared_min"])[0]).sqrt()
        require(pair(core.pack(minimum_surface))[0]>0,"positive whole cached surface Jacobian")
        propagated=majorize(RA,variation)
    except ValueError as error:
        return {**flags,"status":"UNRESOLVED_FRESH_CACHED_ANALYTIC_DOMAIN","reason":str(error)}
    z0,dz=max(center_rows),max(propagated)
    z=z0+dz
    lhs=eta+z*RADIUS
    passed=z<1 and lhs<RADIUS
    return {**flags,"status":"FRESH_CACHED_VARIATION_AWAITS_INDEPENDENT_REVIEW","site_counts":counts,
            "point_minimum_J":core.pack(min_point),"whole_minimum_bulk_J":core.pack(min_j),
            "whole_minimum_surface_J":core.pack(minimum_surface),"whole_analytic_domain":domain,
            "R_inf_exact":str(max(sum(abs(x) for x in rvalues[i:i+N]) for i in range(0,len(rvalues),N))),
            "signed_RFQ_center_rows":signed,"eta_Q_upper_exact":str(eta),"center_Z0_upper_exact":str(z0),
            "center_defect_row_bounds_exact":list(map(str,center_rows)),
            "unmapped_variation_rows_exact":[str(pair(core.pack(variation[i,0]))[1]) for i in range(N)],
            "preconditioned_variation_rows_exact":list(map(str,propagated)),
            "variation_upper_exact":str(dz),"ZQ_upper_exact":str(z),
            "auxiliary_cached_selfmap_lhs_exact":str(lhs),"auxiliary_cached_selfmap_pass":passed,
            "auxiliary_cached_posterior_exact":str(eta/(1-z)) if passed else None}


def bindings():
    m,_=input_phase()
    data={}
    common={"source_sha256":sha(SOURCE),"target_manifest_sha256":T_SHA}
    for stage,root in OUTS.items():
        folder=root/"preparation"
        ph,fh=sha(folder/"protocol.json"),sha(folder/"preflight.json")
        p=read(folder/"protocol.json",ph)
        f=read(folder/"preflight.json",fh)
        require(all(p[k]==f[k]==v for k,v in common.items()) and f["protocol_sha256"]==ph
                and f["status"]=="PASS" and yes(f["checks"]) and p["stage"]==f["stage"]==stage
                and all(zero(f[k]) for k in COLD) and p["case_order"]==[x["case_id"] for x in p["trials"]]
                and sha(folder/"launcher_frozen.py")==common["source_sha256"],
                "current frozen stage preparation")
        require(tuple(p["case_order"])==stage_cases(stage,m["targets"]),"canonical stage coverage")
        data[stage]=(p,ph,fh)
    return m,data,common


def review_document(receipt,expected,orders):
    require(receipt["status"]==STAGE_REVIEW and yes(receipt["checks"])
            and all(receipt[k]==v for k,v in expected.items())
            and receipt["permitted_case_ids"]==orders
            and receipt["exact_spatial_bridge_authorized"] is False,"same cold launcher review and exact phase scope")


def launch_gate(args):
    m,stages,common=bindings()
    target_gate(args.target_review,args.target_review_sha256,m)
    extra={"target_context_review_sha256":args.target_review_sha256}
    for key,(p,ph,fh) in stages.items():
        extra[key+"_protocol_sha256"]=ph
        extra[key+"_preflight_sha256"]=fh
    expected={**common,**extra}
    rev=read(args.review,args.review_sha256)
    review_document(rev,expected,{k:v[0]["case_order"] for k,v in stages.items()})
    return stages[args.stage][0],{**expected,"launcher_review_sha256":args.review_sha256}


def child_admit(run_code,result,terminal,identity,digest):
    require(zero(run_code) and zero(terminal["exit_code"]) and terminal["status"]=="FINISHED"
            and all(positive_int(identity[k]) and positive_int(terminal[k]) and positive_int(result[k])
                    for k in ("pid","supervisor_pid"))
            and all(terminal[k]==result[k]==v for k,v in identity.items())
            and terminal["results_sha256"]==digest and result["root_asserted"] is False,
            "real authenticated child exit and current non-promoted result")


def launch_files(stage,group):
    return OUTS[stage]/"actual"/"launches"/group


def worker(args):
    p,links=launch_gate(args)
    require(re.fullmatch("[0-9a-f]{64}",args.group) is not None and positive_int(args.supervisor_pid),
            "canonical launch group and positive supervisor PID")
    parent=launch_files(args.stage,args.group)
    start=read(parent/"started.json",args.supervisor_started_sha256)
    require(positive_int(start["pid"]) and start["pid"]==os.getppid()==args.supervisor_pid
            and start["bindings"]==links
            and start["stage"]==args.stage and args.case_id[0] in start["case_order"]
            and E._digest(start["case_order"])==args.group,
            "only observed actual authorized supervisor can start worker")
    case=args.case_id[0]
    trial=A.pick(p["trials"],case)
    folder=OUTS[args.stage]/"actual"/"cases"/case
    require(type(start["case_order"]) is list and start["launch_group"]==args.group
            and not (folder/"started.json").exists(),"exclusive case worker")
    mark={"pid":os.getpid(),"supervisor_pid":args.supervisor_pid,"launch_group":args.group,
          "case_id":case,"stage":args.stage,"bindings":links}
    save(folder/"started.json",{"status":"RUNNING",**mark})
    term={"status":"FAILED","exit_code":1,**mark}
    try:
        calc=calc_geometry if args.stage=="geometry" else calc_variation
        result=calc(trial,folder)
        h=save(folder/"results.json",{**result,**mark,"trial":trial})
        term.update(status="FINISHED",exit_code=0,results_sha256=h)
    except BaseException as error:
        save(folder/"failure.json",{**mark,"error":str(error),"traceback":traceback.format_exc()})
        raise
    finally:
        save(folder/"terminal.json",term)


@contextmanager
def stage_lock(folder):
    """One active supervisor per stage. The OS releases the lock on exit."""
    folder.mkdir(parents=True,exist_ok=True)
    with (folder/"supervisor.lock").open("a+b") as lock:
        lock.seek(0)
        try:
            if os.name=="nt":
                import msvcrt
                msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
            else:
                import fcntl
                fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except OSError as error:
            raise ValueError("already an active supervisor for this stage") from error
        yield


def run(args):
    p,links=launch_gate(args)
    with stage_lock(OUTS[args.stage]/"actual"):
        run_locked(args,p,links)


def reap_owned_child(proc):
    """Return only a real code from our process handle, even during interrupt cleanup."""
    terminated=False
    while True:
        try:
            code=proc.poll()
            if type(code) is int:
                return code
            if not terminated:
                proc.terminate()
                terminated=True
            code=proc.wait(timeout=1)
            if type(code) is int:
                return code
        except BaseException:
            # The supervisor still owns the process and its stage lock. Retry
            # until the OS confirms exit, including when terminate itself fails.
            try:
                sleep(.05)
            except BaseException:
                pass


def run_owned_child(cmd,folder,base_identity,stage_workers):
    proc=None
    ident=None
    code=None
    with (folder/"stdout.log").open("x",encoding="utf8") as log:
        try:
            proc=subprocess.Popen(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            ident={**base_identity,"pid":proc.pid}
            require(positive_int(proc.pid),"OS process handle requires a typed positive PID")
            save(folder/"launched.json",{**ident,"command":cmd})
            code=proc.wait()
            require(type(code) is int and code==proc.returncode,"observed child exit code required")
        except BaseException:
            # The finally block reaps before this original exception can leave
            # this function or let the caller release its stage lock.
            code=None
            raise
        finally:
            if proc is not None and code is None:
                real_code=reap_owned_child(proc)
                note={**(ident or {**base_identity,"pid":proc.pid}),
                      "status":"REAPED_OWNED_CHILD_AFTER_SUPERVISOR_ERROR",
                      "observed_exit_code":real_code,"scientific_candidate_admitted":False}
                stage_workers.append(note)
                save(folder/"abnormal_child_reap.json",note)
    return ident,code


def guard_previous_children(folder,stage):
    """An unclosed reservation also blocks restart after abrupt supervisor death."""
    for case in folder.iterdir():
        require(case.is_dir(),"canonical per-case reservation directory")
        normal=case/"observed_os_exit.json"
        abnormal=case/"abnormal_child_reap.json"
        choices=[p for p in (normal,abnormal) if p.is_file()]
        require(len(choices)==1,"unclosed previous child reservation requires separately observed recovery")
        path=choices[0]
        obs=read(path,sha(path))
        require(obs["case_id"]==case.name and obs["stage"]==stage and positive_int(obs["pid"])
                and positive_int(obs["supervisor_pid"]) and type(obs["observed_exit_code"]) is int
                and (path==normal or (
                    obs["status"]=="REAPED_OWNED_CHILD_AFTER_SUPERVISOR_ERROR"
                    and obs["scientific_candidate_admitted"] is False)),
                "typed prior owned-child observation required before next launch")


def run_locked(args,p,links):
    ids=args.case_id or p["case_order"]
    selected=set(ids)
    require(len(selected)==len(ids) and ids==[i for i in p["case_order"] if i in selected],
            "explicit unique canonical stage scope")
    group=E._digest(ids)
    out=OUTS[args.stage]/"actual"
    parent=launch_files(args.stage,group)
    reserve_parent=out/"cases"
    reserve_parent.mkdir(parents=True,exist_ok=True)
    guard_previous_children(reserve_parent,args.stage)
    require(all(not (reserve_parent/i).exists() for i in ids),"retain previous results; reject duplicate case start")
    parent.mkdir(parents=True,exist_ok=False)
    pid=os.getpid()
    marker={"pid":pid,"stage":args.stage,"case_order":ids,"launch_group":group,"bindings":links}
    h=save(parent/"started.json",{"status":"RUNNING",**marker})
    terminal={"status":"FAILED","exit_code":1,"workers":[],**marker}
    try:
        for case in ids:
            folder=reserve_parent/case
            folder.mkdir(exist_ok=False)
            cmd=[sys.executable,"-B",str(SOURCE),"worker","--stage",args.stage,"--case-id",case,
                 "--supervisor-pid",str(pid),"--group",group,"--supervisor-started-sha256",h,
                 "--review",str(args.review),"--review-sha256",args.review_sha256,
                 "--target-review",str(args.target_review),"--target-review-sha256",args.target_review_sha256]
            base_identity={**{k:marker[k] for k in ("stage","launch_group","bindings")},
                           "case_id":case,"supervisor_pid":pid}
            ident,code=run_owned_child(cmd,folder,base_identity,terminal["workers"])
            end=folder/"terminal.json"
            os_result={**ident,"observed_exit_code":code,"terminal_sha256":sha(end) if end.is_file() else None}
            os_h=save(folder/"observed_os_exit.json",os_result)
            closed={**os_result,"os_exit_sha256":os_h}
            terminal["workers"].append(closed)
            require(os_result["terminal_sha256"] is not None,"real worker terminal required")
            t=read(end,os_result["terminal_sha256"])
            result=folder/"results.json"
            rh=sha(result) if result.is_file() else None
            require(rh is not None,"worker did not produce scientific candidate/readback result")
            res=read(result,rh)
            child_admit(code,res,t,ident,rh)
            closed.update(results_sha256=rh,candidate_status=res["status"])
            print(case,res["status"],flush=True)
        candidate={"status":"COMPLETE_ACTUAL_WORK_AWAITS_INDEPENDENT_REPLAY","workers":terminal["workers"],
                   "root_asserted":False,**marker}
        rh=save(parent/"manifest.json",candidate)
        terminal.update(status="FINISHED",exit_code=0,manifest_sha256=rh)
        print({"stage":args.stage,"launch_group":group,"manifest_sha256":rh},flush=True)
    except BaseException as error:
        terminal["error"]=str(error)
        raise
    finally:
        save(parent/"terminal.json",terminal)


class FixtureChild:
    pid=711

    def __init__(self, fail):
        self.fail=fail
        self.returncode=None
        self.trace=[]

    def poll(self):
        self.trace.append("poll")
        return self.returncode

    def wait(self,timeout=None):
        self.trace.append(("wait",timeout))
        if timeout is None and self.fail=="wait":
            raise KeyboardInterrupt()
        self.returncode=23
        return self.returncode

    def terminate(self):
        self.trace.append("terminate")


def lifecycle_fixtures():
    results={}
    initial={"case_id":"test","stage":"geometry","supervisor_pid":313,"launch_group":"f"*64,"bindings":{}}
    real_save=save
    with TemporaryDirectory(dir=GOUT) as tmp:
        scratch=Path(tmp)
        for kind in ("launched_record","wait"):
            root=scratch/kind
            root.mkdir()
            folder=root/"test"
            folder.mkdir()
            child=FixtureChild(kind)
            def save_fixture(p,value):
                if kind=="launched_record" and p.name=="launched.json":
                    raise OSError("synthetic launch-record write failure")
                return real_save(p,value)
            observed=[]
            with stage_lock(root), patch.object(subprocess,"Popen",return_value=child), \
                 patch.object(sys.modules[__name__],"save",side_effect=save_fixture):
                try:
                    run_owned_child(["fixture"],folder,initial,observed)
                except (OSError,KeyboardInterrupt) as error:
                    correct_error=(kind=="wait" and isinstance(error,KeyboardInterrupt)
                                   or kind=="launched_record" and isinstance(error,OSError))
                else:
                    correct_error=False
                require(child.returncode==23 and child.trace[-3:]==["poll","terminate",("wait",1)]
                        and len(observed)==1 and observed[0]["observed_exit_code"]==23
                        and (folder/"abnormal_child_reap.json").is_file(),
                        "exceptional return before actual owned child exited")
            cases=root/"cases"
            cases.mkdir()
            # Guard reads the exact bytes that were saved above under the case reservation.
            folder.rename(cases/"test")
            guard_previous_children(cases,"geometry")
            (cases/"unobserved").mkdir()
            require(A.reject(lambda:guard_previous_children(cases,"geometry")),
                    "abrupt unresolved child reservation must block restart")
            results["reap_on_"+kind+"_failure"]=correct_error
    results["unobserved_previous_launch_blocks_reacquired_supervisor"]=True
    return results


def fixtures():
    A.X.bootstrap()
    w=core.arb_mat(2,2,[1,2,3,4])
    v=core.arb_mat(2,1,[1,2])
    vv=majorize(w,v)
    negative=core.arb_mat(2,1,[-1,2])
    t={"status":"FINISHED","exit_code":0,"pid":7,"supervisor_pid":13,"case_id":"a","results_sha256":"f"*64}
    result={"pid":7,"supervisor_pid":13,"case_id":"a","root_asserted":False}
    ident={"pid":7,"supervisor_pid":13,"case_id":"a"}
    child_admit(0,result,t,ident,"f"*64)
    exp={"source_sha256":"0"*64}
    orders={"geometry":["test_a","test_b"],"variation":["test_b"]}
    review={"status":STAGE_REVIEW,"checks":{"phase":True},**exp,"permitted_case_ids":orders,
            "exact_spatial_bridge_authorized":False}
    review_document(review,exp,orders)
    with TemporaryDirectory() as scratch:
        directory=Path(scratch)
        with stage_lock(directory):
            def duplicate():
                with stage_lock(directory):
                    raise AssertionError("duplicate stage supervisor acquired its lock")
            lock_checked=A.reject(duplicate)
        with stage_lock(directory):
            pass
    return {**lifecycle_fixtures(), "stage_os_lock_exclusive_and_released":lock_checked,
            "reject_pid_float_aliases_and_sibling_case":all(A.reject(lambda:child_admit(0,b,t,ident,"f"*64))
                for b in ({**result,"pid":7.0},{**result,"supervisor_pid":13.0},{**result,"case_id":"b"})),
            "exact_nonnegative_R_times_component_majorant":vv==[5,11] and A.reject(lambda:majorize(w,negative)),
            "reject_fake_numeric_exit_or_foreign_pid":A.reject(lambda:child_admit(False,result,t,ident,"f"*64))
                and A.reject(lambda:child_admit(0,result,t,{**ident,"pid":8},"f"*64)),
            "separate_review_schema_and_no_spatial_promotion":A.reject(lambda:review_document(
               {**review,"exact_spatial_bridge_authorized":True},exp,orders)),
            "reject_mixed_stage_case_bound":A.reject(lambda:review_document(
               {**review,"permitted_case_ids":dict(orders,geometry=["test_b","test_a"])},exp,orders))}


def prepare():
    m,kernels=input_phase()
    require(sha(HOLD)==HOLD_SHA,"retain independent predecessor HOLD")
    checks=fixtures()
    require(yes(checks),"cold launcher synthetic fixture failure")
    for root in OUTS.values():
        (root/"preparation").mkdir(parents=True,exist_ok=False)
    with patch.object(G,"certify_box",side_effect=AssertionError("cold phase forbids geometry")), A.X.no_integrals(), \
         patch.object(core,"variation_rows",side_effect=AssertionError("cold phase forbids cached variation")), \
         patch.object(maps,"terms_from_arrays",side_effect=AssertionError("cold phase forbids cached site build")):
        rows=[pack_trial(x) for x in m["targets"]]
    counts=Counter(x["derivative_policy"] for x in rows)
    require(len(rows)==32 and counts[ZERO_AD]==22 and sum(v for k,v in counts.items() if k in FRESH)==10
            and sum(Q(r["planned_geometry_radius_exact"])==A.R36 for r in rows)==3,
            "all geometry 32, only 10 fresh variation, 3 explicitly enlarged")
    enlarged=[x for x in rows if Q(x["planned_geometry_radius_exact"])==A.R36]
    e=enlarged[0]
    old=read(A.canonical(e["box"]),e["box_sha256"])
    box=read(A.canonical(e["trial_box"]),e["trial_box_sha256"])
    malformed={**box,"coordinate_hulls_exact":[[[0,0],[0,0]],*box["coordinate_hulls_exact"][1:]]}
    require(A.reject(lambda:box_pair(old,malformed)),"reject incorrect enlarged coordinate endpoint")
    checks["three_explicit_same_center_larger_boxes_and_nontransfer_of_small_ZQ"]=all(
        not x["native_cached_bound_usable_at_trial_radius"] and x["trial_box_sha256"]!=x["box_sha256"] for x in enlarged)
    checks["no_variation_for_22_exact_matched_existing_ZQ_true"]=all(
        x["native_cached_bound_usable_at_trial_radius"] for x in rows if x["derivative_policy"]==ZERO_AD)
    checks["malformed_enlarged_exact_cube_rejected"]=True
    links={"source_sha256":sha(SOURCE),"target_manifest_sha256":T_SHA}
    require(yes(checks),"cold stage checks required")
    outputs={}
    for stage,root in OUTS.items():
        folder=root/"preparation"
        order=stage_cases(stage,rows)
        p={"status":"FROZEN_EXPLICIT_REMAINING32_STAGE","stage":stage,**links,**COLD,"scope":SCOPE,
           "geometry_or_variation_actual_requires_independent_review":True,"exact_spatial_bridge_authorized":False,
           "case_order":list(order),"trials":[x for x in rows if x["case_id"] in order],
           "cold_case_payloads_opened":[p for row in enlarged for p in row["opened_preexisting_context_files_at_prep"]],
           "kernel_source_sha256":kernels,"target_protocol_sha256":T_PROTOCOL,"target_preflight_sha256":T_PREFLIGHT,
           "closed_radius_exact":str(RADIUS),"max_concurrent_workers_per_supervisor":1,
           "superseded_cold_hold":{"file":rel(HOLD),"sha256":HOLD_SHA},
           "stage_actual_independent_review_required":True}
        ph=save(folder/"protocol.json",p)
        fh=save(folder/"preflight.json",{"status":"PASS",**links,**COLD,"stage":stage,
            "protocol_sha256":ph,"checks":checks})
        with (folder/"launcher_frozen.py").open("xb") as f:
            f.write(SOURCE.read_bytes())
        outputs[stage]={"protocol_sha256":ph,"preflight_sha256":fh}
    print({"source_sha256":sha(SOURCE),**outputs},flush=True)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    sub=cli.add_subparsers(dest="mode",required=True)
    sub.add_parser("prepare")
    for stage in ("run","worker"):
        p=sub.add_parser(stage)
        p.add_argument("--stage",choices=tuple(OUTS),required=True)
        p.add_argument("--case-id",action="append",required=stage=="worker")
        for flag in ("review","target-review"):
            p.add_argument("--"+flag,type=lambda x:Path(x).resolve(),required=True)
            p.add_argument("--"+flag+"-sha256",required=True)
        if stage=="worker":
            p.add_argument("--group",required=True)
            p.add_argument("--supervisor-pid",required=True,type=int)
            p.add_argument("--supervisor-started-sha256",required=True)
    args=cli.parse_args()
    if args.mode=="prepare":
        prepare()
    elif args.mode=="run":
        run(args)
    else:
        require(len(args.case_id)==1,"exact single-case worker")
        worker(args)
