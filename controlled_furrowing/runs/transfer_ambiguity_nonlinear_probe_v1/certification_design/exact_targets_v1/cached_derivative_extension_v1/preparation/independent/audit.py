"""Independent cold review of the fixed-center cached derivative extension."""
from contextlib import ExitStack
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
PRE=HERE.parent
EXACT=PRE.parents[1]
PRIMARY=EXACT.parents[1]
ROOT=PRIMARY.parents[1]
SRC=ROOT/"scripts/transfer_ambiguity_cached_derivative_extension.py"
GEO=EXACT/"independent_actual/replay/receipt.json"
CACHE=PRIMARY/"certification_design/cached_pilots_v1/independent_actual_preparation/actual_replay/receipt.json"
PINS={
    SRC:"7ad0b7d4d42efe634123989c6092ccf768cedb92a79244cc7aaa5c10a62f370c",
    PRE/"protocol.json":"4af874d6e45836d11b608ae2929ff5120088fdbd549fd3615398c0068b29ebd7",
    PRE/"preflight.json":"899d9c9486aa0412934bf8866d24c68842c0d2085c380089e1852d906f343f4f",
    GEO:"884e0bc76a27465c86b62a202cab491d764aedf030623592238e09adee86d77d",
    CACHE:"4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17"
}
R1,R2=Q(1,2**36),Q(1,2**34)


def read(p):
    return json.loads(p.read_bytes())


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def save(p,a):
    with p.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(a,f,indent=2)
        f.write("\n")
    return sha(p)


def reject(f):
    try:
        f()
    except ValueError:
        return True
    return False


def exact_sizing(m, c):
    import numpy as np
    x=Q(1,7)
    binary_R=[[-2.0,1.25,0.125],[0.5,-0.75,3.0],[1.,0.,-1.125]]
    u=[Q(1,8),Q(1,16),Q(3,8)]
    raw=m.core.arb_mat(3,1,[m.core.A(a) for a in u])
    A=m.abs_matrix(np.asarray(binary_R))
    actual=[sum(abs(Q(b))*a for a,b in zip(u,row,strict=True)) for row in binary_R]
    P,q=m.propagate(A,raw,3)
    assert q==actual and P==max(actual)
    assert reject(lambda:m.propagate(A,-raw,3))
    assert reject(lambda:m.propagate(A,raw,2))
    for z in (Q(0),Q(1,4),Q(1,3),Q(99,100),Q(1),Q(2)):
        observed=m.spend(z,R2,x,c)
        h=1-z
        assert Q(observed["headroom_exact"])==h
        assert observed["strict_positive_headroom"] is (h>0)
        assert observed["below_optional_quarter_reference"] is (z<=Q(1,4))
        assert observed["quarter_reference_is_gate"] is False
        if h<=0:
            assert observed["planning_status"]=="NO_POSITIVE_HEADROOM"
            assert observed["b_star_exact"] is observed["preconditioned_signed_total_force_allowance_exact"] is None
            continue
        assert observed["planning_status"]=="POSITIVE_CONDITIONAL_HEADROOM"
        b=h/2
        F=h*R2/4
        assert Q(observed["b_star_exact"])==b and Q(observed["preconditioned_signed_total_force_allowance_exact"])==F
        assert Q(observed["unspent_strict_selfmap_margin_at_both_planning_caps_exact"])==R2-(z+b)*R2-F==F>0
        for k,bound,ob,old in ((6,F,"force","force"),(72,b,"H","H_entry")):
            a=Q(observed[f"per_component_{ob}_entry_tail_cap_exact"])
            assert a==min(Q(c[f"tail_{old}_budget_exact"]),bound/(10000*x*k))
            assert 0<a and x*k*a<=bound/10000
    return {"synthetic_abs_R_variation_exact":list(map(str,actual)),"synthetic_Z_cases_checked":6}


def main():
    assert all(sha(p)==h for p,h in PINS.items())
    prot,pre,geo,old=map(read,(PRE/"protocol.json",PRE/"preflight.json",GEO,CACHE))
    assert pre["status"]=="PASS" and pre["checks"] and all(v is True for v in pre["checks"].values())
    assert pre["actual_site_AD_assemblies"]==pre["actual_whole_box_variations"]==0
    assert pre["source_sha256"]==prot["source_sha256"]==PINS[SRC]==sha(PRE/"transfer_ambiguity_cached_derivative_extension_frozen.py")
    assert pre["protocol_sha256"]==PINS[PRE/"protocol.json"]
    assert len(prot["input_sha256"])==194 and all(sha(ROOT/f)==h for f,h in prot["input_sha256"].items())
    assert prot["whole_box_site_counts"]=={"bulk":9216,"surface":288} and prot["precision_bits"]==256
    assert prot["threshold_is_gate"] is False and Q(prot["optional_engineering_ZQ_threshold_exact"])==Q(1,4)
    assert Q(prot["old_cached_certificate_radius_exact"])==R1 and Q(prot["larger_radius_exact"])==R2
    assert prot["new_spatial_integrals"]==prot["new_equilibrium_solves"]==prot["new_preconditioners"]==0
    assert not (PRE.parent/"actual").exists()
    geos={x["pilot_id"]:x for x in geo["geometry_records"] if Q(x["radius_exact"])==R2}
    man=read(EXACT/"preparation/target_manifest.json")
    ps={x["pilot_id"]:x for x in man["targets"]}
    cases=prot["ordered_cases"]
    files=[]
    for c,o in zip(cases,old["field_results"],strict=True):
        p=c["pilot_id"]
        entry=ps[p]
        assert p==o["pilot_id"] and p==geos[p]["pilot_id"]
        assert c["whole_box_geometry_review"]==geos[p]["receipt"]
        assert c["whole_box_geometry_review_sha256"]==geos[p]["sha256"]
        assert c["box"]==entry["box_1"] and c["box_sha256"]==entry["box_1_sha256"]==geos[p]["box_sha256"]
        assert c["verified_point_sha256"]==o["receipt_sha256"]
        point=read(ROOT/c["verified_point"])
        model=read(ROOT/entry["model_contract"])
        assert model["larger_radius_cached_derivative_asserted"] is False
        d=model["exact_requested_drive"]
        dose,gain=Q(d["optical_dose_exact"]),Q(d["common_gain_exact"])
        assert Q(c["actual_exact_drive"])==Q(point["exact_strength"])==gain*dose
        assert Q(c["rounded_solver_strength"])==Q(float(gain*dose))!=gain*dose
        assert c["field_sha256"]==point["field_sha256"]==model["source_field_sha256"]
        assert c["cache_packet_sha256"]==point["cache_packet_sha256"]==model["point_cache_sha256"]
        assert c["embedded_cache_sha256"]==point["cache_array_contents_sha256"]==model["frozen_binary_cache_sha256"]
        assert c["preconditioner_binary_sha256"]==model["binary_preconditioner_sha256"]
        assert c["model_contract_sha256"]==entry["model_contract_sha256"]
        eta,z0,norm=map(Q,(c["point_eta_upper_exact"],c["center_Z0_upper_exact"],c["R_inf_exact"]))
        assert eta==Q(point["independent_eta_upper_exact"])
        assert z0==Q(point["independent_Z0_upper_exact"])
        assert norm==Q(point["preconditioner_infinity_norm_exact"])
        signed=point["signed_center_force_intervals"]
        assert len(signed)==611
        assert all(-eta<=Q(h["lower_exact"])<=Q(h["upper_exact"])<=eta for h in signed)
        assert 0<=z0<1 and 0<=eta and norm>0
        assert Q(c["cached_proof_radius_exact"])==Q(point["radius_exact"])==R1
        assert Q(c["new_derivative_radius_exact"])==R2
        assert R1<R2<Q(c["observer_goal_exact"])==Q(model["observer_design_goal_exact"])
        route=read(ROOT/entry["routes"])
        f,h=[0]*611,[0]*611
        for comp in route["bulk"]+route["surface"]:
            inds=comp["global_indices"]
            assert len(inds)==len(set(inds))
            for i in inds:
                f[i]+=1
                h[i]+=len(inds)
        sizing=c["prior_optional_tail_policy"]
        assert max(f)==sizing["max_force_incidence"]==6 and max(h)==sizing["max_H_entry_row_incidence"]==72
        assert Q(sizing["frozen_R_infinity_norm_upper_exact"])>=norm
        assert all(Q(sizing[k])>0 for k in ("tail_force_budget_exact","tail_H_entry_budget_exact"))
        files.append({"pilot_id":p,"own_fixed_cache_sha256":c["embedded_cache_sha256"],"geometry_box_sha256":c["box_sha256"]})
    assert len(files)==4

    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_cached_derivative_extension as m
    def forbid(*args,**kw):
        raise AssertionError("actual calculation forbidden")
    bindings=dict(zip(("source_sha256","protocol_sha256","preflight_sha256",
                       "actual_geometry_review_sha256","cached_root_independent_review_sha256"),
                      PINS.values(),strict=True))
    with ExitStack() as stack:
        for mod,names in ((m,("calculate","worker","launch")),
                          (m.core,("restore_cache","local_terms","assemble_point","variation_rows")),
                          (m.EX,("open_context","geometry_box"))):
            for meth in names:
                stack.enter_context(patch.object(mod,meth,side_effect=forbid))
        stack.enter_context(patch.object(subprocess,"Popen",side_effect=forbid))
        assert prot==m.protocol(m.authority())
        oldprec=m.core.ctx.prec
        m.core.ctx.prec=256
        m.core.A.cache_clear()
        try:
            fixture=exact_sizing(m,cases[0]["prior_optional_tail_policy"])
            case=cases[0]
            proof=read(ROOT/case["verified_point"])
            assert reject(lambda:m.validate_center(proof,cases[-1]))
            assert reject(lambda:m.validate_center(proof,dict(case,actual_exact_drive=case["rounded_solver_strength"])))
            assert reject(lambda:m.validate_center(dict(proof,signed_center_force_intervals=[]),case))
        finally:
            m.core.ctx.prec=oldprec
            m.core.A.cache_clear()
    assert not (PRE.parent/"actual").exists()
    checks={
      "all_194_frozen_pins_and_current_true_terminals_rehashed":True,
      "all_four_eta_Z0_exact_drive_and_R_center_anchors_match":True,
      "new_radius_built_from_same_center_not_inherited_from_geometry_or_old_z":True,
      "all_four_168_route_incidence_counts_recomputed":True,
      "abs_R_times_nonnegative_rows_independently_checked":True,
      "planning_headroom_and_optional_local_caps_exact_and_no_quarter_gate":True,
      "serial_exclusive_worker_exit_and_failure_paths_source_reviewed":True,
      "no_actual_cache_reconstruction_variation_force_or_Hessian":True,
    }
    f=HERE/"support.json"
    h=save(f,{"status":"PASS","checks":checks,"frozen_source_and_gate_pins":bindings,
        "cases":files,"fixture":fixture,"audit_source_sha256":sha(Path(__file__)),
        "actual_recomputation_count":0,
        "math":"At the same fixed requested strength and center, |R| times the whole-row absolute Hessian variation bounds every row of R(H_Q(y)-H_Q(p)). Add prior Z0 only after this same-center match. For h=1-ZQ>0, proposed exact-spatial budgets b=h/2 and total signed force h*r/4 leave h*r/4 strictly unused. Local tails propagate through the verified incidence and exact R norm. Force budget is total, not a defect or bias statement. For h<=0 retain the failed bound and omit budgets.",
        "scope":"Preparation/source review only. Actual larger radius requires a fresh analytic cached-domain and derivative enclosure at all 9504 sites. A valid geometry box alone proves none of those derivative bounds."})
    g=save(PRE/"prelaunch_independent_review.json",{"status":"PASS","checks":checks,**bindings,
             "support_file":f.relative_to(ROOT).as_posix(),"support_sha256":h})
    print(g,h,flush=True)


if __name__=="__main__":
    main()
