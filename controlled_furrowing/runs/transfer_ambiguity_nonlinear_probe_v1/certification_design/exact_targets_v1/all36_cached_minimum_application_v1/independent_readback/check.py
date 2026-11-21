"""Exact independent scalar readback of the 36 existing cached-energy proofs."""
import base64
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import zipfile

import numpy as np
from numpy.lib import format as npfmt

HERE = Path(__file__).resolve().parent
APP = HERE.parent
ROOT = next(p for p in APP.parents if (p / "PLAN.md").is_file())
EXACT = APP.parent
RADIUS = Q(1, 2**34)
EXCEPTION = "max_contrast_change_small_02"
EXPECTED = {
    "apply.py":"613866091ccfe49324d07f5c7ac3cf5f8a8f60c48af6135fa980188168247838",
    "apply_frozen.py":"613866091ccfe49324d07f5c7ac3cf5f8a8f60c48af6135fa980188168247838",
    "protocol.json":"ae5c8423d31192fab7cdb6715a3f791facd7c2c42c0b22f955dd8455a22e3667",
    "preflight.json":"a559ecc59c45e8b0efd0a28bf6bce799185ffe57b8fa21cf465ffaf388a1c692",
    "root_source_review.json":"82af7595bf247a2e496204c9ed2fc63c7e20ee26362dd617dd0030cda371db60",
    "actual/os_exit_observation.json":"560f5b338506d818e086e1473ff9fab3eee639bb41a2c41160ace6f4687e695a",
    "actual/terminal.json":"74722528788760a207457eb77ce380e9c5182a4d9b7cbe2b98c655eba9fda58e",
    "actual/results.json":"7482665cc52e0dea14514c2d539fd470dc5081633305a899cff200ba2b450400",
}
MARGINS = {"bulk_det_min","bulk_hoop_min","bulk_J_min","normalized_weight_min","surface_radius_min",
           "surface_tangent_squared_min","surface_J_min"}


def require(state, message):
    if not state:
        raise ValueError(message)


def q(x):
    require(type(x) in (str,int,Q), "exact scalar required")
    return Q(x)


def select(rows, key, value):
    row, = (r for r in rows if r[key] == value)
    return row


def digest(b):
    return sha256(b).hexdigest()


def raw(path, h):
    b = path.read_bytes()
    require(digest(b) == h, str(path))
    return b


def domain(d):
    require(d["counts"] == {"bulk":9216,"surface":288}, "complete cached rule")
    lower = dict(d["strict_lower_margins"])
    if "bulk_J_raw" in d:
        lower["bulk_J_min"] = d["bulk_J_raw"]
        lower["surface_J_min"] = d["surface_J_lower"]
    require(set(lower) == MARGINS, "all positivity and regularity premises")
    ends = {k: (q(item["lower_exact"]),q(item["upper_exact"])) for k,item in lower.items()}
    require(all(0 < a <= b for a,b in ends.values()) and ends["bulk_J_min"][0] > Q(1,10**6),
            "C2 cached branch needs positive log/hoop/trace/weights and inactive low-J term")
    return {k:str(a) for k,(a,b) in ends.items()}


def native(n, case):
    scheme = n["status"]
    sigma = Q(0)
    if "own_radius_attempts" in n:
        require(scheme == "PASS_INDEPENDENT_CACHED_LOCAL_ROOT", case)
        p = select(n["own_radius_attempts"],"radius_exact",n["chosen_radius_exact"])
        require(p["status"] == scheme and p["whole_cube_AD_recomputed_for_this_radius"] is True, case)
        names = ("eta_upper_exact","ZQ_upper_exact","strict_selfmap_lhs_upper_exact","posterior_error_upper_exact")
        c = "analytic_domain"
        packet = n["cache_packet_sha256"]
        field = n["field_sha256"]
        goal = n["observer_radius_budget_exact"]
    elif scheme == "PASS_INDEPENDENT_RECENTERED_CACHED_LOCAL_ROOT_AND_OLD_OBSERVER_BUDGET":
        p = n
        names = ("eta_upper_exact","ZQ_upper_exact","strict_selfmap_lhs_upper_exact","posterior_error_upper_exact")
        c = "positive_whole_cached_domain"
        packet = n["recentered_cache_packet_sha256"]
        field = n["archived_field_sha256"]
        goal = n["archived_observer_goal_exact"]
        sigma = q(n["actual_binary_center_shift_sigma_exact"])
        require(sigma > 0 and n["fresh_center_and_entire_cube_recomputed"] is True
                and n["original_field_or_scores_replaced"] is False, case)
    elif scheme == "PASS_INDEPENDENT_NEGATIVE_M_SMALLER_CACHED_ROOT":
        p = n
        names = ("eta_upper_exact","independent_ZQ_upper_exact","strict_selfmap_lhs_upper_exact","posterior_error_upper_exact")
        c = "independent_whole_cube_domains"
        require(n["original_unresolved_test_preserved"] is True, case)
        packet,field,goal = (n[k] for k in ("cache_packet_sha256","field_sha256","observer_goal_exact"))
    else:
        p = n
        require(scheme in ("PASS_INDEPENDENT_CACHED_LOCAL_ROOT","PASS_INDEPENDENT_M_CACHED_LOCAL_ROOT"), case)
        if scheme == "PASS_INDEPENDENT_M_CACHED_LOCAL_ROOT":
            names = ("independent_eta_upper_exact","independent_ZQ_upper_exact","strict_selfmap_lhs_upper_exact",
                     "independent_posterior_error_upper_exact")
            goal = n["observer_goal_exact"]
        else:
            names = ("independent_eta_upper_exact","independent_z_upper_exact","independent_selfmap_lhs_upper_exact",
                     "independent_posterior_error_upper_exact")
            goal = n["science_radius_budget_exact"]
        c = "analytic_domain"
        packet,field = (n[k] for k in ("cache_packet_sha256","field_sha256"))
    eta,z,lhs,e = (q(p[k]) for k in names)
    r = q(p["radius_exact"])
    require(n.get("case_id",n.get("pilot_id")) == case and 0 <= z < 1 and eta >= 0
            and eta + r*z == lhs < r and eta == e*(1-z) and 0 <= e < r
            and sigma+e < q(goal) and n.get("rounded_strength_used_in_proof",False) is False, case)
    return eta,z,r,e,sigma,packet,field,q(goal),domain(p[c])


def array_payload(bundle, member, shape):
    with zipfile.ZipFile(BytesIO(bundle)) as zf, zf.open(member+".npy") as f:
        version = npfmt.read_magic(f)
        dims, fort, dtype = { (1,0):npfmt.read_array_header_1_0,(2,0):npfmt.read_array_header_2_0}[version](f)
        require(dims == shape and dtype.str == "<f8", "same binary coordinate payload")
        payload = f.read()
        if fort and len(shape) > 1:
            cells = np.frombuffer(payload, dtype="V8")
            return cells.reshape(shape, order="F").tobytes(order="C")
        return payload


def cube(box,packet,ident):
    b = packet["frozen_binary_cache"]
    cache = base64.b64decode(b["data"],validate=True)
    require(len(cache) == b["bytes"] and digest(cache) == b["sha256"] == ident["embedded_cache_sha256"],
            "frozen NPZ bytes")
    # The preconditioner is authenticated as raw bytes. No numerical matrix is formed.
    rhash = digest(array_payload(cache,"preconditioner",(611,611)))
    uhash = digest(array_payload(cache,"u",(611,)))
    require(rhash == ident["R_binary_sha256"] and uhash == ident.get("center_binary_sha256",uhash),"same R and center")
    with np.load(BytesIO(cache),allow_pickle=False) as data:
        X = data["reference"].tolist()
        u = data["u"].tolist()
        fr = data["free"].tolist()
        mesh = data["triangles"].tolist()
        norm = data["energy_scale"].tolist()
    exact_X = [list(map(Q,row)) for row in X]
    rim = Q(packet["parameters"]["outer_radius"])
    indexes = [2*i+k for i,(x,y) in enumerate(exact_X) for k in range(2) if x != rim and (x != 0 or k == 1)]
    require(len(X)==325 and len(mesh)==144 and fr == indexes == box["free"] and len(fr)==611
            and mesh == box["triangles"] and len(norm)==1 and Q(norm[0]) == q(ident["energy_scale_binary_exact"])
            and q(box["corrector_radius_exact"]) == RADIUS
            and box["root_membership_assumed"] is False and box["any_spatial_derivative_bound_assumed"] is False,
            "same exact finite dimension, reference map, cache quadrature and cube flags")
    require([list(map(q,row)) for row in box["reference_nodes_exact"]] == exact_X
            and len(box["predictor_hulls_exact"]) == 611
            and all(tuple(map(q,h)) == (0,0) for h in box["predictor_hulls_exact"]), "fixed point singleton")
    disp = {i:Q(value) for i,value in zip(fr,u,strict=True)}
    coordinates = list(zip(box["current_centers_exact"],box["coordinate_radii_exact"],box["coordinate_hulls_exact"],exact_X,strict=True))
    require(len(coordinates)==325, "full nodal array")
    for i,(center,width,hulls,ref) in enumerate(coordinates):
        for j,(center_j,w,bounds,x) in enumerate(zip(center,width,hulls,ref,strict=True)):
            k=2*i+j
            c = x+disp.get(k,0)
            r = RADIUS if k in disp else Q(0)
            require(q(center_j)==c and q(w)==r and tuple(map(q,bounds))==(c-r,c+r),"exact X+u ball")
    return rhash,uhash,exact_X,fr,u


def main():
    top = {n:json.loads(raw(APP/n,h)) for n,h in EXPECTED.items() if n.endswith(".json")}
    t = top["actual/terminal.json"]
    o = top["actual/os_exit_observation.json"]
    result = top["actual/results.json"]
    protocol = top["protocol.json"]
    preflight = top["preflight.json"]
    gate = top["root_source_review.json"]
    for s in ("apply.py","apply_frozen.py"):
        raw(APP/s,EXPECTED[s])
    phase = result["source_and_phase"]
    require(type(o["observed_exit_code"]) is int and o["observed_exit_code"] == 0 and o["pid"]==t["pid"]==phase["pid"]==8644
            and o["unified_exec_chunk_id"]=="9835a6" and o["results_sha256"]==t["results_sha256"]==EXPECTED["actual/results.json"]
            and o["terminal_sha256"]==EXPECTED["actual/terminal.json"] and t["exit_code"]==0 and t["status"]=="FINISHED"
            and o["root_source_review_sha256"]==phase["root_source_review_sha256"]==EXPECTED["root_source_review.json"]
            and all(t[k]==v for k,v in phase.items()),"observed same scalar application")
    for key,file in (("source_sha256","apply.py"),("protocol_sha256","protocol.json"),("preflight_sha256","preflight.json")):
        require(phase[key]==gate[key]==EXPECTED[file],"same approved scalar phase")
    require(preflight["status"]=="PASS" and gate["status"]=="PASS_ROOT_ALL36_CACHED_MINIMUM_APPLICATION_SOURCE"
            and all(x is True for a in (preflight,gate) for x in a["checks"].values()) and
            result["status"]=="CANDIDATE_36_CACHED_ENERGY_LOCAL_MINIMA_AWAITING_INDEPENDENT_ACTUAL_READBACK",
            "approved exact scalar application")
    links = result["source_and_input_sha256"]
    require(links == protocol["source_and_input_sha256"] and len(links)==261, "unchanged full input closure")
    docs = {}
    for name,h in links.items():
        pp = (ROOT/name).resolve()
        require(pp.relative_to(ROOT).as_posix() == name, "canonical readback path")
        verified = raw(pp,h)
        if pp.suffix==".json":
            docs[h] = json.loads(verified)
    observer = docs[phase["observer_authority_sha256"]]
    spd4 = docs[phase["four_spd_authority_sha256"]]
    spd32 = docs[phase["remaining32_spd_authority_sha256"]]
    joined = docs[phase["joined_r34_authority_sha256"]]
    derivative4 = docs[phase["four_r34_derivative_authority_sha256"]]
    old = docs["cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3"]["targets"]
    new = docs["dfea134640f99b45e9ca75d1ce718043d67f2917375506997d2465abda8c1242"]["targets"]
    require(observer["status"]=="PASS_INDEPENDENT_ACTUAL_36_CACHED_ROOT_OBSERVERS" and observer["primary_signed_residuals"]==136
            and joined["status"]=="PASS_INDEPENDENT_REMAINING32_SAME_TARGET_BRIDGE_INPUTS"
            and spd4["status"]=="PASS_INDEPENDENT_ACTUAL_FOUR_EXACT_SYMMETRIC_PRECONDITIONER_PARTS_SPD"
            and spd32["status"]=="PASS_INDEPENDENT_ACTUAL_REMAINING32_EXACT_SYMMETRIC_PARTS", "native prior scope")
    require(len(old)==4 and len(new)==32 and len(result["cases"])==36
            and [r["pilot_id"] for r in result["cases"]]==protocol["ordered_case_ids"]
            and len(set(protocol["ordered_case_ids"]))==36, "complete 36-case scope")
    for key in ("four_spd_authority_sha256","remaining32_spd_authority_sha256"):
        proof_path = ROOT/next(name for name,h in links.items() if h==phase[key])
        receipt_name = proof_path.with_name("os_exit_observation.json").relative_to(ROOT).as_posix()
        os_gate = docs[links[receipt_name]]
        term_name = proof_path.with_name("terminal.json").relative_to(ROOT).as_posix()
        term = docs[links[term_name]]
        require(type(os_gate["observed_exit_code"]) is int and os_gate["observed_exit_code"] == term["exit_code"] == 0
                and os_gate["receipt_sha256"] == phase[key] == term["receipt_sha256"]
                and os_gate["terminal_sha256"] == links[term_name] and term["pid"] == os_gate["pid"], "separate actual SPD receipt")
    reviewed=[]
    extra_links={}
    def extra(name,h):
        path=(ROOT/name).resolve()
        require(path.relative_to(ROOT).as_posix()==name,"canonical extra archive evidence")
        require(extra_links.setdefault(name,h)==h,"immutable archive proof digest")
        return json.loads(raw(path,h))
    for row in result["cases"]:
        case = row["pilot_id"]
        b = row["proof_sha256"]
        d,g,box,n,spd = (docs[b[k]] for k in ("cached_derivative_proof_sha256","geometry_sha256","geometry_box_sha256",
                                              "root_proof_sha256","spd_proof_sha256"))
        ident=row["target_identity"]
        pilot=case in {r["pilot_id"] for r in old}
        e=select(old if pilot else new,"pilot_id" if pilot else "case_id",case)
        contract=docs[e["model_contract_sha256"]]
        score=select(observer["observer_transfer_terms"],"case_id",case)
        eta,z_native,r_native,e_native,sigma,cache,field,goal,native_domain=native(n,case)
        packet=docs[cache]
        same_norm=q(spd["exact_R_infinity_norm"])
        qspd=(max(q(spd["fresh_congruence_q_upper_exact"]),q(spd["producer_congruence_q_upper_exact"]))
              if pilot else q(spd["independent_q_infinity_upper_exact"]))
        sm=select(spd4["ordered_cases"] if pilot else spd32["case_results"],"pilot_id" if pilot else "case_id",case)
        require(sm["case_receipt_sha256" if pilot else "sha256"]==b["spd_proof_sha256"]
                and ident==spd["target_identity"]==sm["target_identity"] and q(sm["exact_R_infinity_norm"])==same_norm
                and 0<=qspd<1 and same_norm>0 and (spd["exact_target_symmetric_part_spd"] if pilot else spd["strict_congruence_bound"]) is True,
                "independent same-R congruence")
        drive=q(packet["activation_strength_exact"])
        dose,gain=q(ident["dose_exact"]),q(ident["gain_exact"])
        cell=contract["parameter_cell"]
        norm=q(packet["normalization_binary_exact"])
        require(case==ident["pilot_id"]==g["pilot_id"] and ident["record_id"]==box["record_id"]==g["record_id"]==contract["record_id"]
                and cache==box["point_cache_sha256"]==ident["cache_packet_sha256"]==contract["point_cache_sha256"]
                and field==ident["field_sha256"]==packet["input_state_sha256"]
                and box["model_observation_contract_sha256"]==e["model_contract_sha256"]==ident["model_contract_sha256"]
                and digest(json.dumps(cell,sort_keys=True,separators=(",",":")).encode())==box["parameter_cell_sha256"]==ident["parameter_cell_sha256"]
                and all(tuple(map(q,v))==(gain,gain) if k=="common_gain" else tuple(map(q,v))==(Q(packet["parameters"][k]),)*2 for k,v in cell.items())
                and q(box["dose_exact"])==dose==q(packet["optical_dose_exact"]) and gain==q(packet["common_gain_exact"])
                and drive==dose*gain==q(ident["drive_exact"])==q(n.get("exact_requested_strength",n.get("exact_strength")))
                and drive!=q(packet["activation_strength_binary_exact"]) and drive>0 and norm>0
                and norm==q(ident["energy_scale_binary_exact"])==q(contract["energy_scale_binary_exact"]), "same exact law and drive")
        rhash,uhash,X,fmap,u=cube(box,packet,ident)
        xi=q(score["archived_xi_exact"])
        require(sigma==q(score["sigma_exact"])==q(row["archived_sigma_exact"])
                and xi==q(row["archived_xi_exact"]) and xi>=0 and q(score["posterior_exact"])==e_native
                and q(score["epsilon_um_exact"])==Q(67,2)*(sigma+e_native+xi), "unchanged archived observer geometry/score")
        field_data=extra(packet["input_state"],field)
        archive_u=field_data["free_displacement"]
        require([list(map(Q,v)) for v in field_data["reference_nodes"]]==X and len(archive_u)==611, "frozen archive coordinates")
        archived_steps=dict(zip(fmap,map(Q,archive_u),strict=True))
        observed=[list(map(Q,v)) for v in field_data["current_nodes"]]
        gap=max(abs(x+archived_steps.get(2*i+j,0)-observed[i][j]) for i,ref in enumerate(X) for j,x in enumerate(ref))
        require(gap==xi and sigma==max(abs(Q(a)-Q(b)) for a,b in zip(u,archive_u,strict=True)), "exact archived nodal error and proof shift")
        if sigma:
            old_packet=extra(contract["original_point_cache"],contract["original_point_cache_sha256"])
            older=base64.b64decode(old_packet["frozen_binary_cache"]["data"],validate=True)
            a0=np.frombuffer(array_payload(older,"u",(611,)),dtype="<f8").tolist()
            a1=np.frombuffer(array_payload(base64.b64decode(packet["frozen_binary_cache"]["data"]),"u",(611,)),dtype="<f8").tolist()
            require(a0==archive_u and sigma==max(abs(Q(y)-Q(x)) for x,y in zip(a0,a1,strict=True))
                    and rhash==n["original_binary_R_sha256"]==digest(array_payload(older,"preconditioner",(611,611)))
                    and uhash==n["replayed_u1_binary_sha256"] and xi==q(n["archived_addition_gap_xi0_exact"]), "actual proof recenter")
        require(g["status"]==("PASS_FRESH_SINGLETON_WHOLE_GEOMETRY_ONLY" if pilot else "PASS_INDEPENDENT_REMAINING32_WHOLE_BOX_GEOMETRY_ONLY")
                and g["box_sha256"]==b["geometry_box_sha256"] and q(g["radius_exact"])==RADIUS
                and g["all_144_body_and_24_trace_domains"] is True and g["raw_surface_cache_geometry_lines"]==288
                and 0<q(g["strong_monotonicity_lower_exact"]) and 0<q(g["body_J_bounds_exact"][0])<=q(g["body_J_bounds_exact"][1]),
                "geometry embedding of same cube")
        mode=row["fixed_r34_cached_derivative_source"]
        d_hash=b["cached_derivative_proof_sha256"]
        if pilot:
            mark=select(derivative4["reviewed_cases"],"pilot_id",case)
            require(g["proof_center_field_sha256"]==field and sigma==0 and e["box_1_sha256"]==b["geometry_box_sha256"]
                    and spd["cached_point_proof_sha256"]==b["root_proof_sha256"]
                    and spd["larger_radius_box_sha256"]==b["geometry_box_sha256"]
                    and spd["derivative_case_sha256"]==mark["native_case_sha256"]==d_hash and d["cache_packet_sha256"]==cache
                    and d["box_sha256"]==b["geometry_box_sha256"] and d["status"]=="PASS_INDEPENDENT_LARGER_RADIUS_CACHED_DERIVATIVE"
                    and q(d["radius_exact"])==RADIUS and q(d["prior_point_eta_upper_exact"])==eta
                    and q(d["preconditioner_infinity_norm_exact"])==same_norm and q(d["proof_drive_exact"])==drive, "four older larger boxes")
            z=max(map(q,(d["reviewed_ZQ_upper_exact"],d["producer_ZQ_upper_exact"])))
            margins=domain(d["whole_box_domains"])
            require(mode=="prior_four_independently_extended_r34_cached_derivatives",case)
        else:
            j=select(joined["cases"],"case_id",case)
            require(j["geometry_replay_case_sha256"]==spd["native_geometry_case_sha256"]==b["geometry_sha256"]
                    and spd["trial_box_sha256"]==j["trial_box_sha256"]==g["box_sha256"]
                    and spd["native_cached_derivative_case_sha256"]==j["cached_derivative_case_sha256"]==d_hash
                    and contract["native_case_sha256"]==b["root_proof_sha256"] and q(j["R_inf_exact"])==same_norm
                    and g["point_cache_sha256"]==cache==j["cache_packet_sha256"]==spd["R_file_sha256"]
                    and g["center_binary_sha256"]==uhash and g["R_binary_sha256"]==rhash and q(g["exact_drive"])==drive
                    and g["archived_field_sha256"]==field and q(g["archived_shift_exact"])==sigma
                    and q(g["archived_addition_gap_exact"])==xi, "joined current singletons")
            if j["cached_input_kind"]=="NATIVE_R34":
                require(r_native==RADIUS and d_hash==b["root_proof_sha256"] and q(j["reviewed_ZQ_upper_exact"])==z_native
                        and mode=="native_independent_r34_cached_derivative","native same-radius derivative")
                z,margins=z_native,native_domain
            else:
                require(j["cached_input_kind"]=="FRESH_R34" and mode=="fresh_independent_r34_cached_derivative"
                        and d["status"]=="PASS_REVIEWED_FRESH_R34_CACHED_DERIVATIVE" and d["case_id"]==case
                        and d["cache_packet_sha256"]==cache and d["native_cached_anchor_sha256"]==b["root_proof_sha256"]
                        and d["trial_box_sha256"]==j["trial_box_sha256"] and d["R_binary_sha256"]==rhash
                        and d["center_binary_sha256"]==uhash and d["model_contract_sha256"]==ident["model_contract_sha256"]
                        and q(d["radius_exact"])==RADIUS and q(d["R_inf_exact"])==same_norm
                        and q(d["proof_drive_exact"])==drive and q(d["native_eta_upper_exact"])==eta
                        and q(d["archived_center_shift_exact"])==sigma and q(d["archived_xi_exact"])==xi,
                        "fresh native derivative on same cube and signed point anchor")
                z=max(q(d["independent_ZQ_upper_exact"]),q(d["producer_ZQ_upper_exact"]))
                require(z==q(d["reviewed_ZQ_upper_exact"])==q(j["reviewed_ZQ_upper_exact"])
                        and q(d["independent_ZQ_upper_exact"])==q(d["native_center_Z0_upper_exact"])+q(d["independent_variation_upper_exact"]),case)
                margins=domain(d["analytic_domain"])
        require(0<=z<1 and margins==row["cached_site_analytic_domain_lower_exact"], "whole cube derivative and cached C2")
        if case==EXCEPTION:
            lhs=eta+z*RADIUS
            local=eta/(1-z)
            require(r_native==4*RADIUS and e_native>=RADIUS and lhs<RADIUS and local<RADIUS
                    and q(d["diagnostic_cached_r34_selfmap_lhs_upper_exact"])==lhs
                    and q(row["auxiliary_smaller_cube_selfmap_lhs_upper_exact"])==lhs
                    and q(d["diagnostic_cached_r34_posterior_upper_exact"])==local
                    and row["same_root_localization_method"]=="same_root_via_strict_enclosed_small_cube_selfmap",
                    "same-anchor nested Banach localization is required")
        else:
            local=e_native
            lhs=None
            require(local<RADIUS and row["auxiliary_smaller_cube_selfmap_lhs_upper_exact"] is None
                    and row["same_root_localization_method"]=="native_posterior_inside_r34","native root lies in same cube interior")
        m=(1-z)/same_norm
        interior=RADIUS-local
        require(m>0 and interior>0 and sigma+local<goal and row["archived_observer_score_or_saved_field_replaced"] is False
                and row["conditional_unique_minimum_on_entire_closed_r34_cube"] is True
                and row["exact_spatial_energy_or_loading_claimed"] is False
                and all(q(row[key])==val for key,val in {
                    "native_root_radius_exact":r_native,"native_posterior_exact":e_native,"native_ZQ_upper_exact":z_native,
                    "fixed_r34_cached_ZQ_upper_exact":z,"same_center_interior_error_upper_exact":local,
                    "same_center_interior_margin_lower_exact":interior,"archived_scoring_native_posterior_exact":e_native,
                    "radius_of_geometry_and_curvature_cube_exact":RADIUS,"exact_R_infinity_norm":same_norm,
                    "symmetric_R_congruence_q_upper_exact":qspd,"normalized_cached_H_min_eigenvalue_lower_exact":m,
                    "normalized_cached_energy_gap_coefficient_exact":m/2}.items()),"independent scalar theorem and source equality")
        reviewed.append({"pilot_id":case,"proof_sha256":b,"R_binary_sha256":rhash,"center_binary_sha256":uhash,
            "derivative_kind":mode,"eta_exact":str(eta),"norm_R_exact":str(same_norm),"z_r34_exact":str(z),
            "sigma_exact":str(sigma),"archived_xi_exact":str(xi),"native_scoring_posterior_exact":str(e_native),
            "localized_posterior_exact":str(local),"r34_auxiliary_lhs_exact":None if lhs is None else str(lhs),
            "strict_cube_interior_margin_exact":str(interior),"cached_H_eigenvalue_lower_exact":str(m),"cached_energy_gap_coefficient_exact":str(m/2)})
    require(Counter(r["derivative_kind"] for r in reviewed)=={"prior_four_independently_extended_r34_cached_derivatives":4,
                 "fresh_independent_r34_cached_derivative":10,"native_independent_r34_cached_derivative":22}
            and sum(r["r34_auxiliary_lhs_exact"] is not None for r in reviewed)==1 and sum(q(r["sigma_exact"])>0 for r in reviewed)==2,
            "complete exact actual strata")
    mm=min(reviewed,key=lambda r:q(r["cached_H_eigenvalue_lower_exact"]))
    mint=min(reviewed,key=lambda r:q(r["strict_cube_interior_margin_exact"]))
    doc={"status":"PASS_INDEPENDENT_ACTUAL_ALL36_CACHED_ENERGY_LOCAL_MINIMA_SCALAR_READBACK","producer_bindings":EXPECTED,
         "independent_source_sha256":digest(Path(__file__).read_bytes()),"source_and_input_sha256":links,
         "independently_opened_archived_fields_and_recenter_packets_sha256":extra_links,
         "verified_cached_local_minima":36,"certificate_cases":reviewed,"minimum_curvature_case":mm["pilot_id"],
         "minimum_curvature_exact":mm["cached_H_eigenvalue_lower_exact"],"minimum_interior_case":mint["pilot_id"],
         "minimum_interior_margin_exact":mint["strict_cube_interior_margin_exact"],
         "new_matrix_AD_force_moment_or_solve_jobs":0,"spatial_or_continuum_minima_from_this_receipt":0,
         "archived_scores_or_fields_changed":False,"full_five_gate_goal_achieved":False}
    out=HERE/"receipt.json"
    with out.open("x",encoding="utf8",newline="\n") as f:
        json.dump(doc,f,indent=2)
        f.write("\n")
    print(json.dumps({"sha256":digest(out.read_bytes()),"cases":len(reviewed),"min_m":float(q(doc["minimum_curvature_exact"])),
                     "min_m_case":mm["pilot_id"],"min_inner":float(q(doc["minimum_interior_margin_exact"])),
                     "min_inner_case":mint["pilot_id"]}))


if __name__=="__main__":
    main()
