"""Independent admission checks for the four singleton spatial-FE targets."""
from collections import Counter
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np

HERE = Path(__file__).resolve().parent
PRE = HERE.parent
PRIMARY = PRE.parents[2]
ROOT = PRIMARY.parents[1]
PROOFS = PRIMARY/"certification_design/cached_pilots_v1/independent_actual_preparation/actual_replay"
SRC = ROOT/"scripts/transfer_ambiguity_exact_targets.py"
PINS = {
    SRC:"53ff92a87d164871e04323759b27e015fc8e19d039dc0dd15c2d6080acda7cc8",
    PRE/"protocol.json":"28efd2467bbeae4da32797e0405b9a9595cdbb03ded9ce94221f18592305687c",
    PRE/"target_manifest.json":"cb25fbeaef753cea700c187cd8b377b7dcf4aa1fe2b0b246ed17134a07a89bc3",
    PRE/"preflight.json":"a4dacbab0dafb730db12321892d741979fea4654c436996b964cc9cfde1b0808",
    PROOFS/"receipt.json":"4e891afd79dd2f73508d1f878ddc1113ef107bdff8d9e0f673894f4d2f10ae17",
}
RADII=(Q(1,2**36),Q(1,2**34))
MAX_R=RADII[1]
ORDER=("basal_area_modulus","apical_area_modulus","basal_surface_tension","apical_surface_tension","common_gain")
CHARTS=("off_axis","axis_vertex","axis_edge_physical")
SCOPE=("Only the four frozen singleton target inputs and their eight proposed spatial geometry boxes. "
       "Cached-root force/derivative premises apply only to the separate 2^-36 cached quadrature "
       "root proofs. This gate authorizes geometry alone.")


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def save(path, obj):
    with path.open("x",encoding="utf-8",newline="\n") as out:
        json.dump(obj,out,indent=2)
        out.write("\n")
    return sha(path)


def canon(obj):
    return sha256(json.dumps(obj,sort_keys=True,separators=(",",":")).encode()).hexdigest()


def qarr(arr):
    return [Q(x) for x in arr]


def identity(ctx, row, contract, predinfo, proof, entry, arrays):
    p=arrays["parameters"]
    assert tuple(row["parameter_order"]) == ORDER == tuple(row["parameter_cell"])
    drive=contract["exact_requested_drive"]
    state=read(ROOT/contract["source_field"])
    gain=Q(drive["common_gain_exact"])
    dose=Q(drive["optical_dose_exact"])
    requested=gain*dose
    rounded=Q(state["activation_strength"])
    assert gain>0 and dose>0 and rounded==Q(float(requested))!=requested
    assert Q(drive["activation_strength_exact"])==Q(proof["exact_strength"])==requested
    assert Q(drive["activation_strength_binary_exact"])==rounded
    assert Q(drive["binary_input_minus_exact_strength"])==rounded-requested
    assert requested.denominator & (requested.denominator - 1)
    assert p==state["parameters"]
    expected=dict(zip(ORDER,[Q(p[k]) for k in ORDER[:-1]]+[gain],strict=True))
    assert ctx.base==expected and ctx.dose==dose and ctx.drive==requested
    assert ctx.cell=={k:(v,v) for k,v in expected.items()}
    cell={k:[str(v),str(v)] for k,v in expected.items()}
    assert row["parameter_cell"]==contract["parameter_cell"]==cell
    assert row["parameter_cell_sha256"]==contract["parameter_cell_sha256"]==canon(cell)
    assert {k:Q(v) for k,v in row["base_parameters_exact"].items()}==expected
    assert row["record_id"]==ctx.record_id==contract["record_id"]==f"transfer_exact_{ctx.pilot_id}"
    assert row["source_record_id"]==contract["physical_source_record_id"]==proof["record_id"]
    assert ctx.record_id!=proof["record_id"]
    assert row["point_cache_sha256"]==contract["point_cache_sha256"]==proof["cache_packet_sha256"]
    assert row["frozen_binary_cache_sha256"]==contract["frozen_binary_cache_sha256"]==proof["cache_array_contents_sha256"]
    assert sha256(ctx.R.tobytes(order="C")).hexdigest()==contract["binary_preconditioner_sha256"]
    assert ctx.R.dtype==np.dtype("<f8") and np.isfinite(ctx.R).all() and not ctx.R.flags.writeable
    zeros=bytes(8*611*5)
    assert predinfo["little_endian_zero_binary64_sha256"]==sha256(zeros).hexdigest()
    assert Q(contract["cached_root_radius_exact"])==Q(proof["radius_exact"])==RADII[0]
    assert Q(contract["cached_root_error_exact"])==Q(proof["independent_posterior_error_upper_exact"])<RADII[0]
    assert MAX_R<Q(contract["observer_design_goal_exact"])
    assert contract["larger_radius_cached_derivative_asserted"] is contract["exact_spatial_root_asserted"] is False
    assert predinfo["shape"]==[611,5] and predinfo["identically_zero"] is True
    assert not ctx.S.flags.writeable and ctx.S.tobytes()==zeros
    assert predinfo["parameter_cell_sha256"]==canon(cell)
    assert row["predictor_sha256"]==entry["predictor_sha256"]
    assert row["model_observation_contract_sha256"]==entry["model_contract_sha256"]
    assert Q(contract["energy_scale_binary_exact"])==Q(state["reference_volume"])==Q(float(ctx.target.a["energy_scale"][0]))>0
    return state


def boxes_and_activation(ctx, state, docs):
    a=ctx.target.a
    ref=[qarr(row) for row in a["reference"].tolist()]
    x=[p.copy() for p in ref]
    p=state["parameters"]
    rim,height=Q(p["outer_radius"]),Q(p["height"])
    free=[2*i+k for i,(r,_) in enumerate(ref) for k in (0,1) if r!=rim and (r!=0 or k==1)]
    assert a["free"].tolist()==free and len(free)==611 and len(ref)==325
    assert x==[qarr(row) for row in state["reference_nodes"]]
    for i,v in zip(free,map(Q,state["free_displacement"]),strict=True):
        n,c=divmod(i,2)
        x[n][c]+=v
    stored=[qarr(row) for row in state["current_nodes"]]
    assert all(float(v)==float(w) for a,b in zip(x,stored,strict=True) for v,w in zip(a,b,strict=True))
    gap=max(abs(v-w) for a,b in zip(x,stored,strict=True) for v,w in zip(a,b,strict=True))
    tri=a["triangles"].tolist()
    assert tri==state["triangles"] and len(tri)==144
    patch,lo,hi=(Q(p[k]) for k in ("patch_radius","activation_z_min","activation_z_max"))
    active=[]
    for e in tri:
        corners=[ref[i] for i in e[:3]]
        r,z=(sum(w)/3 for w in zip(*corners,strict=True))
        selected=r<patch and lo<z<hi
        if selected:
            assert all(y<=patch and lo<=v<=hi for y,v in corners)
        active.append(selected)
    assert active==state["activation_elements"]==a["active"][:,0].tolist() and sum(active)==24
    keys=(*ctx.row.keys(),"precision_bits")
    authoritative={"record_id":ctx.record_id,"dose_exact":str(ctx.dose),"precision_bits":640}
    authoritative.update({k:v for k,v in ctx.row.items() if k.endswith("_sha256") and
                          k!="frozen_binary_cache_sha256"})
    for i,r in enumerate(RADII):
        d=docs[f"box_{i}"]
        assert all(d[k]==v for k,v in authoritative.items())
        assert d["root_membership_assumed"] is d["any_spatial_derivative_bound_assumed"] is False
        assert Q(d["corrector_radius_exact"])==r and d["free"]==free and d["triangles"]==tri
        assert [qarr(w) for w in d["reference_nodes_exact"]]==ref
        assert [qarr(w) for w in d["current_centers_exact"]]==x
        assert len(d["predictor_hulls_exact"])==611 and all(qarr(w)==[0,0] for w in d["predictor_hulls_exact"])
        lookup=set(free)
        caps=[[r if 2*j+k in lookup else Q(0) for k in (0,1)] for j in range(325)]
        assert [qarr(w) for w in d["coordinate_radii_exact"]]==caps
        assert [[qarr(w) for w in b] for b in d["coordinate_hulls_exact"]]==[
            [[v-e,v+e] for v,e in zip(now,cap,strict=True)] for now,cap in zip(x,caps,strict=True)]
    return ref,tri,active,{val:i for i,val in enumerate(free)},gap


def check_rules(ctx,p,route,rule, ref, tri, mask, free):
    a=ctx.target.a
    scale=ctx.target.energy_scale
    b=[]
    for v in a["shape"].tolist():
        _,n1,n2,n3,n4,n5=qarr(v)
        s=n1+(n3+n4)/2
        t=n2+(n4+n5)/2
        b.append([str(1-s-t),str(s),str(t)])
        assert 0<s and 0<t and s+t<1
    assert rule["volume_barycentric"]==b and len(b)==64
    vweights=[[str(w/scale) for w in qarr(row)] for row in a["volume_weights"].tolist()]
    assert all(Q(w)>0 for row in vweights for w in row)
    assert vweights==rule["volume_weights"]
    def indices(ids):
        dofs=[2*i+k for i in ids for k in (0,1)]
        loc=[i for i,n in enumerate(dofs) if n in free]
        return loc,[free[dofs[i]] for i in loc]
    counts=Counter()
    for i,(elem, ids) in enumerate(zip(route["bulk"],tri,strict=True)):
        klass=CHARTS[sum(ref[j][0]==0 for j in ids[:3])]
        counts[klass]+=1
        assert elem["element"]==i and elem["nodes"]==ids and elem["chart"]==klass
        assert elem["active"]==str(int(mask[i]))
        assert (elem["free_local"],elem["global_indices"])==indices(ids)
        assert elem["auxiliary_element_rule_sha256"]==canon([b,vweights[i]])

    assert dict(counts)==route["bulk_chart_counts"]=={"off_axis":132,"axis_vertex":6,"axis_edge_physical":6}
    mesh={}
    for i,j,k,l,m,n in tri:
        for x,y,z in ((i,j,l),(j,k,m),(k,i,n)):
            k=tuple(sorted((x,y),key=lambda h:ref[h][0]))
            mesh.setdefault(frozenset((x,y)),[]).append((k[0],z,k[1]))

    for face in (0,1):
        s=[Q(v[1])/2+Q(v[2]) for v in a[f"surface_{face}_N"].tolist()]
        assert len(s)==12 and s==qarr(rule["surface_s"][face]) and all(0<x<1 for x in s)
        traces=sorted((v[0] for v in mesh.values() if len(v)==1 and
                       all(ref[k][1]==face for k in v[0])),key=lambda v:ref[v[0]][0])
        assert len(traces)==12 and list(map(list,traces))==a[f"surface_{face}_edges"].tolist()
        for i,ids in enumerate(traces):
            w=[x/scale for x in qarr(a[f"surface_{face}_weights"][i].tolist())]
            assert w==qarr(rule["surface_weights"][face][i]) and all(x>0 for x in w)
            e=route["surface"][12*face+i]
            assert e["face"]==face and e["edge"]==i and e["nodes"]==list(ids)
            assert (e["free_local"],e["global_indices"])==indices(ids)
            assert Q(e["Rlo_exact"])==ref[ids[0]][0] and Q(e["Rhi_exact"])==ref[ids[-1]][0]
            assert e["trace_rule_sha256"]==canon([list(map(str,s)),list(map(str,w))])


def rejects(action):
    try:
        action()
    except ValueError:
        return True
    return False


def main():
    assert all(sha(p)==h for p,h in PINS.items())
    protocol, pre, man, cached=map(read,(PRE/"protocol.json",PRE/"preflight.json",PRE/"target_manifest.json",PROOFS/"receipt.json"))
    assert all(pre[k]==v for k,v in (
        ("source_sha256",PINS[SRC]),("protocol_sha256",PINS[PRE/"protocol.json"]),
        ("target_manifest_sha256",PINS[PRE/"target_manifest.json"])))
    assert pre["status"]=="PASS" and pre["checks"] and all(v is True for v in pre["checks"].values())
    assert pre["actual_geometry_jobs"]==pre["new_spatial_integrals"]==pre["new_cached_AD_assemblies"]==0
    assert not (PRE.parent/"geometry").exists()
    assert len(protocol["input_hashes"])==130
    assert all(sha(ROOT/k)==v for k,v in protocol["input_hashes"].items())
    assert tuple(map(Q,protocol["radii_exact"]))==tuple(map(Q,man["radii_exact"]))==RADII

    sys.path.insert(0,str(ROOT/"scripts"))
    import transfer_ambiguity_exact_targets as m
    def forbidden(*_a,**_k):
        raise AssertionError("numerical geometry, site AD or exact integral work prohibited")
    assert PINS[SRC]==sha(PRE/"transfer_ambiguity_exact_targets_frozen.py")
    reviewed_pilots=[]
    roots={}
    with (m.no_integrals(),patch.object(m.G,"certify_box",side_effect=forbidden),
          patch.object(m.E.SUR,"uniform_G_bounds",side_effect=forbidden),
          patch.object(m.E.SUR,"local_site_domain",side_effect=forbidden)):
        pro, fresh, data=m.prepared()
        assert pro==protocol and fresh==man
        for entry, pc in zip(man["targets"],cached["field_results"],strict=True):
            pid=entry["pilot_id"]
            assert pid==pc["pilot_id"]
            names=("model_contract","predictor","row","rule","routes","box_0","box_1")
            docs={k:m.pinned(ROOT/entry[k],entry[k+"_sha256"]) for k in names}
            ctx,rule,_=m.open_context(pid,man,data)
            assert ctx.pilot_id==pid
            assert ctx.target.packet["integration_policy"]["volume_duffy_gauss_order"]==8
            assert ctx.target.packet["integration_policy"]["surface_gauss_order"]==12
            assert docs["routes"]==m.routes(ctx,rule)
            contract=docs["model_contract"]
            proof=m.pinned(PROOFS/f"{pid}.json",pc["receipt_sha256"])
            state=identity(ctx,docs["row"],contract,docs["predictor"],proof,entry,ctx.target.packet)
            ref,tri,active,free,gap=boxes_and_activation(ctx,state,docs)
            assert gap==Q(data[pid][0]["xi_dimensionless_exact"])
            check_rules(ctx,state["parameters"],docs["routes"],docs["rule"],ref,tri,active,free)
            assert all(docs[f"box_{i}"]==m.exact_node_box(ctx,r) for i,r in enumerate(RADII))
            assert docs["routes"]["identity"]==m.E._context_identity(ctx)
            assert rejects(lambda:m.check_context(replace(ctx,record_id=proof["record_id"])))
            assert rejects(lambda:m.check_context(replace(ctx,dose=ctx.dose+1)))
            assert rejects(lambda:m.check_context(replace(ctx,drive=Q(contract["exact_requested_drive"]["activation_strength_binary_exact"]))))
            assert rejects(lambda:m.check_context(replace(ctx,cell=dict(ctx.cell,common_gain=(ctx.base["common_gain"],ctx.base["common_gain"]+1)))))
            bad=deepcopy(rule)
            bad["surface_s"][0][0]+=Q(1,1024)
            assert rejects(lambda:m.E._surface_rule(ctx.target,bad,0,0))
            bad=deepcopy(rule)
            bad["volume_weights"][0][0]+=Q(1,1024)
            assert rejects(lambda:m.E._bulk_rule(ctx.target,bad,0))
            assert rejects(lambda:m.exact_node_box(ctx,Q(1,2**35)))
            copied=deepcopy(entry)
            copied["predictor_sha256"]="1"*64
            assert rejects(lambda:m.open_context(pid,dict(man,targets=[copied]),data))
            roots[pid]=ctx
            reviewed_pilots.append({
                "pilot_id":pid,"record_id":ctx.record_id,"source_field_sha256":contract["source_field_sha256"],
                "cache_sha256":contract["frozen_binary_cache_sha256"],"R_binary_sha256":contract["binary_preconditioner_sha256"],
                "parameter_cell_sha256":canon(contract["parameter_cell"]),"exact_binary_gap":str(gap),
                "geometry_box_sha256":[entry["box_0_sha256"],entry["box_1_sha256"]],
                "remaining_original_cube_proofs_reused":False})
        for pid,ctx in roots.items():
            other=next(k for k in roots if k!=pid)
            assert rejects(lambda:m.check_context(replace(ctx,target=roots[other].target)))
            first=data[pid]
            foreign=data[other][1]
            # Even with its own valid digest, the foreign packet fails the target identity.
            assert rejects(lambda:m.objects(pid,{pid:(first[0],foreign,*first[2:])}))
            altered=dict(ctx.target.packet)
            altered["integration_policy"]=dict(altered["integration_policy"],surface_gauss_order=16)
            original=m.pinned
            def changed(f,h):
                return altered if f==ROOT/first[1]["cache_packet"] else original(f,h)
            with patch.object(m,"pinned",side_effect=changed):
                assert rejects(lambda:m.objects(pid,data))

    assert not (PRE.parent/"geometry").exists()
    checks={
        "all_130_dependencies_and_fresh_cached_root_authority_match":True,
        "all_four_exact_requested_drive_and_frozen_current_cache_preconditioners":True,
        "all_four_exact_X_plus_u_and_independent_611_free_144_active_and_24_trace_maps":True,
        "zero_S_and_degenerate_common_gain_exact_cells_digest_verified":True,
        "all_eight_650_coordinate_hulls_and_radii_independently_recomputed":True,
        "exact_frozen_Q8_S12_auxiliary_rules_with_no_integration":True,
        "foreign_record_cache_drive_surface_order_rule_proof_and_domain_rejected":True,
        "source_only_geometry_guard_and_no_new_jobs":True,
    }
    support=HERE/"support.json"
    h=save(support,{"status":"PASS","checks":checks,"scope":SCOPE,"audit_source_sha256":sha(Path(__file__)),
        "source_pins":{f.relative_to(ROOT).as_posix():v for f,v in PINS.items()},
        "fresh_target_count":4,"independently_recomputed_box_count":8,
        "new_selector_tail_force_hessian_or_cached_site_work":0,"actual_geometry_jobs":0,
        "reuse_notes":"Fresh geometry will reuse the reviewed lower geometry and physical-axis kernels. "
        "This review rebuilt exact identities, topology and geometry-free quadrature rules without "
        "certifying any whole-body or trace domain. Public open_context validates serialized cache "
        "identity and freezes arrays. The lower internal helper _context_identity is a metadata check "
        "and by itself is insufficient authority for an altered in-memory target.",
        "targets":reviewed_pilots})
    binding_names=("source_sha256","protocol_sha256","target_manifest_sha256","preflight_sha256",
                   "cached_root_independent_review_sha256")
    gate=PRE/"prelaunch_independent_review.json"
    bindings=dict(zip(binding_names,PINS.values(),strict=True))
    gh=save(gate,{"status":"PASS","checks":checks,**bindings,
                  "support_sha256":h,"support_path":support.relative_to(ROOT).as_posix(),"scope":SCOPE})
    print(gh,h,flush=True)


if __name__=="__main__":
    main()
