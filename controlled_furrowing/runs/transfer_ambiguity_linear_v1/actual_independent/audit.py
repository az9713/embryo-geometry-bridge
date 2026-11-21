"""Independent postrun reconstruction of the saved affine diagnostic."""
import csv
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu

HERE=Path(__file__).resolve().parent
OUT=HERE.parent
ROOT=HERE.parents[2]
ACT=OUT/"actual"
PRE=OUT/"preparation"
PINS={
 ACT/"results.json":"68ddc210523cb40001c99babb251916ac22c62eecd836d976fac1ab78dc0fe3f",
 ACT/"terminal.json":"d49142d9a24f4da88d610acf9979ca942639ba67a513a97245bd35d46af67a9d",
 ACT/"tangents.npz":"b01af5efb7be40684f4baee1fb50e3847fde44775af62502d7833c0b759c0743",
 PRE/"protocol.json":"ece0d458af81906c95faa231633254c3c69d55b3f0026074ccb8e9a0adeed314",
 PRE/"preflight.json":"23e1812067bf2a643ae133a1f7a05bfe1a17a93458f8cf10676ec0ae9ce29bce",
 PRE/"prelaunch_independent_review.json":"13c2289d63fd0838c04c3879de602d342283e059ee9814f45551c8707ade442e",
 ROOT/"scripts/transfer_ambiguity_linear.py":"4e89db4c6e46f219f790408a45a71deb7fba4a641e7b359d7734e09b5082cb8b",
}
OBS=["maximum_shortening_um","region_averaged_shortening_um","apical_max_shortening_um","basal_max_shortening_um"]
DIR=ROOT/"runs/material_cohort_descent_v1"
FIELD=ROOT/"runs/transfer_pair_pilot_v1/jobs/baseline_large/checkpoint_03/field.json"
SMALL=ROOT/"data/derived/calibration_small_18um.csv"
L=Q("33.5")


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def qfloat(value):
    return Q.from_float(float(value))


def dot(x,y):
    return sum((a*b for a,b in zip(x,y,strict=True)),Q(0))


def intervals_product(a,b):
    vals=[x*y for x,y in product(a,b)]
    return min(vals),max(vals)


def scaled(v,x):
    return tuple(sorted((v*x[0],v*x[1])))


def plus(a,b):
    return a[0]+b[0],a[1]+b[1]


def archived_rows():
    allrows=read(DIR/"results.json")
    protocol=read(DIR/"protocol.json")
    (item,)=[x for x in allrows["candidates"] if x["variant"]=="current_disk__spatial_max"]
    af=item["affine_LP"]
    records=protocol["records"]
    assert item["all68_constraint_order"]==[{"record_id":r["record_id"],"observable":k} for r,k in product(records,OBS)]
    with SMALL.open(encoding="utf-8-sig",newline="") as f:
        csv_rows=list(csv.DictReader(f))
    a,J=[],[]
    bindings={SMALL:sha(SMALL)}
    for sr,rec, ar in zip(csv_rows,records,allrows["records"],strict=True):
        assert sr["record_id"]==rec["record_id"]==ar["record_id"]
        p=ROOT/ar["file"]
        assert sha(p)==ar["sha256"]
        bindings[p]=ar["sha256"]
        x=read(p)
        klass=item["variant"]
        for k in OBS:
            assert Q(sr[k])==Q(rec["observed_exact"][k])
            a.append(sum(map(Q,x["baseline_endpoint_intervals"][klass][k]))/2-Q(sr[k]))
            J.append([qfloat(sum(map(Q,d["two_steps"][1]["central_difference_intervals"][klass][k]))/2)
                      for d in x["coordinate_diagnostics"]])
    assert a==list(map(Q,af["affine_residual_intercepts_exact"]))
    assert J==[list(map(Q,row)) for row in af["affine_derivative_coefficients_exact"]]
    return a,J,[tuple(map(Q,x)) for x in af["trust_box_exact"]],bindings


def lp_endpoint(item,c,r,J,delta,bounds):
    """Rebuild exact primal and weak-dual objective without the producer."""
    d=list(map(Q,item["step_exact"]))
    lam=list(map(Q,item["multipliers_exact"]))
    assert len(d)==len(bounds)==5 and len(lam)==2*len(J)==136
    assert all(lo<=x<=hi for x,(lo,hi) in zip(d,bounds,strict=True)) and d[2]>=0
    assert all(x>=0 for x in lam) and 0<=Q(item["repair_alpha_exact"])<=1
    jd=[dot(row,d) for row in J]
    assert all(abs(a+x)<=delta for a,x in zip(r,jd,strict=True))
    left,right=lam[:len(J)],lam[len(J):]
    inner=[a-b for a,b in zip(left,right,strict=True)]
    bmass=delta*sum(lam)-dot(inner,r)
    cg=[c[j]-dot(inner,[row[j] for row in J]) for j in range(5)]
    upper=bmass+sum(max(k*lo,k*hi) for k,(lo,hi) in zip(cg,bounds,strict=True))
    lower=dot(c,d)
    assert Q(item["primal_lower_exact"])==lower>=0
    assert Q(item["dual_upper_exact"])==upper>=lower
    assert Q(item["duality_gap_exact"])==upper-lower
    residuals=[a+x for a,x in zip(r,jd,strict=True)]
    return lower,upper,d,max(map(abs,residuals)),max(map(abs,jd))


def p2controls(state):
    rp=state["reference_nodes"]
    incidence={}
    for p,q,r,pq,qr,rp_mid in state["triangles"]:
        for x,y,m in ((p,q,pq),(q,r,qr),(r,p,rp_mid)):
            incidence.setdefault(tuple(sorted((x,y))),[]).append(m)
    topology=[]
    for (x,y),m in incidence.items():
        if len(m)==1 and rp[x][1]==rp[y][1]==state["parameters"]["height"]:
            if rp[x][0]>rp[y][0]:
                x,y=y,x
            topology.append((x,m[0],y))
    topology.sort(key=lambda ids:rp[ids[0]][0])
    return topology


def basis(t):
    return (2*(t-Q(1,2))*(t-1),4*t*(1-t),2*t*(t-Q(1,2)))


def derivative(t):
    return (4*t-3,4-8*t,4*t-1)


def evaluate(nodes,edge,t,k):
    return dot(basis(t),[nodes[i][k] for i in edge])


def divide(v,w):
    assert w[0]>0
    return intervals_product(v,(1/w[1],1/w[0]))


def coordinate_derivative(nodes,edge,t,k):
    return dot(derivative(t),[nodes[i][k] for i in edge])


def q_probe(nodes,topology,shoulder,chosen):
    assert nodes[topology[0][0]][0]==0
    assert all(coordinate_derivative(nodes,e,t,0)>0 for e,t in product(topology,(0,1)))
    assert [i for i,e in enumerate(topology) if nodes[e[0]][0]<shoulder<nodes[e[-1]][0]]==[chosen]
    edge=topology[chosen]
    lo,hi=Q(0),Q(1)
    for _ in range(112):
        mid=(lo+hi)/2
        val=evaluate(nodes,edge,mid,0)
        if val==shoulder:
            lo=hi=mid
            break
        if val<shoulder:
            lo=mid
        else:
            hi=mid
    # The endpoint height bracket is exact if Z_t does not switch sign.
    dz=[coordinate_derivative(nodes,edge,t,1) for t in (lo,hi)]
    assert min(dz)>0 or max(dz)<0
    heights=sorted(evaluate(nodes,edge,t,1) for t in (lo,hi))
    center=nodes[topology[0][0]][1]
    return (lo,hi),(center-heights[1],center-heights[0])


def independent_scalar_operators(nodes,topology,shoulder,chosen,free,S):
    edge=topology[chosen]
    axis=topology[0][0]
    itr,q0=q_probe(nodes,topology,shoulder,chosen)
    slopes=[]
    for k in (0,1):
        slopes.append(tuple(sorted(coordinate_derivative(nodes,edge,t,k) for t in itr)))
    ratio=divide(slopes[1],slopes[0])
    weight=[]
    for j in range(3):
        values=[basis(t)[j] for t in itr]
        # Derivatives of the quadratic weights can vanish at these exact knots.
        t0=(Q(3,4),Q(1,2),Q(1,4))[j]
        if itr[0]<=t0<=itr[1]:
            values.append(basis(t0)[j])
        weight.append((min(values),max(values)))
    d={}
    for node,w in zip(edge,weight,strict=True):
        d[2*node]=intervals_product(ratio,w)
        d[2*node+1]=(-w[1],-w[0])
    d[2*axis+1]=plus(d.get(2*axis+1,(0,0)),(Q(1),Q(1)))
    dx=[scaled(L,d.get(int(k),(Q(0),Q(0)))) for k in free]
    coefficients=[]
    for j in range(1,6):
        rm=(Q(0),Q(0))
        for x,v in zip(dx,S[:,j],strict=True):
            rm=plus(rm,scaled(qfloat(v),x))
        coefficients.append(rm)
    return dx,coefficients,q0


def probe_directions(nodes,topology,shoulder,chosen,free,S,c,selected):
    directions=[("coordinate_"+str(i),np.eye(5)[i]) for i in range(5)]
    directions.extend(("epsilon_1e-3_"+k,np.array([float(x) for x in step])) for k,step in selected)
    reports=[]
    for name,direction in directions:
        # Pure apical T and any feasible mixed positive-Ta direction use a
        # one-sided step. Other directions preserve Ta=0 for both signs.
        assert direction[2]>=0
        central=direction[2]==0
        nodal=S[:,1:]@direction
        flat=[Q(0)]*(2*len(nodes))
        for index,value in zip(free,nodal,strict=True):
            flat[int(index)]=qfloat(value)*L
        prediction=float(c@direction)
        observed=[]
        for h in (Q(1,512),Q(1,2048),Q(1,8192)):
            y=[]
            for sign in ((-1,1) if central else (0,1)):
                moved=[tuple(x+sign*h*flat[2*i+j] for j,x in enumerate(row)) for i,row in enumerate(nodes)]
                assert moved[topology[0][0]][0]==0
                y.append(q_probe(moved,topology,shoulder,chosen)[1])
            denom=(2 if central else 1)*h
            interval=((y[1][0]-y[0][1])/denom,(y[1][1]-y[0][0])/denom)
            estimate=float(sum(interval)/2)
            observed.append({"step":str(h),"secant":estimate,"absolute_difference":abs(estimate-prediction)})
        assert observed[-1]["absolute_difference"]<max(1e-4,abs(prediction)*1e-3)
        reports.append({"direction":name,"prediction":prediction,"central":bool(central),"finite_differences":observed})
    return reports


def main():
    assert all(sha(p)==h for p,h in PINS.items())
    start,term,result=[read(ACT/f) for f in ("started.json","terminal.json","results.json")]
    pre,p,review=[read(PRE/f) for f in ("preflight.json","protocol.json","prelaunch_independent_review.json")]
    assert start["pid"]==term["pid"]==9264 and type(term["exit_code"]) is int and term["exit_code"]==0
    assert term["status"]=="FINISHED" and term["results_sha256"]==PINS[ACT/"results.json"]
    assert result["prelaunch_review_sha256"]==start["review_sha256"]==PINS[PRE/"prelaunch_independent_review.json"]
    assert result["protocol_sha256"]==PINS[PRE/"protocol.json"] and pre["protocol_sha256"]==result["protocol_sha256"]
    assert review["preflight_sha256"]==PINS[PRE/"preflight.json"]
    assert result["source_sha256"]==PINS[ROOT/"scripts/transfer_ambiguity_linear.py"]
    assert result["tangents_sha256"]==PINS[ACT/"tangents.npz"]
    assert result["heldout_response_tables_read"] is result["new_equilibrium_solves"] is False
    assert p["input_hashes"]==result["input_hashes"]
    assert all(sha(ROOT/k)==v for k,v in p["input_hashes"].items())
    residual,J,bounds,binding=archived_rows()
    assert bounds==[tuple(map(Q,x)) for x in p["trust_box_exact"]]
    assert bounds[2][0]==0 and bounds[2][1]>0
    c=list(map(Q,result["rounded_affine_contrast_gradient_exact"]))
    exported=[tuple(map(Q,x)) for x in result["contrast_gradient_intervals_exact"]]
    assert all(qfloat((a+b)/2)==v for (a,b),v in zip(exported,c,strict=True))
    base=max(map(abs,residual))
    assert Q(result["source_baseline_residual_um_exact"])==base
    expected=[(key,v) for key,bs in (("source_residual_allowance",p["source_discrepancy_offsets_um_exact"]),
                  ("small_prediction_change",p["small_prediction_change_budgets_um_exact"])) for v in bs]
    assert len(expected)==10 and [(r["mode"],r["budget_exact"]) for r in result["rows"]]==expected
    rows,all_gaps,selected=[],[],[]
    for item in result["rows"]:
        mode,budget=item["mode"],Q(item["budget_exact"])
        r,delta=(residual,base+budget) if mode=="source_residual_allowance" else ([Q(0)]*68,budget)
        assert Q(item["allowance_um_exact"])==delta
        on=[]
        for sign,key in ((1,"max_contrast_change"),(-1,"max_negative_contrast_change")):
            q=item[key]
            lower,upper,d,norm,change=lp_endpoint(q,[sign*x for x in c],r,J,delta,bounds)
            all_gaps.append(upper-lower)
            on.append((lower,upper))
            if mode=="small_prediction_change" and budget in (Q(1,1000),Q(1,100)):
                name="maximize_Q" if sign==1 else "minimize_Q"
                if budget==Q(1,1000):
                    selected.append((name,d))
            else:
                name="maximize_Q" if sign==1 else "minimize_Q"
            rows.append({"mode":mode,"budget_exact":str(budget),"objective":name,
                         "signed_objective_exact":[str(lower),str(upper)],"step_exact":list(map(str,d)),
                         "max_calibration_residual_rounded_affine_um":float(max(abs(a+dot(b,d)) for a,b in zip(residual,J,strict=True))),
                         "max_predicted_small_change_um":float(change)})
        width=[sum(x[i] for x in on) for i in (0,1)]
        assert width==list(map(Q,item["contrast_range_width_um_exact"]))
        assert [float(x) for x in width]==item["contrast_range_width_um_display"]
    assert len(rows)==20
    state=read(FIELD)
    with np.load(ACT/"tangents.npz",allow_pickle=False) as archive:
        z={k:archive[k] for k in archive.files}
    free,S,B=z["free"],z["S"],z["B"]
    H=csr_matrix((z["H_data"],z["H_indices"],z["H_indptr"]),shape=tuple(z["H_shape"]))
    assert H.shape==(611,611) and S.shape==B.shape==(611,6) and len(free)==611
    assert np.array_equal(z["u"],state["free_displacement"]) and np.unique(free).size==611
    assert all(np.isfinite(a).all() for a in (S,B,H.data,z["u"]))
    skew=H-H.T
    assert skew.nnz==0 or max(abs(skew.data))<1e-12
    residual_norm=np.linalg.norm(H@S+B,axis=0)/np.linalg.norm(B,axis=0)
    assert max(residual_norm)<1e-9
    topo=p2controls(state)
    nodes=[tuple(qfloat(v)*L for v in row) for row in state["current_nodes"]]
    chosen=result["shoulder"]["edge_index"]
    assert chosen==3 and len(topo)==12
    assert 2*topo[0][0] not in set(free) and 2*topo[0][0]+1 in set(free)
    # Physical axis/rim constraints rebuilt from reference geometry.
    expected_free=[]
    for i,(radial,height) in enumerate(state["reference_nodes"]):
        if radial==state["parameters"]["outer_radius"]:
            continue
        if radial!=0:
            expected_free.append(2*i)
        expected_free.append(2*i+1)
    assert expected_free==free.tolist()
    shoulder=Q(p["shoulder_um_exact"])
    dx,coeff,q0=independent_scalar_operators(nodes,topo,shoulder,chosen,free,S)
    assert all(a<=c<=d<=b for (a,b),(c,d) in zip(exported,coeff,strict=True))
    pQ=tuple(map(Q,result["shoulder"]["Q_um_exact"]))
    assert pQ[0]<=q0[0]<=q0[1]<=pQ[1]
    g=np.array([float(sum(x)/2) for x in dx])
    tangent=g@S
    factor=splu(H.transpose().tocsc())
    lam=factor.solve(g)
    adjoint=-lam@B
    adjoint_residual=np.linalg.norm(H.T@lam-g)/np.linalg.norm(g)
    assert adjoint_residual<=1e-9 and np.max(abs(adjoint-tangent))<=1e-8
    assert np.max(abs(tangent[1:]-np.array(list(map(float,c)))))<=1e-10
    fds=probe_directions(nodes,topo,shoulder,chosen,free,S,np.asarray(list(map(float,c))),selected)
    more={**{ROOT/k:v for k,v in p["input_hashes"].items()},**PINS,**binding,
          ACT/"started.json":sha(ACT/"started.json"),Path(__file__):sha(Path(__file__))}
    assert all(sha(k)==v for k,v in more.items())
    report={
        "status":"PASS_INDEPENDENT_ACTUAL_LOCAL_AFFINE_DIAGNOSTIC","checks":{
            "actual_terminal_and_all_bindings":True,"reconstructed_all_68_rows":True,
            "reconstructed_all_20_LP_primal_dual_and_widths":True,"exact_Ta_one_sided_support_and_inequalities":True,
            "recovered_saved_free_map_HS_B_symmetry_and_adjoint":True,
            "independent_current_radius_derivative_and_feasible_graph_finite_differences":True},
        "scope":"Saved-field floating mechanical linearization and exact rounded-affine optimizations only.",
        "bindings":{k.relative_to(ROOT).as_posix():v for k,v in sorted(more.items())},
        "max_exact_duality_gap_um_display":float(max(all_gaps)),
        "LP_endpoints":rows,"Hs_plus_B_relative_column_residuals":residual_norm.tolist(),
        "adjoint_relative_residual":float(adjoint_residual),"max_adjoint_tangent_absolute_discrepancy":float(np.max(abs(adjoint-tangent))),
        "independent_gradient_display":tangent[1:].tolist(),"finite_difference_rows":fds,
        "source_compatible_witness":False,"validated_equilibrium_or_PDE_claim":False,"novelty_established":False,
        "larger_response_files_read":False,"new_force_or_equilibrium_computation":False,
    }
    o=HERE/"receipt.json"
    with o.open("x",encoding="utf-8") as f:
        json.dump(report,f,indent=2)
        f.write("\n")
    print(sha(o),float(max(all_gaps)),residual_norm.max(),adjoint_residual)
    print("finest graph secant differences",[x["finite_differences"][-1]["absolute_difference"] for x in fds])


if __name__=="__main__":
    main()

