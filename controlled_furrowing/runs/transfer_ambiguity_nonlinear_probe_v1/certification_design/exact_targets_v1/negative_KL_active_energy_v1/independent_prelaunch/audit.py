"""Independent cold tests of exact active-energy and two matched proof-center authorities."""
from fractions import Fraction as Q
from itertools import permutations, product
from pathlib import Path
import json
import sys
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import active_energy as A

SOURCE="072ed83a30565be34bd5fe7dde559526e089a070df0ad95f0e61611db9e9e76d"
PROTOCOL="06dd35434358ad7dcb993349281ed6bf8a84fd0a34f7ca2d4932ec0792842361"
PREFLIGHT="ae11e52676f1e7ecacc5b925de47a553024dec6604d4302ade9ddd2bad43738c"
EDGES=((0,1),(1,2),(2,0))


def singleton(x):
    x=Q(x)
    return (x,x)


def tensor(vertices,f):
    points=vertices+[tuple((a+b)/2 for a,b in zip(vertices[i],vertices[j],strict=True)) for i,j in EDGES]
    return [list(map(singleton,f(r,z))) for r,z in points]


def test_elements():
    oblique=list(map(lambda p:tuple(map(Q,p)),[(2,3),(3,3),(4,6)]))
    # Integral of R is det (2+3+4)/6=9/2; derivative (2,3).
    affine=tensor(oblique,lambda r,z:(r+2*z,3*z))
    assert A.element(oblique,affine)==(Q(117,4),)*2
    # R=2+s+2t, Z=3+3t. The weighted squared derivative is
    # 2s²+8st+80t²+s³+6s²t+48st²+80t³ times det/2.
    quadratic=tensor(oblique,lambda r,z:((z-3)**2,(r-2)*(z-3)))
    assert A.element(oblique,quadratic)==(Q(727,40),)*2
    on_axis=list(map(lambda p:tuple(map(Q,p)),[(0,0),(2,0),(0,3)]))
    # R=2s, Z=3t; weighted squared derivative with det/2 integrates
    # as 216*integral(s t²)+24*integral(s³)=24/5.
    for v in permutations(on_axis):
        v=list(v)
        assert A.element(v,tensor(v,lambda r,z:(z*z,r*z)))==(Q(24,5),)*2
    # Only two midside coordinates are uncertain, with signed crosses.
    lower=[[q[0] for q in row] for row in quadratic]
    box=[list(row) for row in quadratic]
    coords=((3,0),(5,1))
    interval=((Q(-1,5),Q(2,7)),(Q(-3,7),Q(5,11)))
    for (i,j),h in zip(coords,interval,strict=True):
        v=lower[i][j]
        box[i][j]=(v+h[0],v+h[1])
    enclosure=A.element(oblique,box)
    for u in product(*[(a,(a+b)/2,b) for a,b in interval]):
        p=[row.copy() for row in quadratic]
        for (i,j),d in zip(coords,u,strict=True):
            p[i][j]=singleton(lower[i][j]+d)
        m,M=A.element(oblique,p)
        assert enclosure[0]<=m==M<=enclosure[1]
    assert A.classification((1,2),(2,3))[0]=="UNRESOLVED_NECESSARY_ENERGY_ORDER"
    assert A.classification((2,3),(1,2))[0]=="NECESSARY_ENERGY_ORDER_VERIFIED_ONLY"
    assert A.classification((1,2),(Q(9,4),3))[0]=="EXCLUDES_SMOOTH_STABLE_SAME_OPERATOR_CONNECTION"


def verify_scope(cases,pins,p):
    K,L=cases
    assert tuple(x["case_id"] for x in cases)==(A.K,A.L)
    k=json.loads((A.KPHASE/(A.K+".json")).read_bytes())
    l=json.loads((A.LPHASE/"actual_replay"/A.L/"receipt.json").read_bytes())
    assert k["target"]["archived_center_shift_exact"]!="0"
    for basis,a,proof in ((k["arithmetic"],K,k),(l["independent_bounds"],L,l)):
        eta,z=map(Q,(basis["eta_upper_exact"],basis["z_upper_exact"]))
        r=Q(proof["radius_exact"])
        e=eta/(1-z)
        assert r==A.RADIUS and eta+z*r<r and Q(basis["selfmap_lhs_upper_exact"])==eta+z*r
        assert e==Q(basis["posterior_error_upper_exact"])==a["posterior"] and e<r
        ref=[q for row in a["reference"] for q in row]
        fixed=set(range(650))-set(map(int,a["a"]["free"]))
        vectors={int(i):Q(float(v)) for i,v in zip(a["a"]["free"],a["a"]["u"],strict=True)}
        assert all(x+vectors.get(i,0)==a["center"][i] and a["radii"][i]==(0 if i in fixed else r)
                   for i,x in enumerate(ref))
        assert a["exact_drive"]==a["gain"]*a["dose"]
        assert all(a["reference"][int(i)][0]>=0 for ii,flag in enumerate(a["a"]["active"]) if int(flag[0])
                   for i in a["a"]["triangles"][ii,:3])
    assert K["dose"]<L["dose"] and K["gain"]==L["gain"]>0 and K["normalization"]==L["normalization"]>0
    assert K["material"]==L["material"] and len(pins)==len(p["inputs_sha256"])==22
    for path in (A.KPHASE,A.LPHASE/"actual_replay"/A.L):
        native=json.loads((path/"os_exit_observation.json").read_bytes())
        assert type(native["unified_exec_session"]) is int and native["unified_exec_session"]>0
        assert type(native["observed_exit_code"]) is int and native["observed_exit_code"]==0


def main():
    assert A.digest(A.HERE/"active_energy.py")==SOURCE
    p=A.load(A.HERE/"protocol.json",PROTOCOL,{})
    pf=A.load(A.HERE/"preflight.json",PREFLIGHT,{})
    assert p["source_sha256"]==pf["source_sha256"]==SOURCE and pf["protocol_sha256"]==PROTOCOL
    assert pf["actual_energy_evaluations"]==p["cold_actual_energy_evaluations"]==0
    # Confirm all trusted geometry and OS artifacts using ordinary source input API
    # before exercising only synthetic elements. Real evaluate() is barred.
    with patch.object(A,"evaluate",side_effect=AssertionError("actual endpoint energy forbidden in prep")):
        cases,shared,pins=A.inputs()
        assert p==A.protocol(cases,shared,pins)
        verify_scope(cases,pins,p)
        test_elements()
    print(A.save(HERE/"synthetic_and_authority_readback.json",{"status":"PASS","source_sha256":SOURCE,
        "protocol_sha256":PROTOCOL,"preflight_sha256":PREFLIGHT,"fixture_source_sha256":A.digest(Path(__file__)),
        "checks":{"exact_oblique_affine_and_quadratic_scalar_oracles":True,
        "axis_triangle_all_six_vertex_and_midside_permutations":True,
        "signed_coordinate_intervals_include_nine_independent_values":True,
        "same_certified_K_L_armed_proofs_centers_e_no_archive_sigma":True,
        "exact_nonzero_dose_gap_fixed_support_and_common_positive_scale":True,
        "correct_sign_conclusions_at_ordering_boundary":True},
        "new_actual_endpoint_energy_evaluations":0,"other_mechanical_jobs":0}))


if __name__=="__main__":
    main()
