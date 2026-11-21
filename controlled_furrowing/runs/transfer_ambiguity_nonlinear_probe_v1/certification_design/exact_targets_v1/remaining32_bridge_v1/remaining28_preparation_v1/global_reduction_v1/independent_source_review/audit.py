"""Independent cold wire and exact signed endpoint oracle for the later-wave reducer."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import reduce as A

SOURCE="3051ac6d35ee18407e29e8b82e96d789e685073dd3f067b01b50b459c0b970e5"
PROTOCOL="1198124fe22265e81113a3f383a82c6b816334e5d3a54df837bd16eb491aadea"
PREFLIGHT="b3310545e2ab4e15e5988c90d4011565eebd5797f0c9449c75bcfad52a5f8b88"


def rejected(f):
    try:
        f()
    except (KeyError,ValueError):
        return True
    return False


def point_data(ids,s,f,h):
    route={"kind":"bulk","slot":s,"global_indices":ids}
    upper=[[i,j,*v] for (i,j),v in zip(((0,0),(0,1),(1,1)),h,strict=True)]
    return route,{**route,"point_force":{"scope":A.S.POINT,"F_I_minus_F_cache_intervals_exact":f},
                  "whole_box_hessian":{"scope":A.S.BOX,"radius_exact":str(A.RADIUS),"H_I_minus_H_cache_upper_exact":upper}}


def oracle(I):
    p0,x=point_data([2,0],0,[["2","3"],["-1","1"]],[["-1","2"],["3","4"],["2","3"]])
    p1,y=point_data([1,0],1,[["0","2"],["-2","0"]],[["0","1"],["-2","-1"],["-3","-2"]])
    f,h=A.K.assemble([x,y],[p0,p1],3)
    independentF=[(Q(0),Q(0)) for i in range(3)]
    independentH={}
    for route,packet in ((p0,x),(p1,y)):
        I.scatter_add(independentF,independentH,route["global_indices"],
                list(map(lambda v:tuple(map(Q,v)),packet["point_force"]["F_I_minus_F_cache_intervals_exact"])),
                [tuple(map(Q,v[2:])) for v in packet["whole_box_hessian"]["H_I_minus_H_cache_upper_exact"]])
    assert independentF==f==[(-3,1),(0,2),(2,3)]
    assert independentH==h=={(0,0):(-1,1),(0,1):(-2,-1),(0,2):(3,4),(1,1):(0,1),(2,2):(-1,2)}
    R=((Q(1,7),Q(-1,5),Q(3,7)),(Q(-5,6),Q(2,5),Q(-1,7)),(Q(0),Q(-4,9),Q(-1,3)))
    anchor=((Q(-7,19),Q(-2,7)),(Q(2,5),Q(3,7)),(Q(-3,11),Q(4,13)))
    A.K.ctx.prec=640
    one=A.K.contract(R,f,h,anchor,Q(1,4),radius=A.RADIUS)
    two=I.classify(f,h,anchor,R,Q(1,4),A.RADIUS)
    for i,row in enumerate(R):
        corners=[sum(x*y for x,y in zip(row,v,strict=True))+a for v in product(*f) for a in anchor[i]]
        for data,key in ((one,"signed_RF_I_rows_exact"),(two,"signed_R_F_I_rows_exact")):
            lo,hi=map(Q,data[key][i])
            assert lo<=min(corners)<=max(corners)<=hi
    # Independent per-entry exact rational R*H extrema on the signed symmetric box.
    rows=[]
    for r in R:
        entries=[]
        for j in range(3):
            values=[tuple(q*end for end in h.get(tuple(sorted((i,j))),(0,0))) for i,q in enumerate(r)]
            entries.append(max(abs(sum(min(v) for v in values)),abs(sum(max(v) for v in values))))
        rows.append(sum(entries))
    producer=Q(one["b_upper_exact"])
    independent=Q(two["b_H_upper_exact"])
    exact=max(rows)
    assert producer>=exact and independent>=exact and independent-exact<Q(1,2**700)
    assert Q(two["z_upper_exact"])==Q(1,4)+independent
    radius=A.RADIUS
    case=lambda c:I.classify([(c,c)],{},[(Q(0),Q(0))],[[Q(1)]],Q(1,2),radius)
    small=case(radius/4)
    assert small["strict_selfmap"] and Q(small["posterior_error_upper_exact"])==radius/2
    assert not case(radius/2)["strict_selfmap"] and not case(radius/2)["strict_cube_exclusion"]
    assert not case(3*radius/2)["strict_cube_exclusion"] and case(2*radius)["strict_cube_exclusion"]
    fixture={"archived_center_shift_exact":str(radius/16),"original_addition_gap_exact":str(radius/128),
             "observer_goal_exact":str(3*radius/4)}
    t=I.observer_minimum(small,fixture,True,Q(3))
    assert t["strict_minimum_proved"] and t["archived_state_budget_proved"]
    assert Q(t["physical_nodal_error_um_upper_exact"])==Q(67,2)*(radius/2+radius/16+radius/128)
    assert Q(t["normalized_H_euclidean_lower_exact"])==Q(1,6)
    assert I.observer_minimum(case(2*radius),fixture,True,Q(3))["strict_minimum_proved"] is False


def wire(I):
    name="point_same_Y_cache_correction"
    row={"index":0,"independent_exact":["-1/3","2/7"],"producer_exact":["0","2/3"],
         "hull_exact":["-1/3","2/3"],"producer_contains_independent":False,"disjoint":False}
    c={"rows":[row],"producer_containment_misses":[0],"disjoint_indices":[]}
    def check(ch):
        x,n=A.D.hull_channel(ch,[0],[["0","2/3"]])
        y,m=I.validated_channel({"signed_channels":{name:ch}},name,[0],[["0","2/3"]])
        assert list(map(lambda v:tuple(map(Q,v)),x))==y and n==m==1
    check(c)
    for mutation in ({"index":False},{"hull_exact":["0","2/7"]},{"producer_exact":["0","1"]},
                     {"independent_exact":["1","2"],"hull_exact":["0","2"],"disjoint":True}):
        assert rejected(lambda:check({**c,"rows":[{**row,**mutation}]}))


def observed_producer_mock(profile):
    cases=profile["case_order"]
    with TemporaryDirectory(dir=HERE) as tmp:
        d=Path(tmp)
        p=d/"producer"
        p.mkdir()
        a=SimpleNamespace(wave=profile["wave_id"],hull_admission_sha256="a"*64,raw_wave_admission_sha256="b"*64,
                 raw_root_observed_exit_sha256="c"*64,protocol_sha256=PROTOCOL,preflight_sha256=PREFLIGHT,review_sha256="d"*64)
        phase={"mode":"producer","wave_id":a.wave,"pid":144,
               "hull_admission_sha256":a.hull_admission_sha256,"raw_wave_admission_sha256":a.raw_wave_admission_sha256,
               "raw_root_observed_exit_sha256":a.raw_root_observed_exit_sha256,"source_sha256":SOURCE,
               "protocol_sha256":PROTOCOL,"preflight_sha256":PREFLIGHT,"source_review_sha256":a.review_sha256}
        rows=[]
        for case in cases:
            folder=p/case
            folder.mkdir()
            rows.append({"case_id":case,"operands_sha256":A.save(folder/"operands.json",{"inert":True}),
                         "case_results_sha256":A.save(folder/"results.json",{"inert":True})})
        sr=A.save(p/"started.json",{"status":"RUNNING",**phase})
        rr=A.save(p/"results.json",{"status":A.PROD_STATUS,"cases":rows,**phase})
        tr=A.save(p/"terminal.json",{"status":"FINISHED","results_sha256":rr,"exit_code":0,**phase})
        obs={"status":A.PROD_OS,"wave_id":a.wave,"pid":144,"hull_admission_sha256":a.hull_admission_sha256,
             "observed_exit_code":0,"unified_exec_session":44,"started_sha256":sr,"terminal_sha256":tr,"results_sha256":rr}
        a.producer_os_sha256=A.save(p/"root_observed_exit.json",obs)
        assert tuple(A.producer_observed(a,d,cases))==tuple(cases)
        assert rejected(lambda:A.producer_observed(a,d,list(reversed(cases))))
        p0=p/cases[-1]/"operands.json"
        p0.write_text('{"corrupted":true}',encoding="utf-8")
        assert rejected(lambda:A.producer_observed(a,d,cases))


def main():
    p,profiles=A.static()
    pf=A.read(A.HERE/"preflight.json",PREFLIGHT)
    assert p==A.read(A.HERE/"protocol.json",PROTOCOL)
    assert p["source_sha256"]==pf["source_sha256"]==A.sha(A.SOURCE)==SOURCE
    assert A.V.yes(pf["checks"]) and A.V.yes(A.fixtures(p))
    first=profiles[A.X.WAVES[0]]
    obs={"status":A.HULL_OS,"wave_id":first["wave_id"],
         "cases":[{"case_id":c,"pid":20+i,"unified_exec_chunk_id":str(23+i),"observed_exit_code":0}
                  for i,c in enumerate(first["case_order"])]}
    assert len(A.hull_observation(obs,first["wave_id"],first["case_order"],{}))==4
    for wrong in ({"unified_exec_chunk_id":" "},{"observed_exit_code":False},{"pid":20.0},{"pid":21}):
        assert rejected(lambda:A.hull_observation({**obs,"cases":[{**obs["cases"][0],**wrong},*obs["cases"][1:]]},
                                                  first["wave_id"],first["case_order"],{}))
    I=A.endpoint_kernel()
    oracle(I)
    wire(I)
    observed_producer_mock(first)
    print(A.save(HERE/"synthetic_readback.json",{"status":"PASS","source_sha256":SOURCE,"protocol_sha256":PROTOCOL,
         "preflight_sha256":PREFLIGHT,"fixture_source_sha256":A.sha(Path(__file__)),
         "checks":{"seven_cold_profiles_and_32_positive_sym_R_inputs":True,
         "exact_3d_signed_force_and_symmetric_H_oracle_two_independent_kernels":True,
         "no_interval_shrinking_or_disjoint_same_Y_wire_in_either_mode":True,
         "strict_boundary_and_cube_only_exclusion_and_original_observer_terms":True,
         "typed_four_native_hull_exits_and_separate_sealed_producer_metadata":True},
         "scientific_hull_packets_opened":0,"actual_mechanical_or_global_jobs":0}))


if __name__=="__main__":
    main()
