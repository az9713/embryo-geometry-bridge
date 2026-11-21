"""Cold independent checks for typed three-channel hull replay authority."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
import hull_replay as H


SOURCE="d5cb7ccaba8526e9ca3f9d413ba0f6ea2ae24328f43cc41825af5314044adbfc"
PROTOCOL="d7c4d611814350dff83aba7559ebfaaaff852bb232c2a0199b903a4740d8d2b3"
PREFLIGHT="f51a9789357238b6e39f6fe1806ffcafd5a66f693f6f3bc9919bba5516757dce"


def reject(call,reason=None):
    try:
        call()
    except ValueError as e:
        return reason is None or reason in str(e)
    return False


def test_hull(local):
    # Separate point force, signed upper Hessian, and same-Y correction.
    indices=[7,2]
    s={"kind":"surface","slot":5,"global_indices":indices}
    F0={0:(F(-1,3),F(1,7)),1:(F(-3),F(-1))}
    F1={0:(F(0),F(2,7)),1:(F(-4),F(-1))}
    H0={(0,0):(F(2),F(5)), (0,1):(F(-3),F(-2)),(1,1):(F(6),F(7))}
    H1={(0,0):(F(2),F(3)), (0,1):(F(-2),F(0)),(1,1):(F(5),F(7))}
    C0={0:(F(-1,8),F(0)),1:(F(0),F(0))}
    C1={0:(F(0),F(1,8)),1:(F(0),F(0))}
    good={**s,"checked_cache_sites":12,"root_asserted":False,"only_independently_valid_enclosures_used_for_hull_truth":True,
          "overlap_not_used_as_proof":True,"point_force_scope":H.S.POINT,"whole_H_scope":H.S.BOX,"radius_exact":str(H.S.RADIUS),
          "status":"PASS_INDEPENDENT_COMPONENT_WITH_EXACT_HULLS","signed_channels":{
          "point_force":local.labelled_channel(F0,F1),"whole_box_H":local.labelled_channel(H0,H1),
          "point_same_Y_cache_correction":local.labelled_channel(C0,C1)}}
    assert H.require_hull_shape(s,good,local,"surface",5) is None
    assert good["signed_channels"]["point_force"]["rows"][0]["hull_exact"]==["-1/3","2/7"]
    assert good["signed_channels"]["whole_box_H"]["rows"][1]["hull_exact"]==["-3","0"]
    assert good["signed_channels"]["point_same_Y_cache_correction"]["rows"][0]["hull_exact"]==["-1/8","1/8"]
    for key in H.CHANNELS:
        missing=deepcopy(good)
        missing["signed_channels"].pop(key)
        assert reject(lambda:H.require_hull_shape(s,missing,local,"surface",5))
    typed_upper=deepcopy(good)
    typed_upper["signed_channels"]["whole_box_H"]["rows"][1]["index"]=[False,1]
    assert reject(lambda:H.require_hull_shape(s,typed_upper,local,"surface",5),"typed ordered")
    for field,spoof in (("hull_exact",["0","1/8"]),("producer_contains_independent",True),("disjoint",0)):
        bad=deepcopy(good)
        bad["signed_channels"]["point_same_Y_cache_correction"]["rows"][0][field]=spoof
        assert reject(lambda:H.require_hull_shape(s,bad,local,"surface",5),"exact full union")
    below={**C1,0:(F(1,32),F(1,8))}
    bad=deepcopy(good)
    bad["signed_channels"]["point_same_Y_cache_correction"]=local.labelled_channel(C0,below)
    assert reject(lambda:H.require_hull_shape(s,bad,local,"surface",5),"forces HOLD")
    bad["status"]="HOLD_DISJOINT_SIGNED_ENCLOSURES"
    assert H.require_hull_shape(s,bad,local,"surface",5) is None


def test_phase(profiles,p,f,b):
    current={"status":H.SOURCE_STATUS,"checks":{"dummy":True},**b,"permitted_waves":list(H.X.WAVES),
             "requires_external_raw_observed_exit_for_each_wave":True,"new_raw_spatial_workers_authorized":False}
    assert H.review_fields(p,current,f,b) is None
    for fields in ({"checks":{"dummy":1}},{"requires_external_raw_observed_exit_for_each_wave":False},
                   {"new_raw_spatial_workers_authorized":True},{"permitted_waves":list(H.X.WAVES[::-1])},
                   {"wave_execution_source_sha256":"f"*64}):
        assert reject(lambda:H.review_fields(p,{**current,**fields},f,b))
    first_id=profiles[H.X.WAVES[0]]["case_order"][0]
    assert reject(lambda:H.wave_scope(profiles[H.X.WAVES[1]],first_id))
    first=profiles[H.X.WAVES[0]]
    t=H.wave_scope(first,first_id)
    assert t["original_observer_field_sha256"]==first["targets"][first_id]["original_observer_field_sha256"]
    # Full exterior metadata admission is tested in the upstream adapter audit.
    # Here test that its referenced bytes are read before any local arithmetic.
    args=SimpleNamespace(wave=H.X.WAVES[0],case_id=first_id,raw_wave_admission_sha256="c"*64,
                         raw_root_observed_exit_sha256="d"*64)
    with TemporaryDirectory(dir=HERE) as scratch:
        stub=Path(scratch)/"not_a_real_packet.txt"
        stub.write_text("stub",encoding="utf-8")
        mark=H.rel(stub)
        with patch.object(H,"gate",return_value=(first,b)), patch.object(H.X,"open_raw_phase",
                 return_value=(first,{}, {}, {"sha256":{mark:"0"*64}})), patch.object(H.X,"exact_context",
                 return_value=({"target":t},{},(F(1),F(1)),dict.fromkeys(H.S.OPS))), patch.object(H,"math",
                 side_effect=AssertionError("corrupt raw bytes advanced to local arithmetic")) as kernel:
            assert reject(lambda:H.run(args),"exact root-observed complete later-wave raw bytes")
            assert not kernel.called


def main():
    p=H.read(H.PREP/"protocol.json",PROTOCOL)
    f=H.read(H.PREP/"preflight.json",PREFLIGHT)
    assert SOURCE==H.sha(H.SOURCE)==H.sha(H.PREP/"hull_replay_frozen.py")==p["source_sha256"]==f["source_sha256"]
    profiles,release=H.static()
    assert list(profiles)==list(H.X.WAVES) and len({c for v in profiles.values() for c in v["case_order"]})==28
    assert H.V.yes(H.fixtures(profiles))
    b={"source_sha256":SOURCE,"protocol_sha256":PROTOCOL,"preflight_sha256":PREFLIGHT,
       "later_waves_release_sha256":release["later_waves_release_sha256"],
       "wave_execution_source_sha256":H.PINS[H.X.SOURCE],"local_bound_and_hull_source_sha256":H.PINS[H.HULL_SOURCE]}
    test_phase(profiles,p,f,b)
    test_hull(H.math())
    out={"status":"PASS","source_sha256":SOURCE,"protocol_sha256":PROTOCOL,"preflight_sha256":PREFLIGHT,
         "independent_fixture_source_sha256":H.sha(Path(__file__)),"checks":{
             "all_seven_frozen_wave_profiles_and_original_scope":True,
             "same_active_review_release_phase_guards":True,"raw_native_stub_changed_bytes_fail_before_local_math":True,
             "full_signed_three_channels_exact_own_union_replayed":True,
             "same_Y_correction_alone_disjoint_forces_hold":True,"typed_order_fake_diagnostics_and_channel_narrowing_rejected":True},
         "actual_numerical_packets_opened":0,"scientific_jobs":0}
    print(H.save(HERE/"synthetic_readback.json",out),flush=True)


if __name__=="__main__":
    main()
