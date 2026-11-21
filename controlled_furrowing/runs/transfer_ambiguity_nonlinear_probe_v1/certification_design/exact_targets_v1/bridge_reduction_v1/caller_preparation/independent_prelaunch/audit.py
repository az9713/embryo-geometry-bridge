"""Cold authority and exact arithmetic checks for the frozen bridge reducer caller."""
from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import json
import sys
from unittest.mock import patch

SOURCE = Path(__file__).resolve()
ROOT = SOURCE.parents[7]
sys.path.insert(0,str(ROOT/"scripts"))
import transfer_ambiguity_exact_bridge_reduce_driver as R

SRC = "6eafff943ae116d4f9550f480f5a3c6e75648e5428291aa263238a762555a998"
PRO = "bdb7d8576aad46407a2a6caef8dfc40b79b54034feb0c01803d77e470df94a40"
PRE = "7b61454b090aa18cfb4fcf6188a0a3b79047b8ea323e3b678d7dd8c67cce8213"


def main():
    assert R.sha(R.SOURCE)==R.sha(R.PREP/"transfer_ambiguity_exact_bridge_reduce_driver_frozen.py")==SRC
    protocol=R.pinned(R.PREP/"protocol.json",PRO)
    pre=R.pinned(R.PREP/"preflight.json",PRE)
    assert pre["status"]=="PASS" and len(pre["checks"])==8 and R.yes(pre["checks"])
    assert pre["actual_packet_payloads_opened"]==pre["actual_matrix_work"]==0
    assert protocol["radius_exact"]==str(R.RADIUS) and protocol["precision_bits"]==640
    deny=AssertionError("cold caller attempted actual admission or payload processing")
    with patch.object(R,"authorize",side_effect=deny),patch.object(R,"context",side_effect=deny), \
            patch.object(R,"run",side_effect=deny),patch.object(R,"loader",side_effect=deny):
        assert R.authority_protocol()==protocol and R.fixtures()==pre["checks"]
    kernel=R.loader(R.KERNEL_SOURCE,"_cold_tested_pure_bridge_kernel")
    kernel.ctx.prec=640
    assert kernel.RADIUS==R.RADIUS and kernel.POINT==R.POINT and kernel.BOX==R.BOX
    width=R.RADIUS
    # Nonidentity R distinguishes an already-preconditioned anchor from raw gradient.
    # H is zero so the contraction check is force-dominated, and each exact row cancels.
    actual=kernel.contract([[Q(2),Q(1)],[Q(0),Q(3)]],
                            [(width,width),(-width,-width)],{},
                            [(-width,-width),(3*width,3*width)],Q(1,4))
    assert actual["classification"]=="SUFFICIENT_UNIQUE_LOCAL_ZERO"
    assert actual["signed_RF_defect_rows_exact"]==[[str(width)]*2,[str(-3*width)]*2]
    assert actual["signed_RF_I_rows_exact"]==[["0","0"],["0","0"]]
    assert actual["b_upper_exact"]==actual["posterior_error_upper_exact"]=="0"
    assert actual["minimum_or_continuum_claimed"] is False
    z=Q(1,4)
    assert kernel.classify([((1-z)*width,(1-z)*width)],z,width)["classification"]=="UNRESOLVED"
    assert kernel.classify([((1+z)*width,(1+z)*width)],z,width)["classification"]=="UNRESOLVED"
    assert kernel.classify([(2*width,3*width)],z,width)["classification"]=="SUFFICIENT_NO_ZERO_IN_THIS_CUBE"

    # Additional two-column signed containment and typed ordering fixtures.
    f="a"*64
    case={"kind":"surface","slot":4,"global_indices":[9,2],"point_force_scope":R.POINT,
          "whole_box_H_scope":R.BOX,"plan_file_sha256":f,"packet_file":"unused","packet_sha256":f}
    indep={**case,"checked_force_rows":2,"checked_upper_H_entries":3,"checked_cache_sites":12,
           "containment_failures":dict.fromkeys(R.FAILURE_KEYS,[]),
           "independent_point_F_exact":[["-1","1"],["2","3"]],
           "independent_whole_box_H_exact":[[0,0,"1","2"],[0,1,"-1","1"],[1,1,"2","3"]]}
    comp={k:case[k] for k in ("kind","slot","global_indices")}
    packet={**comp,"bindings":{"case":1},"record_id":"fixture","root_asserted":False,
            "precision_bits":640,"plan_file_sha256":f,"proof_sha256":f,
            "point_force":{"scope":R.POINT,"F_I_minus_F_cache_intervals_exact":[["-2","2"],["1","3"]],
                           "Qstar_minus_cache_correction_exact":[["0","0"],["0","0"]]},
            "whole_box_hessian":{"scope":R.BOX,"radius_exact":str(R.RADIUS),
                                 "H_I_minus_H_cache_upper_exact":[[0,0,"0","2"],[0,1,"-2","2"],[1,1,"2","4"]]}}
    def apply(m=comp,p=packet,s=indep):
        return R.component_operands(case,s,p,comp,m,"fixture",{"case":1},f)
    assert apply()==(comp,packet)
    assert R.rejects(lambda:apply(m={**comp,"global_indices":[2,9]}))
    altered=deepcopy(packet)
    altered["whole_box_hessian"]["H_I_minus_H_cache_upper_exact"][1]=[0,1,"0","2"]
    assert R.rejects(lambda:apply(p=altered))
    assert R.rejects(lambda:apply(s={**indep,"checked_cache_sites":True}))
    assert R.rejects(lambda:apply(s={**indep,"checked_upper_H_entries":2}))
    assert R.rejects(lambda:apply(p={**packet,"root_asserted":0}))
    assert not R.ACTUAL.exists()

    g={"status":R.REVIEW_STATUS,"source_sha256":SRC,"protocol_sha256":PRO,"preflight_sha256":PRE,
       "pure_kernel_review_sha256":R.PINS[R.KERNEL_REVIEW],
       "native_replay_source_sha256":R.PINS[R.REPLAY_SOURCE],
       "native_replay_protocol_sha256":R.PINS[R.REPLAY_PROTOCOL],
       "native_replay_preflight_sha256":R.PINS[R.REPLAY_PREFLIGHT],
       "review_source_sha256":R.sha(SOURCE),
       "checks":{"unchanged_frozen_source_and_current_kernel_native_protocol_pins":True,
                 "native_replay_requires_true_OS_and_exact_complete_168_catalogues":True,
                 "same_exact_singleton_drive_cube_and_binary_R_reverified_before_payload_reduction":True,
                 "signed_cached_rows_preconditioned_once_plus_signed_R_defect":True,
                 "all_producer_point_and_H_local_endpoints_contain_independent_replay":True,
                 "complete_typed_and_ordered_projection_and_hessian_fixtures":True,
                 "strict_classification_and_scope_on_synthetic_balls":True,
                 "no_actual_payload_matrix_or_mechanics_entered_during_review":True},
       "scope":"Cold caller approval only. Runtime requires native 168 replay containment, producer and native process evidence, and a separate independent readback of any actual arithmetic candidate. Signed cached anchor is the independent Arb enclosure of R g, not raw g.",
       "direct_cache_diagnostic":"Caller reads native pinned per-component no-failure flags from direct cache correction. It independently checks exact endpoint inclusion for the total forces and Hessians it actually consumes.",
       "pins_rebound":len(protocol["source_and_prerequisite_sha256"])}
    path=SOURCE.with_name("review.json")
    h=R.save(path,g)
    p,bindings=R.review_gate(path,h)
    assert p==protocol and bindings["prelaunch_review_sha256"]==h
    print(json.dumps({"caller_prelaunch_gate_sha256":h,"pins_rebound":g["pins_rebound"]}))


if __name__=="__main__":
    main()
