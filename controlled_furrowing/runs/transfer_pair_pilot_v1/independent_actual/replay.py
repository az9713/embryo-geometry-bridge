"""Independently review four frozen saved P2 fields, with no solve or response score."""

from collections import Counter
from dataclasses import asdict
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import math
import platform
import sys
import time

HERE = Path(__file__).resolve().parent
PILOT = HERE.parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import numpy as np
import scipy
import torch

from quadratic_endpoint_observer import Curve, positive_thickness
from quadratic_quadrature_study import AxisymmetricModel as CrossModel, exact_geometry, evaluate
from robust_quadratic import Parameters, SolvePolicy
from sparse_quadratic import AxisymmetricModel as MainModel, IntegrationPolicy, SparsePolicy

BINDINGS = {
    "driver_sha256": ("scripts/transfer_pair_pilot.py", "a05d677ac7d919c2f012b20040c663e6f181563bc008fd4c216a81e6a30b7d18"),
    "protocol_sha256": ("runs/transfer_pair_pilot_v1/preparation/protocol.json", "a2745a7e50e9a64ee6e5e597f04b30981c7910ad5e3ca2d3a16e6d4d0f3de927"),
    "preflight_sha256": ("runs/transfer_pair_pilot_v1/preparation/preflight.json", "c6defad69457de5d28cffb1c1a657ddf36b71575b75a7078a2dff7df875a816d"),
    "prelaunch_review_sha256": ("runs/transfer_pair_pilot_v1/prelaunch_independent_review.json", "07b15616876d068a48fc0e0177ffaf2e63085db6a8bb46573cdb407214c0a21f"),
}
RESULTS = "6a2fcb2a0c4b7a51301d025a7725247ec81a999764b7e85a0e97b5f6e139be05"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def bound(path, digest, bindings):
    path = path.resolve()
    key = path.relative_to(ROOT).as_posix()
    assert sha(path) == digest == bindings.setdefault(key, digest), key
    return read(path)


def write(path, data):
    with path.open("x", encoding="utf8") as output:
        json.dump(data, output, indent=2)
        output.write("\n")


def yes(checks):
    return bool(checks) and all(x is True for x in checks.values())


def signed_difference(a, b):
    return a[0]-b[1], a[1]-b[0]


def contains(big, small):
    return big[0] <= small[0] <= small[1] <= big[1]


def packed(a):
    return {"exact": list(map(str, a)), "display_um": list(map(float, a)),
            "interval_width_um_display": float(a[1]-a[0])}


def parsed(row):
    assert row["lower_closed"] is True and row["upper_closed"] is True
    return tuple(Q(row[s]["numerator"], row[s]["denominator"]) for s in ("lower","upper"))


def geometric_controls(state, length):
    """Independently count boundary edges and read the P2 midpoint of each."""
    reference = [tuple(map(Q, xy)) for xy in state["reference_nodes"]]
    current = [tuple(map(Q, xy)) for xy in state["current_nodes"]]
    assert len(reference) == len(current) == 325
    triangle_edges = Counter()
    midpoints = {}
    for a,b,c,ab,bc,ca in state["triangles"]:
        for i,j,m in ((a,b,ab),(b,c,bc),(c,a,ca)):
            edge = tuple(sorted((i,j)))
            triangle_edges[edge] += 1
            assert midpoints.setdefault(edge, m) == m
    assert set(triangle_edges.values()) == {1,2}
    out = {}
    for face, saved, z in (("apical", "top", Q(state["parameters"]["height"])),
                           ("basal", "bottom", Q(0))):
        triples = []
        for (i,j), n in triangle_edges.items():
            if n == 1 and reference[i][1] == reference[j][1] == z:
                a,b = sorted((i,j), key=lambda k: reference[k][0])
                triples.append((reference[a][0], [current[a], current[midpoints[i,j]], current[b]]))
        controls = [points for _,points in sorted(triples)]
        assert len(controls) == 12
        original = [[tuple(map(Q,point)) for point in edge]
                    for edge in state[f"boundary_{saved}_quadratic_controls"]]
        assert controls == original
        physical = [[(length*r,length*z) for r,z in edge] for edge in controls]
        assert all(a[-1]==b[0] for a,b in zip(physical,physical[1:]))
        assert physical[0][0][0] == 0 and physical[-1][-1][0] == length*Q(state["parameters"]["outer_radius"])
        out[face] = physical
    return out


def independent_height(controls, radius, bisections):
    """Exact current-radius inversion using Lagrange nodal shape functions."""
    found = []
    for edge in controls:
        if not edge[0][0] <= radius <= edge[-1][0]:
            continue
        def eval_axis(t,k):
            return sum(v * node[k] for v,node in
                       zip(((1-t)*(1-2*t),4*t*(1-t),t*(2*t-1)),edge))
        lo, hi = Q(0),Q(1)
        if radius==edge[0][0]:
            hi=lo
        elif radius==edge[-1][0]:
            lo=hi
        else:
            for _ in range(bisections):
                m=(lo+hi)/2
                r=eval_axis(m,0)
                if r==radius:
                    lo=hi=m
                    break
                if r<radius:
                    lo=m
                else:
                    hi=m
        assert eval_axis(lo,0)<=radius<=eval_axis(hi,0)
        z0,zm,z1=(v[1] for v in edge)
        nodes=[lo,hi]
        curvature=2*(z0-2*zm+z1)
        if curvature:
            p=(3*z0-4*zm+z1)/(2*curvature)
            if lo<=p<=hi:
                nodes.append(p)
        found.extend(eval_axis(t,1) for t in nodes)
    assert found
    return min(found),max(found)


def replay_surface(state, stored, length, shoulder, bisections, is_final):
    controls = geometric_controls(state,length)
    derivatives = {}
    for face, edges in controls.items():
        expected = [[tuple(map(Q,p)) for p in edge] for edge in stored["controls_um"][face]]
        assert edges == expected
        end_derivatives = []
        for (ra,_),(rm,_),(rb,_) in edges:
            end_derivatives.extend((4*rm-3*ra-rb,ra-4*rm+3*rb))
        lower = min(end_derivatives)
        assert lower == Q(stored["minimum_radius_derivative_um_per_local_parameter_exact"][face]) > 0
        derivatives[face]=str(lower)
    top,bottom = Curve(controls["apical"],bisections), Curve(controls["basal"],bisections)
    assert top.domain==bottom.domain
    thickness = positive_thickness(top,bottom,top.domain,Q(1,1000),8192)
    assert thickness == Q(stored["minimum_current_thickness_um_lower_exact"]) > 0
    center = independent_height(controls["apical"],Q(0),bisections)
    outer = independent_height(controls["apical"],shoulder,bisections)
    result = signed_difference(center,outer)
    assert center == parsed(stored["apical_center_um"]) and outer==parsed(stored["apical_shoulder_um"])
    assert result==tuple(map(Q,stored["Q_um_exact"]))==parsed(stored["Q_um"])
    profiles = stored.get("sampled_full_current_profiles",[])
    assert bool(profiles)==is_final
    if profiles:
        assert len(profiles)>=101
        radii = [Q(row["radius_um_exact"]) for row in profiles]
        assert radii==sorted(set(radii)) and radii[0]==top.domain[0] and radii[-1]==top.domain[1]
        assert shoulder in radii
        for row, radius in zip(profiles,radii):
            for face in ("apical","basal"):
                assert parsed(row[f"{face}_height_um"]) == independent_height(controls[face],radius,bisections)
    return result, thickness, derivatives, len(profiles)


def assert_mechanical_close(new, old):
    for key in ("scaled_energy","energy","scaled_force_linf","sampled_J_min","sampled_J_max"):
        assert math.isclose(new[key],old[key],rel_tol=5e-12,abs_tol=3e-14), (key,new[key],old[key])
    for key,v in new["energy_components"].items():
        assert math.isclose(v,old["energy_components"][key],rel_tol=5e-12,abs_tol=3e-14), key
    return abs(new["scaled_force_linf"]-old["scaled_force_linf"])


def replay_case(case, job, result, protocol, hashes):
    label = case["case_id"]
    directory = PILOT/"jobs"/label
    local = read(directory/"job.json")
    assert local==job and local["case_id"]==label
    assert job["started"] is True and job["terminal"] is True and job["exit_code"]==0
    assert job["status"]=="PROCESS_EXITED" and job["started_at_utc"]<job["ended_at_utc"]
    assert job["bindings"]==result["bindings"] and job["worker_status"]=="COMPLETE_CONDITIONAL_STORED_FIELD"
    assert (directory/"stderr.log").stat().st_size==0
    hashes[(directory/"job.json").relative_to(ROOT).as_posix()]=sha(directory/"job.json")
    worker = bound(ROOT/job["worker_receipt"],job["worker_receipt_sha256"],hashes)
    assert worker["worker_pid"]==job["pid"] and worker["case"]==case
    assert worker["status"]==job["worker_status"] and worker["bindings"]==result["bindings"]
    assert worker["source_hashes"]==protocol["source_hashes"] and worker["input_hashes"]==protocol["input_hashes"]
    assert worker["root_validated"] is False and len(worker["checkpoints"])==4
    assert job["started_at_utc"] <= worker["started_at_utc"] <= worker["worker_finished_at_utc"] <= job["ended_at_utc"]
    assert result["diagnostic_fields_accepted"][label]==worker["checkpoints"][-1]==worker["last_retained_field"]

    params = Parameters(**case["parameters"])
    options = protocol["numerical_policy"]
    primal = IntegrationPolicy(**options["solve_quadrature"])
    alternate = IntegrationPolicy(**options["cross_force_quadrature"])
    solve = SolvePolicy(**options["solve"])
    main = MainModel(params,solve,primal,SparsePolicy(**options["sparse"]))
    cross = CrossModel(params,solve,alternate)
    for attr in ("reference","triangles","free","fixed","activation_elements"):
        assert np.array_equal(getattr(main,attr),getattr(cross,attr))
    mesh = worker["mesh"]
    assert (mesh["reference_node_count"],mesh["triangle_count"],mesh["free_dof_count"])==(325,144,611)
    assert all(mesh[k] is True for k in ("patch_radius_is_vertex_aligned","axial_activation_faces_aligned",
                                        "has_axis_and_outer_rim_fixed_components","positive_reference_element_area"))
    gain = Q(case["gain_exact"])
    candidate = protocol["calibration_candidates"][case["candidate"]]
    assert gain==Q(candidate["gain_exact"])
    assert all(case["parameters"][k]==v for k,v in candidate["parameters"].items() if k!="patch_radius")
    length = Q(protocol["observable"]["length_unit_um_exact"])
    assert length==Q("33.5") and Q(case["diameter_um_exact"])/(2*length)==Q(case["patch_radius_model_exact"])
    assert Q(case["parameters"]["patch_radius"])==Q(case["patch_radius_model_binary64_exact"])==Q(float(Q(case["patch_radius_model_exact"])))
    shoulder=Q(protocol["observable"]["shoulder_radius_um_exact"])
    bisections=protocol["observable"]["root_bisections"]
    final_dose=Q(protocol["input_selection"]["matched_dose_exact"])
    rows=[]
    stdout=[json.loads(line) for line in (directory/"stdout.log").read_text().splitlines()]
    assert [(v["case_id"],v["checkpoint"],v["status"]) for v in stdout]==[(label,i,"PASS") for i in range(4)]
    hashes[(directory/"stdout.log").relative_to(ROOT).as_posix()]=sha(directory/"stdout.log")
    hashes[(directory/"stderr.log").relative_to(ROOT).as_posix()]=sha(directory/"stderr.log")
    for i,(step,row) in enumerate(zip(case["schedule"],worker["checkpoints"])):
        stage=directory/f"checkpoint_{i:02d}"
        assert row["index"]==i and all(row[k]==v for k,v in step.items())
        dose=i*final_dose/3
        target=dose*gain
        rounded=Q(float(target))
        assert dose==Q(step["dose_exact"]) and target==Q(step["requested_strength_exact"])
        assert rounded==Q(step["solver_strength_binary64_exact"])==Q(row["actual_strength_binary64_exact"])
        assert rounded-target==Q(step["binary64_rounding_difference_exact"]) and rounded/gain==Q(row["actual_dose_implied_by_binary_strength_exact"])
        assert ROOT/row["field"]==stage/"field.json" and ROOT/row["verification"]==stage/"verification.json"
        state=bound(stage/"field.json",row["field_sha256"],hashes)
        verification=bound(stage/"verification.json",row["verification_sha256"],hashes)
        assert verification["status"]==row["numerical_diagnostic_status"]=="PASS"
        assert yes(verification["checks"]) and verification["checks"]==row["checks"] and verification["root_validated"] is False
        assert Q(state["activation_strength"])==rounded and state["reached_target"] is True and state["converged"] is True
        assert state["parameters"]==asdict(params)==case["parameters"]
        assert state["integration_policy"]==asdict(primal)
        assert np.array_equal(state["reference_nodes"],main.reference) and np.array_equal(state["triangles"],main.triangles)
        assert np.array_equal(state["activation_elements"],main.activation_elements)
        current=main.coordinates(torch.tensor(state["free_displacement"],dtype=torch.float64))
        assert np.array_equal(state["current_nodes"],current.detach().cpu().numpy())
        assert read(stage/"solver/final.json")==state
        hashes[(stage/"solver/final.json").relative_to(ROOT).as_posix()]=sha(stage/"solver/final.json")
        ledger=state["robust_continuation"]
        assert ledger["failed_state_reuse"] is False and ledger["stop_reason"]==row["continuation_stop_reason"]
        assert len(ledger["attempts"])==row["attempts"] and sum(not a["accepted"] for a in ledger["attempts"])==row["rejected_attempts"]
        if i:
            assert Q(ledger["accepted_increments"][0]["from_strength"])==Q(rows[-1]["strength_binary64_exact"])
            assert Q(ledger["accepted_increments"][-1]["to_strength"])==rounded
        assert [(str(Q(x["from_strength"])/gain),str(Q(x["to_strength"])/gain))
                for x in ledger["accepted_increments"]]==[(v["from_dose_exact"],v["to_dose_exact"])
                                                            for v in row["accepted_path_implied_optical_doses"]]
        mr, cr=evaluate(main,state),evaluate(cross,state)
        drift=(assert_mechanical_close(mr,verification["Q8_mechanical_evaluation"]),
               assert_mechanical_close(cr,verification["Q12_cross_evaluation"]))
        assert max(mr["scaled_force_linf"],cr["scaled_force_linf"])<solve.force_tolerance
        assert row["Q8_residual"]==verification["Q8_mechanical_evaluation"]["scaled_force_linf"]
        assert row["Q12_residual"]==verification["Q12_cross_evaluation"]["scaled_force_linf"]
        geom=exact_geometry(state)
        assert geom==verification["geometry"] and len(geom["triangle_certificates"])==144
        assert geom["global_embedding_certified"] and geom["all_exact_J_positive"] and Q(geom["J_lower_exact"])>0
        assert geom["J_lower_exact"]==row["exact_J_lower"] and geom["J_upper_exact"]==row["exact_J_upper"]
        surface=verification["surface"]
        q,thickness,derivatives,profile_count=replay_surface(state,surface,length,shoulder,bisections,i==3)
        assert tuple(map(Q,row["Q_um_exact"]))==q
        if i==3:
            sampled=bound(ROOT/worker["profiles"],worker["profiles_sha256"],hashes)
            assert sampled=={"case_id":label,"planned_final_checkpoint":True,
                             "diagnostic_gate_passed":True,"field_sha256":row["field_sha256"],**surface}
        rows.append({"index":i,"dose_exact":str(dose),"strength_binary64_exact":str(rounded),
                     "Q8_force_linf":mr["scaled_force_linf"],"Q12_force_linf":cr["scaled_force_linf"],
                     "same_backend_force_residual_abs_drift":list(drift),
                     "J_bound_exact":[geom["J_lower_exact"],geom["J_upper_exact"]],
                     "current_thickness_lower_um_exact":str(thickness),"radial_derivative_lower":derivatives,
                     "Q_um":packed(q),"profile_sample_count":profile_count,
                     "accepted_attempts":len(ledger["accepted_increments"]),"rejected_attempts":row["rejected_attempts"]})
        print(f"Replayed {label} checkpoint {i}",flush=True)
    return rows


def review():
    start=time.perf_counter()
    hashes={}
    expected={k:h for k,(_,h) in BINDINGS.items()}
    paths={k:ROOT/v[0] for k,v in BINDINGS.items()}
    assert all(sha(paths[k])==v for k,v in expected.items())
    for k,v in expected.items():
        hashes[paths[k].relative_to(ROOT).as_posix()]=v
    protocol=read(paths["protocol_sha256"])
    preflight=read(paths["preflight_sha256"])
    prelaunch=read(paths["prelaunch_review_sha256"])
    assert preflight["status"]==prelaunch["status"]=="PASS" and yes(preflight["checks"]) and yes(prelaunch["checks"])
    assert all(prelaunch[k]==v for k,v in expected.items() if k!="prelaunch_review_sha256")
    assert preflight["driver_sha256"]==expected["driver_sha256"] and preflight["protocol_sha256"]==expected["protocol_sha256"]
    bound(ROOT/prelaunch["review_receipt_file"],prelaunch["review_receipt_sha256"],hashes)
    inputs={**prelaunch["input_hashes"],**protocol["input_hashes"],**protocol["source_hashes"]}
    assert all(not any(s in name for s in ("evaluation_large_","evaluation_larger_bending","all_bending_records")) for name in inputs)
    for name, digest in inputs.items():
        assert sha(ROOT/name)==digest==hashes.setdefault(name,digest),name
    results=bound(PILOT/"results.json",RESULTS,hashes)
    job=read(PILOT/"job.json")
    hashes[(PILOT/"job.json").relative_to(ROOT).as_posix()]=sha(PILOT/"job.json")
    assert job["results_sha256"]==RESULTS and results["bindings"]==job["bindings"]==expected
    assert results["source_hashes"]==protocol["source_hashes"] and results["input_hashes"]==protocol["input_hashes"]
    assert job["status"]==results["status"]=="COMPLETE_FOUR_CONDITIONAL_FIELDS"
    assert job["terminal_cases"]==results["terminal_case_receipts"] and len(job["terminal_cases"])==4
    assert all(results[k] is False for k in ("source_outcomes_read_or_scored","root_validated",
                                            "source_compatible_witness","novelty_established"))
    assert {k for k in protocol["claims"] if protocol["claims"][k]}==set()
    runtime=protocol["runtime"]
    torch.set_num_threads(runtime["torch_threads"])
    assert (platform.python_version(),np.__version__,scipy.__version__,torch.__version__,torch.get_num_threads()) == (
        runtime["python"],runtime["numpy"],runtime["scipy"],runtime["torch"],runtime["torch_threads"])
    cases={r["case_id"]:r for r in protocol["cases"]}
    assert set(cases)=={"baseline_small","baseline_large","joint_small","joint_large"}
    assert len({j["pid"] for j in job["terminal_cases"]})==4
    replays={}
    for j in job["terminal_cases"]:
        label=j["case_id"]
        assert job["started_at_utc"]<=j["started_at_utc"]<=j["ended_at_utc"]<=job["ended_at_utc"]
        replays[label]=replay_case(cases[label],j,results,protocol,hashes)
    Qs={label:tuple(map(Q,rows[-1]["Q_um"]["exact"])) for label,rows in replays.items()}
    assert all(x==tuple(map(Q,results["signed_Q_um"][k]["exact"])) for k,x in Qs.items())
    delta={name:signed_difference(Qs[name+"_large"],Qs[name+"_small"]) for name in ("baseline","joint")}
    assert all(x==tuple(map(Q,results["signed_diameter_change_large_minus_small_Q_um"][k]["exact"]))
               for k,x in delta.items())
    gap=signed_difference(Qs["joint_large"],Qs["baseline_large"])
    assert gap==tuple(map(Q,results["large_diameter_gap_test"]["joint_minus_baseline_large_Q_um"]["exact"]))
    assert gap[1]<0 and results["large_diameter_gap_test"]["status"]=="DISJOINT_STORED_CURVE_INTERVALS"
    interaction=signed_difference(delta["joint"],delta["baseline"])
    assert interaction[1]<0
    assert all(sha(ROOT/key)==value for key,value in hashes.items())
    return {"status":"PASS_ACTUAL_FOUR_SAVED_FIELD_REPLAY","replay_sha256":sha(Path(__file__)),
            "producer_results_sha256":RESULTS,"all_replayed_hashes":hashes,
            "os_terminal_receipts":4,"checkpoints":sum(map(len,replays.values())),
            "positive_stored_element_geometry_replays":144*sum(map(len,replays.values())),
            "numeric_force_samples":32,"cases":replays,"Q_um":{k:packed(v) for k,v in Qs.items()},
            "joint_minus_baseline_large_Q_um":packed(gap),
            "large_minus_small_Q_um":{k:packed(v) for k,v in delta.items()},
            "joint_minus_baseline_diameter_change_um":packed(interaction),
            "elapsed_seconds":time.perf_counter()-start,
            "scope":"Independent control flow, topology and physical-radius Lagrange inversion, and Fraction comparison. "
                    "Reuses reviewed Q8/Q12 evaluation, exact geometry, and positive-thickness kernels. Enclosures "
                    "belong to stored-coordinate P2 curves. No optimizer, held-out response read, validation of "
                    "equilibrium or integration error, continuum/stability/source claim, or public scoring."}


if __name__ == "__main__":
    packet=review()
    out=HERE/"actual_review.json"
    write(out,packet)
    print(json.dumps({"status":packet["status"],"receipt_sha256":sha(out), "inputs":len(packet["all_replayed_hashes"]),
                      "elapsed_seconds":packet["elapsed_seconds"],
                      "joint_minus_baseline_large_Q_um":packet["joint_minus_baseline_large_Q_um"]["display_um"],
                      "joint_minus_baseline_diameter_change_um":packet["joint_minus_baseline_diameter_change_um"]["display_um"]},
                     indent=2),flush=True)
