"""Exact additive attribution of the final L preconditioned force interval."""
from fractions import Fraction as Q
import json
from pathlib import Path
import sys

from replay import ROOT, STAGE, OUT, ACTUAL_PINS, N, ZERO, canonical, inter, pack, record_inputs, write_once

PIN = "54dee8659601c387712c16732ab6aed90fbf313c1e3b0f77930b3e12b19220e6"
ORIGINAL = "runs/exact_integration_cohort_v1/records/optoGEF2_small_col_L/original_cube_v1"
OLD_ADMISSION = "runs/exact_integration_cohort_v1/reduction_preparation/optoGEF2_small_col_L_input_admission.json"


def sum_ranges(ranges):
    return tuple(sum((v[i] for v in ranges),Q(0)) for i in (0,1))


def width(pair):
    return pair[1]-pair[0]


def magnitude(pair):
    return max(map(abs,pair))


def run():
    if not __debug__:
        raise RuntimeError("Replay requires enabled assertions")
    hashes,read = record_inputs()
    final = read(STAGE/"reduction_preparation/actual_independent_review.json",PIN)
    assert final["status"] == "PASS_ACTUAL_SCALED_ORIGINAL_CUBE_REDUCTION_REVIEW"
    for path,digest in final["input_hashes"].items():
        read(path,digest,False)
    admission = read(OLD_ADMISSION,hashes[OLD_ADMISSION])
    manifest = read(STAGE/"actual/results.json",ACTUAL_PINS["actual/results.json"])
    scaled_by_elem = {p["element"]:p for p in manifest["files"]}
    groups = {g:[ZERO]*N for g in ("scaled_offaxis","retained_axis","retained_surface")}
    locals_by_group = {k:[] for k in groups}
    for part in ("bulk","surface"):
        for row in admission["expected_packets"][part]:
            old = read(row["file"],row["sha256"])
            new = None
            if part == "bulk":
                i = old["element"]
                owner = "scaled_offaxis" if i in scaled_by_elem else "retained_axis"
                if i in scaled_by_elem:
                    selected = scaled_by_elem[i]
                    new = read(selected["file"],selected["sha256"])
            else:
                owner = "retained_surface"
            force = list(map(pack,old["force_integral_minus_cache"]))
            if new is not None:
                force = [inter(a,b) for a,b in zip(force,map(pack,new["force_integral_minus_cache"]))]
            ids = old["global_indices"]
            assert len(force) == len(ids)
            for i,value in zip(ids,force):
                groups[owner][i] = sum_ranges((groups[owner][i],value))
            locals_by_group[owner].append((row["file"],ids,force))
    sys.path.insert(0,str(ROOT/"scripts"))
    import review_cohort_exact_record as D
    context = D.C.load_context("optoGEF2_small_col_L")
    for path,digest in context.input_hashes.items():
        read(path,digest,False)
    assert [len(v) for v in locals_by_group.values()] == [132,12,24]
    R = context.R
    projected = {name:D._dyadic_reduce(R,g,[Q(0)]*N)[0] for name,g in groups.items()}
    RF = [sum_ranges(rows) for rows in zip(*projected.values())]
    assert canonical(RF) == final["result"]["RF_endpoint_sha256"]
    assert max(map(magnitude,RF)) == Q(final["result"]["RF_uniform_upper_exact"])
    chosen = {"maximum_absolute_upper": max(range(N),key=lambda i:magnitude(RF[i])),
              "maximum_width":max(range(N),key=lambda i:width(RF[i]))}
    rows = {}
    for key,i in chosen.items():
        weights = [Q.from_float(float(v)) for v in R[i]]
        summaries = []
        for group,values in projected.items():
            cells = []
            for path,ids,force in locals_by_group[group]:
                factors = [((weights[j]*a,weights[j]*b) if weights[j] >= 0 else (weights[j]*b,weights[j]*a))
                           for j,(a,b) in zip(ids,force)]
                bound = sum_ranges(factors)
                cells.append((path,bound))
            assert sum_ranges([v for _,v in cells]) == values[i]
            leaders = sorted(cells,key=lambda item:(-width(item[1]),item[0]))[:5]
            summaries.append({"group":group,"signed_enclosure":values[i],"width_exact":width(values[i]),
                              "width_share_exact":width(values[i])/width(RF[i]),
                              "largest_width_packets":[{"file":p,"signed_interval":v,"width_exact":width(v)}
                                                      for p,v in leaders]})
        assert sum(Q(s["width_exact"]) for s in summaries) == width(RF[i])
        rows[key] = {"free_index_zero_based":i,"global_coordinate":int(context.free[i]),
                     "signed_RF_interval":RF[i],"groups":summaries}
    doc = {"status":"PASS_EXACT_SIGNED_FORCE_ATTRIBUTION","final_reduction_review_sha256":PIN,
           "scope":"Exact width attribution of the retained preconditioned force enclosures. Interval width measures enclosure uncertainty.",
           "rows":rows,"all_RF_intervals_contain_zero":all(a <= 0 <= b for a,b in RF),
           "input_hashes":dict(sorted(hashes.items()))}
    digest = write_once(OUT/"force_attribution.json",doc)
    print(json.dumps({"status":doc["status"],"sha256":digest,"rows":{
         key:{"index":r["free_index_zero_based"],"width":float(width(r["signed_RF_interval"])),
         "shares":{s["group"]:float(s["width_share_exact"]) for s in r["groups"]},
         "largest":{s["group"]: (Path(s["largest_width_packets"][0]["file"]).name,
                           float(s["largest_width_packets"][0]["width_exact"])) for s in r["groups"]}}
         for key,r in rows.items()}}),flush=True)


if __name__ == "__main__":
    run()
