"""Exact force and Hessian width attribution for independently reviewed C/D endpoints."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys

import replay as T

REPLAY_SHA = "d60b0a952ea55938f70f9c67c85207a3208eb888a904e5b063f05ad3c4d975bd"
RECEIPTS = {
    "C":"b66b821666c74ee68d2c0105df3087c4bec9e3a408d92f14cfa88231c9fc7d08",
    "D":"03017ddf88803327d2bf24cb4d22c8b4957ae4969d999ed062c66924d1e5b411",
}
OUT = Path(__file__).resolve().parent
KEYS = ("off_axis","axis_vertex","axis_edge_physical","surface")


def width(pair):
    return pair[1] - pair[0]


def prep():
    T.need(sha256(Path(T.__file__).read_bytes()).hexdigest() == REPLAY_SHA, "current reviewed replay source")
    for key, digest in RECEIPTS.items():
        data = T.read(f"{OUT.relative_to(T.ROOT).as_posix()}/optoGEF2_small_col_{key}.json", digest)
        T.need(data["status"] == "PASS_DIRECT_REDUCTION_UNRESOLVED" and data["root_certified"] is False,
               "already independently reviewed direct endpoints")
    h = T.write(OUT/"CD_width_protocol.json", {"status":"FROZEN_ATTRIBUTION_PREP_ONLY",
        "source_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),"independent_review_pins":RECEIPTS,
        "method":"Partition local interval vectors and symmetric Hessian contributions by reviewed chart/part; add each group exactly. Signed preconditioned interval widths add at every aggregate force row. Absolute Hessian row envelopes of the separated groups are attribution diagnostics.",
        "component_payloads_opened_for_attribution":False})
    print(h)


def run(protocol_sha):
    proto = T.read(f"{OUT.relative_to(T.ROOT).as_posix()}/CD_width_protocol.json",protocol_sha)
    T.need(proto["source_sha256"] == sha256(Path(__file__).read_bytes()).hexdigest()
           and proto["independent_review_pins"] == RECEIPTS
           and sha256(Path(T.__file__).read_bytes()).hexdigest() == REPLAY_SHA,"frozen attribution")
    sys.path.insert(0,str(T.ROOT/"scripts"))
    import exact_cohort_context as C
    final = {}
    for key, digest in RECEIPTS.items():
        record = f"optoGEF2_small_col_{key}"
        data = T.read(f"{OUT.relative_to(T.ROOT).as_posix()}/{record}.json",digest)
        external = data["producer_pins"]
        admission = T.read(T.PREP+f"/{record}_input_admission.json",external["admission_sha256"])
        prod = T.read(T.PREP+f"/{record}/results.json",external["results_sha256"])
        f = {name:[T.ZERO]*611 for name in KEYS}
        H = {name:{} for name in KEYS}
        for part, manifest in admission["expected_packets"].items():
            for meta in manifest:
                vp = T.read(meta["file"],meta["sha256"],payload=True)
                group = "surface" if part == "surface" else vp["chart"]
                T.scatter(f[group],H[group],vp["global_indices"],
                          vp["force_integral_minus_cache"],vp["H_integral_minus_cache_upper"])
        ctx = C.load_context(record)
        exact_R = [[Q.from_float(float(x)) for x in row] for row in ctx.R]
        rhs = {name:[width(v) for v in f[name]] for name in KEYS}
        prewidth = {} 
        for name in KEYS:
            prewidth[name] = [sum(abs(r)*v for r,v in zip(row,rhs[name])) for row in exact_R]
        original = [tuple(map(Q,v)) for v in prod["preconditioned_RF_intervals_exact"]]
        width_i = max(range(611),key=lambda i:width(original[i]))
        mag_i = max(range(611),key=lambda i:T.magnitude(original[i]))
        T.need(all(sum(prewidth[name][i] for name in KEYS) == width(v) for i,v in enumerate(original)),
               "exact RF endpoint width partition")
        Hwidths = {}
        for name in KEYS:
            hs = [Q(0)]*611
            for (i,j),v in H[name].items():
                for k in {i,j}:
                    hs[k] += width(v)
            Hwidths[name] = max(sum(abs(r)*v for r,v in zip(row,hs)) for row in exact_R)
        report = {"independent_review_sha256":digest,
            "RF_rows_excluding_zero":sum(not lo<=0<=hi for lo,hi in original),
            "RF_widest_row":width_i,"RF_largest_magnitude_row":mag_i,
            "largest_group_by_max_RF_width":max(KEYS,key=lambda g:max(prewidth[g])),
            "max_width_by_group":{k:max(v) for k,v in prewidth.items()},
            "width_by_group_at_widest_RF_row":{k:v[width_i] for k,v in prewidth.items()},
            "width_share_by_group_at_widest_RF_row":{k:v[width_i]/width(original[width_i]) for k,v in prewidth.items()},
            "width_by_group_at_largest_magnitude_RF_row":{k:v[mag_i] for k,v in prewidth.items()},
            "max_preconditioned_symmetric_H_row_width_by_group":Hwidths,
            "force_defect_allowance_with_retained_H":(1-Q(prod["total_Z_upper_exact"]))*Q(prod["old_radius_exact"])-Q(prod["old_cached_eta_upper_exact"]),
            "signed_residual_allowance_with_retained_H":(1-Q(prod["total_Z_upper_exact"]))*Q(prod["old_radius_exact"])}
        final[record] = report
        print(key, "excludes zero",report["RF_rows_excluding_zero"], "max RF width by group",
              {k:float(v) for k,v in report["max_width_by_group"].items()}, "shares",
              {k:float(v) for k,v in report["width_share_by_group_at_widest_RF_row"].items()}, "H group widths",
              {k:float(v) for k,v in Hwidths.items()},flush=True)
    h = T.write(OUT/"CD_width_results.json",{"status":"PASS_REVIEWED_DIRECT_ENDPOINT_WIDTH_ATTRIBUTION",
        "protocol_sha256":protocol_sha,"records":final,
        "scope":"Widths and proof margins only, from existing component interval endpoints. Numerical root status remains unresolved by the retained direct sufficient test."})
    print(h)


if __name__ == "__main__":
    prep() if len(sys.argv) == 1 else run(sys.argv[1])

