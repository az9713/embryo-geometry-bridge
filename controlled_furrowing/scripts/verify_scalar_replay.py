"""Replay the frozen predictor with evaluation outcomes physically absent."""
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import uuid

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    replay=ROOT/"runs/independent_scalar_replay"/uuid.uuid4().hex
    replay.mkdir(parents=True,exist_ok=True)
    files=["PROTOCOL.json","scripts/scalar_reproduction.py",
           "data/derived/calibration_small_18um.csv","data/derived/large_inputs_only.json"]
    for relative in files:
        destination=replay/relative
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT/relative,destination)
    responses_absent=not (replay/"data/derived/evaluation_large_50um.csv").exists()
    frozen=replay/"runs/frozen_scalar_predictions.json"
    if frozen.exists():
        # A replay already frozen under the same source is read, never overwritten.
        completed=None
    else:
        completed=subprocess.run([sys.executable,str(replay/"scripts/scalar_reproduction.py"),"predict"],capture_output=True,text=True,timeout=120)
        (replay/"execution.log").write_text(completed.stdout+completed.stderr,encoding="utf-8")
        if completed.returncode:
            raise RuntimeError("Independent replay process failed")
    original=ROOT/"runs/frozen_scalar_predictions.json"
    same=sha(frozen)==sha(original)
    prediction=json.loads(original.read_text(encoding="utf-8"))
    with (ROOT/"data/derived/evaluation_large_50um.csv").open(encoding="utf-8",newline="") as stream:
        observed=[float(row["maximum_shortening_um"]) for row in csv.DictReader(stream)]
    score=json.loads((ROOT/"evidence/scalar_scores.json").read_text(encoding="utf-8"))
    metric_checks={}
    for name in ("primary_mechanical","OLS","zero_shortening"):
        residual=[p-y for p,y in zip(prediction[name]["evaluation_prediction_um"],observed)]
        calculated={"RMSE_um":math.sqrt(math.fsum(x*x for x in residual)/len(residual)),
                    "MAE_um":math.fsum(abs(x) for x in residual)/len(residual),
                    "signed_bias_um":math.fsum(residual)/len(residual)}
        metric_checks[name]=all(abs(calculated[k]-score["scores"][name][k])<1e-12 for k in calculated)
    checks={"evaluation_responses_physically_absent_from_prediction_replay":responses_absent,
            "prediction_bytes_identical":same,
            "all_metrics_independently_match":all(metric_checks.values()),
            "no_geometric_output_claim":prediction["geometric_prediction"] is None}
    result={"status":"PASS" if all(checks.values()) else "FAIL", "checks":checks,
            "metric_checks":metric_checks,"prediction_sha256":sha(original),
            "replay_input_files":{f:sha(replay/f) for f in files},
            "limitations":"Same software environment, isolated data inputs; not analyst blinding, external expert review or biological validation",
            "source_sha256":sha(Path(__file__))}
    (ROOT/"evidence/scalar_independent_replay.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__=="__main__":
    main()
