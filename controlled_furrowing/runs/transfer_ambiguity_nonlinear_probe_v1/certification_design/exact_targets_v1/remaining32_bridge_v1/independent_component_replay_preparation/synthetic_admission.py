"""Synthetic metadata-only test of the frozen first-wave admission."""
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import replay as R


def attempt(root):
    with TemporaryDirectory(prefix="dummy_seals_", dir=root) as tmp, patch.object(R, "RAW", Path(tmp)):
        phase = {"dummy_bindings": "synthetic-only"}
        supervisor = 1234
        plans = [{"case_id": case} for case in R.FIRST]
        targets = {case: {"case_id": case} for case in R.FIRST}
        auth = {"targets": targets, "plans": plans}
        parent = {"pid": supervisor, "stage": R.STAGE, "jobs": [list(job) for job in R.JOBS], "bindings": phase}
        pins = {}
        cases = []
        for idx, case in enumerate(R.FIRST):
            parts = []
            for j, part in enumerate(R.COUNTS):
                folder = R.RAW / "cases" / case / part
                folder.mkdir(parents=True)
                meta = {"pid": 4050 + 2 * idx + j, "supervisor_pid": supervisor, "stage": R.STAGE,
                        "case_id": case, "part": part, "bindings": phase}
                rows = []
                for slot in range(R.COUNTS[part]):
                    path = folder / f"{part}_{slot:03d}.json"
                    digest = R.save(path, {"dummy_numeric_fields": False})
                    pins[R.name(path)] = digest
                    rows.append({"kind": part, "slot": slot, "file": path.name, "sha256": digest})
                counts = {"components": len(rows)}
                m_sha = R.save(folder / "manifest.json", {"status": "COMPLETE_RAW_REMAINING32_SPATIAL_PART_AWAITS_INDEPENDENT_REPLAY",
                               **meta, "target": targets[case], "plan": plans[idx], "components": rows,
                               "counts": counts, "point_force_scope": R.POINT, "whole_H_scope": R.BOX,
                               "radius_exact": str(R.RADIUS), "root_asserted": False,
                               "packet_independent_replay_required": True})
                t_sha = R.save(folder / "terminal.json", {"status": "FINISHED", "exit_code": 0, "results_sha256": m_sha, **meta})
                obs = {**meta, "observed_exit_code": 0, "terminal_sha256": t_sha}
                o_sha = R.save(folder / "observed_os_exit.json", obs)
                parts.append({**obs, "manifest_sha256": m_sha, "observed_os_exit_sha256": o_sha})
                pins.update({R.name(folder / f"{k}.json"): h for k, h in
                             (("manifest", m_sha), ("terminal", t_sha), ("observed_os_exit", o_sha))})
                for k, doc in (("launched", meta), ("started", {"status": "RUNNING", **meta})):
                    path = folder / f"{k}.json"
                    pins[R.name(path)] = R.save(path, doc)
            cases.append({"case_id": case, "target": targets[case], "parts": parts})
        top = R.save(R.RAW / "manifest.json", {"status": "COMPLETE_FOUR_RAW_SPATIAL_BRIDGES_AWAITING_INDEPENDENT_PACKET_REPLAY",
                     **parent, "cases": cases, "raw_components": 672, "root_asserted": False,
                     "later_waves_launched": False, "independent_packet_replay_required": True})
        pins.update({R.name(R.RAW / name): digest for name, digest in
                     (("manifest.json", top), ("started.json", R.save(R.RAW / "started.json", {"status": "RUNNING", **parent})),
                      ("terminal.json", R.save(R.RAW / "terminal.json", {"status": "FINISHED", "exit_code": 0,
                                                                        "manifest_sha256": top, **parent})))})
        p = R.RAW / "root_terminal_admission.json"
        ext = {"status": "PASS_ROOT_OBSERVED_REMAINING32_FIRST_WAVE_SPATIAL_TERMINALS", "pid": supervisor,
               "unified_exec_session": 4567, "observed_exit_code": 0, "sha256": pins}
        def check(value):
            if p.exists():
                p.unlink()
            checked = {}
            rows, observed = R.full_actual(R.save(p, value), phase, auth, checked)
            return len(rows) == 8 and len(checked) == 716 and observed == value
        def reject(value, message):
            try:
                check(value)
            except ValueError as error:
                return message in str(error)
            return False
        good = check(ext)
        path = R.RAW / "cases" / R.FIRST[0] / "bulk" / "bulk_000.json"
        miss = dict(pins)
        miss.pop(R.name(path))
        short = reject({**ext, "sha256": miss}, "672 raw packet SHA pins required")
        fake_zero = reject({**ext, "observed_exit_code": False}, "external real complete")
        bad_parent = reject({**ext, "pid": supervisor+1}, "supervisor identity")
        path.write_text('{"corrupted_inert_fixture": true}', encoding="utf-8")
        bad_bytes = reject(ext, "compiled evidence changed")
        return {"fully_typed_715_pins_admit": good, "a_missing_packet_holds": short,
                "boolean_native_exit_rejected": fake_zero, "wrong_observed_supervisor_rejected": bad_parent,
                "post_admission_changed_bytes_rejected": bad_bytes}


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    checks = attempt(here)
    R.require(R.yes(checks), "synthetic complete raw-admission closure test")
    print(R.save(here / "synthetic_admission_readback.json", {"status": "PASS",
          "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
          "replay_source_sha256": sha256(R.SRC.read_bytes()).hexdigest(), "checks": checks,
          "actual_raw_packets_opened": 0, "mechanical_jobs": 0}))
