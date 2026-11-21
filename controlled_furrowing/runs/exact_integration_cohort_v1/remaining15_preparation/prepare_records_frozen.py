"""Prepare the 15 unproved cohort records without opening a numerical job."""
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np


OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import cohort_exact_pilot as P
import exact_cohort_context as C
import exact_cohort_integrator as E
import review_cohort_exact_record as D


BASE = D.BASE
OUTPUTS = ("protocol.json", "preflight.json", "coordinate_box.json")
EXCLUDED = {"optoGEF2_small_col_B", "optoGEF2_small_col_L"}
EXPECTED = {"optoGEF2_small_col_" + label for label in "CDEFGHIJKMNOPQR"}
PINS = {
    "scripts/cohort_exact_pilot.py": "a21b680117f236cc38efc46e086be1d9aff80cfb09ae613691b1642db8bed25e",
    "scripts/exact_cohort_context.py": "c22612a1de1ccfd34e8b0f0859f6413c4da09c765ba4945d060340e93e928f9f",
    "scripts/exact_cohort_integrator.py": "092dc51161e708cdba3791ec47a7f8e843dd0e1afc866b6f65f286b8657d1745",
    "scripts/review_cohort_exact_record.py": "67db748a5d40f31471ddd3d6c84fae280ff7ead650e6eb4e950a4766e54e0285",
    "evidence/exact_cohort_context_independent_review.json": "898fb6828c16fc7f0db410c9d677373cf1696db5d0d418c08ef0004a30e215c5",
    "evidence/exact_cohort_integrator_independent_review.json": "1d233c61e5b7283ddfeb856e206a064d28cb2731246f13bff716e9a4063dbfa8",
    "evidence/exact_integral_root_B_safe_independent_review.json": "4500e5b34039eafe2e8af742cad739b1133ab414521fd6ed604cf4787d80d561",
    f"{BASE}/L_defect_diagnostic/signed_predictor_L_v1/actual_independent_review.json":
        "ea1ebcb2f1788eaa3fd491cca6d519e8f40568e9cdc870919a9b32f9f4adc0a9",
    C.AUTHORITY: C.AUTHORITY_SHA256,
}


@dataclass(frozen=True)
class Layout:
    kind: str
    present: tuple
    other_files: tuple


def _relative(path):
    return Path(path).relative_to(ROOT).as_posix()


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _record_ids(rows):
    ids = [row["record_id"] for row in rows]
    D.require(len(ids) == len(set(ids)) == 17 and set(ids) == EXPECTED | EXCLUDED,
              "exact original seventeen record IDs and known B/L exclusions required")
    selected = sorted(x for x in ids if x not in EXCLUDED)
    D.require(len(selected) == len(set(selected)) == 15 and set(selected) == EXPECTED,
              "all and only the fifteen remaining records required")
    return selected


def _classify(present, other, regular=True):
    if regular and tuple(present) == OUTPUTS:
        return Layout("existing_complete_read_only", tuple(present), tuple(other))
    if regular and not present and not other:
        return Layout("absent_safe_to_prepare", (), ())
    return Layout("partial_or_foreign_preserve", tuple(present), tuple(other))


def _inspect(record):
    D.require(record in EXPECTED, "completed B/L or unrecognized IDs cannot be prepared")
    folder = P.outdir(record)
    D.require(folder.resolve().is_relative_to(ROOT) and not folder.is_symlink(), "canonical local record path required")
    entries = set(p.name for p in folder.iterdir()) if folder.is_dir() else set()
    present = tuple(name for name in OUTPUTS if name in entries)
    other = tuple(sorted(entries - set(OUTPUTS)))
    regular = (not folder.exists() or folder.is_dir()) and all(
        (folder / p).is_file() and not (folder / p).is_symlink() for p in present
    )
    return _classify(present, other, regular)


def _independent_Rnorm(R):
    high = Q(0)
    for row in R:
        entries = [float(x).as_integer_ratio() for x in row]
        denominator = max(b for _, b in entries)
        absolute = sum(abs(a) * (denominator // b) for a, b in entries)
        high = max(high, Q(absolute, denominator))
    return high


def _independent_incidence(ctx):
    force, hessian = [0] * 611, [0] * 611
    groups = [*ctx.target.a["triangles"], *ctx.target.a["surface_0_edges"], *ctx.target.a["surface_1_edges"]]
    for group in groups:
        free, _ = D._mesh_map(ctx, group)
        D.require(bool(free), "nonempty finite-element free projection required")
        for index in free:
            force[index] += 1
            hessian[index] += len(free)
    return force, hessian


def _tube(ctx, box):
    d = [(ctx.cell[key][0] - ctx.base[key], ctx.cell[key][1] - ctx.base[key]) for key in C.ORDER]
    mapping = {int(dof): j for j, dof in enumerate(ctx.free)}
    n_radii = box["coordinate_radii_exact"]
    centers = box["current_centers_exact"]
    reference = box["reference_nodes_exact"]
    hulls = box["predictor_coordinate_hulls_exact"]
    D.require(box["root_membership_assumed"] is False and Q(box["corrector_radius_exact"]) == ctx.old_radius
              and all(Q(x) == 0 for x in box["shift_exact"]) and box["free"] == list(map(int, ctx.free))
              and box["triangles"] == ctx.target.a["triangles"].tolist()
              and box["parameter_cell_sha256"] == ctx.row["parameter_cell_sha256"],
              "same unshifted explicit selected original tube required")
    fixed_axis = 0
    for k in range(650):
        n, axis = divmod(k, 2)
        x = Q.from_float(float(ctx.target.a["reference"][n, axis]))
        y, width = Q(centers[n][axis]), Q(n_radii[n][axis])
        p_lo, p_hi = map(Q, hulls[n][axis])
        D.require(Q(reference[n][axis]) == x, "changed exact reference coordinate")
        if k in mapping:
            j = mapping[k]
            direction = [Q.from_float(float(s)) for s in ctx.S[j]]
            amin = sum((s * bounds[0 if s >= 0 else 1] for s, bounds in zip(direction, d)), Q(0))
            amax = sum((s * bounds[1 if s >= 0 else 0] for s, bounds in zip(direction, d)), Q(0))
            D.require(y == x + Q.from_float(float(ctx.u[j])) and [amin, amax] == list(map(Q, box["predictor_hulls_exact"][j]))
                      and (p_lo, p_hi) == (y + amin, y + amax)
                      and y - width <= p_lo - ctx.old_radius <= p_hi + ctx.old_radius <= y + width
                      and width == max(-amin, amax) + ctx.old_radius,
                      "exact predictor extrema and full old-radius coordinate enclosure required")
        else:
            D.require(x == y == p_lo == p_hi and width == 0, "constrained coordinate changed")
            fixed_axis += int(axis == 0 and x == 0)
    return fixed_axis


def _profile(ctx):
    a = ctx.target.a
    scale = ctx.target.energy_scale
    layouts = {"shape": (64, 6), "shape_grad": (144, 64, 6, 2), "volume_weights": (144, 64),
               "surface_0_N": (12, 3), "surface_1_N": (12, 3),
               "surface_0_weights": (12, 12), "surface_1_weights": (12, 12)}
    policy = ctx.target.packet["parameters"]
    D.require(scale > 0 and all(a[k].shape == shape for k, shape in layouts.items())
              and all(np.min(a[k]) > 0 for k in ("volume_weights", "surface_0_weights", "surface_1_weights"))
              and policy["mu"] == 1 and policy["bulk"] == 5000
              and np.count_nonzero(a["active"]) == 24 and set(np.unique(a["active"])) == {0.0, 1.0},
              "retained positive normalized bulk/surface weights and physical model required")
    keys = ("reference", "free", "triangles", "shape", "shape_grad", "active", "qr", "energy_scale", "volume_weights")
    keys += tuple(sorted(k for k in a if k.startswith("surface_")))
    geometry = {k: sha256(a[k].tobytes(order="C")).hexdigest() for k in keys}
    return {"bulk_elements": 144, "surface_edges": 24, "bulk_sites": 9216, "surface_sites": 288,
            "free_dofs": 611, "active_bulk_elements": 24, "energy_scale_exact": scale,
            "normalization": "actual binary cached weights divided by positive exact binary energy scale; physical measure 2*pi/scale",
            "geometry_and_rule_sha256": D._digest(geometry)}


def _validate(record, before, hashes):
    folder = P.outdir(record)
    D.require(_inspect(record).kind == "existing_complete_read_only", "all expected immutable divisor/protocol/box outputs required")
    verified = {name: {"file": _relative(folder / name), "sha256": _sha(folder / name)} for name in OUTPUTS}
    values = {name: D._bound(item["file"], item["sha256"], hashes) for name, item in verified.items()}
    p, pre, box = (values[n] for n in OUTPUTS)
    for name, digest in p["source_and_input_hashes"].items():
        D._bound(name, digest, hashes, parse=False)
    ctx = C.load_context(record)
    D.require(all(p["source_and_input_hashes"][name] == h for name, h in ctx.input_hashes.items())
              and p["record_id"] == record and p["row"] == ctx.row and p["precision_bits"] == 640
              and p["parameter_cell_sha256"] == ctx.row["parameter_cell_sha256"]
              and Q(p["normal_radius_exact"]) == ctx.old_radius and p["root_asserted"] is False
              and p["coordinate_box_file"] == verified["coordinate_box.json"]["file"]
              and p["coordinate_box_sha256"] == verified["coordinate_box.json"]["sha256"]
              and pre["source_sha256"] == PINS["scripts/cohort_exact_pilot.py"]
              and pre["protocol_sha256"] == verified["protocol.json"]["sha256"]
              and pre["status"] == "PASS_PREPARATION_ONLY" and D._yes(pre["checks"])
              and pre["actual_body_or_force_job_run"] is False, "record preparation bindings are inconsistent")
    D.require(ctx.dose > 0 and ctx.drive == ctx.base["common_gain"] * ctx.dose
              and all(ctx.cell[key][0] <= ctx.base[key] <= ctx.cell[key][1] for key in C.ORDER),
              "exact source dose and model base in closed C required")
    eta, z, r, bound = D._old_premises(ctx, p, hashes)
    force_count, hessian_count = _independent_incidence(ctx)
    Rnorm = _independent_Rnorm(ctx.R)
    D.require(Rnorm == Q(p["R_infinity_norm_exact"]) > 0 and force_count == p["force_incidence"]
              and hessian_count == p["H_row_incidence"] and max(force_count) == p["max_force_incidence"]
              and max(hessian_count) == p["max_H_row_incidence"], "exact topology counts or binary R infinity norm differ")
    fbudget, hbudget = (Q(p[k]) for k in ("force_tail_budget_exact", "H_entry_tail_budget_exact"))
    divisor = Q(p["margin_divisor"])
    D.require(divisor == 10_000 and fbudget > 0 and hbudget > 0
              and Rnorm * max(force_count) * fbudget == Q(p["global_force_tail_budget_exact"]) == (1 - z) * r / divisor
              and Rnorm * max(hessian_count) * hbudget == Q(p["global_derivative_tail_budget_exact"]) == (1 - z) / divisor,
              "unchanged per-component exact tail allocations required")
    fixed_axis = _tube(ctx, box)
    row = ctx.row
    return {
        "record_id": record, "incoming_layout": asdict(before),
        "preparation_status": "PREPARED_VIA_FROZEN_DRIVER" if before.kind == "absent_safe_to_prepare"
                              else "VALIDATED_EXISTING_WITHOUT_REWRITE",
        "selected_context": row, "exact_base": ctx.base, "exact_dose": ctx.dose, "exact_drive": ctx.drive,
        "selected_corrector_radius_exact": r, "audited_radius_and_cached_eta_Z_provider": bound,
        "exact_cached_eta_upper": eta, "exact_cached_Z_upper": z,
        "exact_R_infinity_norm": Rnorm, "R_binary_sha256": sha256(ctx.R.tobytes(order="C")).hexdigest(),
        "S_binary_sha256": sha256(ctx.S.tobytes(order="C")).hexdigest(),
        "u_binary_sha256": sha256(ctx.u.tobytes(order="C")).hexdigest(),
        "force_tail_budget_exact": fbudget, "H_entry_tail_budget_exact": hbudget,
        "normalized_profile": _profile(ctx), "exact_box_coordinate_checks": 650, "fixed_axis_coordinates": fixed_axis,
        "files": verified, "independent_domain_prelaunch_review_required": True,
    }


def tests(authority):
    checks = {
        "blank_layout_only_is_new": _classify((), ()).kind == "absent_safe_to_prepare",
        "complete_layout_is_read_only": _classify(OUTPUTS, ("old_note.json",)).kind == "existing_complete_read_only",
        "partial_layout_is_preserved": _classify(OUTPUTS[:1], ()).kind == "partial_or_foreign_preserve",
        "foreign_layout_is_preserved": _classify((), ("domains.json",)).kind == "partial_or_foreign_preserve",
        "nonregular_summary_is_rejected": _classify(OUTPUTS, (), False).kind == "partial_or_foreign_preserve",
        "binary_Rnorm_toy_signed_entries": _independent_Rnorm([[0.5, -0.25], [0.125, -1.0]]) == Q(9, 8),
    }
    for name, action in (
        ("duplicate_authority_ID_rejected", lambda: _record_ids(authority["records"] + authority["records"][:1])),
        ("excluded_record_cannot_enter_prepare", lambda: _inspect("optoGEF2_small_col_L")),
        ("changed_driver_digest_rejected", lambda: D._bound("scripts/cohort_exact_pilot.py", "0" * 64, {}, parse=False)),
    ):
        try:
            action()
        except ValueError:
            checks[name] = True
        else:
            checks[name] = False
    D.require(D._yes(checks), "write-once batch and exact-context guards must pass")
    return checks


def run():
    hashes = {}
    authority = D._bound(C.AUTHORITY, C.AUTHORITY_SHA256, hashes)
    records = _record_ids(authority["records"])
    for name, digest in PINS.items():
        D._bound(name, digest, hashes, parse=False)
    for module, name in ((P, "scripts/cohort_exact_pilot.py"), (C, "scripts/exact_cohort_context.py"),
                         (E, "scripts/exact_cohort_integrator.py"), (D, "scripts/review_cohort_exact_record.py")):
        D.require(_sha(module.__file__) == PINS[name], "imported runtime source differs from authority")
    checks = tests(authority)
    incoming = {record: _inspect(record) for record in records}
    source = Path(__file__).resolve()
    protocol, preflight, manifest = (OUT / n for n in ("protocol.json", "preflight.json", "manifest.json"))
    snapshot = OUT / "prepare_records_frozen.py"
    D.require(not any(p.exists() for p in (protocol, preflight, manifest, snapshot)), "preserve write-once batch preparation")
    with snapshot.open("xb") as f:
        f.write(source.read_bytes())
    psha = D._write(protocol, {
        "status": "FROZEN_REMAINING15_ORIGINAL_CUBE_PREPARATION_ONLY", "source_file": _relative(source),
        "source_sha256": _sha(source), "source_and_authority_hashes": hashes,
        "selection": records, "excluded_completed": sorted(EXCLUDED),
        "incoming_layouts": {k: asdict(row) for k, row in incoming.items()},
        "write_policy": "Call the frozen driver.prepare exactly once for each empty record layout. Revalidate a complete existing "
                        "layout as read-only. Preserve every partial, foreign, or failed record as-is.",
        "gate": "No B or L preparation and no domains, force, cache derivative, global matrix or root job. "
                "A separate record-specific independent review must bind each protocol before domains run.",
    })
    receipts, rejected = [], {}
    forbidden = AssertionError("a numerical job was attempted during write-once record preparation")
    with (
        patch.object(P, "domains_work", side_effect=forbidden),
        patch.object(P.G, "certify_box", side_effect=forbidden),
        patch.object(E, "bulk_component", side_effect=forbidden),
        patch.object(E, "surface_component", side_effect=forbidden),
        patch.object(E.U, "uniform_cache", side_effect=forbidden),
        patch.object(E.U.C, "exact_bulk_maps", side_effect=forbidden),
        patch.object(E.core, "assemble_point", side_effect=forbidden),
        patch.object(E.core, "local_terms", side_effect=forbidden),
    ):
        for record in records:
            try:
                before = incoming[record]
                D.require(_inspect(record) == before, "record directory layout changed since inventory")
                D.require(before.kind != "partial_or_foreign_preserve", "partial or foreign output preserved")
                if before.kind == "absent_safe_to_prepare":
                    P.prepare(record)
                receipts.append(_validate(record, before, hashes))
            except Exception as error:
                rejected[record] = {"type": type(error).__name__, "reason": str(error),
                                    "incoming_layout": asdict(incoming[record]), "terminal_layout": asdict(_inspect(record))}
                print(json.dumps({"record_id": record, **rejected[record]}), flush=True)
    checks.update({
        "all15_unique_records_prepared_and_attested": [row["record_id"] for row in receipts] == records,
        "original_9750_exact_coordinate_checks": sum(r["exact_box_coordinate_checks"] for r in receipts) == 15 * 650,
        "same_normalized_mesh_and_rule_geometry_across_15": len({
            r["normalized_profile"]["geometry_and_rule_sha256"] for r in receipts}) == 1,
        "selected_F_K_recentered_and_actual_M_nested_review": all(
            any(row["record_id"].endswith("_" + letter)
                and row["selected_context"]["proof_center_state_file"].endswith("/proof_center_state.json")
                for row in receipts) for letter in ("F", "K")) and any(
                r["record_id"].endswith("_M") and "parametric_actual_M_independent_review" in
                r["audited_radius_and_cached_eta_Z_provider"] for r in receipts),
        "no_numerical_phase_called": True,
    })
    ok = not rejected and D._yes(checks)
    manifest_sha = D._write(manifest, {
        "status": "READY_FOR_INDEPENDENT_RECORD_DOMAIN_PRELAUNCH_REVIEWS" if ok else "PARTIAL_PRESERVED_PREPARATION",
        "protocol_sha256": psha, "frozen_source_sha256": _sha(source), "expected_record_ids": records,
        "excluded_completed_records": sorted(EXCLUDED), "records": receipts, "rejected_attempts": rejected,
        "utc": datetime.now(timezone.utc).isoformat(), "actual_new_body_or_integral_or_cache_jobs_run": False,
        "phase_authorization": "Record-specific domain launch review required. The L prelaunch gates do not apply.",
    })
    pre_sha = D._write(preflight, {"status": "PASS_REMAINING15_PREPARATION_ONLY" if ok else "PARTIAL_PREPARATION_RETAINED",
                                   "protocol_sha256": psha, "manifest_sha256": manifest_sha, "source_sha256": _sha(source),
                                   "checks": checks, "source_input_and_record_file_hashes": hashes,
                                   "actual_new_body_or_integral_or_cache_jobs_run": False})
    print(json.dumps({"status": "PASS" if ok else "PARTIAL", "selected": len(records), "prepared": len(receipts),
                      "source_sha256": _sha(source), "protocol_sha256": psha, "preflight_sha256": pre_sha,
                      "manifest_sha256": manifest_sha}), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(run())
