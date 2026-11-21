"""Read-only published Movie3 geometry/calibration audit.

Run from controlled_furrowing or embryo-geometry-bridge. No model code is
imported, no fit is performed, and no source files are modified. Source-image
diagnostics in this run directory are local-only and are not project figures.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import sys
import xml.etree.ElementTree as ET

import cv2
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "runs/movie_geometry"
OUT = ROOT / "evidence/movie_geometry_audit.json"
MOVIE = ROOT / "data/raw/41467_2025_62483_MOESM7_ESM.avi"
FIGURE = RUN / "source_fig4_unmodified.jpg"
TIMES = [0, 21, 43, 64, 86, 108, 129, 151, 172, 194, 216, 237,
         259, 280, 302, 324, 345, 367, 388, 410, 432, 453, 475, 496, 518, 540]


def digest(path):
    data = path.read_bytes()
    return {"path": path.relative_to(ROOT).as_posix(), "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "md5": hashlib.md5(data).hexdigest()}


def correlation(a, b):
    return float(np.corrcoef(a.ravel(), b.ravel())[0, 1])


def main():
    sources = [MOVIE, FIGURE, ROOT / "data/raw/primary.xml",
               ROOT / "data/raw/PMC12361397.1.json",
               ROOT / "data/derived/additional_files_description.txt",
               ROOT / "data/derived/evaluation_large_50um.csv"]
    hashes_before = [digest(p) for p in sources]
    assert hashes_before[0]["sha256"] == "260555e647c56a74b157ef624fe263568cc47fbda553d3e8aa7fe0ad05e14be7"
    assert hashes_before[1]["sha256"] == "c587150d86e7326b9774643410b1242b61053326f194cab8e0b91ce61dfe6473"
    official = json.loads((ROOT / "data/raw/PMC12361397.1.json").read_text())
    checksum_checks = []
    for path in [MOVIE, FIGURE]:
        original_name = "41467_2025_62483_Fig4_HTML.jpg" if path == FIGURE else path.name
        candidates = [url for url in official["media_urls"] if original_name in url]
        assert len(candidates) == 1
        official_md5 = candidates[0].split("md5=")[1]
        actual = digest(path)["md5"]
        assert actual == official_md5
        checksum_checks.append({"file": path.name, "official_md5": official_md5,
                                "actual_md5": actual, "pass": True})

    capture = cv2.VideoCapture(str(MOVIE))
    fps = float(capture.get(cv2.CAP_PROP_FPS))
    frames = []
    while True:
        ok, frame = capture.read()
        if not ok:
            break
        frames.append(frame)
    capture.release()
    assert len(frames) == 26 and frames[0].shape == (826, 1908, 3)
    gray = [cv2.cvtColor(f, cv2.COLOR_BGR2GRAY) for f in frames]
    assert TIMES == [int(np.floor(21.6 * i)) for i in range(26)]
    frame_records = [{"index": i, "burned_in_seconds_manually_transcribed": TIMES[i],
                      "decoded_bgr_shape": list(f.shape),
                      "decoded_bgr_sha256": hashlib.sha256(f.tobytes()).hexdigest()}
                     for i, f in enumerate(frames)]

    # Preselected display region containing only the white 20-um scale bar.
    # Pixel geometry here is display geometry, not microscope metadata.
    bars = []
    for threshold in [128, 200, 230, 250]:
        binary = (gray[0][775:810, 1180:1380] >= threshold).astype(np.uint8)
        n, labels, stats, centroids = cv2.connectedComponentsWithStats(binary, 8)
        idx = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
        x, y, w, h, area = [int(v) for v in stats[idx]]
        bars.append({"threshold": threshold, "x": x + 1180, "y": y + 775,
                     "width": w, "height": h, "area": area})
    assert all(b["width"] == 160 and b["height"] == 15 for b in bars)

    row_diagnostics = [{"source_y": y,
                        "mean_intensity": float(gray[0][y, 20:576].mean()),
                        "max_intensity": int(gray[0][y, 20:576].max())}
                       for y in [168, 169, 170, 171, 172, 535, 536, 537, 538]]
    figure_gray = cv2.imread(str(FIGURE), cv2.IMREAD_GRAYSCALE)
    crops = {"c0": [326, 14, 460, 95], "c2": [326, 97, 460, 178],
             "c8": [326, 180, 460, 261], "d0": [509, 25, 627, 143],
             "d2": [509, 145, 627, 261]}
    matches = {}
    for label, (x0, y0, x1, y1) in crops.items():
        template = figure_gray[y0:y1, x0:x1].astype(float)
        roi = [20, 170, 576, 538] if label.startswith("c") else [675, 170, 1232, 728]
        u0, v0, u1, v1 = roi
        values = [correlation(template, cv2.resize(g[v0:v1, u0:u1],
                  (template.shape[1], template.shape[0])).astype(float)) for g in gray]
        order = sorted(range(26), key=lambda i: -values[i])
        matches[label] = {"figure_crop_xyxy": [x0, y0, x1, y1],
                          "movie_roi_xyxy": roi, "pearson_ncc_all_frames": values,
                          "best_three": [{"frame": i, "ncc": values[i]} for i in order[:3]]}

    x0, y0, x1, y1 = crops["d2"]
    template = figure_gray[y0:y1, x0:x1].astype(float)
    high_template = template - cv2.GaussianBlur(template, (0, 0), 2)
    hp = {"middle": [], "right_panel_negative_control": []}
    for g in gray:
        for key, roi in [("middle", [675, 170, 1232, 728]),
                         ("right_panel_negative_control", [1330, 170, 1887, 728])]:
            u0, v0, u1, v1 = roi
            a = cv2.resize(g[v0:v1, u0:u1], (template.shape[1], template.shape[0])).astype(float)
            hp[key].append(correlation(high_template, a - cv2.GaussianBlur(a, (0, 0), 2)))
    assert int(np.argmax(hp["middle"])) == 5
    assert hp["middle"][5] > 0.6
    assert max(hp["right_panel_negative_control"]) < 0.1

    root = ET.parse(ROOT / "data/raw/primary.xml").getroot()
    source_locators = {}
    for sec in root.iter("sec"):
        if sec.get("id") in ["Sec7", "Sec14", "Sec15", "Sec17"]:
            source_locators[sec.get("id")] = " ".join(sec.itertext())
    fig4 = next(fig for fig in root.iter("fig") if fig.get("id") == "Fig4")
    source_locators["Fig4_caption"] = " ".join(fig4.itertext())
    description = (ROOT / "data/derived/additional_files_description.txt").read_text(encoding="utf8")
    start = description.index("File Name: Supplementary  Movie 3")
    stop = description.index("File Name: Supplementary  Movie 4")
    source_locators["Supplementary_Movie3_description"] = description[start:stop]

    with (ROOT / "data/derived/evaluation_large_50um.csv").open(encoding="utf8", newline="") as handle:
        records = list(csv.DictReader(handle))
    # Find the actual canonical dose-column name without assuming a specimen ID.
    dose_key = next(k for k in records[0] if k in ["dose", "activation_magnitude"])
    rounded_matches = [r for r in records if round(float(r[dose_key]), 2) == 1.02]
    assert len(rounded_matches) == 1
    assert [digest(p) for p in sources] == hashes_before

    result = {
        "audit_version": 1,
        "status": "PASS provenance/display checks; INSUFFICIENT certified anatomical geometry",
        "scope": "One published illustrative Movie3; no cohort substitution, no model import or fitting",
        "software": {"python": platform.python_version(), "opencv": cv2.__version__, "numpy": np.__version__},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_hashes": hashes_before, "official_checksum_checks": checksum_checks,
        "source_locator_text_hashes": {
            locator: {"utf8_bytes": len(text.encode("utf8")),
                      "sha256_of_extracted_text": hashlib.sha256(text.encode("utf8")).hexdigest()}
            for locator, text in source_locators.items()},
        "container": {"decoded_frames": 26, "width": 1908, "height": 826,
                      "playback_fps": fps, "biological_time_is_not_playback_time": True,
                      "ffmpeg_log": "runs/movie_geometry/ffmpeg_decode_local.txt",
                      "display_sample_aspect": "1:1 container metadata; not proof of biological reslice calibration"},
        "frames": frame_records,
        "timing": {"method": "Manual transcription of all burned-in timestamps; not raw acquisition timestamps",
                   "formula_check": "All 26 equal floor(21.6*frame_index)",
                   "source_nominal_interval_seconds": 22,
                   "precise_activation_origin_and_subsecond_timestamps": "not provided"},
        "scale_bar": {"caption_micrometers": 20, "threshold_checks": bars,
                      "nominal_horizontal_um_per_display_pixel": 20 / 160,
                      "assumed_endpoint_allowance_pixels_each": 1,
                      "conditional_width_interval_pixels": [158, 162],
                      "conditional_horizontal_um_per_pixel_interval": [20 / 162, 20 / 158],
                      "interval_scope": "Raster endpoint allowance only; physical calibration and panel transform uncertainty not bounded",
                      "per_panel_xy_and_transverse_xz_isotropy": "consistent with source display, not explicitly documented"},
        "crop_diagnostics": {"transverse_roi_xyxy_half_open": [20, 170, 576, 538],
                             "frame0_rows": row_diagnostics,
                             "visual_inspection_frames": [0, 13],
                             "finding": "Lateral membrane fluorescence reaches upper display crop; continuous apical boundary not visible. Crop edge is not anatomy."},
        "figure_correspondence": {"direct_resize_is_not_physical_calibration": True,
                                  "unfiltered": matches, "highpass_sigma_resized_pixels": 2,
                                  "highpass_d2_all_frames": hp,
                                  "finding": "Image correspondence supports common published field. No acquisition/specimen ID crosswalk supplied."},
        "cohort_identity": {"caption_dose": 1.02, "rounded_dose_candidates": rounded_matches,
                            "classification": "unique rounded-dose inference, not verified specimen ID"},
        "plane_convention": {"caption": "2 micrometers below apical surface, middle; 30 micrometers, right",
                             "initial_vs_current_global_reference": "not specified",
                             "pointwise_surface_following_plane": "not assumed; not a single horizontal plane on curved tissue",
                             "Sec17_tracking_is_separate_myosin_z_spreading_procedure": True},
        "geometry_decision": {
            "full_initial_apical_basal_height": "NO-GO: apical censoring and missing z calibration/registration",
            "peak_height_or_shortening_profile": "NO-GO: no uncensored time-matched surface pair or independently known peak",
            "rigorous_biological_censored_micrometer_bound": "NO-GO: no validated basal annotation/z calibration error bound",
            "basal_display_proxy": "Potential exploratory expert annotation; not extracted as certified anatomy",
            "annulus_Ai_Ao_with_certified_uncertainty": "NO-GO from this raster alone: no resolved unique inner contour, PSF or semantic occupancy validation",
            "annulus_image_proxy_future_work": "Possible conditional fixed annotation/segmentation sensitivity study, explicitly not anatomical ground truth",
            "source_qualitative_bending": "Supported by source Sec7 and visible central plug/dark moat; fluorescence alone does not certify solid occupancy",
            "no_model_parameters_fitted": True,
            "no_source_frames_for_public_site": True,
            "one_movie_does_not_validate_full_cohort_or_identified_set": True},
        "sources_unchanged_before_after": True,
    }
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf8")
    print(json.dumps({"status": result["status"], "frames": len(frames),
                      "bar_width": bars[0]["width"], "d2_highpass_frame5_ncc": hp["middle"][5],
                      "receipt": str(OUT.relative_to(ROOT)), "source_unchanged": True}))


if __name__ == "__main__":
    main()
