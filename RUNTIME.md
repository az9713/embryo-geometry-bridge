# Active image runtime and historical reproducibility

Active requirements are in `controlled_furrowing/requirements.txt`. They now pin **Pillow 13.0.0**. The image privacy checker rejects Pillow older than 13 before opening an image. Use an isolated environment for image processing and publication.

For a fresh Windows environment:

```powershell
python -m venv .venv
.venv/Scripts/python.exe -m pip install -r controlled_furrowing/requirements.txt
.venv/Scripts/python.exe scripts/check_public_privacy.py --output privacy-check.json
```

The requirements include the project's CUDA-specific Torch pin and its additional package index. CPU-only or different CUDA deployments need a separately validated environment specification. Python 3.13 was used for the present local comparison. Do not assume a change to Torch or other numerical packages is bit-for-bit compatible without checking it.

The actual local security comparison used a virtual environment with access to the unchanged numerical site packages and a **local Pillow 13.0.0 override**. This avoided changing the historical interpreter or reinstalling the numerical stack. Loading the override from the selected environment was verified. It is an image-runtime isolation/comparison, not a claim that every inherited package is newly audited.

The exact previous requirements bytes are retained in `reproducibility/environment-before-pillow13.record`, with SHA256 and version provenance in `reproducibility/ENVIRONMENT-PROVENANCE.json`. The record is historical documentation, not an active install manifest. The earlier interpreter with Pillow 11.0.0 was preserved for historical reproduction and is rejected by the active image privacy checker.

All 45 previously published raster files have identical decoded RGBA pixel hashes under Pillow 11 and 13. Matplotlib canvas pixels are identical; generated PNG/TIFF/JPEG/GIF checks pass. Scoped mathematical publication replay also passed under the new interpreter before migration. Retained checks are linked from the [executed research report](controlled_furrowing/research_directions/relaxed_slab_execution_v1/REPORT.html).

The separate low-severity Torch advisory was not addressed by this Pillow migration. Security advisories describe potential vulnerability conditions, not proof of exploitation. The site serves static HTML and assets; these Python libraries run in research/publication tools rather than visitors' browsers. No malicious exploit files were executed during compatibility verification.
