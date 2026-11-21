"""Archive the prior publication checkpoint, then link the new verified phase.

This updates reporting metadata only, never numerical certificates or PLAN.md.
"""
from pathlib import Path
import hashlib,json,shutil

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent
ARCHIVE=ROOT/'runs/checkpoint_history/point_root_phase'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def change(path,old,new):
    text=path.read_text(encoding='utf-8')
    if new in text:return
    if old not in text:raise RuntimeError('expected text absent:'+str(path))
    path.write_text(text.replace(old,new,1),encoding='utf-8')

if not ARCHIVE.exists():
    files=['docs/CURRENT-FIVE-GATE-AUDIT.md','evidence/current_five_gate_independent_review.json',
           'evidence/gates.json','evidence/final_component_binding_check.json',
           'evidence/readiness_snapshot_verification.json','scripts/build_report.py']
    manifest={}
    for rel in files:
        source=ROOT/rel
        target=ARCHIVE/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target);manifest[rel]=sha(target)
    for rel in ['README.md','scripts/build_dossier.py']:
        target=ARCHIVE/'parent'/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(PARENT/rel,target);manifest['parent/'+rel]=sha(target)
    (ARCHIVE/'archive_manifest.json').write_text(json.dumps({'status':'ARCHIVED_PRIOR_CHECKPOINT_BYTES',
        'files':manifest,'scope':'Historical exact checkpoint, not current publication verification.'},indent=2)+'\n',encoding='utf-8')
else:
    manifest=json.loads((ARCHIVE/'archive_manifest.json').read_text(encoding='utf-8'))['files']
    if any(sha(ARCHIVE/rel)!=h for rel,h in manifest.items()):raise RuntimeError('prior archive changed')

path=ROOT/'evidence/gates.json';g=json.loads(path.read_text(encoding='utf-8'))
g['snapshot']='Uniform discrete-root and body-geometry families plus actual-M conditional cell exclusion are independently reviewed. All17 shared-cell uniform root/observer extension is being prepared. Full material/geometry inference, source-matched temporal observables, continuum accuracy and revised larger-region shape predictions remain incomplete.'
gate=g['gates'][2]
addition=' A new old-material weak-load gain family has uniform energy-norm root/local-minimum and entire-body geometry certificates. A distinct actual-M five-coordinate box has a uniform stationary root/apical observer certificate and is conditionally excluded below approximately5.18638714micrometres apical discrepancy. This is one selected branch, not all68 coverage or full-domain exclusion. The reviewed inference layer retains all constraints and unresolved finite leaves, other branches, zero faces and unbounded tails.'
if addition not in gate['scientific_outcome']:gate['scientific_outcome']+=addition
for rel in ['docs/PARAMETER-INFERENCE-CHECKPOINT.md','docs/ROOT-BODY-GEOMETRY.md',
            'docs/ENERGY-GAIN-CELL.md','docs/ENERGY-GAIN-GEOMETRY.md','docs/PARAMETRIC-ROOT-CELLS.md',
            'docs/ACTUAL-M-CELL-EXCLUSION.md','docs/VALIDATED-IDENTIFIED-SET-REVIEW.md',
            'evidence/parametric_actual_M_independent_review.json','evidence/validated_identified_set_actual_M_cover.json']:
    if rel not in gate['evidence']:gate['evidence'].append(rel)
g['latest_uniform_parameter_phase']={
    'status':'PARTIAL_CONDITIONAL_PARAMETER_COVERAGE',
    'old_material_gain_family':'Relative active increment [-1/1000,1/1000]; strict local discrete energy minima, globally injective bodies, J safely in[0.992464257,1.006576818].',
    'actual_M_cell':'Shared five-coordinate computational box with relativehalfwidth1e-8; apicaltension[0,2e-10]. ActualM selectedbranch staticapical contradiction only; not all68 fit or model rejection.',
    'all17_uniform_cohort':'PREPARATION_ONLY_AT_THIS_CHECKPOINT',
    'full_identified_set_computed':False,
    'continuum_accuracy_certified':False,
    'confirmed_novelty':False}
g['next_executable_inference_step']='Certify the SAME five-coordinate cell at all17 source doses and compose all68 static output enclosures, retaining failed tests as unresolved. Report the conditional discrepancy-family threshold without relabeling it measurement noise. This local stationary-tuple component does not replace full parameter/geometry/branch/tail inference, uniform integration/continuum bounds, source-time/reference/ROI resolution, or revised larger-region predictions. A witnesspair search requires a defensible nonempty calibration contract; the alreadyexcluded M family cannot be IN at5.1micrometres.'
path.write_text(json.dumps(g,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')

old='## Findings\n'
new='## Findings\n\n* **Uniform parameter certificates, with a narrow scope:** [current checkpoint](controlled_furrowing/docs/PARAMETER-INFERENCE-CHECKPOINT.html). One old weak-load model now has certified local energy minima and globally injective geometry throughout a relative gain interval of ±0.1%. Separately, an extremely small five-coordinate box at actual calibration record M predicts only about6.81361µm static apical motion against a reported12µm. Its selected root family is excluded below approximately5.18638714µm apical discrepancy; this is not a calibrated noise level, all-observation fit, full-model rejection or continuum result. The shared-cell all17 extension is being prepared. [Actual-M proof and cover visual](controlled_furrowing/docs/ACTUAL-M-CELL-EXCLUSION.html).\n'
change(PARENT/'README.md',old,new)

path=ROOT/'scripts/build_report.py'
old="    body+='<h2>Latest: completed candidate replay and material inference</h2>"
new="    body+='<h2>Latest: uniform parameter families and a conditional exclusion</h2><p><a href=\"PARAMETER-INFERENCE-CHECKPOINT.html\">Read the current checkpoint, proofs and remaining original requirements</a>. One old weak-load model has uniform local energy minima and globally injective bodies throughout a relative gain interval of ±0.1%. A different, extremely narrow five-coordinate cell at actual calibration record M has a certified static-apical response near6.81361µm against the reported12µm. Its selected nearby stationary family is excluded below approximately5.18638714µm declared apical discrepancy. This does not reject all branches or establish a noise estimate, all68 fit, continuum accuracy or novelty. The shared-cell all17 extension is being prepared. <a href=\"ACTUAL-M-CELL-EXCLUSION.html\">Actual-M proof and finite-cover visual</a>; <a href=\"ENERGY-GAIN-CELL.html\">energy-norm interval proof</a>; <a href=\"ENERGY-GAIN-GEOMETRY-REVIEW.html\">independent geometry review</a>.</p>'\n    body+='<h2>Completed candidate replay and material inference</h2>"
change(path,old,new)

path=PARENT/'scripts/build_dossier.py'
old="content+='<h2>Latest: controlled-furrowing five-gate study</h2>"
new="content+='<h2>Latest: uniform finite-model certificates and remaining scientific gaps</h2><p>The controlled-furrowing study now verifies an old weak-load gain family with local energy minima and whole-body geometry, and a distinct actual-M narrow parameter-cell exclusion under a static observer. These are conditional finite-model results. All68 compatible calibration, source-matched time summaries, full parameter/geometry inference and continuum shape prediction remain incomplete. <a href=\"../controlled_furrowing/docs/PARAMETER-INFERENCE-CHECKPOINT.html\">Current proof-by-proof checkpoint</a>; <a href=\"../controlled_furrowing/docs/ACTUAL-M-CELL-EXCLUSION.html\">actual-M result and cover visual</a>.</p>'\ncontent+='<h2>Controlled-furrowing five-gate study</h2>"
change(path,old,new)

print(json.dumps({'status':'REPORTING_CHECKPOINT_UPDATED','historical_archive':str(ARCHIVE),
    'gates_sha256':sha(ROOT/'evidence/gates.json'),'plan_unchanged':True}))
