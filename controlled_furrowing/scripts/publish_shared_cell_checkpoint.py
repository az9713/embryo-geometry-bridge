"""Update local reporting only; archive prior control bytes and preserve PLAN."""
from pathlib import Path
import hashlib,json,shutil

ROOT=Path(__file__).resolve().parents[1];PARENT=ROOT.parent
ARCHIVE=ROOT/'runs/checkpoint_history/parameter_family_phase'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def replace(path,old,new):
    text=path.read_text(encoding='utf-8')
    if new in text:return
    if old not in text:raise RuntimeError('reporting insertion absent '+str(path))
    path.write_text(text.replace(old,new,1),encoding='utf-8')

def main():
    inference=ROOT/'evidence/shared_cell_cohort_inference.json';v=json.loads(inference.read_text(encoding='utf-8'))
    if v['records']!=17 or v['constraints']!=68 or v['full_five_gate_goal_complete'] is not False:raise RuntimeError('actual scoped inference required')
    if not ARCHIVE.exists():
        paths={
            'evidence/gates.json':ROOT/'evidence/gates.json',
            'evidence/parameter_inference_checkpoint_review.json':ROOT/'evidence/parameter_inference_checkpoint_review.json',
            'docs/CURRENT-FIVE-GATE-AUDIT.md':ROOT/'docs/CURRENT-FIVE-GATE-AUDIT.md',
            'scripts/build_report.py':ROOT/'scripts/build_report.py',
            'docs/index.html':ROOT/'docs/index.html',
            'parent/README.md':PARENT/'README.md','parent/scripts/build_dossier.py':PARENT/'scripts/build_dossier.py',
            'parent/docs/index.html':PARENT/'docs/index.html'}
        hashes={}
        for name,source in paths.items():
            target=ARCHIVE/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target);hashes[name]=sha(target)
        (ARCHIVE/'archive_manifest.json').write_text(json.dumps({'status':'ARCHIVED_EXACT_PREVIOUS_REPORTING_SNAPSHOT','files':hashes},indent=2)+'\n',encoding='utf-8')
    else:
        manifest=json.loads((ARCHIVE/'archive_manifest.json').read_text(encoding='utf-8'))
        if any(sha(ARCHIVE/k)!=h for k,h in manifest['files'].items()):raise RuntimeError('previous reporting archive changed')
    path=ROOT/'evidence/gates.json';g=json.loads(path.read_text(encoding='utf-8'))
    g['snapshot']='Seventeen local stationary-root families and all68 conditional static outputs now have uniform bounds on ONE tiny shared five-parameter cell. F/K required one explicit numerical-center refinement, with physical potential/R/S/C unchanged. Full physical compatible-set inference, source-matched time/geometry observations, integration/continuum accuracy and revised larger-region shape prediction remain incomplete.'
    g['goal_achieved']=False
    if 'latest_uniform_parameter_phase' in g:
        historical=g.pop('latest_uniform_parameter_phase')
        historical.update(historical_snapshot=True,superseded_by='latest_shared_cell_inference')
        g['historical_uniform_parameter_phase']=historical
    gate=g['gates'][2]
    gate['scientific_outcome']='Finite-height P2 mechanics, integration corrections and sparse Newton methods are implemented. Strong-load inter-mesh spatial discrepancies remain approximately0.67–0.68micrometres, without a continuum error theorem. One common five-coordinate computational cell now supports17 independently verified nearby stationary roots and68 uniform static output intervals. Two old-center sufficient tests failed; larger radii also failed; one numerical-center correction percase passed fresh independent proofs without changing potential/material/gain/dose/C/R/S. Exact geometry and explicit center-to-original-graph transfers are reviewed. The selected tuple worst-case residual is safely enclosed by[5.18638714,5.18638737]micrometres throughoutC. Illustrative5.1/5.2 allowances exclude/admit this conditional stationary family, not measured noise or stable biological trajectories. The extremely small C is not the full compatible set or a physiological prior. Source temporal/reference/ROI/geometry uncertainty and uniform integration/continuum accuracy remain unresolved. The separate old weak-load gain family has local discrete energy minima; that stability result is not transferred to the actual cohort.'
    for name in ['docs/SHARED-CELL-COHORT.md','docs/SHARED-CELL-COHORT-REVIEW.md',
        'docs/COHORT-UNIFORM-ROOTS-REVIEW.md','docs/COHORT-RECENTERED-ROOTS-REVIEW.md',
        'docs/COHORT-STATIC-APPLICATION-REVIEW.md','docs/RECENTERED-OBSERVERS-REVIEW.md',
        'evidence/shared_cell_cohort_inference.json','evidence/cohort_recentered_roots_independent_review.json',
        'evidence/cohort_static_observers_independent_review.json','evidence/recentered_observers_independent_review.json']:
        if name not in gate['evidence']:gate['evidence'].append(name)
    g['latest_shared_cell_inference']={'status':v['status'],'file':'evidence/shared_cell_cohort_inference.json','sha256':sha(inference),
        'records':17,'constraints':68,'parameter_cell_sha256':v['parameter_cell_sha256'],'model_observation_contract_sha256':v['model_observation_contract_sha256'],
        'worst_case_residual_lower_exact':v['maximum_residual_lower_exact'],'worst_case_residual_upper_exact':v['maximum_residual_upper_exact'],
        'interpretation':'Pointwise enclosure of D(theta) on a tiny computational cell, not a proved attained sharp range or a global optimum.',
        'estimated_noise':False,'full_identified_set':False,'continuum_accuracy':False,'source_temporal_validation':False,'novelty_confirmed':False}
    g['remaining_original_requirements']['uniform_validated_forward_bound']='PARTIAL: uniform algebraic-root/static-output bounds on one tiny fixed-geometry parameter cell; integral/continuum accuracy and full admitted domain remain incomplete'
    g['next_executable_inference_step']='Resolve a defensible source observation/time/reference/geometry and discrepancy contract, then continue full conditional identified-set coverage including remaining finite cells, zero faces, other branches and unbounded tails. The common-cell stationary witness is now complete; it is not a substitute for the full domain, a stable activation path, continuum error control or revised source-matched larger-region predictions. A substantive transfer-ambiguity result additionally requires distinct compatible explanations and a separating feasible observation; neither has been established.'
    path.write_text(json.dumps(g,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')
    readme=PARENT/'README.md';text=readme.read_text(encoding='utf-8');start=text.index('* **Uniform parameter certificates, with a narrow scope:**');end=text.index('\n\n',start)
    text=text[:start]+'* **All68 conditional static outputs now have uniform certificates:** [proof, full-cohort visual and receipts](controlled_furrowing/docs/SHARED-CELL-COHORT.html). The SAME tiny five-parameter cell supports17 independently verified local stationary-root families. Two cases required a single numerical-center correction; their physical target and parameters stayed fixed. Throughout that cell, the maximum absolute residual is safely enclosed by[5.18638714,5.18638737]µm. Illustrative allowances5.1/5.2µm exclude/admit the selected stationary tuple; these are not estimates of measurement noise. Full physical parameter/geometry inference, stable time paths, source temporal matching, continuum accuracy and novelty remain unresolved.'+text[end:]
    readme.write_text(text,encoding='utf-8')
    replace(ROOT/'scripts/build_report.py',"    body+='<h2>Latest: uniform parameter families and a conditional exclusion</h2>",
        "    body+='<h2>Latest: all68 static outputs on one shared parameter cell</h2><p><a href=\"SHARED-CELL-COHORT.html\">Read the theorem, full-cohort visual and verification receipts</a>. Seventeen independently checked local stationary-root families now support68 uniform output intervals on the SAME extremely small five-parameter box. Two cases required a single numerical-center correction; the physical potential, parameters, preconditioner and predictor stayed fixed. The worst-case residual D(theta) is safely enclosed by[5.18638714,5.18638737]µm throughout that box. Illustrative5.1/5.2µm allowances exclude/admit this selected stationary tuple. These are not measurement-noise estimates, a stable activation history, full identified-set inference or continuum accuracy. Source temporal/geometry matching and novelty remain unresolved.</p>'\n    body+='<h2>Earlier: parameter families and one conditional exclusion</h2>")
    replace(ROOT/'scripts/build_report.py','The shared-cell all17 extension is being prepared.','At that earlier checkpoint the shared-cell all17 extension was being prepared; its current result is linked above.')
    replace(PARENT/'scripts/build_dossier.py',"content+='<h2>Latest: uniform finite-model certificates and remaining scientific gaps</h2>",
        "content+='<h2>Latest: a verified common-cell cohort, with scientific limits</h2><p>The controlled-furrowing model now has17 local stationary-root proofs and68 uniform static output bounds on one tiny common parameter box. Its maximum residual lies in the safe enclosure[5.18638714,5.18638737]µm. A single proof-center refinement resolved two numerical cases without changing physics or parameters. Source temporal matching, full material/geometry inference, stable histories, continuum accuracy and novelty remain unresolved. <a href=\"../controlled_furrowing/docs/SHARED-CELL-COHORT.html\">Full theorem, comparison visual and receipts</a>.</p>'\ncontent+='<h2>Earlier finite-model certificates and scientific gaps</h2>")
    print(json.dumps({'status':'LOCAL_REPORTING_UPDATED','gates_sha256':sha(ROOT/'evidence/gates.json'),'README_sha256':sha(readme),'original_goal_complete':False}))

if __name__=='__main__':main()
