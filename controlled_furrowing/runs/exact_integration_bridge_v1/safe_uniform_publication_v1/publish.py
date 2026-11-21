"""Publish locally only after the mathematical reviews; preserve old receipts."""
from pathlib import Path
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];PROJECT=ROOT.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
def once(text,old,new):
 if text.count(old)!=1:raise ValueError('unique publication replacement required: '+old[:80])
 return text.replace(old,new)
checks={}
def ck(k,v):
 checks[k]=bool(v)
 if not v:raise AssertionError(k)
certificate=ROOT/'runs/exact_integration_bridge_v1/integral_root_B_safe_v2/results.json'
review=ROOT/'evidence/exact_integral_root_B_safe_independent_review.json'
ck('literal_mathematical_certificate',sha(certificate)=='42d97691311c057d8152f54a6ea31a7218763e7b451df63358bdf423010e18d4')
ck('literal_final_independent_review',sha(review)=='4500e5b34039eafe2e8af742cad739b1133ab414521fd6ed604cf4787d80d561')
ck('actual_root_and_minimum_certified',read(certificate)['root_certified']is True and read(certificate)['strict_local_minimum_certified']is True)
ck('all_independent_checks_true',read(review)['checks']and all(v is True for v in read(review)['checks'].values()))
manifest=read(HERE/'prepublication_snapshot/manifest.json')['files']
controls=[PROJECT/'README.md',ROOT/'scripts/build_report.py',PROJECT/'scripts/build_dossier.py',ROOT/'evidence/gates.json',ROOT/'docs/EXACT-INTEGRATION-BRIDGE.md',ROOT/'docs/CURRENT-FIVE-GATE-AUDIT.md',ROOT/'docs/COMPLETION-AUDIT.md']
for p in controls:ck('unchanged_archived_control_'+p.relative_to(PROJECT).as_posix(),sha(p)==manifest[p.relative_to(PROJECT).as_posix()])
pointdraft=ROOT/'runs/exact_integration_bridge_v1/safe_v2_publication_preparation/SAFE-POINT-RESTORATION.md'
rootdraft=HERE/'preparation/EXACT-INTEGRAL-ROOT.md'
point=pointdraft.read_text(encoding='utf8').replace('(safe_point_restoration.png)','(../runs/exact_integration_bridge_v1/safe_v2_publication_preparation/safe_point_restoration.png)')
point=point[:point.index('This draft and original')].rstrip()+'\n\nThis corrected point result is now independently reviewed. The subsequent [uniform integral-root theorem](EXACT-INTEGRAL-ROOT.md) adds full-domain derivative and residual bounds; it proves more than the point certificate alone. Historical unsafe receipts remain withdrawn.\n'
roottext=rootdraft.read_text(encoding='utf8')
ck('new_draft_has_final_review_and_two_figures','1,500' in roottext and 'root_margins.png' in roottext and 'root_location.png' in roottext)
updates={ROOT/'docs/SAFE-INTEGRATION-REPLAY.md':point,ROOT/'docs/EXACT-INTEGRAL-ROOT.md':roottext}
old=controls[4];text=old.read_text(encoding='utf8');lines=text.splitlines();warning=next(v for v in lines if v.startswith('**WITHDRAWN PENDING SAFE REPLAY.**'))
updates[old]=once(text,warning,'**HISTORICAL UNSAFE RECEIPTS REMAIN WITHDRAWN.** The sparse interval-polynomial failure invalidated the earlier point certificates retained below. New, separately reviewed safe computations now establish the [corrected point result](SAFE-INTEGRATION-REPLAY.md) and a [uniform exactly integrated finite-element root for record B](EXACT-INTEGRAL-ROOT.md). The original unsafe receipts are not restored. Read the newer chapters for current proof status; the rest of this chapter is the retained historical development record.')
readme=PROJECT/'README.md';lines=readme.read_text(encoding='utf8').splitlines()
for i,line in enumerate(lines):
 if line.startswith('* **Correction: integration certificates withdrawn'):
  lines[i]='* **Arithmetic correction completed:** the unsafe interval-polynomial receipts remain withdrawn. New exact-zero-only calculations and independent reviews support the [corrected point proof](controlled_furrowing/docs/SAFE-INTEGRATION-REPLAY.html). The older dense quadrature-root proofs were unaffected.'
 if line.startswith('* **WITHDRAWN pending replay — proposed integration'):
  lines[i]='* **Exactly integrated finite-element equilibrium certified:** for one B-dose case and every parameter in the same small cell, the [new proof, two visuals and receipts](controlled_furrowing/docs/EXACT-INTEGRAL-ROOT.html) establish a unique nearby equilibrium, a strict local energy minimum, and exclusion from the parameter-matched old quadrature cube. Root error from the translated predictor is below 2.557e-12 in normalized coordinates. A second residual route also passes. This concerns a fixed 611-coordinate model; it does not establish all17 integral roots, PDE accuracy, biological validation or novelty. The full five-gate study remains incomplete.'
updates[readme]='\n'.join(lines)+'\n'
child=ROOT/'scripts/build_report.py';text=child.read_text(encoding='utf8');line=next(v for v in text.splitlines()if 'Correction: integration certificates withdrawn pending replay</h2>'in v)
updates[child]=once(text,line,'    body+=\'<h2>Corrected arithmetic and an exactly integrated equilibrium</h2><p>New independent reviews establish an exact-spatial-integral finite-element equilibrium and strict local energy minimum for record B throughout the same small parameter cell. The new root is outside its parameter-matched old quadrature cube; common-region uniqueness excludes integral roots in that old cube. A second residual composition also passes. This is a fixed 611-coordinate model, not a continuum or biological validation. <a href="EXACT-INTEGRAL-ROOT.html">Proof, visuals and receipts</a>; <a href="SAFE-INTEGRATION-REPLAY.html">arithmetic correction and safe replay</a>; <a href="EXACT-INTEGRATION-BRIDGE.html">historical withdrawn receipts</a>. The original five gates remain incomplete.</p>\'')
parent=PROJECT/'scripts/build_dossier.py';text=parent.read_text(encoding='utf8');line=next(v for v in text.splitlines()if 'Correction: integration certificates withdrawn pending replay</h2>'in v)
updates[parent]=once(text,line,'content+=\'<h2>Corrected arithmetic and an exact finite-element equilibrium</h2><p>The controlled-furrowing study now certifies one B-dose integral equilibrium and strict local energy minimum throughout the same small parameter cell, with a second successful residual route. This does not establish continuum accuracy, biological prediction or novelty. <a href="../controlled_furrowing/docs/EXACT-INTEGRAL-ROOT.html">New proof, visuals and retained verification</a>; <a href="../controlled_furrowing/docs/SAFE-INTEGRATION-REPLAY.html">safe replay and arithmetic correction</a>. Historical unsafe receipts remain withdrawn; the five-gate study remains incomplete.</p>\'')
gates=ROOT/'evidence/gates.json';d=read(gates);d['snapshot']='Corrected SAFE arithmetic independently replayed. ONE B-dose EXACT-SPATIAL-INTEGRAL finite-element root is uniformly certified on the same small C, with strict local energy minimum and parameter-matched old-cube exclusion. Two residual routes pass. Old unsafe receipts remain withdrawn. All17 integral roots, PDE/biology/full identified set and full five-gate completion remain unproved.';d['goal_achieved']=False
g=d['gates'][2];mark=' CORRECTION: proposed integration and center-escape result is WITHDRAWN pending safe sparse-Arb cleanup and full replay. Exact geometry boxes remain certified independently of any root. Old numerical proposal is retained; its integral-defect accuracy must be re-established.'
ck('old_gate_withdrawal_marker_found',mark in g['scientific_outcome'])
g['scientific_outcome']=g['scientific_outcome'].replace(mark,' New SAFE replay supersedes the withdrawn point proof with independent evidence. Uniform derivative/residual estimates now certify ONE B-dose exactly integrated finite-element equilibrium, a strict local minimum and exclusion from the old parameter-matched cube; a second direct residual route also passes. This does not certify the other16 integral roots or close continuum, physical-domain, observation or larger-prediction gaps.')
for name in ['docs/SAFE-INTEGRATION-REPLAY.md','docs/EXACT-INTEGRAL-ROOT.md','evidence/exact_integral_root_B_safe_independent_review.json','evidence/exact_uniform_direct_root_B_crosscheck.json']:
 if name not in g['evidence']:g['evidence'].append(name)
updates[gates]=json.dumps(d,indent=2)+'\n'
checkpoint='\n**Current integration checkpoint.** A new [exactly integrated finite-element root theorem](EXACT-INTEGRAL-ROOT.md) covers ONE B dose uniformly on the existing small parameter cell, including strict local minimality and old-cube exclusion. Both reviewed residual routes pass. This supersedes blanket statements that no integral-root transfer has been established; it does not supply the remaining16 integral roots, PDE error, full physical identified set or source-compatible larger-region prediction. Old unsafe receipts remain withdrawn; the five-gate objective remains incomplete.\n'
for name in ['CURRENT-FIVE-GATE-AUDIT.md','COMPLETION-AUDIT.md']:
 p=ROOT/'docs'/name;text=p.read_text(encoding='utf8');a,b=text.split('\n',1);updates[p]=a+'\n'+checkpoint+b
for p,text in updates.items():p.write_text(text,encoding='utf8')
for command,cwd in [([sys.executable,str(ROOT/'scripts/build_report.py')],ROOT),([sys.executable,str(PROJECT/'scripts/build_dossier.py')],PROJECT)]:
 result=subprocess.run(command,cwd=cwd,capture_output=True,text=True)
 if result.returncode:raise RuntimeError(result.stdout+result.stderr)
 print(result.stdout.strip(),flush=True)
out={'status':'LOCAL_DOCUMENTS_BUILT_BROWSER_VERIFICATION_PENDING','checks':checks,'source_and_output_hashes':{p.relative_to(PROJECT).as_posix():sha(p)for p in list(updates)+[ROOT/'docs/EXACT-INTEGRAL-ROOT.html',ROOT/'docs/SAFE-INTEGRATION-REPLAY.html',ROOT/'docs/index.html',PROJECT/'docs/index.html',Path(__file__)]},'github_pushed':False,'public_site_deployed':False,'full_five_gate_goal_complete':False}
(HERE/'local_publication.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({'status':out['status'],'changed_source_files':len(updates),'local_receipt_sha256':sha(HERE/'local_publication.json')}),flush=True)
