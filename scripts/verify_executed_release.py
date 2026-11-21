"""Check release identities and exact scalar conclusions for the executed steps."""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import argparse
import json

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'controlled_furrowing/research_directions/relaxed_slab_execution_v1'


def read(p):return json.loads(p.read_text(encoding='utf8'))
def digest(p):return sha256(p.read_bytes()).hexdigest()


def main():
    cli=argparse.ArgumentParser();cli.add_argument('--output',type=Path,required=True);a=cli.parse_args()
    provenance=read(ROOT/'reproducibility/ENVIRONMENT-PROVENANCE.json')
    old=(ROOT/provenance['historical_requirements_record']).read_bytes()
    assert sha256(old).hexdigest()==provenance['historical_requirements_sha256']
    assert old.count(b'Pillow==11.0.0')==1
    assert (ROOT/'controlled_furrowing/requirements.txt').read_bytes()==old.replace(b'Pillow==11.0.0',b'Pillow==13.0.0')
    producer=read(RUN/'virtual_work_radial_2048.json');replay=read(RUN/'virtual_work_integer_recentered_4096.json')
    assert producer['source_sha256']==digest(RUN/'virtual_work.py')
    assert replay['source_sha256']==digest(RUN/'replay_virtual_work.py')
    checks=0
    for r in [producer,replay]:
        for p,h in r['input_sha256'].items():assert digest(ROOT/'controlled_furrowing'/p)==h;checks+=1
        assert len(r['rows'])==2
        for row in r['rows']:
            lo,hi=map(Q,row['reduced_virtual_work_exact']);assert lo<=hi<0
            assert row['node']==323 and row['component']=='radial'
    readable=[(Q('-0.001129'),Q('-0.000978')),(Q('-0.001180'),Q('-0.001028'))]
    for row,(a0,b0) in zip(producer['rows'],readable,strict=True):
        lo,hi=map(Q,row['reduced_virtual_work_exact']);assert a0<=lo<=hi<=b0
    slab=read(RUN/'slab_verification.json');independent=read(RUN/'slab_independent_constants.json')
    assert slab['source_sha256']==digest(RUN/'check_slab.py')
    assert slab['written_theorem_sha256']==digest(RUN/'THEOREM.md')
    assert independent['source_sha256']==digest(RUN/'replay_slab.py')
    assert independent['producer_receipt_sha256']==digest(RUN/'slab_verification.json')
    assert len(independent['angle_checks'])==4 and len(independent['bound_checks'])==7
    for row in read(RUN/'document_build.json')['rows']:
        assert digest(RUN/row['source'])==row['source_sha256']
        assert digest(RUN/row['html'])==row['html_sha256']
    assert read(RUN/'pillow13_regression.json')['status']=='PASS_PILLOW13_ISOLATION_AND_PIXEL_REGRESSION'
    assert read(RUN/'old_runtime_guard.json')['status']=='PASS_UNPATCHED_IMAGE_RUNTIME_REJECTED'
    r={'status':'PASS_EXECUTED_RELEASE_IDENTITIES_AND_STRICT_INEQUALITIES',
       'source_bindings_checked':checks,'strict_virtual_work_cases':2,'arithmetic_implementations':2,
       'slab_angle_checks':4,'slab_bound_rows':7,'only_dependency_pin_changed':'Pillow 11.0.0 to 13.0.0',
       'scope':'Evidence identities and exact signed scalar conclusions; written proof assessment and external peer review are distinct.',
       'full_cylinder_or_biological_theorem_proved':False,'sharp_rate_or_novelty_claim':False,
       'verifier_sha256':digest(Path(__file__))}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps(r))


if __name__=='__main__':main()
