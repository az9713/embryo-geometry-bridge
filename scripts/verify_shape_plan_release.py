"""Verify saved authoring evidence and scientific-input identities for this update."""
from pathlib import Path
from hashlib import sha256
import argparse
import json

ROOT=Path(__file__).resolve().parents[1]


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--output',type=Path,required=True)
    args=cli.parse_args()
    plan=ROOT/'controlled_furrowing/research_directions/regularization_shape_v1'
    manifest=json.loads((plan/'MANIFEST.json').read_text())
    for relative,expected in manifest['files'].items():
        assert digest(plan/relative)==expected,relative
    guide=json.loads((ROOT/'docs/SHAPE-DESCRIPTION-BUILD.json').read_text())
    for relative,expected in guide['sha256'].items():
        assert digest(ROOT/relative)==expected,relative
    source=ROOT/'controlled_furrowing/runs/nonaxisymmetric_surface_sign_preparation_v1/independent_actual_readback.json'
    prior=json.loads(source.read_text())
    for relative,expected in prior['verified_source_and_actual_sha256'].items():
        assert digest(ROOT/'controlled_furrowing'/relative)==expected,relative
    receipt={'status':'PASS_SHAPE_PLAN_RELEASE_IDENTITIES',
             'plan_bundle_files':len(manifest['files']),
             'guide_build_files':len(guide['sha256']),
             'retained_scientific_inputs':len(prior['verified_source_and_actual_sha256']),
             'main_theorem_proved':False,'novelty_confirmed':False,
             'scope':'Saved bundle identities and retained input hashes. Does not rerun historical integrations or prove the proposed continuum theorem.',
             'verifier_sha256':digest(Path(__file__))}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print(json.dumps(receipt))


if __name__=='__main__':
    main()
