"""Freeze all 44 primary checkpoints after independently successful replay."""
from hashlib import sha256
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
REVIEW = BASE/'independent_actual_preparation/actual_replay/receipt.json'
REVIEW_SHA = 'c795ae78ba0b53a01c3d65f466d5fa2f1d233c2a68c4275c04260c3cb5120ebd'


def sha(p):
    return sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_bytes())


def main():
    assert sha(REVIEW) == REVIEW_SHA
    review = read(REVIEW)
    assert review['status'] == 'PASS_INDEPENDENT_ACTUAL_NONLINEAR_STORED_FIELD_REPLAY'
    done = read(REVIEW.parent/'terminal.json')
    assert done['status'] == 'FINISHED' and done['exit_code'] == 0 and done['receipt_sha256'] == REVIEW_SHA
    result_path, terminal_path = BASE/'actual/results.json', BASE/'actual/terminal.json'
    result, terminal = read(result_path), read(terminal_path)
    assert sha(result_path) == review['actual_results_sha256'] == terminal['results_sha256']
    assert review['all_hashes'][terminal_path.relative_to(ROOT).as_posix()] == sha(terminal_path)
    fields = []
    for endpoint in result['endpoints']:
        for row in endpoint['checkpoints']:
            for path_key, digest_key in (('field','field_sha256'), ('verification','verification_sha256')):
                assert review['all_hashes'][row[path_key]] == row[digest_key] == sha(ROOT/row[path_key])
            fields.append({'endpoint_id': endpoint['endpoint_id'], **row})
    assert len(fields) == 44 and len({x['field'] for x in fields}) == 44
    output = {'status':'COMPLETE_44_PRIMARY_CHECKPOINTS', 'primary_results_sha256':sha(result_path),
              'primary_terminal_sha256':sha(terminal_path), 'independent_primary_review_sha256':REVIEW_SHA,
              'independent_replay_unified_exec_session':62471, 'independent_replay_observed_os_exit_code':0,
              'fields':fields, 'source_sha256':sha(Path(__file__)),
              'scope':'All saved fields and verification links, not an equilibrium error certificate.'}
    p = BASE/'field_manifest.json'
    with p.open('x',encoding='utf-8') as f:
        json.dump(output,f,indent=2)
        f.write('\n')
    print(json.dumps({'field_manifest_sha256':sha(p),'fields':len(fields)}))


if __name__=='__main__':
    main()
