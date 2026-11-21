"""Single-run wrapper retaining actual process and terminal evidence."""
from pathlib import Path
import datetime, hashlib, json, os, subprocess, sys, time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / 'scripts/exact_uniform_bulk_defect_B_safe.py'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def save(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2) + '\n', encoding='utf8')

expected = {
    SOURCE: '7c43896c220090025e3bc98969340f98c2b145099905649d12b16bb57d6fac9c',
    HERE / 'protocol.json': 'bbc559fab371051ab272db511dc3b2a6fffa8fde7797a419ef3a33a1f050d266',
    HERE / 'preflight.json': '627b4aec091ce3568015a58654ea09e458596955e08acd544721a33ba722f030',
    HERE / 'prelaunch_independent_review.json': '49f69791f98ecc4e5fc1f335194f71cd7b0a2f1c512e0b73b23d6f4641033efe',
}
assert all(sha(p) == h for p, h in expected.items())
review = json.loads((HERE / 'prelaunch_independent_review.json').read_text())
assert review['status'] == 'PASS' and review['checks'] and all(v is True for v in review['checks'].values())
assert not any((HERE / n).exists() for n in ('job.json', 'started.json', 'results.json', 'terminal.json'))
start = time.monotonic()
child = subprocess.Popen([sys.executable, '-u', str(SOURCE), 'run'], cwd=ROOT,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, encoding='utf8')
save('job.json', {'status': 'RUNNING_AT_LAUNCH', 'pid': child.pid, 'wrapper_pid': os.getpid(),
                 'utc': stamp(), 'command': [str(SOURCE.relative_to(ROOT)), 'run'],
                 'bindings': {str(p.relative_to(ROOT)): h for p, h in expected.items()},
                 'prelaunch_checks_passed': len(review['checks'])})
print(json.dumps({'phase': 'launched', 'pid': child.pid}), flush=True)
with (HERE / 'stdout.log').open('x', encoding='utf8') as log:
    for line in child.stdout:
        log.write(line)
        log.flush()
        print(line, end='', flush=True)
code = child.wait()
terminal = {'status': 'FINISHED' if code == 0 else 'FAILED', 'exit_code': code,
            'pid': child.pid, 'utc': stamp(), 'seconds': time.monotonic() - start,
            'all_launch_bindings_current': all(sha(p) == h for p, h in expected.items()),
            'stdout_sha256': sha(HERE / 'stdout.log')}
if (HERE / 'results.json').exists():
    terminal['results_sha256'] = sha(HERE / 'results.json')
save('terminal.json', terminal)
print(json.dumps(terminal), flush=True)
sys.exit(code)
