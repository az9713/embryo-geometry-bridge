"""Single-run wrapper retaining actual process and terminal evidence."""
from pathlib import Path
import datetime, hashlib, json, os, subprocess, sys, time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / 'scripts/exact_integral_root_B_safe.py'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def save(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2) + '\n', encoding='utf8')

expected = {
    SOURCE: '8404a27fa339f4990c9a3fb784a84f5576cfa1fea8eaa20f1f58f8cf6ff6633c',
    HERE / 'protocol.json': 'ce2a816ba037725395cc4dc2b6a8f419edd1c645c071d542893615b228fefd41',
    HERE / 'preflight.json': '547f58b221ee0ef1d9857b908edd5f7d5543fefe7cf650755691badc4d8dc1df',
    HERE / 'prelaunch_independent_review.json': 'd3cafecf9a3152f6cbae39d10cf7bfcedd73123b4215cdde5e5777c212f08478',
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
