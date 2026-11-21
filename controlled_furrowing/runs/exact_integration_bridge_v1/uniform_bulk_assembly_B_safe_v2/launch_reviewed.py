"""Single-run wrapper retaining actual process and terminal evidence."""
from pathlib import Path
import datetime, hashlib, json, os, subprocess, sys, time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / 'scripts/exact_uniform_bulk_assembly_B_safe.py'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def save(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2) + '\n', encoding='utf8')

expected = {
    SOURCE: '3fd3fa4a2b70ddf3a390f42289bb95ac85c59129daafa3e4e647f1f28c51086a',
    HERE / 'protocol.json': '9ab7c568d5cff63a53fdcec16fcdac34d8b004ab735190ad176ccc06ae978828',
    HERE / 'preflight.json': 'ee80aa589fb56d03cd476f62889ca87d7b0afe0563d51ed2579d49769ec763da',
    HERE / 'prelaunch_independent_review.json': '99a8fa459fd331fe78392703d4d71b9cb0c9a9c08a07fc8be6c124c397648690',
}
assert all(sha(p) == h for p, h in expected.items())
review = json.loads((HERE / 'prelaunch_independent_review.json').read_text())
assert review['status'] == 'PASS' and review['checks'] and all(v is True for v in review['checks'].values())
assert not any((HERE / n).exists() for n in ('launch_job.json', 'started.json', 'results.json', 'launch_terminal.json'))
start = time.monotonic()
child = subprocess.Popen([sys.executable, '-u', str(SOURCE), 'run'], cwd=ROOT,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, encoding='utf8')
save('launch_job.json', {'status': 'RUNNING_AT_LAUNCH', 'pid': child.pid, 'wrapper_pid': os.getpid(),
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
save('launch_terminal.json', terminal)
print(json.dumps(terminal), flush=True)
sys.exit(code)
