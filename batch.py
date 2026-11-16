#!/usr/bin/env python
"""batch.py T n1 n2 ... : run prove.py for each n with a time limit of T seconds; append results to batch_results.jsonl"""
import sys, subprocess, os, signal, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
PY = 'python3'
T = int(sys.argv[1])
for a in " ".join(sys.argv[2:]).split():
    n = int(a)
    t0 = time.time()
    logf = open(os.path.join(HERE, 'logs', f'prove_n{n}.log'), 'w')
    p = subprocess.Popen(['nice', '-n', '10', PY, os.path.join(HERE, 'prove.py'), str(n), '-', '1e9'], stdout=logf, stderr=subprocess.STDOUT,
                         cwd=HERE, start_new_session=True)
    try:
        rc = p.wait(timeout=T); status = 'closed' if rc == 0 else 'failed'
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL); p.wait(); status = 'timeout'
    res = {'n': n, 'status': status, 'seconds': round(time.time() - t0)}
    with open(os.path.join(HERE, 'batch_results.jsonl'), 'a') as f: f.write(json.dumps(res) + '\n')
    print(res, flush=True)
