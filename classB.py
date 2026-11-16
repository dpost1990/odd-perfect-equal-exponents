#!/usr/bin/env python
"""classB.py r T : class version.  Hypothesis: N odd perfect, r prime, r == 3 (mod 4), and r | 2*beta_i + 1 for EVERY
non-special prime q_i of N.  Then r | N and (same forcing rules as in prove.py) we try to exhibit >= T distinct
non-special primes == 1 (mod r) dividing N in every scenario.  Since each of them contributes at least one factor r
to sigma(N), this proves v_r(N) >= T.   Output: certs/class_r<r>.json (checked by verify.py)."""
import sys as _sys; _sys.set_int_max_str_digits(0)   # certificates contain integers with > 4300 digits
import sys, json, os, time
from closure import Engine
from prove import run_scenario, HERE
from ecmtools import factor_partial
r = int(sys.argv[1]); T = int(sys.argv[2])
assert r % 4 == 3
eng = Engine(r, 10**9)
t0 = time.time()
out = {'class': True, 'n': None, 'r': r, 'need': T, 'scenarios': []}
G = run_scenario(eng, r, [r], None, T)
if G is None: print('generic FAILED'); sys.exit(1)
out['scenarios'].append(G)
specials = [int(q) for q in G['E'] if int(q) % 4 == 1]
print(f'class r={r}: generic ok ({G["count"]}/{T}); specials: {len(specials)}', flush=True)
for Qs in specials:
    primes, comps = factor_partial((Qs + 1) // 2, ((2000, 40), (11000, 100), (50000, 300), (250000, 600)))
    sc = run_scenario(eng, r, [r] + [q for q in primes if q != r], Qs, T, log=lambda *a: None)
    if sc is None: print(f'scenario p={Qs} FAILED'); sys.exit(1)
    sc['half_factors'] = [str(q) for q in primes]
    out['scenarios'].append(sc)
fn = os.path.join(HERE, 'certs', f'class_r{r}.json')
json.dump(out, open(fn + '.tmp', 'w')); os.replace(fn + '.tmp', fn)
print(f'class r={r}: v_r(N) >= {T} certified ({len(out["scenarios"])} scenarios, {time.time()-t0:.0f}s)', flush=True)
