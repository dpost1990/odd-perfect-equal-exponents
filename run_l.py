#!/usr/bin/env python
"""Driver: for prime l, try to show the generic closure from {l} has count >= target, using cycfac + ECM."""
import sys, json, gmpy2
from closure import Engine, closure, phi_l
from ecmtools import factor_partial
l = int(sys.argv[1]); X = int(float(sys.argv[2])); target = int(sys.argv[3]) if len(sys.argv) > 3 else l + 6
rounds = int(sys.argv[4]) if len(sys.argv) > 4 else 6
eng = Engine(l, X)
tried = set()
for rnd in range(rounds):
    log = []
    res = closure(eng, [l], need=target, log=log)
    print(f'round {rnd}: E={len(res["E"])} Uprimes={len(res["Uprimes"])} cof={len(res["cof"])} count={res["count"]} coprime={res["coprime_ok"]}', flush=True)
    if res['count'] >= target: break
    # ECM on the smallest untried cofactors
    cofs = sorted(res['cof'].items(), key=lambda kv: kv[1])
    progress = False
    for Q, V in cofs:
        if (Q, V) in tried: continue
        tried.add((Q, V))
        d = len(str(V))
        plan = ((2000, 30), (11000, 80), (50000, 200), (250000, 300)) if d <= 70 else ((2000, 30), (11000, 60), (50000, 60)) if d <= 200 else ((2000, 20), (11000, 30))
        primes, comps = factor_partial(V, plan)
        print(f'  ecm on cof of {Q} ({d} digits): primes {primes} comps {[len(str(c)) for c in comps]}', flush=True)
        for q in primes:
            if len(comps) == 0 and q == max(primes) and len(primes) == 1: continue
            eng.add_extra(Q, q); progress = True
        if progress and d > 70: break
    if not progress: print('no progress'); break
for Q, facs, V in log:
    print(Q, '->', facs, 'cof digits', len(str(V)) if V > 1 else 0)
print('E', res['E']); print('U', res['Uprimes'])
