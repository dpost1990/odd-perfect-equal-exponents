#!/usr/bin/env python
"""prove.py n [r] [X]: try to prove 'no odd perfect number N = p^a * M^(n-1)' (all non-special exponents
equal to 2*beta = n-1), n = 2*beta+1 odd, using a prime r | n (default: n itself if prime).

Facts used (proved in report.md):
 (F1) r | N, r is not the special prime (r == 3 mod 4, or n a power of r), so r^(n-1) || N.
 (F2) For a non-special prime q | N, q != r: Phi_r(q) | 2N; its prime factors are r (iff q == 1 mod r) or
      primes == 1 (mod r) dividing N.  Phi_r(r) | N likewise.
 (F3) S = {non-special primes == 1 mod r dividing N}: each contributes exactly e = v_r(n) factors r to sigma(N),
      so |S| <= (n-1)/e.
 (F4) (p+1)/2 | N.
Scenario G: p is not one of the expanded primes.  Scenario p=Q*: one per expanded prime Q* == 1 (mod 4) of G.
In each scenario we exhibit >= need = floor((n-1)/e)+1 distinct members of S.
Counting: E = expanded primes == 1 mod r (proven prime, not special); items = pairwise coprime integers > 1,
coprime to r, to E and to the special prime, each dividing some Phi_r(Q) with Q expanded; an item has weight 2
if it is composite and not a perfect power (>= 2 distinct primes), else 1.  count = |E| + sum(weights) - [G].
Output: cert_n<n>.json
"""
import sys as _sys; _sys.set_int_max_str_digits(0)   # certificates contain integers with > 4300 digits
import sys, json, math, os, gmpy2, time
from closure import Engine, phi_l
from ecmtools import factor_partial

HERE = os.path.dirname(os.path.abspath(__file__))

def coprime_base(items):
    items = [int(a) for a in items if a > 1]
    changed = True
    while changed:
        changed = False
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                g = math.gcd(items[i], items[j])
                if g > 1:
                    a, b = items[i], items[j]
                    while a % g == 0: a //= g
                    while b % g == 0: b //= g
                    rest = [x for k, x in enumerate(items) if k not in (i, j)]
                    items = [x for x in rest + [a, b, g] if x > 1]
                    changed = True; break
            if changed: break
    return items

def weight(a):
    if gmpy2.is_prime(a): return 1
    if gmpy2.is_power(a): return 1
    return 2

def run_scenario(eng, r, roots, special, need, ecm_rounds=80, maxexp=10**50, log=print, ecm_effort=1):
    tried = set()
    tag = 'p=' + str(special) if special else 'G'
    for rnd in range(ecm_rounds + 1):
        known = {int(q): 'root' for q in roots}
        expanded = []; steps = []; cof = {}
        def tally():
            E = [q for q in expanded if q % r == 1]
            U = [q for q in known if q not in expanded and q % r == 1 and q != special]
            cofs = []
            for Q, V in cof.items():
                for q in list(known) + ([special] if special else []):
                    while V % q == 0: V //= q
                if V > 1: cofs.append(V)
            items = coprime_base(U + cofs)
            cnt = len(E) + sum(weight(a) for a in items) - (0 if special else 1)
            return E, items, cnt
        while True:
            E, items, cnt = tally()
            if cnt >= need: break
            cand = [q for q in known if q not in expanded and q != special and q < maxexp and gmpy2.is_prime(q)]
            if not cand: break
            Q = min(cand)
            facs, V = eng.expand(Q)
            expanded.append(Q)
            for q, m in facs:
                if q != special and q not in known: known[q] = Q
            for q in list(known) + ([special] if special else []):
                while V > 1 and V % q == 0: V //= q
            fl = [q for q, m in facs]
            if V > 1 and gmpy2.is_prime(V):
                if V not in known: known[int(V)] = Q
                fl.append(int(V)); V = 1
            steps.append({'Q': str(Q), 'factors': [str(q) for q in fl]})
            if V > 1: cof[Q] = V
        log(f'   [{tag}] round {rnd}: expanded={len(expanded)} E={len(E)} items={len(items)} count={cnt}/{need}')
        if cnt >= need:
            return {'special': str(special) if special else None, 'roots': [str(q) for q in roots], 'steps': steps,
                    'E': [str(q) for q in E], 'items': [str(a) for a in items], 'count': cnt}
        progress = False
        for level in range(3):
            for Q, V in sorted(cof.items(), key=lambda kv: kv[1]):
                if (Q, V, level) in tried or eng.ecm_level(V) > level: continue
                tried.add((Q, V, level)); d = len(str(V))
                if level == 0:
                    if d <= 75: plan = ((2000, 30), (11000, 80), (50000, 200), (250000, 400))
                    elif d <= 220: plan = ((2000, 30), (11000, 60), (50000, 40 * ecm_effort))
                    elif d <= 600: plan = ((2000, 20), (11000, 30 * ecm_effort))
                    else: plan = ((2000, 10 * ecm_effort),)
                elif level == 1:
                    if d <= 75: plan = ((1000000, 500),)
                    elif d <= 220: plan = ((50000, 150), (250000, 100))
                    elif d <= 600: plan = ((11000, 60), (50000, 60))
                    else: plan = ((2000, 30), (11000, 30))
                else:
                    if d <= 220: plan = ((250000, 300), (1000000, 200))
                    elif d <= 600: plan = ((50000, 200), (250000, 80))
                    else: plan = ((11000, 80), (50000, 40))
                primes, comps = factor_partial(V, plan)
                eng.set_ecm_level(V, level + 1)
                log(f'     ecm[L{level}] cof of {Q} ({d} digits): primes {primes} leftover digits {[len(str(c)) for c in comps]}')
                for q in primes:
                    eng.add_extra(Q, q); progress = True
                if progress: break
            if progress: break
        if not progress: return None
    return None

def prime_factors(n):
    f = []; m = n; d = 3
    while d * d <= m:
        if m % d == 0:
            f.append(d)
            while m % d == 0: m //= d
        d += 2
    if m > 1: f.append(m)
    return f

def vr(n, r):
    e = 0
    while n % r == 0: n //= r; e += 1
    return e

def main():
    n = int(sys.argv[1])
    X = int(float(sys.argv[3])) if len(sys.argv) > 3 else 10**9
    assert n % 2 == 1 and n >= 3
    pf = prime_factors(n)
    good = [q for q in pf if q % 4 == 3]
    if len(sys.argv) > 2 and sys.argv[2] != '-': r = int(sys.argv[2])
    elif good: r = min(good)
    else: r = min(pf)
    assert r in pf
    partB = not (r % 4 == 3 or len(pf) == 1)      # r could be the special prime: needs a second prime factor
    e = vr(n, r); need = (n - 1) // e + 1
    eng = Engine(r, X)
    out = {'n': n, 'beta': (n - 1) // 2, 'r': r, 'e': e, 'need': need, 'X': X, 'scenarios': []}
    t0 = time.time()
    G = run_scenario(eng, r, [r], None, need)
    if G is None:
        print(f'n={n}: generic scenario FAILED', flush=True); return 1
    out['scenarios'].append(G)
    specials = [int(q) for q in G['E'] if int(q) % 4 == 1]
    print(f'n={n}: generic ok, count {G["count"]}/{need}; special candidates: {len(specials)}', flush=True)
    for Qs in specials:
        primes, comps = factor_partial((Qs + 1) // 2, ((2000, 40), (11000, 100), (50000, 300), (250000, 600)))
        roots = [r] + [q for q in primes if q != r]
        sc = run_scenario(eng, r, roots, Qs, need)
        if sc is None:
            print(f'n={n}: scenario p={Qs} FAILED (roots {roots})', flush=True); return 1
        sc['half_factors'] = [str(q) for q in primes]
        out['scenarios'].append(sc)
    if partB:
        # case p = r: use another prime factor r2 of n (r2 is then non-special and divides N)
        r2 = min(q for q in pf if q != r)
        e2 = vr(n, r2); need2 = (n - 1) // e2 + 1
        eng2 = Engine(r2, X)
        primes, comps = factor_partial((r + 1) // 2)
        roots = [r2] + [q for q in primes if q != r2]
        sc = run_scenario(eng2, r2, roots, r, need2)
        if sc is None:
            print(f'n={n}: part B (p={r}, r2={r2}) FAILED', flush=True); return 1
        sc['half_factors'] = [str(q) for q in primes]
        out['partB'] = {'r': r2, 'e': e2, 'need': need2, 'scenario': sc}
    fn = os.path.join(HERE, 'certs', f'cert_n{n}.json')
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump(out, open(fn + '.tmp', 'w'), indent=1); os.replace(fn + '.tmp', fn)
    print(f'n={n}: ALL SCENARIOS CLOSED ({len(out["scenarios"])} scenarios{" + part B" if partB else ""}, {time.time()-t0:.0f}s) -> {fn}', flush=True)
    return 0

if __name__ == '__main__':
    sys.exit(main())
