#!/usr/bin/env python
"""Equal-exponent closure engine.

Setting: N = p^alpha * prod q_i^(l-1) odd perfect, l = 2*beta+1 prime, l | N as a non-special prime.
For a non-special prime q | N:  Phi_l(q) | 2N and every prime factor of Phi_l(q) other than l is == 1 (mod l)
and divides N, so it is the special prime p or lies in S = {non-special primes == 1 mod l}; |S| <= l-1.

expand(Q, X): prime factors q <= X of Phi_l(Q) (found by ./cycfac), plus the remaining cofactor.
Results cached in cache_l<l>.json.
"""
import sys as _sys; _sys.set_int_max_str_digits(0)   # certificates contain integers with > 4300 digits
import json, os, subprocess, sys, math
import gmpy2
from gmpy2 import mpz

HERE = os.path.dirname(os.path.abspath(__file__))

def phi_l(l, Q):
    Q = mpz(Q)
    return (Q**l - 1) // (Q - 1)

class Engine:
    def __init__(self, l, X=10**10, verbose=True):
        self.l = l; self.X = X; self.verbose = verbose
        self.cache_file = os.path.join(HERE, f'cache_l{l}.json')
        self.cache = json.load(open(self.cache_file)) if os.path.exists(self.cache_file) else {}

    def save(self):
        tmp = self.cache_file + '.tmp'
        json.dump(self.cache, open(tmp, 'w'), indent=0)
        os.replace(tmp, self.cache_file)

    def expand(self, Q, X=None):
        """returns (factors: list of (prime, mult)), cofactor) of Phi_l(Q) with l removed."""
        X = X or self.X
        key = str(Q)
        ent = self.cache.get(key)
        if ent is None and phi_l(self.l, Q) < 10**48:
            # fast path: small value, factor it completely (rho / ECM); leftovers stay in the cofactor
            from ecmtools import factor_partial
            V0 = phi_l(self.l, Q)
            while V0 % self.l == 0: V0 //= self.l
            primes, comps = factor_partial(V0, ((2000, 40), (11000, 100), (50000, 300), (250000, 600)))
            ent = {'X': 10**30, 'found': [], 'extra': [str(q) for q in primes]}
            self.cache[key] = ent
            self.save()
        if ent is None or ent['X'] < X:
            lo = 1 if ent is None else ent['X'] + 1
            out = subprocess.run(['nice', '-n', '10', os.path.join(HERE, 'cycfac'), str(self.l), str(Q), str(lo), str(X)],
                                 capture_output=True, text=True, check=True).stdout.split()
            found = [] if ent is None else list(ent['found'])
            extra = [] if ent is None else list(ent.get('extra', []))
            for s in out:
                q = int(s)
                if gmpy2.is_prime(q) and q not in found:
                    found.append(q)
            ent = {'X': X, 'found': found, 'extra': extra}
            self.cache[key] = ent
            self.save()
        V = phi_l(self.l, Q)
        if V % self.l == 0:
            V //= self.l
            assert V % self.l != 0
        facs = []
        for q in sorted(set(ent['found']) | set(int(x) for x in ent.get('extra', []))):
            q = mpz(q); m = 0
            while V % q == 0:
                V //= q; m += 1
            if m: facs.append((int(q), m))
        return facs, int(V)

    def ecm_level(self, V):
        return self.cache.get('__ecm__', {}).get(f'{len(str(V))}:{V % 10**40}', 0)

    def set_ecm_level(self, V, lev):
        self.cache.setdefault('__ecm__', {})[f'{len(str(V))}:{V % 10**40}'] = lev
        self.save()

    def add_extra(self, Q, q):
        """register an extra (e.g. ECM-found) factor q of Phi_l(Q)."""
        assert phi_l(self.l, Q) % q == 0
        ent = self.cache[str(Q)]
        ent.setdefault('extra', [])
        if str(q) not in ent['extra']:
            ent['extra'].append(str(q)); self.save()

def closure(eng, roots, special=None, need=None, maxexp=10**40, maxnodes=400, log=None):
    """Breadth-first closure, smallest unexpanded prime first.
    roots: non-special primes known to divide N.  special: the special prime (never expanded) or None.
    Returns dict with expanded S-primes E, unexpanded coprime items U, count."""
    l = eng.l
    known = {}            # prime -> parent
    for r in roots: known[int(r)] = 'root'
    expanded = []         # in order
    cof = {}              # Q -> cofactor (>1)
    while True:
        cand = [q for q in known if q not in expanded and q != special and q < maxexp and gmpy2.is_prime(q)]
        E = [q for q in expanded if q % l == 1]
        if not cand or len(expanded) >= maxnodes: break
        # count so far
        Uprimes = [q for q in known if q not in expanded and q != special and q % l == 1]
        cnt = len(E) + len(Uprimes) + len(cof) - (0 if special else 1)
        if need and cnt >= need: break
        Q = min(cand)
        facs, V = eng.expand(Q)
        expanded.append(Q)
        if log is not None: log.append((Q, facs, V))
        for q, m in facs:
            if q != special and q not in known: known[q] = Q
        if V > 1:
            if gmpy2.is_prime(V):
                if V != special and V not in known: known[V] = Q
            else:
                cof[Q] = V
    # sanitize unexpanded items: pairwise coprime, coprime to known primes
    E = [q for q in expanded if q % l == 1]
    Uprimes = [q for q in known if q not in expanded and q != special and q % l == 1]
    items = list(Uprimes) + list(cof.values())
    allknown = set(known) | ({special} if special else set())
    ok = True
    for i, a in enumerate(items):
        if i >= len(Uprimes):
            for q in allknown:
                if a % q == 0: ok = False
        for b in items[i+1:]:
            if math.gcd(a, b) != 1: ok = False
    cnt = len(E) + len(items) - (0 if special else 1)
    return {'E': E, 'Uprimes': Uprimes, 'cof': cof, 'count': cnt, 'coprime_ok': ok, 'known': known, 'expanded': expanded}

if __name__ == '__main__':
    l = int(sys.argv[1]); X = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**9
    eng = Engine(l, X)
    log = []
    res = closure(eng, [l], need=l + 3, log=log)
    for Q, facs, V in log:
        print(Q, '->', facs, 'cof digits', len(str(V)) if V > 1 else 0, 'prp' if V > 1 and gmpy2.is_prime(V) else '')
    print('E', len(res['E']), 'Uprimes', len(res['Uprimes']), 'cof', len(res['cof']), 'count', res['count'], 'coprime', res['coprime_ok'])
