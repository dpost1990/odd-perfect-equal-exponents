#!/usr/bin/env python
"""Independent referee checker for the equal-exponent OPN certificates (certs/).

Written from the mathematics, NOT from the author's verify.py (which was not imported or copied).

Claim being checked, for one certificate with budget prime r | n (e = v_r(n)):
  Let N = p^a * M^(n-1) be an OPN (M squarefree, p special).  Then
   (K) r | N  (Kanold 1941: t | 2b_i+1 for all i  =>  t^4 | N),
   (NS) r != p  (r == 3 mod 4; or n a power of r [Kanold]; or else part B handles p = r),
   so v_r(N) = n-1, and the set S of non-special primes Q == 1 (mod r) dividing N has |S| <= (n-1)/e.
  Scenario G (p not among G's expanded primes) and one scenario per G-expanded prime Q* == 1 (mod 4), Q* != r.
  In each scenario we must exhibit > (n-1)/e distinct elements of S.

Rules re-derived here:
  * expanded prime Q: proven prime, != special, and either a root or a factor listed in an EARLIER step;
    then Phi_r(Q) | N.
  * roots: r (non-special, divides N), and in scenario p = Q*: divisors of (Q*+1)/2 (they divide N, are != p).
  * every listed step factor f must divide Phi_r(Q).
  * E := expanded primes == 1 (mod r) (never r, never the special prime) -> members of S.
  * items: integers > 1, each dividing Phi_r(Q) for some expanded Q, coprime to r, to the special prime
    (if fixed), to every element of E, and pairwise coprime.  All their prime factors are == 1 (mod r) and divide N,
    so each prime factor is p or in S.  weight 1 always; weight 2 only if provably composite (failed SPRP/Lucas)
    and not a perfect power (=> >= 2 distinct primes).
  * count = |E| + sum(weights) - [scenario G]  (in G the special prime may hide in one item).
"""
import sys, json, math, os, time
sys.set_int_max_str_digits(0)
import gmpy2
from gmpy2 import mpz

def vr(n, r):
    e = 0
    while n % r == 0:
        n //= r; e += 1
    return e

def prime_factors_small(n):
    out = []; d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: out.append(n)
    return out

def phi_r(r, Q):
    Q = mpz(Q)
    return (Q ** r - 1) // (Q - 1)

def is_prime_small(x):
    return x > 1 and all(x % d for d in range(2, int(math.isqrt(x)) + 1))

class Fail(Exception):
    pass

def check_scenario(sc, r, need, special, roots_allowed, prime_obligations, tag):
    """Returns (count, details).  special: int or None (scenario G)."""
    roots = [int(x) for x in sc['roots']]
    for x in roots:
        if x not in roots_allowed:
            raise Fail(f'{tag}: root {x} not justified')
    if special is not None and special in roots:
        raise Fail(f'{tag}: special prime used as root')
    known = set(roots)           # integers known to divide N (roots + listed step factors)
    expanded = []
    phis = {}
    for st in sc['steps']:
        Q = int(st['Q'])
        if Q not in known:
            raise Fail(f'{tag}: expanded {Q} not derived from an earlier step/root')
        if special is not None and Q == special:
            raise Fail(f'{tag}: special prime expanded')
        if Q in phis:
            raise Fail(f'{tag}: {Q} expanded twice')
        prime_obligations.add(Q)
        F = phi_r(r, Q)
        phis[Q] = F
        for fs in st['factors']:
            f = int(fs)
            if f <= 1 or F % f != 0:
                raise Fail(f'{tag}: factor {f} does not divide Phi_{r}({Q})')
            known.add(f)
        expanded.append(Q)
    E = [Q for Q in expanded if Q % r == 1]
    assert r not in E
    if special is not None:
        assert special not in E
    certE = set(int(x) for x in sc.get('E', []))
    if not certE <= set(E):
        raise Fail(f'{tag}: certificate E contains non-expanded/invalid members {certE - set(E)}')
    items = [int(x) for x in sc['items']]
    # each item divides some Phi_r(Q), Q expanded
    for a in items:
        if a <= 1:
            raise Fail(f'{tag}: item <= 1')
        if not any(F % a == 0 for F in phis.values()):
            # allowed alternative: an unexpanded root (divisor of (p+1)/2) that is PRIME and == 1 (mod r):
            # it divides N, is not p, so it lies in S.  Primality becomes a proof obligation.
            if a in roots and a % r == 1:
                prime_obligations.add(a)
            else:
                raise Fail(f'{tag}: item {str(a)[:30]} divides no Phi_r(Q) of an expanded Q')
        if math.gcd(a, r) != 1:
            raise Fail(f'{tag}: item not coprime to r')
        if special is not None and math.gcd(a, special) != 1:
            raise Fail(f'{tag}: item not coprime to special prime')
        for q in E:
            if a % q == 0:      # q prime, so coprime <=> not divisible
                raise Fail(f'{tag}: item divisible by E member {q}')
    # pairwise coprime: product-tree free O(k^2) gcds
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if gmpy2.gcd(items[i], items[j]) != 1:
                raise Fail(f'{tag}: items {i},{j} not coprime')
    w = 0; nw2 = 0
    for a in items:
        A = mpz(a)
        composite = not gmpy2.is_strong_prp(A, 2) or not gmpy2.is_strong_prp(A, 3)
        if composite and not gmpy2.is_power(A):
            w += 2; nw2 += 1
        else:
            w += 1
    cnt = len(E) + w - (1 if special is None else 0)
    if cnt < need:
        raise Fail(f'{tag}: count {cnt} < need {need}')
    return cnt, {'E': len(E), 'items': len(items), 'w2': nw2, 'expanded': expanded}

def check_cert(c):
    obligations = set()
    out = {}
    cls = bool(c.get('class'))
    r = int(c['r'])
    if not is_prime_small(r):
        raise Fail('r not prime')
    if cls:
        if r % 4 != 3:
            raise Fail('class cert requires r == 3 mod 4')
        need = int(c['need'])
        n = None
    else:
        n = int(c['n'])
        if n % 2 == 0 or n < 3 or n % r:
            raise Fail('bad n / r')
        e = vr(n, r)
        need = (n - 1) // e + 1
        if int(c['need']) != need or int(c['e']) != e:
            raise Fail('need/e mismatch')
        pf = prime_factors_small(n)
        partB_needed = not (r % 4 == 3 or len(pf) == 1)
    scs = c['scenarios']
    G = [s for s in scs if s['special'] is None]
    if len(G) != 1:
        raise Fail('need exactly one generic scenario')
    G = G[0]
    cntG, dG = check_scenario(G, r, need, None, {r}, obligations, 'G')
    if [int(x) for x in G['roots']] != [r]:
        raise Fail('G roots must be [r]')
    specials_needed = sorted(Q for Q in dG['expanded'] if Q % 4 == 1 and Q != r)
    have = {}
    for s in scs:
        if s['special'] is None: continue
        Qs = int(s['special'])
        half = (Qs + 1) // 2
        hf = [int(x) for x in s.get('half_factors', [])]
        for x in hf:
            if x <= 1 or half % x:
                raise Fail(f'p={Qs}: half factor {x} does not divide (p+1)/2')
        allowed = {r} | set(hf)
        cnt, d = check_scenario(s, r, need, Qs, allowed, obligations, f'p={Qs}')
        have[Qs] = cnt
    missing = [Q for Q in specials_needed if Q not in have]
    if missing:
        raise Fail(f'uncovered special candidates {missing[:5]}')
    if r % 4 == 1 and not cls:
        if partB_needed:
            if 'partB' not in c:
                raise Fail('part B missing')
            pb = c['partB']; r2 = int(pb['r'])
            if r2 == r or n % r2 or not is_prime_small(r2):
                raise Fail('bad r2')
            e2 = vr(n, r2); need2 = (n - 1) // e2 + 1
            if int(pb['need']) != need2:
                raise Fail('need2 mismatch')
            s = pb['scenario']
            if int(s['special']) != r:
                raise Fail('part B special must be r')
            hf = [int(x) for x in s.get('half_factors', [])]
            for x in hf:
                if x <= 1 or ((r + 1) // 2) % x:
                    raise Fail('part B half factor')
            cntB, _ = check_scenario(s, r2, need2, r, {r2} | set(hf), obligations, 'partB')
            out['partB'] = cntB
        # else: n = r^k -> p != r by Kanold's lemma (all 2b_i+1 powers of r => r non-special)
    out.update({'n': n, 'r': r, 'need': need, 'G': cntG, 'scenarios': len(scs), 'specials_needed': len(specials_needed),
                'min_count': min([cntG] + list(have.values()))})
    return out, obligations

if __name__ == '__main__':
    import glob
    certdir = sys.argv[1]
    files = sys.argv[2:] or sorted(glob.glob(os.path.join(certdir, '*.json')))
    allob = set()
    res = {}
    for f in files:
        t0 = time.time()
        c = json.load(open(f))
        try:
            o, ob = check_cert(c)
            allob |= ob
            res[os.path.basename(f)] = {'ok': True, **o}
            print(f'{os.path.basename(f)}: STRUCTURE OK {o} ({time.time()-t0:.1f}s)', flush=True)
        except Fail as ex:
            res[os.path.basename(f)] = {'ok': False, 'why': str(ex)}
            print(f'{os.path.basename(f)}: REJECTED {ex}', flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    tag = os.environ.get('TAG', 'all')
    json.dump(res, open(os.path.join(here, f'rcheck_{tag}.json'), 'w'), indent=1)
    with open(os.path.join(here, f'prime_obligations_{tag}.txt'), 'w') as fh:
        for q in sorted(allob):
            fh.write(f'{q}\n')
    print('prime obligations:', len(allob))
