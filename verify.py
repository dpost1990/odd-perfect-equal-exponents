#!/usr/bin/env python
"""verify.py certs/cert_n<n>.json ...   Independent checker for the equal-exponent certificates.

Second code path: pure Python integers only in the trusted part (no gmpy2): own Miller-Rabin, own Pocklington
(N-1) primality proofs, own evaluation of Phi_r by Horner's rule.  Untrusted helpers (ecmtools) are used only
to FIND factors of Q-1 for the Pocklington proofs; every factor found is re-checked here.

What is checked, per certificate (n, r):
  V0  n odd, r a prime factor of n; e = v_r(n), need = floor((n-1)/e) + 1.  r == 3 (mod 4) or n a power of r;
      otherwise the certificate must contain 'partB': the case p = r, treated with another prime factor r2 of n
      (one scenario with special prime r, roots r2 and the primes of (r+1)/2, need2 = floor((n-1)/v_r2(n)) + 1).
  V1  scenario 0 is generic; every expanded prime Q == 1 (mod 4) of scenario 0 has its own scenario p = Q.
  per scenario:
  V2  every root is r, or (special scenario) a prime dividing (p+1)/2.
  V3  every expanded Q: is a root or a listed factor of an earlier expansion; Q != p; Q is PROVEN prime.
  V4  every listed factor f of an expansion of Q divides Phi_r(Q) and gcd(f, r) = 1.
  V5  every item a: a > 1, divides Phi_r(Q) for some expanded Q (or: a is a root that is a proven prime
      == 1 mod r, i.e. a member of S coming from (p+1)/2), gcd(a, r) = 1, gcd(a, p) = 1,
      gcd(a, Q') = 1 for all expanded Q', items pairwise coprime.
  V6  count = #(expanded Q == 1 mod r) + sum of weights - [generic] >= need, where weight(a) = 2 if a fails a
      strong pseudoprime test (hence composite) and is not a perfect power, else 1.
"""
import sys as _sys; _sys.set_int_max_str_digits(0)   # certificates contain integers with > 4300 digits
import sys, json, math, os
sys.setrecursionlimit(10000)
MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
MR_LIMIT = 3317044064679887385961981   # Sorenson-Webster 2017: the 13 bases above are a proof below this

def strong_prp(n, a):
    if n % a == 0: return n == a
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    x = pow(a, d, n)
    if x == 1 or x == n - 1: return True
    for _ in range(s - 1):
        x = x * x % n
        if x == n - 1: return True
    return False

def mr_all(n):
    if n < 2: return False
    for a in MR_BASES:
        if n == a: return True
        if n % a == 0: return False
    return all(strong_prp(n, a) for a in MR_BASES)

_proven = {}
def find_factors(m):
    """untrusted: return list of probable primes whose product (with multiplicity ignored) covers as much of m as possible"""
    from ecmtools import factor_partial
    primes, comps = factor_partial(m, ((2000, 40), (11000, 100), (50000, 300), (250000, 700), (1000000, 1200)))
    return [int(p) for p in primes]

def prove_prime(n, depth=0):
    """True only if n is proven prime (deterministic MR below MR_LIMIT, Pocklington above)."""
    n = int(n)
    if n in _proven: return _proven[n]
    if not mr_all(n):
        _proven[n] = False; return False
    if n < MR_LIMIT:
        _proven[n] = True; return True
    # Pocklington: n-1 = F*R, F^2 > n, all primes f | F proven, witnesses a_f.
    fs = sorted(set(find_factors(n - 1)), reverse=True)
    F = 1; used = []
    m = n - 1
    for f in sorted(fs):          # small factors first (cheap), stop when F^2 > n
        if m % f: continue
        if not prove_prime(f, depth + 1): continue
        k = 0
        while m % f == 0: m //= f; k += 1
        F *= f ** k; used.append(f)
        if F * F > n: break
    ok = F * F > n
    if ok:
        for f in used:
            good = False
            for a in range(2, 200):
                if pow(a, n - 1, n) != 1: ok = False; break
                if math.gcd(pow(a, (n - 1) // f, n) - 1, n) == 1: good = True; break
            if not good: ok = False
            if not ok: break
    _proven[n] = ok
    return ok

_phi_cache = {}
def phi(r, Q, mod=None):
    """Phi_r(Q) = 1 + Q + ... + Q^(r-1) by Horner; the exact value is cached, then reduced modulo mod."""
    key = (r, Q)
    if key not in _phi_cache:
        acc = 0
        for _ in range(r):
            acc = acc * Q + 1
        _phi_cache[key] = acc
    val = _phi_cache[key]
    return val % mod if mod else val

def iroot(a, k):
    if a < 2: return a
    x = 1 << ((a.bit_length() + k - 1) // k)
    while True:
        y = ((k - 1) * x + a // x ** (k - 1)) // k
        if y >= x: return x
        x = y

def is_perfect_power(a):
    """a = b^k for some k >= 2?  It is enough to test prime exponents k <= log2(a)."""
    L = a.bit_length()
    sieve = [True] * (L + 1)
    for k in range(2, L + 1):
        if not sieve[k]: continue
        for j in range(k * k, L + 1, k): sieve[j] = False
        rt = iroot(a, k)
        if rt < 2: break
        if rt ** k == a: return True
    return False

_weight_cache = {}
def weight(a):
    """2 if a is proven composite (fails a strong pseudoprime test) and is not a perfect power, else 1 (memoised)."""
    if a in _weight_cache: return _weight_cache[a]
    if a < MR_LIMIT and mr_all(a): w = 1
    elif strong_prp(a, 2) and strong_prp(a, 3): w = 1        # possibly prime: count 1
    else: w = 1 if is_perfect_power(a) else 2
    _weight_cache[a] = w
    return w

def is_prime_small(m):
    return m >= 2 and all(m % d for d in range(2, int(m ** 0.5) + 1))

def v_of(n, r):
    e = 0
    while n % r == 0: n //= r; e += 1
    return e, n

def check_scenario(sc, r, need, tag, req, stats):
    """checks V2..V6 for one scenario; returns the list of expanded primes"""
    p = int(sc['special']) if sc['special'] else None
    roots = [int(x) for x in sc['roots']]
    if p is not None:
        req(p % 4 == 1 and prove_prime(p), f'{tag} special prime is prime and 1 mod 4')
    for q in roots:
        if q == r: continue
        req(p is not None and ((p + 1) // 2) % q == 0, f'V2 {tag} root {q}')
    justified = set(roots); expanded = []
    for st in sc['steps']:
        Q = int(st['Q'])
        req(Q in justified, f'V3 {tag} {Q} justified')
        req(Q != p, f'V3 {tag} {Q} is not special')
        req(Q not in expanded, f'V3 {tag} {Q} expanded once')
        req(prove_prime(Q), f'V3 {tag} {Q} proven prime')
        stats['expanded'] += 1; stats['maxQ'] = max(stats['maxQ'], Q)
        if Q >= MR_LIMIT: stats['pock'] += 1
        for fs in st['factors']:
            f = int(fs)
            req(f > 1 and phi(r, Q, f) == 0 and math.gcd(f, r) == 1, f'V4 {tag} {f} | Phi({Q})')
            justified.add(f)
        expanded.append(Q)
    E = [Q for Q in expanded if Q % r == 1]
    req(sorted(E) == sorted(int(x) for x in sc['E']), f'{tag} E list')
    items = [int(x) for x in sc['items']]
    for i, a in enumerate(items):
        req(a > 1 and math.gcd(a, r) == 1, f'V5 {tag} item {i} > 1, coprime to r')
        is_root_member = (a in roots and a != r and a % r == 1 and prove_prime(a))   # prime == 1 mod r dividing (p+1)/2
        req(is_root_member or any(phi(r, Q, a) == 0 for Q in expanded), f'V5 {tag} item {i} divides some Phi(Q)')
        if p: req(math.gcd(a, p) == 1, f'V5 {tag} item {i} coprime to p')
        req(all(math.gcd(a, Q) == 1 for Q in expanded), f'V5 {tag} item {i} coprime to expanded')
        req(all(math.gcd(a, b) == 1 for b in items[i + 1:]), f'V5 {tag} item {i} pairwise coprime')
    cnt = len(set(E)) + sum(weight(a) for a in items) - (0 if p else 1)
    req(cnt >= need, f'V6 {tag} count {cnt} >= need {need}')
    req(cnt == sc['count'], f'V6 {tag} count matches certificate ({cnt} vs {sc["count"]})')
    return expanded

def check(fn):
    c = json.load(open(fn)); n = c['n']; r = c['r']; errs = []
    def req(cond, msg):
        if not cond: errs.append(msg)
    if c.get('class'):
        # class certificate: hypothesis r | 2*beta_i+1 for all i, r == 3 (mod 4); claim: >= need members of S
        req(r % 4 == 3 and is_prime_small(r), 'V0 class: r prime, r == 3 (mod 4)')
        need = c['need']; stats = {'expanded': 0, 'pock': 0, 'maxQ': 0}
        scs = c['scenarios']
        req(scs[0]['special'] is None, 'V1 scenario 0 generic')
        have = {int(s['special']) for s in scs[1:]}
        for si, sc in enumerate(scs):
            if si > 0: req(sc['special'] is not None, f'V1 S{si} has a special prime')
            expanded = check_scenario(sc, r, need, f'S{si}', req, stats)
            if si == 0:
                for Q in expanded:
                    if Q != r and Q % 4 == 1: req(Q in have, f'V1 scenario for p={Q} present')
        return f'class(r={r})', r, need, len(scs), stats, errs
    req(n % 2 == 1 and n >= 3, 'V0 n odd')
    req(r >= 3 and n % r == 0 and is_prime_small(r), 'V0 r is a prime factor of n')
    e, m = v_of(n, r)
    need = (n - 1) // e + 1
    req(c['need'] == need, 'V0 need')
    stats = {'expanded': 0, 'pock': 0, 'maxQ': 0}
    if not (r % 4 == 3 or m == 1):
        # r might be the special prime: part B must treat p = r with another prime factor r2 of n
        B = c.get('partB')
        req(B is not None, 'V0 r == 3 (mod 4), or n a power of r, or part B present')
        if B:
            r2 = B['r']
            req(r2 != r and n % r2 == 0 and is_prime_small(r2), 'V0 part B: r2 is another prime factor of n')
            e2, _ = v_of(n, r2); need2 = (n - 1) // e2 + 1
            req(B['need'] == need2, 'V0 part B need')
            req(B['scenario']['special'] is not None and int(B['scenario']['special']) == r, 'V0 part B: special prime is r')
            check_scenario(B['scenario'], r2, need2, 'B', req, stats)
    scs = c['scenarios']
    req(scs[0]['special'] is None, 'V1 scenario 0 generic')
    specials_have = {int(s['special']) for s in scs[1:]}
    for si, sc in enumerate(scs):
        if si > 0: req(sc['special'] is not None, f'V1 S{si} has a special prime')
        expanded = check_scenario(sc, r, need, f'S{si}', req, stats)
        if si == 0:
            for Q in expanded:
                if Q != r and Q % 4 == 1:
                    req(Q in specials_have, f'V1 scenario for p={Q} present')
    return n, r, need, len(scs), stats, errs

if __name__ == '__main__':
    bad = 0
    for fn in sys.argv[1:]:
        _phi_cache.clear()
        n, r, need, ns, stats, errs = check(fn)
        status = 'VERIFIED' if not errs else 'REJECTED'
        if errs: bad += 1
        print(f'n={n} beta={(n-1)//2 if isinstance(n, int) else "-"} r={r} need={need} scenarios={ns} expansions={stats["expanded"]} pocklington={stats["pock"]} '
              f'maxQ~1e{len(str(stats["maxQ"]))-1} : {status}', flush=True)
        for e_ in errs[:8]: print('    ', e_)
    sys.exit(1 if bad else 0)
