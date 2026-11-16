#!/usr/bin/env python
"""Numerical sanity checks of the lemmas used in report.md (not a proof; guards against mis-stated lemmas)."""
import random, math, json
from sympy import factorint, isprime, primerange
random.seed(1)
def sigma_pp(q, a): return (q**(a + 1) - 1) // (q - 1)
def v(n, r):
    e = 0
    while n % r == 0: n //= r; e += 1
    return e
ok = True
# L1: for an odd prime r | n and a prime Q == 1 (mod r), Q != r:  v_r(sigma(Q^(n-1))) = v_r(n).
#     for a prime q != r, q != 1 (mod r): v_r(sigma(q^(n-1))) > 0 iff 1 < ord_r(q) | n   (only >= 0 is used)
cnt = 0
for r in (3, 5, 7, 11, 13, 19):
    for n in (r, r * r, 3 * r, 5 * r, 7 * r, r * r * 5, r**3):
        if n % 2 == 0: continue
        qs = [q for q in primerange(3, 4000) if q % r == 1][:25]
        for Q in qs:
            cnt += 1
            if v(sigma_pp(Q, n - 1), r) != v(n, r): ok = False; print('L1 FAIL', r, n, Q)
print('L1 checked', cnt, 'cases:', ok)
# L2: prime factors of Phi_r(q), q prime (q may equal r): r (iff q == 1 mod r, exactly once) or == 1 (mod r)
cnt = 0
for r in (3, 5, 7, 11, 13):
    for q in primerange(2, 400):
        val = (q**r - 1) // (q - 1)
        f = factorint(val)
        cnt += 1
        for u, m in f.items():
            if u == r:
                if not (q % r == 1 and m == 1): ok = False; print('L2 FAIL r-part', r, q)
            elif u % r != 1: ok = False; print('L2 FAIL', r, q, u)
        if q % r == 1 and r not in f: ok = False; print('L2 FAIL missing r', r, q)
print('L2 checked', cnt, 'cases:', ok)
# L3 (special prime != r when n = r^e): for primes q != r:
#     q^(n-1) == 1 (mod r);  sigma(q^(n-1)) / r^(e*[q==1 mod r]) == 1 (mod r);  sigma(r^a) == 1 (mod r)
cnt = 0
for r in (3, 5, 7, 13, 17):
    for e in (1, 2, 3):
        n = r**e
        for q in primerange(2, 300):
            if q == r: continue
            cnt += 1
            s = sigma_pp(q, n - 1)
            if q % r == 1:
                if s % r**e or (s // r**e) % r != 1: ok = False; print('L3 FAIL a', r, e, q)
            elif s % r != 1: ok = False; print('L3 FAIL b', r, e, q)
            if pow(q, n - 1, r) != 1: ok = False; print('L3 FAIL c', r, e, q)
print('L3 checked', cnt, 'cases:', ok)
# L4 (Kanold gap): if q1 < q2 are primes with Phi_r(q1) = p^c1 and Phi_r(q2) = p^c2 for the same prime p, then q2 > Phi_r(q1).
#     Search r in {3,5,7}, q < 2*10^5; report pairs sharing p.
from collections import defaultdict
for r in (3, 5, 7):
    byp = defaultdict(list)
    for q in primerange(2, 200000):
        val = (q**r - 1) // (q - 1)
        # prime power test
        for c in range(1, 60):
            rt = round(val ** (1.0 / c))
            hit = None
            for cand in (rt - 1, rt, rt + 1):
                if cand > 1 and cand**c == val and isprime(cand): hit = cand
            if hit: byp[hit].append((q, c)); break
            if rt < 3: break
    shared = {p: L for p, L in byp.items() if len(L) > 1}
    for p, L in shared.items():
        L.sort()
        for (q1, c1), (q2, c2) in zip(L, L[1:]):
            if not q2 > (q1**r - 1) // (q1 - 1): ok = False; print('L4 FAIL', r, p, L)
    print(f'L4 r={r}: {len(byp)} primes q<2e5 with Phi_r(q) a prime power; primes p hit more than once: {len(shared)}')
# L5: (p+1)/2 | sigma(p^a)/2 for a == 1 (mod 4), p == 1 (mod 4)
for p in (5, 13, 17, 29, 37):
    for a in (1, 5, 9):
        if (sigma_pp(p, a) // 2) % ((p + 1) // 2): ok = False; print('L5 FAIL', p, a)
print('L5 ok')
print('ALL LEMMA CHECKS PASSED' if ok else 'SOME LEMMA CHECK FAILED')
