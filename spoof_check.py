#!/usr/bin/env python
"""Barrier test: replay the n = 3 certificate (Steuerwald's theorem: no OPN p^a * M^2) against Descartes' spoof
D = 3^2 * 7^2 * 11^2 * 13^2 * 22021^1 (22021 = 19^2 * 61 treated as a 'prime').  D satisfies the sigma-formal
equation and all exponent hypotheses of the theorem, so a proof that never used the primality of the bases would
also 'prove' that D does not exist.  We locate the exact inference of the certificate that is false for D."""
import json
bases = {3: 2, 7: 2, 11: 2, 13: 2, 22021: 1}
N = 1; sig = 1
for b, e in bases.items():
    N *= b**e; sig *= (b**(e + 1) - 1) // (b - 1)
print('D =', N, ' spoof sigma = 2D:', sig == 2 * N)
r = 3
S = [b for b, e in bases.items() if e == 2 and b % r == 1]
print('bases == 1 mod 3 with exponent 2:', S, '(bound l-1 = 2 is respected: v_3(D) = 2)')
c = json.load(open('certs/cert_n3.json'))
sc = c['scenarios'][0]       # generic scenario: special prime not among the expanded primes (22021 is not)
print('replaying generic scenario of cert_n3 against D:')
for st in sc['steps']:
    Q = int(st['Q'])
    for f in st['factors']:
        f = int(f)
        divides = N % f == 0
        is_base = f in bases
        status = 'OK (a base of D)' if is_base else ('divides D but is NOT a base: hidden inside the fake prime 22021' if divides else 'does not divide D')
        print(f'  step {Q} -> {f}: certificate infers "{f} is a prime factor of N, hence a base (in S or = p)"; for D: {status}')
        if not is_base:
            print('  => this inference (a prime dividing N is one of the bases) is the place where primality of the bases is used.')
            raise SystemExit
