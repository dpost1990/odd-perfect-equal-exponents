#!/usr/bin/env python
"""Prove the primality obligations collected by rcheck.py, with methods different from the author's
(author: deterministic Miller-Rabin below 3.3e24, own recursive Pocklington above).

  screen : gmpy2 BPSW (is_bpsw_prp) AND sympy.isprime (both must say prime)
  proof  : q < 2^64           -> BPSW is a proof (no BPSW pseudoprime below 2^64: Feitsma/Galway table)
           q >= 2^64          -> PARI/GP isprime(q) (BPSW + Pocklington-Lehmer/APR-CL, a rigorous proof),
                                  run in a separate gp process.
Usage: rprime.py obligations.txt [more.txt ...]
"""
import sys, subprocess, os, json
sys.set_int_max_str_digits(0)
import gmpy2, sympy

GP = __import__('shutil').which('gp') or 'gp'
qs = set()
for fn in sys.argv[1:]:
    for line in open(fn):
        line = line.strip()
        if line: qs.add(int(line))
qs = sorted(qs)
bad = []
big = []
for q in qs:
    if not (gmpy2.is_bpsw_prp(q) and sympy.isprime(q)):
        bad.append(q); continue
    if q >= int(__import__("os").environ.get("PARI_FROM", 2**64)):
        big.append(q)
print(f'{len(qs)} obligations, screen failures: {len(bad)}, need PARI proof: {len(big)}', flush=True)
here = os.path.dirname(os.path.abspath(__file__))
inp = os.path.join(here, 'gp_in.txt')
with open(inp, 'w') as fh:
    fh.write('default(parisize, 200000000);\n')
    for q in big:
        fh.write(f'print(isprime({q}));\n')
    fh.write('quit;\n')
out = subprocess.run(['nice', '-n', '10', GP, '-q', '-f', inp], capture_output=True, text=True).stdout.split()
assert len(out) == len(big), (len(out), len(big))
notproven = [q for q, o in zip(big, out) if o != '1']
print('PARI isprime failures:', len(notproven), flush=True)
res = {'obligations': len(qs), 'screen_fail': [str(q) for q in bad], 'pari_checked': len(big),
       'pari_fail': [str(q) for q in notproven], 'max_digits': max(len(str(q)) for q in qs)}
json.dump(res, open(os.path.join(here, 'rprime_result.json'), 'w'), indent=1)
print(res)
