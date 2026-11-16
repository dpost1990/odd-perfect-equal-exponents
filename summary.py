#!/usr/bin/env python
"""summary.py: table of all certificates (reads certs/*.json and logs/verify_all.out)."""
import json, glob, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
ver = {}; cls = {}
for fn in glob.glob(os.path.join(HERE, 'logs', 'verify_final*.out')):
    for line in open(fn):
        m = re.match(r'n=(\d+) .*: (VERIFIED|REJECTED)', line)
        if m: ver[int(m.group(1))] = m.group(2)
        m = re.match(r'n=class\(r=(\d+)\) .* need=(\d+) scenarios=(\d+) expansions=(\d+) .*: (VERIFIED|REJECTED)', line)
        if m: cls[int(m.group(1))] = (int(m.group(2)), int(m.group(3)), int(m.group(4)), m.group(5))
KNOWN = {3: 'Steuerwald 1937', 5: 'Kanold 1941', 7: 'Hagis-McDaniel 1972', 11: 'McDaniel-Hagis 1975', 25: 'McDaniel-Hagis 1975',
         35: 'McDaniel-Hagis 1975', 49: 'McDaniel-Hagis 1975', 125: 'McDaniel-Hagis 1975', 13: 'Cohen-Williams 1985',
         17: 'Cohen-Williams 1985', 23: 'Cohen-Williams 1985', 29: 'Cohen-Williams 1985', 37: 'Cohen-Williams 1985'}
def prior(n):
    if n in KNOWN: return KNOWN[n]
    if n % 3 == 0: return 'McDaniel 1970 (3 | 2b+1)'
    if n % 5 == 0: return 'Brown 2015, arXiv:1512.01270 (5 | 2b+1)'
    return None
rows = []
for fn in glob.glob(os.path.join(HERE, 'certs', 'cert_n*.json')):
    c = json.load(open(fn)); n = c['n']
    exp = sum(len(s['steps']) for s in c['scenarios']) + (len(c['partB']['scenario']['steps']) if 'partB' in c else 0)
    rows.append((n, c['beta'], c['r'], c['need'], len(c['scenarios']) + (1 if 'partB' in c else 0), exp, ver.get(n, 'not yet'), prior(n)))
rows.sort()
if __name__ == '__main__':
    print('| n = 2b+1 | b | r | S-members needed | scenarios | expansions | verify.py | prior result |')
    print('|---|---|---|---|---|---|---|---|')
    for n, b, r, need, ns, ex, v, pr in rows:
        print(f'| {n} | {b} | {r} | {need} | {ns} | {ex} | {v} | {pr or "**none found (new)**"} |')
    new = [b for n, b, r, need, ns, ex, v, pr in rows if pr is None and v == 'VERIFIED']
    allv = sorted(n for n, b, r, need, ns, ex, v, pr in rows if v == 'VERIFIED')
    print('\nnew beta (verified):', new)
    print('\n| r | certified: v_r(N) >= | hence v_r(N) >= (v == r-1 mod 2r) | scenarios | expansions | verify.py |')
    print('|---|---|---|---|---|---|')
    for r in sorted(cls):
        T, ns, ex, v = cls[r]
        low = T + ((r - 1 - T) % (2 * r))
        print(f'| {r} | {T} | {low} | {ns} | {ex} | {v} |')
    k = 3
    while k in allv: k += 2
    print('all odd n verified up to', k - 2, '=> all beta <=', (k - 3) // 2)
