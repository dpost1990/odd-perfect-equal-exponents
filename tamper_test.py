#!/usr/bin/env python
"""Negative controls for verify.py: every tampered certificate must be REJECTED."""
import sys as _sys; _sys.set_int_max_str_digits(0)   # certificates contain integers with > 4300 digits
import json, copy, sys, os, tempfile
import verify
src = sys.argv[1] if len(sys.argv) > 1 else 'certs/cert_n19.json'
c0 = json.load(open(src))
def run(c, name):
    fd, fn = tempfile.mkstemp(suffix='.json', dir=os.path.dirname(os.path.abspath(src))); os.close(fd)
    json.dump(c, open(fn, 'w')); res = verify.check(fn); os.remove(fn)
    print(f'{name:55s} -> {"REJECTED" if res[-1] else "ACCEPTED (BAD!)"}  {res[-1][:1]}')
    return bool(res[-1])
tests = []
c = copy.deepcopy(c0); c['scenarios'].pop(); tests.append((c, 'drop last special scenario'))
c = copy.deepcopy(c0); c['scenarios'][0]['steps'][1]['factors'][0] = str(int(c['scenarios'][0]['steps'][1]['factors'][0]) + 38); tests.append((c, 'wrong factor in a step'))
c = copy.deepcopy(c0); c['scenarios'][0]['items'].append(c['scenarios'][0]['items'][0]); c['scenarios'][0]['count'] += 1; tests.append((c, 'duplicate item'))
c = copy.deepcopy(c0); c['scenarios'][0]['items'].append('191'); c['scenarios'][0]['count'] += 1; tests.append((c, 'item equal to an expanded prime'))
c = copy.deepcopy(c0); c['scenarios'][0]['items'].append(str(2**61 - 1)); c['scenarios'][0]['count'] += 1; tests.append((c, 'unrelated prime as item'))
c = copy.deepcopy(c0); c['scenarios'][0]['count'] += 3; tests.append((c, 'inflated count'))
c = copy.deepcopy(c0); c['need'] -= 1; tests.append((c, 'need lowered'))
c = copy.deepcopy(c0); s = c['scenarios'][1]; s['roots'].append('229'); s['steps'].append({'Q': '229', 'factors': []}); s['E'].append('229'); s['count'] += 1; tests.append((c, 'unjustified root'))
c = copy.deepcopy(c0); s = c['scenarios'][0]; Q = int(s['steps'][-1]['Q']); s['steps'].append({'Q': str(191 * 229), 'factors': []}); tests.append((c, 'unjustified composite expansion'))
c = copy.deepcopy(c0); s = c['scenarios'][1]; s['steps'].append({'Q': s['special'], 'factors': []}); tests.append((c, 'expanding the special prime'))
c = copy.deepcopy(c0); s = c['scenarios'][1]; s['special'] = str(int(s['special']) + 76); tests.append((c, 'special prime changed'))
# item with weight 2 replaced by its square (perfect power) -> weight drops
c = copy.deepcopy(c0)
big = max(range(len(c['scenarios'][0]['items'])), key=lambda i: len(c['scenarios'][0]['items'][i]))
c['scenarios'][0]['items'][big] = str(int(c['scenarios'][0]['items'][big]) * 191); tests.append((c, 'item multiplied by an expanded prime'))
ok = all([run(c, name) for c, name in tests])
n, r, need, ns, stats, errs = verify.check(src)
print('original:', 'VERIFIED' if not errs else errs[:3])
print('ALL TAMPERS REJECTED' if ok and not errs else 'PROBLEM')
