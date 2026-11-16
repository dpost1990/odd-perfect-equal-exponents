#!/usr/bin/env python
"""Sensitivity tests for rcheck.py: every mutation below must be REJECTED."""
import sys, json, copy
sys.set_int_max_str_digits(0)
sys.path.insert(0, __import__('os').path.dirname(__file__))
from rcheck import check_cert, Fail

base = json.load(open(sys.argv[1]))
def t(name, mut):
    c = copy.deepcopy(base); mut(c)
    try:
        check_cert(c); print(f'{name}: ACCEPTED  <-- checker too weak'); return False
    except Fail as ex:
        print(f'{name}: rejected ({ex})'); return True
    except Exception as ex:
        print(f'{name}: rejected by exception ({type(ex).__name__}: {ex})'); return True

def drop_last_scenario(c): c['scenarios'].pop()
def drop_item(c):
    s = c['scenarios'][0]
    # remove items until count would fall below need by exactly one weight
    s['items'] = s['items'][:-3]
def dup_item(c):
    s = c['scenarios'][1]; s['items'].append(s['items'][0])
def fake_factor(c):
    c['scenarios'][0]['steps'][1]['factors'].append('1000003')
def orphan_step(c):
    st = c['scenarios'][0]['steps']; st.insert(1, {'Q': '1000033', 'factors': []})
def special_expanded(c):
    s = c['scenarios'][1]; s['steps'].append({'Q': s['special'], 'factors': []})
def bad_root(c):
    s = c['scenarios'][1]; s['roots'].append('3'); s['half_factors'].append('3')
def item_shares_E(c):
    s = c['scenarios'][0]; q = int(s['E'][-1]); s['items'][0] = str(int(s['items'][0]) * q)
def raise_need(c):
    c['need'] += 1
def item_with_r(c):
    s = c['scenarios'][0]; s['items'][0] = str(int(s['items'][0]) * c['r'])
ok = all([t('drop scenario', drop_last_scenario), t('drop items', drop_item), t('duplicate item', dup_item),
          t('fake factor', fake_factor), t('orphan expansion', orphan_step), t('special expanded', special_expanded),
          t('unjustified root', bad_root), t('item shares E prime', item_shares_E), t('need+1', raise_need),
          t('item divisible by r', item_with_r)])
print('ALL TAMPERS REJECTED' if ok else 'SOME TAMPER ACCEPTED')
# targeted extras
def item_is_E(c):
    s = c['scenarios'][0]; s['items'][0] = s['E'][-1]
def special_item(c):
    s = c['scenarios'][2]; s['items'][0] = s['special']
def drop_partB(c):
    c.pop('partB', None); c['_force'] = 1
t('item equals an E prime', item_is_E); t('item equals special', special_item)
