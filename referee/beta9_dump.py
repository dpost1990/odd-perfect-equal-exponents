import json, sys, gmpy2
sys.set_int_max_str_digits(0)
c = json.load(open('../certs/cert_n19.json'))
r = 19
def short(x):
    s = str(x); return s if len(s) <= 30 else f'{s[:8]}...{s[-4:]} ({len(s)} digits)'
for s in c['scenarios']:
    print('### special =', s['special'], ' roots', s['roots'], ' half', s.get('half_factors'))
    for st in s['steps']:
        Q = int(st['Q']); F = (gmpy2.mpz(Q)**r - 1)//(Q-1)
        print(f'- Q={short(Q)} (Q mod 19={Q%19}, Q mod 4={Q%4}, digits Phi={len(str(F))}): ' + ', '.join(short(f) for f in st['factors']))
    E = [int(q) for q in s['E']]
    print('E =', E, len(E))
    tot = 0
    for a in s['items']:
        A = gmpy2.mpz(a); comp = not gmpy2.is_strong_prp(A, 2); pw = gmpy2.is_power(A)
        w = 2 if comp and not pw else 1; tot += w
        src = [st['Q'] for st in s['steps'] if ((gmpy2.mpz(int(st['Q']))**r-1)//(int(st['Q'])-1)) % A == 0]
        print(f'  item {short(a)} w={w} {"composite" if comp else "prp"} from Q={[short(q) for q in src][:2]}')
    print('count =', len(E), '+', tot, '-', 1 if s['special'] is None else 0, '=', len(E) + tot - (1 if s['special'] is None else 0))
