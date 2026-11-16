import json, sys, glob, gmpy2
sys.set_int_max_str_digits(0)
items = set()
for f in glob.glob('../certs/*.json'):
    c = json.load(open(f))
    scs = c['scenarios'] + ([c['partB']['scenario']] if 'partB' in c else [])
    for s in scs:
        for a in s['items']:
            A = gmpy2.mpz(a)
            if (not gmpy2.is_strong_prp(A, 2) or not gmpy2.is_strong_prp(A, 3)) and not gmpy2.is_power(A):
                items.add(int(a))
items = sorted(items)
with open('gp_w2b.txt', 'w') as fh:
    fh.write('default(parisize, 400000000);\nbad=0;\n')
    for a in items:
        fh.write(f'a={a}; if(Mod(2,a)^(a-1)==1 || ispower(a), bad++; print("BAD ",#Str(a)));\n')
    fh.write('print("W2 checked ", %d, " bad ", bad);\nquit;\n' % len(items))
print(len(items), 'weight-2 items, max digits', max(len(str(a)) for a in items))
