import sys as _sys; _sys.set_int_max_str_digits(0)   # certificates contain integers with > 4300 digits
import subprocess, gmpy2, math
_SMALL = [p for p in range(3, 20000, 2) if gmpy2.is_prime(p)]

def ecm_once(n, B1, curves, timeout=900):
    """run GMP-ECM; return a nontrivial factor of n or None."""
    try:
        out = subprocess.run(['nice', '-n', '10', 'ecm', '-q', '-c', str(curves), str(B1)], input=str(n) + '\n',
                             capture_output=True, text=True, timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        return None
    for tok in out.split():
        if tok.isdigit():
            f = int(tok)
            if 1 < f < n and n % f == 0:
                return f
    return None

def rho(n, limit=200000):
    n = int(n)
    for c in (1, 3, 5, 7):
        x = y = 2; d = 1; it = 0
        while d == 1 and it < limit:
            x = (x * x + c) % n; y = (y * y + c) % n; y = (y * y + c) % n
            d = math.gcd(abs(x - y), n); it += 1
        if 1 < d < n: return d
    return None

def factor_partial(n, plan=((2000, 30), (11000, 60), (50000, 100)), timeout=900):
    """returns (sorted list of distinct prp factors, list of composite leftovers)"""
    n = int(n); primes = set()
    for p in _SMALL:
        if p * p > n: break
        while n % p == 0:
            primes.add(p); n //= p
    todo = [n]; comps = []
    while todo:
        m = todo.pop()
        if m == 1: continue
        if gmpy2.is_prime(m): primes.add(m); continue
        r, exact = gmpy2.iroot(m, 2)
        if exact: todo += [int(r)]; continue
        if gmpy2.is_power(m):
            for k in range(3, m.bit_length()):
                r, exact = gmpy2.iroot(m, k)
                if exact: todo += [int(r)]; break
            continue
        f = rho(m) if m < 10**30 else None
        if f is None:
            for B1, c in plan:
                f = ecm_once(m, B1, c, timeout)
                if f: break
        if f is None: comps.append(m)
        else: todo += [int(f), m // int(f)]
    return sorted(primes), comps
