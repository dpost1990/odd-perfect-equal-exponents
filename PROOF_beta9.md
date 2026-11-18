# β = 9: no odd perfect number N = p^α · M^18 (M squarefree)

This is a human-readable version of certificate `certs/cert_n19.json`. It was written out from the certificate, and every number in it was checked with the second, independently written checker `referee/rcheck.py`.
- Every divisibility below was recomputed by `rcheck.py`, which shares no code with the author's `verify.py`.
- Every prime used as an expansion base was proven prime by PARI/GP `isprime` (≥ 2^64) or by BPSW (< 2^64), after a BPSW plus sympy screen.
- Every "composite, not a perfect power" claim was checked by both gmpy2 and PARI (`ispseudoprime = 0`, `ispower = 0`).

Notation: Φ(x) := Φ_19(x) = (x^19 − 1)/(x − 1).
- "Prime items" below are SPRP probable primes. They count with weight 1, and weight 1 needs only "integer > 1", so their primality is never used.
- Primality is used only for the expanded bases Q and the (p+1)/2 roots. P_k denotes a prime with k digits, and C_k a composite with k digits that is not a perfect power. Full values are in the certificate.

## 1. Setting and lemmas

Let N be an odd perfect number with N = p^α · M^18, where:
- M is squarefree and gcd(p, M) = 1;
- p ≡ α ≡ 1 (mod 4) (Euler).

Call p the *special* prime. The primes q | M are the *non-special* primes; each has exponent exactly 18.

**L1.** For every non-special prime q: Φ(q) | σ(q^18) | N.
- σ(q^18) is a sum of 19 odd terms, so it is odd.
- σ(q^18) divides σ(N) = 2N, hence it divides N.
- Φ_19 divides (x^19 − 1)/(x − 1) as polynomials.

**L2.** (p+1)/2 | N, and every prime factor of (p+1)/2 is non-special.
- α is odd, so p + 1 | σ(p^α), and σ(p^α) | 2N.
- Since p ≡ 1 (mod 4), (p+1)/2 is odd.
- p does not divide p + 1.

Nothing here needs α = 1.

**L3.** Let x ≥ 2 be an integer. Every prime factor ℓ of Φ(x) is either:
- 19, which happens iff x ≡ 1 (mod 19), and then 19 ∥ Φ(x); or
- ≡ 1 (mod 19).

This is classical (Nagell, Thms 94–95).

**L4 (Kanold 1941).** 19 | N. Kanold proves t^4 | N whenever t divides every 2β_i + 1; here t = 19. Yamada (Indian J. Pure Appl. Math. 2023) improves this to t^5.

Self-contained check:
1. Suppose 19 ∤ N.
2. Then no non-special prime q ≡ 1 (mod 19) exists, since otherwise 19 | Φ(q) | N.
3. So every Φ(q) has only prime factors ≡ 1 (mod 19) that divide N. These must all be p, so Φ(q) = p^c with c ≥ 1.
4. Take a prime q | (p+1)/2 (L2). Then p ≡ −1 (mod q) and Φ(q) ≡ 1 (mod q), so (−1)^c ≡ 1 (mod q) and c is even.
5. So Φ_19(q) is a square. Ljunggren (1943) shows (x^n − 1)/(x − 1) = y² has only the solutions x = 3, n = 5 and x = 7, n = 4. Neither has n = 19. Contradiction.

**L5.** 19 ≠ p, because 19 ≡ 3 (mod 4).
- By L4 and L5, v_19(N) = 18.
- Since 19 is non-special, Φ(19) | N (L1).

**L6 (budget; McDaniel–Hagis 1975, Lemma).** Let S be the set of non-special primes Q ≡ 1 (mod 19) dividing N. Then |S| ≤ 18.
- For Q ∈ S, lifting the exponent gives v_19(σ(Q^18)) = v_19(19) = 1.
- Every other factor of σ(N) contributes ≥ 0.
- Hence |S| ≤ v_19(σ(N)) = v_19(2N) = 18.

**Counting rule.** It suffices to exhibit 19 distinct members of S.
- Every prime factor of an integer d > 1 with d | Φ(Q) and 19 ∤ d is ≡ 1 (mod 19) (L3).
- If Q is a non-special prime of N, every such prime factor divides N, so it is p or lies in S.
- A composite d that is not a perfect power has ≥ 2 distinct prime factors.
- Pairwise coprime integers have disjoint sets of prime factors.

## 2. Case split on the special prime

Compute the forcing chain G starting from 19, assuming each prime we expand is non-special:

| expanded Q | Φ(Q), completely split (the product is exact) |
|---|---|
| 19 | **P24** = 109912203092239643840221 (prime) |
| P24 | 19 · **682621057** · C405 |
| 682621057 | 19 · **1124041** · C152 |
| 1124041 | 19 · **191** · **7069** · 32407098973 · 26013368986349 · C78 |
| 191 | 19 · 1277351533 · P31 |
| 7069 | 19 · **827360473** · 1341646577593559 · 14307862907610435557 · 6439376875503349484048791 |
| 827360473 | 19 · C160 |

The expanded primes other than 19 are P24, 682621057, 1124041, 191, 7069 and 827360473. All are ≡ 1 (mod 19), and all except 191 are ≡ 1 (mod 4).

Either p is one of the five candidates {P24, 682621057, 1124041, 7069, 827360473}, or it is not. The case p = 19 is impossible (L5). A prime ≡ 3 (mod 4) such as 191 cannot be p.

### Case G: p is none of the expanded primes

Then every expanded prime is non-special, so all the Φ(Q) above divide N.
- **S-members from expansions:** 6 (P24, 682621057, 1124041, 191, 7069, 827360473).
- **Further prime items:** 7 (32407098973, 26013368986349, 1277351533, P31, 1341646577593559, 14307862907610435557, 6439376875503349484048791).
- **Composite items:** C405, C152, C78 and C160, each with ≥ 2 distinct primes, giving 8.

These 11 items are pairwise coprime and coprime to the six expanded primes.
- Total distinct primes ≡ 1 (mod 19) dividing N: 6 + 7 + 8 = 21.
- At most one of them is p, so |S| ≥ 20 > 18. **Contradiction.**

### Case p = P24 = 109912203092239643840221

- (p+1)/2 = 3262949 · 16842464147039939. Both are prime, non-special and divide N.
- Φ(19) = p contributes nothing new.

| Q | Φ(Q) split |
|---|---|
| 3262949 | **365257** · 2013013 · P106 |
| 365257 | 19 · **457** · **99371** · C92 |
| 457 | 19 · 47564981 · P39 |
| 99371 | 19 · **229** · C87 |
| 229 | 19 · **1103** · 23599 · 54493 · P30 |
| 1103 | 19 · 9463 · C50 |

- Expanded members of S: 365257, 457, 99371, 229 and 1103, giving 5.
- Prime items: 2013013, P106, 47564981, P39, 23599, 54493, P30 and 9463, giving 8.
- Composite items: C92, C87 and C50, giving 6.
- All are pairwise coprime and coprime to p, so |S| ≥ 19 > 18. **Contradiction.**

### Case p = 682621057

- (p+1)/2 = 7² · 67 · 103963.
- New S-primes come from expanding 7, 19, 67, 419, 22573, 103963 and 251941. Their Φ values are split exactly in the certificate:
  - Φ(7) = 419 · 4534166740403
  - Φ(67) = P33
  - Φ(419) = 19 · 22573 · 251941 · 64267121 · C29
  - and so on.
- Expanded members of S: 419, 22573 and 251941, giving 3.
- Prime items: 9, namely 4534166740403, P24, P33, 64267121, 5375177, 2502719, 34570387, 4409 and 37963.
- Composite items: C29, C71, C77 and C88, giving 8.
- Total 20 > 18. **Contradiction.**

### Case p = 1124041

- (p+1)/2 = 562021 (prime).
- Chain: 562021 → 1559 → 40357 and 122209; 40357 → 229 → 1103 → 9463.
- Expanded members of S: 562021, 1559, 40357, 229 and 1103, giving 5.
- Prime items: P24, 122209, 23599, 54493, P30 and 9463, giving 6.
- Composite items: C100, C47, C80 and C50, giving 8.
- Total 19 > 18. **Contradiction.**

### Case p = 7069

- (p+1)/2 = 5 · 7 · 101.
- Φ(5) = 191 · 6271 · 3981071.
- Φ(101) = 7069 · C33, and 7069 = p is excluded.
- The chain continues through 191, 419, 6271 and 22573.
- Expanded members of S: 191, 419, 6271 and 22573, giving 4.
- Prime items: 3981071, 4534166740403, P24, 1277351533, P31, 251941, 64267121 and 5375177, giving 8.
- Composite items: C33, C29, C68 and C71, giving 8.
- Total 20 > 18. **Contradiction.**

### Case p = 827360473

- (p+1)/2 = 547 · 756271.
- Φ(547) = P50.
- Φ(756271) = 457 · 96749 · 753617 · C93.
- The chain continues through 457, 96749, 419 and 22573.
- Expanded members of S: 457, 96749, 419 and 22573, giving 4.
- Prime items: P24, P50, 753617, 47564981, P39, 251941, 64267121 and 5375177, giving 8.
- Composite items: C93, C86, C29 and C71, giving 8.
- Total 20 > 18. **Contradiction.**

Every case leads to a contradiction, so **no odd perfect number p^α·M^18 exists.** ∎

## 3. Notes on what this uses

- **Classical steps:** L1–L6, all standard. L6 is exactly the McDaniel–Hagis (1975) lemma, and Cohen–Williams (1985) use the same scheme for β = 6, 8, 11, 14 and 18.
- **New computational facts:**
  - Φ_19(19) is the prime P24.
  - Φ_19(P24) contains 682621057.
  - The full splits in the tables above.
  - Cohen–Williams (1985) searched only factors < 5·10^6 and report β = 9 "without further success".
- **Spoof sanity check:** the argument uses that p is prime (a prime factor of N that is ≡ 1 mod 19 is either p or in S). A Descartes-type spoof with a composite "p" breaks exactly there. No step would also rule out the spoof structure without using primality.
- **What a specialist should re-check:**
  - the primality proofs of the expanded bases (all ≤ 24 digits here; PARI/GP isprime);
  - the compositeness of the C-items (witness: any failed strong-probable-prime test);
  - the fact that the C-items are not perfect powers.
