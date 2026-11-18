# No odd perfect number p^α M^(2β) for β ≤ 100 (candidate result)

**Claim (computer-assisted, awaiting expert review).** There is no odd perfect number of the form
N = p^α · M^(2β), with M squarefree and gcd(p, M) = 1, for any β ≤ 100.

Here every non-special prime has the same exponent 2β. Hagis and McDaniel conjectured that no such N exists, for
any β.

**Status: candidate result, not yet peer reviewed.**
- Of these 100 values of β, 45 were not found excluded in the literature, starting with **β = 9**. Cohen and
  Williams (1985) attempted β = 9 "without further success", and recent papers still list it as open.
- The 45 values are: 9, 15, 20, 21, 23, 26, 29, 30, 33, 35, 36, 38, 39, 41, 44, 45, 48, 50, 51, 53, 54, 56, 59,
  60, 63, 65, 66, 68, 69, 71, 74, 75, 78, 80, 81, 83, 84, 86, 89, 90, 93, 95, 96, 98, 99.
- The other 55 values reproduce known results. They serve as a consistency check.
- Expert feedback is very welcome; please open an issue.

**This is not a proof that odd perfect numbers do not exist.** It settles cases of one classical thin family.

**Disclosure.** I am not a professional number theorist. The certificates and both checkers were produced with
substantial help from AI agents (Anthropic's Claude). The safeguards are:
- two checkers written independently of each other;
- primality proofs done with PARI/GP;
- tamper tests;
- a literature review.

— Dave Post, November 2026

---

## Start here

- **[PROOF_beta9.md](PROOF_beta9.md)** is a human-readable proof for β = 9, about 3 pages. Every number in it can
  be checked against `certs/cert_n19.json`.
- `certs/cert_n<2β+1>.json` holds one machine-readable certificate per β, for n = 2β+1 = 3 … 201.
- `certs/class_r<r>.json` holds the certificates for Theorem B below.
- `logs/summary_final.txt` is a table of all 100 cases, with the prior literature result where one exists.

## Method

The method is the classical one of Kanold, McDaniel–Hagis and Cohen–Williams. Write n = 2β+1 and take a prime r | n
with e = v_r(n).

1. **r | N** (Kanold 1941), and r is not the special prime.
2. **Budget** (McDaniel–Hagis 1975): at most ⌊(n−1)/e⌋ non-special primes Q ≡ 1 (mod r) divide N.
3. **Forcing chain.**
   - For a non-special prime Q, Φ_r(Q) | N, and its prime factors are r or ≡ 1 (mod r).
   - Also (p+1)/2 | N.
   - Starting from r, the chain forces more primes ≡ 1 (mod r) into N than the budget allows.
4. **Scenarios.** The special prime p could be one of the primes that the chain expands, so each such candidate gets
   its own scenario.
5. **Counting composites.** A composite cofactor that is not a perfect power counts as ≥ 2 distinct primes.

For β = 9 the new ingredient is computational:
- Φ_19(19) = 109912203092239643840221 is prime.
- Φ_19 of that number (415 digits) has the factor 682621057. This is far beyond the trial-division range used in 1985.
- Every scenario then forces 19 or 20 members, against a budget of 18.

**Barrier check.** `spoof_check.py` replays the n = 3 certificate against Descartes' spoof
3²·7²·11²·13²·22021. The argument fails exactly at the step that uses primality of the bases, as it must.
Run it from the repository root.

## Theorem B (class version)

Let r ≡ 3 (mod 4) be prime, and suppose r | 2β_i + 1 for every non-special exponent of an odd perfect number N. Then:

| r | 7 | 11 | 19 | 23 | 31 | 43 | 47 |
|---|---|---|---|---|---|---|---|
| v_r(N) ≥ | 412 | 318 | 208 | 206 | 154 | 214 | 234 |

## Verification

There are two independent checkers. The second was written without reading the first.

**Checker 1: `verify.py`** (pure Python integers).
- Primality is proven by deterministic Miller–Rabin below 3.3·10^24 (Sorenson–Webster), and by recursive
  Pocklington above.
- It re-checks every divisibility, every count and the coverage of scenarios.
- Result: 107 of 107 certificates VERIFIED, in `logs/verify_final_a.out` and `logs/verify_final_b.out`.
- `tamper_test.py`: 16 tampered certificates, all rejected.

**Checker 2: `referee/rcheck.py`** (gmpy2). It was written separately, without reading `verify.py`. "Referee" here means a second program, not a human referee.
- Result: 107 of 107 accepted, in `referee/runA.log`, `referee/runB.log` and `referee/rcheck_*.json`.
- `referee/tamper.py`: 12 tampered certificates plus a dropped part, all rejected.

**Primality and compositeness.**
- `referee/rprime.py` proves all 6 714 primality obligations: BPSW below 2^64, PARI/GP `isprime` above.
  Result in `referee/rprime_result.json`.
- `referee/gp_w2b.txt` checks, in PARI/GP, that all 5 976 items counted as composite fail a base-2 Fermat test
  and are not perfect powers. Result: `W2 checked 5976 bad 0` (`referee/gp_w2b.out`).
- A base-3 Fermat test is unsuitable here. Items coming from Φ_r(3) are base-3 Fermat pseudoprimes by construction.

### Re-run the checks

You need Python 3 with `gmpy2` and `sympy`. PARI/GP is optional.

```sh
python3 verify.py certs/cert_n19.json                 # β = 9, under a second
python3 verify.py certs/*.json                        # all 107
cd referee && python3 rcheck.py ../certs              # independent checker, all certificates
python3 rprime.py prime_obligations_all.txt           # primality proofs (PARI/GP)
gp -q gp_w2b.txt                                      # compositeness / not-a-power (PARI/GP)
```

### Re-generate the certificates

This is optional and needs GMP-ECM (`ecm`) plus a C compiler for `cycfac.c`.

```sh
cc -O2 -o cycfac cycfac.c
python3 prove.py 19            # writes certs/cert_n19.json
```

## Literature checked

- Steuerwald 1937; Kanold 1941; McDaniel 1970; Hagis–McDaniel 1972; McDaniel–Hagis 1975.
- Cohen–Williams 1985 (β = 9, 15, 20, 21, 23 attempted without success).
- Fletcher–Nielsen–Ochem 2012 (β ≡ 2 mod 5).
- Brown 2015 (arXiv:1512.01270).
- Yamada 2019–2024 (arXiv:2103.16936, arXiv:2312.17059: "β ≥ 9" still open).
- Gallardo–Rahavandrainy 2023.
- Pearlman 2005, which excludes β = 38 only under 3 | N.
- Ochem–Zelinsky 2026.

No published exclusion was found for any of the 45 values above. Hagis–McDaniel 1972 and the full text of FNO 2012
could not be read directly.

## Caveats

- Neither a human referee nor a proof assistant has checked this yet.
- β = 89 uses a 33-digit factor taken from factordb. It is verified inside the certificate.

## License

Code: MIT (see `LICENSE`). Data and certificates: CC BY 4.0.
