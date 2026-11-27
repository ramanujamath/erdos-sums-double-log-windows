"""Check C: the large-prime part (Lemmas 4.1 and 5.1).

* a_R = (1/R!) int_0^inf E_1(u)^R du: a_1 = 1, a_2 = log 2, agreement with the
  (R-1)-fold integral (4.7), and |a_R - e^{-gamma}| 2^R bounded (estimate (4.3));
* the generating function H(v) = int_0^inf (e^{v E_1(u)} - 1) du = sum_R a_R v^R;
* the prime-sum approximation |P_z(1+u/L) - E_1(u)| <= eps_z (E_1(u) + 2e^{-u}/L);
* (log z) F_R(z) close to a_R for R = 1, 2, 3 (estimate (4.2)), at z = 1000;
* the downward and upward transition probabilities of the large-prime chain:
  rows of D sum to 1, the adjoint relation, the marginal recursion, and the
  row-sum bound sum_p U(m, mp) <= 1 + 3 eps_z + 1/z.
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from eslib import *

EG = math.exp(-EULER_GAMMA)

print("C1. The integrals a_R")
aR = {R: a_R(R) for R in range(1, 31)}
check(abs(aR[1] - 1) < 1e-12, f"a_1 = int_0^inf E_1 = 1 (computed {aR[1]:.15f})")
check(abs(aR[2] - math.log(2)) < 1e-12, f"a_2 = log 2 (computed {aR[2]:.15f})")


def I3():
    # int_0^1 int_0^1 dx dy / (1 + x + x y) ; inner integral in y exactly
    n = 20000
    tot = 0.0
    for i in range(n):
        x = (i + 0.5) / n
        tot += math.log((1 + 2 * x) / (1 + x)) / x
    return tot / n


def I4(n=500):
    # int dx dy dw / (1 + x + x y + x y w); inner integral in w exactly
    tot = 0.0
    for i in range(n):
        x = (i + 0.5) / n
        for j in range(n):
            y = (j + 0.5) / n
            c = 1 + x + x * y
            tot += math.log((c + x * y) / c) / (x * y)
    return tot / (n * n)


check(abs(aR[3] - I3()) < 1e-7, f"a_3 = int_[0,1]^2 dx dy/(1+x+xy) = {I3():.10f}")
check(abs(aR[4] - I4()) < 1e-5, f"a_4 = int_[0,1]^3 1/(1+x1+x1x2+x1x2x3) = {I4():.7f}")
ratios = [(aR[R] - EG) * 2 ** R for R in range(1, 26)]
check(all(0 < r < 1 for r in ratios), f"0 < (a_R - e^-gamma) 2^R < 1 for R <= 25 (limit about {ratios[-1]:.4f})")
check(all(aR[R] > aR[R + 1] for R in range(1, 25)), "a_R decreases to e^{-gamma}")


def Hv(v, smin=-80, smax=6.5, n=6000):
    hstep = (smax - smin) / n
    tot = 0.0
    for i in range(n + 1):
        s = smin + i * hstep
        u = math.exp(s)
        e = E1(u)
        val = math.expm1(v * e) * u
        tot += (0.5 if i in (0, n) else 1.0) * val
    return tot * hstep


for v in [0.3, -0.4, 0.5]:
    series = sum(aR[R] * v ** R for R in range(1, 31))
    check(abs(Hv(v) - series) < 1e-9, f"H({v}) = int (e^(v E_1) - 1) du = sum a_R v^R = {series:.10f}")

print("C2. Primes beyond z")
N = 300000
primes = sieve(N)
# theta(t)/t - 1 on [z, N]
theta = 0.0
eps_by_z = {}
zs = [1000, 3000]
sup_after = {z: 0.0 for z in zs}
pi_idx = 0
for i, p in enumerate(primes):
    # just before p: theta(t) for t in [prev p, p) is constant
    prev = primes[i - 1] if i else 2
    for z in zs:
        if p > z:
            lo = max(prev, z)
            sup_after[z] = max(sup_after[z], abs(theta / lo - 1), abs(theta / p - 1))
    theta += math.log(p)
    for z in zs:
        if p >= z:
            sup_after[z] = max(sup_after[z], abs(theta / p - 1))
for z in zs:
    check(sup_after[z] < 0.1, f"z={z}: eps_z ~ sup_(z<=t<=3e5) |theta(t)/t - 1| = {sup_after[z]:.4f}")


def Pz(z, s):
    """sum_{p > z} p^{-s}, with the primes beyond N replaced by int_N^inf t^{-s}/log t dt."""
    tot = 0.0
    for p in primes:
        if p > z:
            tot += math.exp(-s * math.log(p))
    if s > 1:
        tot += E1((s - 1) * math.log(N))
    return tot


z = 1000
L = math.log(z)
eps = sup_after[z]
ok = True
for u in [0.05, 0.2, 0.5, 1.0, 2.0, 4.0]:
    P = Pz(z, 1 + u / L)
    bound = eps * (E1(u) + 2 * math.exp(-u) / L)
    ok &= abs(P - E1(u)) <= bound
check(ok, "z=1000: |P_z(1+u/L) - E_1(u)| <= eps_z (E_1(u) + 2e^{-u}/L) at u = 0.05,...,4")

# (log z) F_R(z) via F_R = (1/L) int_0^inf Z_R(z, 1+u/L) du, Z from Newton's identities
logp = [math.log(p) for p in primes if p > z]
smin, smax, n = -30.0, 4.0, 220
hstep = (smax - smin) / n
acc = {1: 0.0, 2: 0.0, 3: 0.0}
for i in range(n + 1):
    s = smin + i * hstep
    u = math.exp(s)
    sig = 1 + u / L
    P1 = sum(math.exp(-sig * lp) for lp in logp) + E1(u * math.log(N) / L)
    P2 = sum(math.exp(-2 * sig * lp) for lp in logp)
    P3 = sum(math.exp(-3 * sig * lp) for lp in logp)
    Z1 = P1
    Z2 = (P1 * P1 + P2) / 2
    Z3 = (P1 ** 3 + 3 * P1 * P2 + 2 * P3) / 6
    w = (0.5 if i in (0, n) else 1.0) * u * hstep
    acc[1] += w * Z1
    acc[2] += w * Z2
    acc[3] += w * Z3
for R in (1, 2, 3):
    # tail of the u-integral below e^{-30}: int_0^{e^-30} Z_R du ~ (30+...)^R e^{-30}/R!, negligible
    diff = abs(acc[R] - aR[R])
    check(diff <= 2 ** R * (eps + 1 / z), f"z=1000, R={R}: (log z) F_R(z) = {acc[R]:.5f}, a_R = {aR[R]:.5f}, "
          f"difference {diff:.5f} <= 2^R (eps_z + 1/z) = {2**R*(eps+1/z):.4f}")

print("C3. Transition probabilities of the large-prime chain")


def nu(n):
    return 1.0 / (n * math.log(n))


def Ddown(n):
    f = factor(n)
    ln = math.log(n)
    return {n // p: e * math.log(p) / ln for p, e in f.items()}


def Uup(m, p):
    f = factor(m)
    v = f.get(p, 0)
    return math.log(m) * (v + 1) * math.log(p) / (p * math.log(m * p) ** 2)


ok_rows = ok_adj = ok_rec = True
zz = 3
rough = [n for n in range(5, 60000) if all(p > zz for p in factor(n))]
for n in rough:
    d = Ddown(n)
    ok_rows &= abs(sum(d.values()) - 1) < 1e-12
    if Omega(n) >= 2:
        tot = 0.0
        for p in factor(n):
            m = n // p
            ok_adj &= abs(nu(n) * d[m] - nu(m) * Uup(m, p)) < 1e-15 * nu(m) * 10
            tot += nu(m) * Uup(m, p)
        ok_rec &= abs(tot - nu(n)) < 1e-12 * nu(n)
check(ok_rows, f"sum_p D(n, n/p) = 1 for all {len(rough)} 3-rough n < 60000")
check(ok_adj, "nu(m) U(m, mp) = nu(mp) D(mp, m) (adjoint relation)")
check(ok_rec, "sum_{p | n} nu(n/p) U(n/p, n) = nu(n): the marginal recursion behind (5.1)")

z = 1000
eps = sup_after[z]
ok = True
worst = 0.0
for m in [1009, 1013 * 1019, 1009 ** 2, 1009 * 1013 * 1019, 1031 ** 3, 104729]:
    Lm = math.log(m)
    f = factor(m)
    tot = 0.0
    for p in primes:
        if p > z:
            tot += Lm * (f.get(p, 0) + 1) * math.log(p) / (p * (Lm + math.log(p)) ** 2)
    tot += Lm / (Lm + math.log(N))   # primes beyond N, by the prime number theorem
    worst = max(worst, tot)
    ok &= tot <= 1 + 3 * eps + 1 / z
check(ok, f"z=1000: sum_(p>z) U(m, mp) <= 1 + 3 eps_z + 1/z for sample m (largest row sum {worst:.4f})")

finish("check_C_rough")
