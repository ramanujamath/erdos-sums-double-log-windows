"""Check D: combinatorial and elementary steps.

* antichains in {0,1,...} x {a,...,b}: distinct coordinates, |H| <= ell, and Lemma 8.2
  (|S \\ I| <= m and |k + j - (K + b)| <= m), exhaustively on small grids;
* the antichain property of the random set H (product-chain bound) on a toy example;
* r^{J_z(n)} = sum_{d | n} g_r(d) (Lemma 6.1, upper tail);
* E(V* - B)^4 <= B + 3B^2 for independent Bernoulli sums (exact rationals);
* the comparison between (1/X) sum_{n<=X} (V(n)-B)^4 and E(V*-B)^4;
* the Chernoff algebra in Lemma 6.1;
* Proposition 8.3: f = Omega + v_2 is completely additive and its level sets are primitive;
  numerical illustration of the halving by powers of 2.
"""

import itertools
import math
import os
import random
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(__file__))
from eslib import *

random.seed(7)

print("D1. Antichains in the grid (Lemma 8.2)")
total = 0
ok = True
for ell in range(1, 5):
    a, b = 3, 3 + ell - 1
    cells = [(k, j) for k in range(0, 7) for j in range(a, b + 1)]
    for r in range(0, ell + 2):
        for H in itertools.combinations(cells, r):
            anti = all(not (p != q and p[0] <= q[0] and p[1] <= q[1]) for p in H for q in H)
            if not anti:
                continue
            total += 1
            ks = [p[0] for p in H]
            js = [p[1] for p in H]
            ok &= len(set(ks)) == len(ks) and len(set(js)) == len(js) and len(H) <= ell
            S = set(ks)
            for K in range(0, 7):
                I = set(range(K, K + ell))
                m = ell - len(S & I)
                ok &= len(S - I) <= m
                for (k, j) in H:
                    if k in I:
                        ok &= abs(k + j - (K + b)) <= m
check(ok, f"all {total} antichains (ell <= 4, k <= 6): distinct coordinates, |H| <= ell, "
          f"|S \\ I| <= m and |k + j - (K+b)| <= m")

print("D2. Product chains give antichains (toy example)")
# A primitive set A of integers; for chains D_0 | D_1 | ... (smooth) and M_a | ... | M_b (rough),
# the set {(k,j) : D_k M_j in A} is an antichain.
z = 5
A = [n for n in range(2, 4000) if Omega(n) == 5]  # primitive
Aset = set(A)
ok = True
for trial in range(300):
    D = [1]
    for _ in range(8):
        D.append(D[-1] * random.choice([2, 3, 5]))
    M = [random.choice([7, 11, 13])]
    for _ in range(3):
        M.append(M[-1] * random.choice([7, 11, 13, 17]))
    H = [(k, j) for k in range(len(D)) for j in range(len(M)) if D[k] * M[j] in Aset]
    ok &= all(not (p != q and p[0] <= q[0] and p[1] <= q[1]) for p in H for q in H)
check(ok, "300 random chain pairs: {(k, j) : D_k M_j in A} is an antichain for primitive A")

print("D3. The identity r^{J_z(n)} = sum_{d | n} g_r(d)")
for z in [3, 7]:
    for r in [Fraction(3, 2), Fraction(2)]:
        ok = True
        for n in range(1, 3000):
            f = factor(n) if n > 1 else {}
            J = sum(e for p, e in f.items() if p > z)
            s = Fraction(0)
            divisors = [1]
            for p, e in f.items():
                divisors = [d * p ** i for d in divisors for i in range(e + 1)]
            for d in divisors:
                g = Fraction(1)
                for p, e in (factor(d).items() if d > 1 else []):
                    g *= 0 if p <= z else (r - 1) * r ** (e - 1)
                s += g
            ok &= s == r ** J
        check(ok, f"z={z}, r={r}: r^J_z(n) = sum_(d|n) g_r(d) for n < 3000")
# Euler product bound: sum_{d} g_r(d)/d over primes in (z, X] equals prod (1 + (r-1)/(p-r))
r = 1.5
ps = [p for p in sieve(200) if p > 7]
prod = 1.0
for p in ps:
    prod *= 1 + (r - 1) / (p - r)
series = 1.0
for p in ps:
    series *= 1 + sum((r - 1) * r ** (e - 1) / p ** e for e in range(1, 60))
check(abs(prod / series - 1) < 1e-12, "sum_e g_r(p^e)/p^e = (r-1)/(p-r) for p > r")

print("D4. Fourth moments")
ok = True
for trial in range(200):
    ps = [Fraction(1, random.randint(3, 200)) for _ in range(random.randint(1, 12))]
    B = sum(ps)
    # exact fourth central moment of an independent Bernoulli sum via cumulants
    k2 = sum(p * (1 - p) for p in ps)
    k4 = sum(p * (1 - p) * (1 - 6 * p * (1 - p)) for p in ps)
    m4 = k4 + 3 * k2 ** 2
    ok &= m4 <= B + 3 * B ** 2
check(ok, "E(V* - B)^4 <= B + 3B^2 for 200 random independent Bernoulli families (exact)")
# brute-force distribution check for one family
ps = [Fraction(1, 5), Fraction(1, 7), Fraction(1, 11)]
B = sum(ps)
m4 = Fraction(0)
for bits in itertools.product([0, 1], repeat=3):
    pr = Fraction(1)
    for bit, p in zip(bits, ps):
        pr *= p if bit else 1 - p
    m4 += pr * (sum(bits) - B) ** 4
k2 = sum(p * (1 - p) for p in ps)
k4 = sum(p * (1 - p) * (1 - 6 * p * (1 - p)) for p in ps)
check(m4 == k4 + 3 * k2 ** 2, "cumulant formula matches brute force")
# integers versus the independent model
X = 200000
z, Y = 5, 40
P = [p for p in sieve(Y) if p > z]
B = sum(1 / p for p in P)
emp = sum((sum(1 for p in P if n % p == 0) - B) ** 4 for n in range(1, X + 1)) / X
pf = [Fraction(1, p) for p in P]
Bf = sum(pf)
k2 = sum(p * (1 - p) for p in pf)
k4 = sum(p * (1 - p) * (1 - 6 * p * (1 - p)) for p in pf)
model = float(k4 + 3 * k2 ** 2)
check(abs(emp - model) <= (1 + B + Y) ** 4 / X, f"(1/X) sum (V(n)-B)^4 = {emp:.5f}, model {model:.5f}, "
      f"difference within (1+B+Y)^4/X = {(1+B+Y)**4/X:.3f}")

print("D5. Chernoff step")
ok = True
for t in [i / 1000 * math.log(2) for i in range(1001)]:
    ok &= math.expm1(t) <= t + t * t + 1e-15
check(ok, "e^theta - 1 <= theta + theta^2 for 0 <= theta <= log 2")
ok = True
for T in [1e4, 1e8, 1e12]:
    for lam in [0.0, 0.5, 2.0, 6.0]:
        for C0 in [1.0, 5.0]:
            D = T ** 0.4
            h = lam * math.sqrt(T)
            Q = h + 2 * D + C0
            th = D / (4 * Q)
            b = math.ceil(h + 3 * D)
            if D / 2 < C0:
                continue
            ok &= th <= math.log(2) and b >= Q + D / 2
            ok &= -th * b + math.expm1(th) * Q <= -D * D / (16 * Q) + 1e-9
check(ok, "-theta b + (e^theta - 1) Q <= -D^2/(16 Q) with theta = D/(4Q), b >= Q + D/2")

print("D6. Proposition 8.3: the level sets of Omega(n) + v_2(n)")


def v2(n):
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c


Nmax = 20000
f = [0] * (Nmax + 1)
for n in range(2, Nmax + 1):
    f[n] = Omega(n) + v2(n)
ok = all(f[u * v] == f[u] + f[v] for u in range(1, 120) for v in range(1, Nmax // 120))
check(ok, "f(uv) = f(u) + f(v) (complete additivity), checked for u < 120")
ok = True
for n in range(2, Nmax + 1):
    for mult in range(2 * n, Nmax + 1, n):
        if f[mult] == f[n]:
            ok = False
            break
    if not ok:
        break
check(ok, f"no n properly divides a multiple with the same value of f (n <= {Nmax})")

finish("check_D_combinatorics")
