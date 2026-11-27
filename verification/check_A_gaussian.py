"""Check A: Gaussian ingredients and exponent bookkeeping.

* G(lambda) = 2 Phi(lambda/2) - 1, its derivatives, and the second-difference bound of Lemma 8.1;
* the curvature of H_x(t) = Phi(t + x/2) - Phi(t - x/2) at t = 0;
* the local limit step in Lemma 2.1 (Poisson weights versus phi), with error O(T^{-1/2});
* the change of variables s = log log t in the partial summation of Lemma 2.1;
* the bounds used for the theorem's error term (ell/sigma, T/D^4, D^2/(16Q), K+b);
* the arithmetic behind the upper bound min{1, lambda/sqrt(2 pi)} and Corollary 7.2.
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from eslib import *

print("A1. The function G")
for lam in [0.1, 0.5, 1, 2, 3, 5]:
    hstep = 1e-4
    d1 = (G(lam + hstep) - G(lam - hstep)) / (2 * hstep)
    d2 = (G(lam + hstep) - 2 * G(lam) + G(lam - hstep)) / hstep ** 2
    check(abs(d1 - phi(lam / 2)) < 1e-7 and abs(d2 + lam * phi(lam / 2) / 4) < 1e-5,
          f"G'({lam}) = phi(lam/2), G''({lam}) = -lam phi(lam/2)/4")
check(abs(G(1e-9)) < 1e-9 and abs(G(40) - 1) < 1e-15, "G(0) = 0 and G(lambda) -> 1")
check(abs(G(2.0) - (2 * Phi(1.0) - 1)) < 1e-15 and abs(G(2.0) - math.erf(1 / math.sqrt(2))) < 1e-15,
      "G(lambda) = erf(lambda/(2 sqrt 2)), the normal mass of an interval of length lambda centred at 0")
# normal mass of an interval of length x is maximal when centred
ok = True
for x in [0.3, 1.0, 2.5]:
    for t in [0.01, 0.1, 0.5, 1.0]:
        ok &= Phi(t + x / 2) - Phi(t - x / 2) < G(x)
check(ok, "H_x(t) = Phi(t+x/2) - Phi(t-x/2) < G(x) for t != 0")
for x in [0.3, 1.0, 2.5]:
    hstep = 1e-4
    H = lambda t: Phi(t + x / 2) - Phi(t - x / 2)
    d2 = (H(hstep) - 2 * H(0) + H(-hstep)) / hstep ** 2
    check(abs(d2 + x * phi(x / 2)) < 1e-5, f"H_x''(0) = -x phi(x/2) at x = {x}")

print("A2. Second-difference bound 2G(x) - G(x-y) - G(x+y) >= c y^2 (0 <= y <= x)")
for lam in [0.25, 1.0, 3.0]:
    cmin = float("inf")
    for i in range(1, 41):
        x = lam * (0.5 + i / 40)          # x in a compact neighbourhood of lam
        for j in range(1, 201):
            y = x * j / 200
            val = (2 * G(x) - G(x - y) - G(x + y)) / (y * y)
            cmin = min(cmin, val)
    check(cmin > 0, f"lambda = {lam}: min of the second difference / y^2 on the compact set = {cmin:.4f} > 0")

print("A3. Local limit step in the lower bound: sqrt(T) e^{-s} s^m/m! vs phi(u - c)")
for T in [1e4, 1e6, 1e8]:
    worst = 0.0
    for c in [-1.5, 0.0, 0.7, 2.0]:
        k = math.floor(T + c * math.sqrt(T))
        m = k - 1
        for u in [-2.0, -0.5, 0.0, 0.3, 1.0, 2.5]:
            s = T + u * math.sqrt(T)
            logpois = -s + m * math.log(s) - math.lgamma(m + 1)
            ratio = math.sqrt(T) * math.exp(logpois) / phi(u - c)
            worst = max(worst, abs(ratio - 1) * math.sqrt(T))
    check(worst < 20, f"T = {T:.0e}: |ratio - 1| <= {worst:.2f} T^(-1/2)")

print("A4. Partial summation: the substitution s = log log t")
# integrand N(t) (log t + 1)/(t^2 (log t)^2) dt with N(t) = t s^m/(m! log t)
# equals e^{-s} s^m/m! (1 + e^{-s}) ds, since ds = dt/(t log t)
ok = True
for s in [2.0, 3.0, 5.5]:
    for m in [3, 7]:
        logt = math.exp(s)
        lhs = (s ** m / (math.factorial(m) * logt)) * (logt + 1) / (logt ** 2) * logt  # per ds, times t
        rhs = math.exp(-s) * s ** m / math.factorial(m) * (1 + math.exp(-s))
        ok &= abs(lhs / rhs - 1) < 1e-12
check(ok, "N(t)(log t + 1)/(t^2 log^2 t) dt = e^{-s} s^m/m! (1 + e^{-s}) ds")

print("A5. Exponent bookkeeping in the proof of the main theorem")
for L in [1.0, 4.0]:
    for T in [1e6, 1e10, 1e16]:
        worst = 0.0
        for lam in [0.0, 0.3, 1.0, L]:
            D = T ** 0.4
            h = lam * math.sqrt(T)
            a = math.floor(D)
            b = math.ceil(h + 3 * D)
            ell = b - a + 1
            sigma = math.sqrt(T - 2 * D)          # sigma_z^2 = T - 2T^{2/5} + O(1)
            worst = max(worst, abs(ell / sigma - lam) * T ** 0.1)
        check(worst < 3 + L, f"L={L}, T={T:.0e}: |ell/sigma - lambda| <= {worst:.2f} T^(-1/10)")
for T in [1e6, 1e10]:
    D = T ** 0.4
    check(abs(T / D ** 4 - T ** -0.6) < 1e-12 * T ** -0.6, f"T={T:.0e}: T/D^4 = T^(-3/5)")
    for L in [0.0, 1.0, 5.0]:
        h = L * math.sqrt(T)
        Q = h + 2 * D + 10
        check(D ** 2 / (16 * Q) >= T ** 0.3 / (16 * (L + 3)),
              f"T={T:.0e}, lambda={L}: D^2/(16Q) >= c_L T^(3/10)")
        check(h * T ** -0.6 <= (L + 1e-12) * T ** -0.1, f"T={T:.0e}, lambda={L}: h T^(-3/5) <= L T^(-1/10)")
# K + b = mu + (a+b)/2 + O(sigma^{1/2}) = T + lambda sqrt T / 2 + O(T^{1/4})
for T in [1e8, 1e12]:
    D = T ** 0.4
    lam = 1.3
    h = lam * math.sqrt(T)
    a, b = math.floor(D), math.ceil(h + 3 * D)
    mu = T - 2 * D
    val = mu + (a + b) / 2
    check(abs(val - T - h / 2) < 5, f"T={T:.0e}: mu + (a+b)/2 = T + lambda sqrt(T)/2 + O(1)")

print("A6. Harmonic bounds")
for T in [1e4, 1e8]:
    for lam in [0.5, 2.0]:
        h = lam * math.sqrt(T)
        val = (1 / math.sqrt(T + h) + 2 * (math.sqrt(T + h) - math.sqrt(T))) / math.sqrt(2 * math.pi)
        check(abs(val - lam / math.sqrt(2 * math.pi)) < 2 / math.sqrt(T) + lam ** 2 / (4 * math.sqrt(T)),
              f"T={T:.0e}, lambda={lam}: (1/sqrt(T+h) + 2(sqrt(T+h)-sqrt T))/sqrt(2pi) -> lambda/sqrt(2pi)")
for T in [1e8, 1e12]:
    h = T ** 0.25
    val = 1 / math.sqrt(T + h) + 2 * (math.sqrt(T + h) - math.sqrt(T))
    check(abs(val / (h / math.sqrt(T)) - 1) < 2 / h, f"T={T:.0e}, h=T^(1/4): bound ~ h/sqrt(T) (Corollary 7.2)")
    b = h / math.sqrt(T)
    check(abs(G(b) / (b * phi(0)) - 1) < b * b, f"T={T:.0e}: G(b) = b/sqrt(2 pi) (1 + O(b^2))")

finish("check_A_gaussian")
