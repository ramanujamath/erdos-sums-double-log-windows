"""Shared helpers for the verification scripts (Python standard library only)."""

import math
import sys

EULER_GAMMA = 0.57721566490153286060651209008240243

_results = {"pass": 0, "fail": 0}


def check(cond, msg):
    if cond:
        _results["pass"] += 1
        print("  [ok]  " + msg)
    else:
        _results["fail"] += 1
        print("  [FAIL] " + msg)


def finish(name):
    p, f = _results["pass"], _results["fail"]
    print(f"{name}: {p} passed, {f} failed")
    sys.exit(0 if f == 0 else 1)


# ----------------------------------------------------------------------
# Gaussian functions
# ----------------------------------------------------------------------

def phi(u):
    return math.exp(-u * u / 2) / math.sqrt(2 * math.pi)


def Phi(v):
    return 0.5 * (1 + math.erf(v / math.sqrt(2)))


def G(lam):
    """G(lambda) = 2 Phi(lambda/2) - 1."""
    return 2 * Phi(lam / 2) - 1


# ----------------------------------------------------------------------
# primes
# ----------------------------------------------------------------------

def sieve(N):
    s = bytearray([1]) * (N + 1)
    s[0] = s[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i:: i] = bytearray(len(range(i * i, N + 1, i)))
    return [i for i in range(N + 1) if s[i]]


def factor(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def Omega(n):
    return sum(factor(n).values())


# ----------------------------------------------------------------------
# exponential integral E_1
# ----------------------------------------------------------------------

def E1(u):
    """E_1(u) = int_u^inf e^{-v}/v dv for u > 0 (double precision)."""
    if u <= 0:
        raise ValueError
    if u < 1.5:
        s = 0.0
        term = 1.0
        k = 1
        while True:
            term *= -u / k
            add = term / k
            s += add
            if abs(add) < 1e-18 * max(1.0, abs(s)):
                break
            k += 1
        return -EULER_GAMMA - math.log(u) - s
    if u > 700:
        return 0.0
    # modified Lentz continued fraction (Numerical Recipes expint, n = 1)
    b = u + 1.0
    c = 1e300
    d = 1.0 / b
    h = d
    i = 1
    while True:
        an = -i * i
        b += 2.0
        d = 1.0 / (an * d + b)
        c = b + an / c
        de = c * d
        h *= de
        if abs(de - 1.0) < 1e-16:
            break
        i += 1
    return h * math.exp(-u)


def a_R(R, smin=None, smax=6.5, n=None):
    """a_R = (1/R!) int_0^inf E_1(u)^R du, via u = e^s and the trapezoid rule."""
    if smin is None:
        smin = -(3 * R + 60)
    if n is None:
        n = 40 * int(smax - smin) + 2000
    hstep = (smax - smin) / n
    tot = 0.0
    logfact = math.lgamma(R + 1)
    for i in range(n + 1):
        s = smin + i * hstep
        u = math.exp(s)
        e = E1(u)
        if e <= 0:
            continue
        val = math.exp(R * math.log(e) + s - logfact)
        w = 0.5 if i in (0, n) else 1.0
        tot += w * val
    return tot * hstep
