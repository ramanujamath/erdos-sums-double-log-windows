"""Check B: the chain on small primes and the rank weights p_k(z).

* exact normalized flows (couplings supported on d -> dp) between consecutive ranks of
  z-smooth integers for small z, found by an exact max-flow computation, and the
  resulting Markov chain with the prescribed marginals (Lemma 3.1);
* log-concavity of W_k(z) (exact rationals);
* the mean/variance of K_z, the Berry-Esseen scale, max_k p_k(z) << 1/sigma, and
  F*_ell(z) = G(ell/sigma) + O(1/sigma) (Lemma 3.2) for z up to 10^6;
* the block of ell largest weights is consecutive and centred at mu_z + O(sigma^{1/2});
* the rearrangement inequality (8.3) for random index sets S.
"""

import math
import os
import random
import sys
from fractions import Fraction
from collections import deque

sys.path.insert(0, os.path.dirname(__file__))
from eslib import *

random.seed(20261127)

print("B1. Normalized flows between consecutive ranks (exact), and the chain")


def smooth_rank(primes, k):
    """all products of k primes (with repetition) from `primes`, as integers."""
    out = set()

    def rec(i, k, cur):
        if k == 0:
            out.add(cur)
            return
        if i == len(primes):
            return
        p = primes[i]
        for e in range(k + 1):
            rec(i + 1, k - e, cur * p ** e)
    rec(0, k, 1)
    return sorted(out)


def max_flow_exact(lower, upper, primes, wl, wu):
    """max flow source -> lower -> upper -> sink with exact Fraction capacities;
    lower->upper edges (d -> dp) have infinite capacity."""
    nodes = ["s", "t"] + [("L", d) for d in lower] + [("U", d) for d in upper]
    cap = {}
    adj = {v: [] for v in nodes}

    def add(u, v, c):
        if (u, v) not in cap:
            adj[u].append(v)
            adj[v].append(u)
            cap[(u, v)] = Fraction(0)
            cap.setdefault((v, u), Fraction(0))
        cap[(u, v)] += c
    INF = Fraction(10 ** 9)
    upset = set(upper)
    for d in lower:
        add("s", ("L", d), wl[d])
        for p in primes:
            if d * p in upset:
                add(("L", d), ("U", d * p), INF)
    for e in upper:
        add(("U", e), "t", wu[e])
    flow = Fraction(0)
    while True:
        prev = {"s": None}
        q = deque(["s"])
        while q and "t" not in prev:
            u = q.popleft()
            for v in adj[u]:
                if v not in prev and cap[(u, v)] > 0:
                    prev[v] = u
                    q.append(v)
        if "t" not in prev:
            break
        # bottleneck
        v = "t"
        b = INF
        while prev[v] is not None:
            u = prev[v]
            b = min(b, cap[(u, v)])
            v = u
        v = "t"
        while prev[v] is not None:
            u = prev[v]
            cap[(u, v)] -= b
            cap[(v, u)] += b
            v = u
        flow += b
    # extract the flow on middle edges: f(u,v) = cap(v,u) (reverse residual) for INF edges
    f = {}
    for d in lower:
        for p in primes:
            if d * p in upset:
                val = cap[(("U", d * p), ("L", d))]
                if val:
                    f[(d, d * p)] = val
    return flow, f


for primes in ([2, 3, 5], [2, 3, 5, 7]):
    kmax = 6 if len(primes) == 4 else 8
    W = {}
    ranks = {}
    for k in range(kmax + 2):
        ranks[k] = smooth_rank(primes, k)
        W[k] = sum(Fraction(1, d) for d in ranks[k])
    ok_flow = True
    dist = {1: Fraction(1)}  # law of D_0
    ok_marg = True
    for k in range(kmax + 1):
        wl = {d: Fraction(1, d) / W[k] for d in ranks[k]}
        wu = {e: Fraction(1, e) / W[k + 1] for e in ranks[k + 1]}
        val, f = max_flow_exact(ranks[k], ranks[k + 1], primes, wl, wu)
        ok_flow &= (val == 1)
        # transitions from the coupling: P(d -> dp) = f(d,dp)/wl[d]
        new = {}
        for (d, e), fv in f.items():
            new[e] = new.get(e, Fraction(0)) + dist.get(d, Fraction(0)) * fv / wl[d]
        ok_marg &= all(new.get(e, 0) == wu[e] for e in ranks[k + 1])
        dist = new
    check(ok_flow, f"primes {primes}: a normalized flow exists between ranks k, k+1 for k <= {kmax}")
    check(ok_marg, f"primes {primes}: the resulting chain has P(D_k = d) = 1/(d W_k) for k <= {kmax+1}")
    check(all(W[k] ** 2 >= W[k - 1] * W[k + 1] for k in range(1, kmax + 1)),
          f"primes {primes}: W_k is log-concave (exact)")

print("B2. Rank weights p_k(z) = W_k(z)/W(z): the law of a sum of geometric variables")


def rank_weights(z, kmax=120):
    ps = sieve(z)
    p = [0.0] * (kmax + 1)
    p[0] = 1.0
    mu = var = 0.0
    for q in ps:
        r = 1.0 / q
        new = [0.0] * (kmax + 1)
        new[0] = (1 - r) * p[0]
        for k in range(1, kmax + 1):
            new[k] = (1 - r) * p[k] + r * new[k - 1]
        p = new
        mu += 1 / (q - 1)
        var += q / (q - 1) ** 2
    return p, mu, math.sqrt(var), ps


diffs = []
for z in [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6]:
    p, mu, sigma, ps = rank_weights(z)
    tot = sum(p)
    mean = sum(k * x for k, x in enumerate(p))
    var = sum((k - mean) ** 2 * x for k, x in enumerate(p))
    check(abs(tot - 1) < 1e-12 and abs(mean - mu) < 1e-9 and abs(var - sigma ** 2) < 1e-8,
          f"z={z:.0e}: sum p_k = 1, mean = mu_z = {mu:.4f}, variance = sigma_z^2 = {sigma**2:.4f}")
    lls = math.log(math.log(z))
    diffs.append((mu - lls, sigma ** 2 - lls))
    check(abs(mu - lls) < 3 and abs(sigma ** 2 - lls) < 3,
          f"z={z:.0e}: mu_z - log log z = {mu - lls:.4f}, sigma_z^2 - log log z = {sigma**2 - lls:.4f} (bounded)")
    check(all(p[k] ** 2 >= p[k - 1] * p[k + 1] * (1 - 1e-12) for k in range(1, 60)), f"z={z:.0e}: p_k log-concave")
    cdf = 0.0
    be = 0.0
    for k in range(60):
        cdf += p[k]
        be = max(be, abs(cdf - Phi((k - mu) / sigma)), abs(cdf - p[k] - Phi((k - mu) / sigma)))
    check(be * sigma < 1.0, f"z={z:.0e}: sup |P(K_z <= t) - Phi((t-mu)/sigma)| = {be:.3f} = {be*sigma:.2f}/sigma")
    check(max(p) * sigma < 1.0, f"z={z:.0e}: max_k p_k(z) = {max(p):.3f} = {max(p)*sigma:.2f}/sigma")
    srt = sorted(p, reverse=True)
    worst = 0.0
    for ell in range(0, int(3 * sigma) + 1):
        Fs = sum(srt[:ell])
        worst = max(worst, abs(Fs - G(ell / sigma)))
        # block of ell largest weights is consecutive
        if ell >= 1:
            idx = sorted(sorted(range(len(p)), key=lambda k: -p[k])[:ell])
            if idx != list(range(idx[0], idx[0] + ell)):
                check(False, f"z={z:.0e}: block of {ell} largest weights is not consecutive")
    check(worst * sigma < 1.5, f"z={z:.0e}: |F*_ell - G(ell/sigma)| <= {worst*sigma:.2f}/sigma for ell <= 3 sigma")

check(max(d[0] for d in diffs) - min(d[0] for d in diffs) < 0.01 and
      max(d[1] for d in diffs) - min(d[1] for d in diffs) < 0.01,
      "the O(1) terms in mu_z and sigma_z^2 have stabilised (spread < 0.01 for z = 10^3..10^6)")

print("B3. Centre of the maximizing block and the rearrangement inequality (8.3)")
p, mu, sigma, _ = rank_weights(10 ** 6)
order = sorted(range(len(p)), key=lambda k: -p[k])
srt = [p[k] for k in order]
Fstar = lambda j: sum(srt[:max(0, j)])
for ell in [1, 2, 3, 4, 5]:
    block = sorted(order[:ell])
    K = block[0]
    check(abs(K + (ell - 1) / 2 - mu) <= 1.0 + math.sqrt(sigma),
          f"ell={ell}: block {{{K},...,{K+ell-1}}}, K + (ell-1)/2 - mu = {K + (ell-1)/2 - mu:+.3f}")
    I = set(block)
    ok = True
    for _ in range(3000):
        size = random.randint(0, ell)
        S = random.sample(range(0, 25), size)
        m = ell - len(I & set(S))
        lhs = Fstar(ell) - sum(p[k] for k in S)
        rhs = 2 * Fstar(ell) - Fstar(ell - m) - Fstar(ell + m)
        ok &= lhs >= rhs - 1e-15 and len(set(S) - I) <= m
    check(ok, f"ell={ell}: F*_ell - sum_S p_k >= 2F*_ell - F*_(ell-m) - F*_(ell+m) for 3000 random S")

finish("check_B_smooth")
