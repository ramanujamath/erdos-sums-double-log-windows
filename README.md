# Gaussian concentration of Erdős sums in double-logarithmic windows

**Sumit Kumar Jha**, Greater Noida, Uttar Pradesh, India (sumitkumarjha75@gmail.com)

This repository holds the manuscript (LaTeX source and PDF), an audit of the earlier version, and pure-Python scripts that check the elementary and numerical ingredients of the proofs.

## The paper in brief

A set of integers > 1 is *primitive* if no element divides another; its *Erdős sum* is Σ 1/(n log n). Let 𝓜(T, h) be the largest Erdős sum of a primitive set inside the window T < log log n ≤ T + h.

* **Concentration (Theorem 1.1).** For every fixed λ₀ > 0, 𝓜(T, λ√T) = 2Φ(λ/2) − 1 + O(T^{−1/10}) uniformly for 0 ≤ λ ≤ λ₀. The integers in the window with exactly ⌊T + λ√T/2⌋ prime factors attain the limit.
* **Stability (Theorem 1.2).** Near-maximizers put almost all their mass on integers whose number of prime factors is within η√T of T + λ√T/2. The leftover is ≪ (deficit + T^{−1/10})^{1/2}.
* **Sharpness (Proposition 8.3).** The level sets of Ω(n) + v₂(n) are asymptotically maximizing, yet no single value of Ω carries more than half of their mass.
* **Corollaries.** 𝓜(T, h) − G(h/√T) → 0 uniformly in h, and 𝓜(T, h) ~ h/√(2πT) when h → ∞ with h = o(√T).

The upper bound runs two independent random divisibility chains, one through the small and one through the large prime factors (after Harper and after Alexeev–Barreto–Li–Lichtman–Price–Shah–Tang–Tao). It then bounds the resulting antichain by a normal approximation.

## Repository layout

| Path | Contents |
|---|---|
| `paper/main.tex`, `paper/main.pdf` | The rewritten manuscript (16 pages) |
| `AUDIT.md` | Audit of the earlier version: clarified steps, citation fixes, notation changes, independent confirmations |
| `original/` | The earlier version that was audited |
| `verification/` | Verification scripts (Python standard library only) |

## Running the checks

```bash
python3 verification/run_all.py
```

This runs 123 checks in about 5 seconds and needs no third-party packages.

| Script | What it checks |
|---|---|
| `check_A_gaussian.py` | derivatives and concavity of G, the curvature used for the block centre, the Poisson local limit step of the lower bound, and the exponent bookkeeping behind T^{−1/10} |
| `check_B_smooth.py` | exact normalized flows and the small-prime chain for small prime sets, log-concavity, and the normal approximation of the rank weights for z up to 10⁶ |
| `check_C_rough.py` | the integrals a_R and their 2^{−R} convergence to e^{−γ}, the prime-sum approximation, (log z)F_R(z) ≈ a_R, and the large-prime chain's transition identities and row sums |
| `check_D_combinatorics.py` | the antichain lemma (exhaustively on small grids), the divisor identity and moment and Chernoff bounds for the large-prime rank, and primitivity in the sharpness example |

The checks are numerical sanity tests of the ingredients. The proofs in the paper are complete without them.

## Building the PDF

```bash
cd paper && latexmk -pdf main.tex
```
