# Audit of the original manuscript

**Original:** *Gaussian concentration of Erdős sums in double-logarithmic windows*, by Sumit Jha, dated 11 November 2026, 16 pages (`original/Erdos_Primitive_Concentration.pdf`).
**Rewrite:** `paper/main.tex` (compiled to `paper/main.pdf`, 16 pages).
**Checks:** `verification/`. There are 123 automated checks and all pass; run `python3 verification/run_all.py`.

Every proof was read line by line, and every elementary inequality, identity, exponent and constant that can be tested numerically was tested. **No mathematical error was found.** Both main theorems, the corollaries and the example are correct as stated, with the proofs given. The issues were heavy reuse of symbols, displays that were garbled in the PDF, a few steps stated without justification, and one citation number. Lemma and equation numbers below refer to the original.

## A. Steps that needed justification or precision

| # | Location | Issue | Resolution |
|---|----------|-------|------------|
| A1 | Lemma 6, proof of (14) | "∫ against dt is E₁(u)" stated without the substitution. | Substitution t = e^{Lw/u} given. |
| A2 | Lemma 6, convergence of H_z, H | "Converge locally uniformly on a neighbourhood of \|v\| ≤ 1/2" justified only by "E₁(u) = −γ − log u + O(u)". | Added that the integrands are O(u^{−1/2−o(1)}) near 0, and that this makes Ψ_z, Ψ analytic, which Cauchy's estimate needs. |
| A3 | Lemma 7, row-sum bound | "Multiplying by L and using L ≥ log z > 1" hides two facts. | Spelled out: L_m/(L_m + log z) < 1 and L_m + log z > 2. Both are needed for the constant 3 in δ_z = 3ε_z + 1/z. |
| A4 | Lemma 9, upper tail | θ ≤ log 2 was asserted. | Noted that Q ≥ 2κ gives θ ≤ 1/8, uniformly down to λ = 0. |
| A5 | Proof of Thm 2 | "The computation in Prop. 8" was reused for subsets B without saying why it applies. | Restated as a displayed identity valid for every subset 𝓑 ⊂ A_g. The bound V ≤ U ≤ F* ≤ 1 used in the indicator term is now stated. |
| A6 | Lemma 5, (7) | The upper bound "any block has mass ≤ G(ℓ/σ) + O(1/σ)" and the centring argument were merged into one sentence. | Separated: the normal mass of an interval of length x is at most G(x), with equality only when centred. |
| A7 | Intro, choice of window length | "Windows of length √T" was not motivated. | One sentence: log log n is the mean and the variance of Ω(n). A heuristic before Lemma 9 explains the interval [a, b]. |

## B. Citations

| # | Location | Issue | Resolution |
|---|----------|-------|------------|
| B1 | Lemma 6, "[9, Theorem 6]" | In Gorodetsky–Lichtman–Wong (arXiv:2303.08277) the statement I_k = e^{−γ} + O(2^{−k}) is **Theorem 1.6**, with I_k = ∫_{[0,1]^k} dx/(1 + x₁(1 + x₂(⋯(1 + x_k)))) as in (16). | Cited as Theorem 1.6; the arXiv number added. The numbering in the published C. R. Math. version should be confirmed before submission. |
| B2 | Reference [5] | Confirmed: arXiv:2605.00301, 1 May 2026; Theorem 1.1 is the bound E(B) ≤ 1 + O(1/log x). | Kept. |
| B3 | Sathe–Selberg, "[4, Theorem 2.3]" | Not checked against the source. | Kept; worth confirming the theorem number. |

## C. Presentation and notation

* **Overloaded symbols.** Many letters had two to four meanings. The rewrite gives each object its own symbol:

  | Object | Original | Rewrite |
  |---|---|---|
  | Erdős sum / expectation | E(A) / E | 𝓔(A) / 𝔼 |
  | bound on λ | L | λ₀ |
  | log z, log m | L, L | L, L_m |
  | window offset T^{2/5} | D | κ |
  | downward/upward transitions | D(n, n/p), U(m, mp) | π↓, π↑ |
  | small-/large-prime chains | D_k, M_j | **s**_k, **r**_j (M is the maximum 𝓜) |
  | antichain | H | 𝓗 |
  | generating integrals | H_z(v), H(v) | Ψ_z, Ψ |
  | Φ(t+x/2) − Φ(t−x/2) | H_x(t) | Φ_x(t) |
  | harmonic sums | H_A(t) | 𝓛_A(t) |
  | Sathe–Selberg factor | F(w) | 𝔖(w) |
  | sum of geometric variables | K_z | 𝒦_z |
  | exponent cap / block start / k_T | K, K, K_T | E, k₀, k_T |
  | missing indices ℓ − \|S ∩ I\| | m (also the rough integer) | 𝔪 |
  | Σ 1/p in Lemma 9 | B | β |
  | variables in (16) | P, B | Π, Σ |
  | subsets in Thm 2 proof | B | 𝓑 |
  | j ≥ 2 remainder | B_z(v, u) | 𝓑_z(v, u) |
  | sums in Thm 2 proof | U, V | 𝓤, 𝓥 |
  | normalization factor | Q_T | 𝓠_T |
  | example set / counting function | C_T, C(t) | 𝓒_T, N^exc(t) |
  | rank interval in Lemma 3 | a, b, c | u₁, u₂, c |
  | maximizing block | I_z | 𝓘 |

* **Definitions added.** Erdős sum, rank, normalized flow, covering relation, ranked poset, coupling, absorbing state, large-prime rank J_z, exponential integral, and Harper's product theorem stated in words.
* **Garbled displays** in the original PDF (abstract, Theorem 1, Lemma 3, the definition of δ_z, Section 7, Prop. 14) are typeset correctly.
* **Structure.** Results, context and method now have separate subsections in the introduction. The method subsection explains the two-chain argument and where the exponent 1/10 comes from (T/κ⁴ = T^{−3/5}, times h ≍ √T).
* **Author information.** Full name, postal address and e-mail; an author-contributions statement; the AI-assistance statement is kept.

## D. Independent confirmations (see `verification/`)

* **Gaussian facts:** G′, G″, the curvature Φ_x″(0) = −xφ(x/2), and positivity of the second-difference constant on compact sets.
* **Lemma 3 (lower bound):** the Poisson local limit step holds with \|ratio − 1\| ≈ 6/√T at T = 10⁴, 10⁶, 10⁸. The s = log log t substitution was verified.
* **Lemma 4 (small-prime chain):** for the primes {2,3,5} and {2,3,5,7}, exact normalized flows between consecutive ranks were found by exact rational max-flow. The resulting Markov chain was checked to have exactly the marginals 1/(d W_k). W_k is log-concave.
* **Lemma 5 (normal approximation):** for z = 10³ to 10⁶, μ_z − log log z → 1.0347 and σ_z² − log log z → 2.4098. Berry–Esseen distance ≈ 0.32/σ, max p_k ≈ 0.43/σ, and \|F*_ℓ − G(ℓ/σ)\| ≤ 0.11/σ. The top-ℓ weights are always consecutive.
* **Lemma 6:**
  * a₁ = 1 and a₂ = log 2.
  * a₃ and a₄ agree with the cube integral (16).
  * **(a_R − e^{−γ})·2^R → 0.31523…**, monotonically, for R ≤ 30. This is an independent numerical confirmation of (12).
  * Ψ(v) = Σ a_R v^R at three values of v.
  * Inequality (14) holds at z = 1000.
  * (log z)F_R(z) is within 0.003 of a_R for R = 1, 2, 3 at z = 1000.
* **Lemma 7 (large-prime chain):** the downward rows sum to 1, the adjoint relation holds, and the marginal recursion Σ ν(n/p)U(n/p, n) = ν(n) holds for all 3-rough n < 60000. Sample row sums at z = 1000 are about 0.75, below 1 + δ_z.
* **Prop. 8 and Lemma 13 (antichains):** checked on all 494 antichains with ℓ ≤ 4; the antichain property also holds for random chain pairs.
* **Lemma 9:**
  * r^{J_z(n)} = Σ_{d\|n} g_r(d) exactly.
  * E(V* − β)⁴ ≤ β + 3β² exactly.
  * The integer fourth moment matches the independent model.
  * The Chernoff algebra holds over a parameter grid.
  * Exponent bookkeeping: ℓ/σ = λ + 2T^{−1/10} + …
* **Prop. 14:** Ω + v₂ is completely additive, and its level sets are primitive up to 20000.
