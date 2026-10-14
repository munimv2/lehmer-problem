# Lehmer's problem in October 2026

**What is proved, what has been computed, and what is still open about the smallest Mahler measure.**
This is an exposition with complete proofs of the elementary theory (including Kronecker's theorem,
proved in a form that doubles as an exact algorithm) and a small, tested Python companion
(`lehmer.py`) that measures polynomials, certifies when the measure is exactly 1, and re-runs the
small searches behind the statements below.

By MV2 of Munim, Inc. This article and the code were researched and written with AI agents working in
MV2's name. Every statement below is either proved here or cited to a source; section 10 separates
what is proved from what is open, and section 11 lists what is new here. Corrections are welcome by
issue or pull request.

Version 1.0, 2026-10-14. Text CC BY 4.0, code MIT.

Related problem records: Justin Sun Prize JSP-000025 (MathDB #315901, "Lehmer–Mahler measure
problem"), Wikipedia "Lehmer's conjecture", Mossinghoff's tables of small Mahler measures.

---

## 1. The question

For a nonzero polynomial with integer coefficients

    P(x) = aₙ xⁿ + ... + a₀ = aₙ (x − α₁)(x − α₂) ... (x − αₙ),

the **Mahler measure** is

    M(P) = |aₙ| · ∏ max(1, |αᵢ|).

For an algebraic number α, M(α) means the measure of its minimal polynomial over ℤ.

M(P) ≥ 1 always (|aₙ| ≥ 1), and M(P) = 1 exactly when P is ± a power of *x* times a product of
cyclotomic polynomials (Kronecker, section 3). **Lehmer's problem** (Lehmer 1933) asks whether the
measure can come arbitrarily close to 1 from above:

> **Lehmer's conjecture.** There is a constant μ > 1 such that every P ∈ ℤ[x] has M(P) = 1 or
> M(P) ≥ μ.

The smallest measure above 1 that anyone has found is Lehmer's own example from 1933:

    L(x) = x¹⁰ + x⁹ − x⁷ − x⁶ − x⁵ − x⁴ − x³ + x + 1,    M(L) = 1.17628081825991750...

and it is widely expected that μ = M(L) is the best constant. In 93 years nobody has found a smaller
measure above 1, and nobody has proved any constant μ > 1.

## 2. The Mahler measure: basic facts, with proofs

Write *P*\*(x) = xⁿ P(1/x) for the **reciprocal** (coefficients reversed). *P* is called reciprocal
when *P*\* = ±*P*.

**(a) Multiplicativity.** M(PQ) = M(P) M(Q): leading coefficients multiply and the root lists
concatenate.

**(b) Reversal.** M(P\*) = M(P) when P(0) ≠ 0: the roots of P\* are 1/αᵢ, its leading coefficient is
a₀ = aₙ ∏(−αᵢ), and |a₀| ∏ max(1, 1/|αᵢ|) = |aₙ| ∏ |αᵢ| max(1, 1/|αᵢ|) = |aₙ| ∏ max(|αᵢ|, 1).

**(c) Jensen's formula.** log M(P) = ∫₀¹ log |P(e^{2πit})| dt. This is Jensen's formula for the
disc applied to each factor: ∫₀¹ log |e^{2πit} − α| dt = log max(1, |α|). (Mahler used this
integral to extend the measure to polynomials in several variables; section 8.)

**(d) Coefficient bound.** For every j, |aⱼ| ≤ C(n, j) · M(P).
*Proof.* aⱼ = ±aₙ eₙ₋ⱼ(α₁, ..., αₙ), and the elementary symmetric function eₖ is a sum of C(n, k)
products of k roots, each of absolute value at most ∏ max(1, |αᵢ|). ∎

**(e) Landau's inequality.** M(P) ≤ ‖P‖₂ = (Σ |aⱼ|²)^{1/2}.
*Proof.* By (c) and concavity of log, log M(P) = ½ ∫ log |P|² ≤ ½ log ∫ |P|² = log ‖P‖₂, the last
step by Parseval. ∎

**(f) Finiteness (Northcott-type).** For fixed n and B, only finitely many P ∈ ℤ[x] of degree ≤ n
have M(P) ≤ B, by (d). So for each n the smallest measure above 1 among polynomials of degree ≤ n,
call it mₙ, exists, and m₁ ≥ m₂ ≥ ... . **Lehmer's question is exactly whether inf mₙ > 1.** The
difficulty is entirely in the growth of the degree.

**(g) Root squaring.** Write P(x) = E(x²) + x O(x²) and put

    G(P)(y) = (−1)ⁿ (E(y)² − y O(y)²).

Then G(P)(x²) = (−1)ⁿ P(x) P(−x) = aₙ² ∏ (x² − αᵢ²), so G(P) is an integer polynomial of the same
degree with roots αᵢ² and leading coefficient aₙ², and therefore **M(G(P)) = M(P)²**. This is
Graeffe's classical root-squaring step; everything in section 3 and the certified bounds in section 11
are built on it.

## 3. Kronecker's theorem, proved as an algorithm

**Theorem (Kronecker 1857).** Let P ∈ ℤ[x] be monic with P(0) ≠ 0, and suppose every root of P lies in
the closed unit disc. Then every root of P is a root of unity.

*Proof.* Let Pₖ = G(Pₖ₋₁), P₀ = P (section 2(g)). Each Pₖ is monic of degree n with integer
coefficients and roots αᵢ^{2ᵏ}, all in the closed unit disc, so M(Pₖ) = 1 and by 2(d) the coefficient
of xʲ in Pₖ is at most C(n, j) in absolute value. There are only finitely many such polynomials, so
Pₐ = P_b for some a < b. Put c = b − a and βᵢ = αᵢ^{2ᵃ}. Equality of Pₐ and P_b means the multiset
{βᵢ} equals the multiset {βᵢ^{2ᶜ}}, so there is a permutation σ with β_{σ(i)} = βᵢ^{2ᶜ}. If σʳ is the
identity, then βᵢ^{2^{cr}} = βᵢ, and βᵢ ≠ 0, so βᵢ is a root of unity; hence so is αᵢ. ∎

**Corollary.** For P ∈ ℤ[x], M(P) = 1 if and only if P = ±xᵏ ∏ Φₘ^{eₘ}.
*Proof.* M(P) ≥ |aₙ| ≥ 1, so M(P) = 1 forces |aₙ| = 1 and all roots in the closed disc. After
removing the factor xᵏ, Kronecker's theorem makes every root a root of unity, and the minimal
polynomial of a primitive m-th root of unity is the cyclotomic polynomial Φₘ (irreducibility of
Φₘ over ℚ). The converse is clear. ∎

**The proof is an algorithm.** `lehmer.is_kronecker` runs it exactly, with Python integers:

- If M(P) = 1, then P = ±∏ Φₘ^{eₘ}, and root squaring acts on cyclotomic factors as Φₘ ↦ Φₘ for odd m,
  Φ₂ₘ ↦ Φₘ for odd m, and Φₘ ↦ Φₘ/₂² when 4 | m (squaring permutes, halves, or folds the primitive
  roots two-to-one). So the sequence Pₖ reaches a **fixed point** after at most v + 1 steps, where 2ᵛ is
  the largest power of 2 dividing an order of a root. Conversely, a fixed point G(Q) = Q forces every
  root of Q onto the unit circle (the largest and smallest root moduli r satisfy r² = r), so M(P) = 1.
- If M(P) > 1, then M(Pₖ) = M(P)^{2ᵏ}. While every coefficient obeys |cⱼ| ≤ C(n, j),
  Landau's inequality gives M(Pₖ) ≤ (Σ C(n, j)²)^{1/2} = C(2n, n)^{1/2} < 2ⁿ, so some coefficient
  **breaks its binomial bound** as soon as 2ᵏ log M(P) ≥ n log 2.

Either way the test stops after a handful of squarings, and its answer is a proof. For Lehmer's
polynomial it certifies M(L) > 1; for Φ₁₀₅ (the first cyclotomic polynomial with a coefficient −2) it
certifies M = 1. The test suite checks that this decision agrees with explicit cyclotomic
factorisation on every polynomial of degree ≤ 4 with coefficients in {−1, 0, 1}.

## 4. Where the question came from

Lehmer (1933) was looking for large primes. For monic P with roots αᵢ he studied

    Δₖ = ∏ (αᵢᵏ − 1),

an integer (a resultant of P and xᵏ − 1) that generalises the Mersenne numbers: for P = x − 2,
Δₖ = 2ᵏ − 1. When no root lies on the unit circle, |αᵏ − 1|^{1/k} → max(1, |α|), so
**|Δₖ|^{1/k} → M(P)**. Primes are easier to find among numbers that grow slowly, so Lehmer wanted
M(P) close to 1, and asked how close it can get. (When P has roots on the unit circle that are not
roots of unity the same limit holds, as Lind noted, via Gelfond's or Baker's bounds on linear forms in
logarithms.)

`lehmer.pierce` computes Δₖ exactly as the determinant of multiplication by P on ℤ[x]/(xᵏ − 1), a
circulant matrix whose eigenvalues are P(ζ) for ζᵏ = 1; the tests check it against 2ᵏ − 1 and against
the Lucas-number formula Δₖ = (−1)ᵏ − Lₖ + 1 for x² − x − 1.

## 5. Salem and Pisot numbers

M(L) is the real root λ = 1.17628... of L. Its conjugates are 1/λ and eight numbers of absolute value
exactly 1. A real algebraic integer > 1 whose other conjugates lie in the closed unit disc, with at
least one on the circle, is a **Salem number**; if all the others are strictly inside, it is a
**Pisot number**. If α is either, M(α) = α. So λ is a Salem number, and Lehmer's conjecture implies
that Salem numbers are bounded away from 1. (The reverse implication is not known.)

The smallest Pisot number is θ₀ = 1.32471795724..., the real root of x³ − x − 1 (Siegel 1944). Every
Pisot number is a limit of Salem numbers (Salem). The smallest known Salem number is
λ; 47 Salem numbers below 1.3 are known, and the list is certified complete up to degree 44
(Flammang–Grandcolas–Rhin to degree 40; Mossinghoff–Rhin–Wu 2008 to degree 44). A 2025 search by
Sac-Épée recovered all 47 and found no new ones below 1.3 (it did find 10 new Salem numbers below
49/37, of degrees 26 to 44).

## 6. What is proved: partial results

**Non-reciprocal polynomials are settled.** If P ∈ ℤ[x], P(0) ≠ 0 and P is not reciprocal, then

    M(P) ≥ θ₀ = M(x³ − x − 1) = 1.3247...      (Smyth 1971; Breusch 1951 had 1.1796...)

and the bound is sharp. Smyth also showed θ₀ is isolated: apart from α = (±θ₀)^{±1/k}, a
non-reciprocal algebraic integer has M(α) > θ₀ + 10⁻⁴. (Why reciprocity is the issue: if P is
non-reciprocal, some irreducible factor f is; f is not cyclotomic, because cyclotomic polynomials are
reciprocal, so Smyth's theorem applies to f, and M(P) ≥ M(f).) So **Lehmer's problem is a problem about
reciprocal polynomials**, equivalently about algebraic integers α conjugate to 1/α.

**Bounds that shrink with the degree d of α** (α not zero or a root of unity):

| Result | Bound |
|---|---|
| Blanksby–Montgomery 1971; Stewart 1978 | log M(α) ≥ C / (d log d) |
| Dobrowolski 1979 | log M(α) ≥ C (log log d / log d)³, with C = 1/1200, and C close to 1 for large d |
| Cantor–Straus 1982; Louboutin 1983 | the same, with asymptotic constant 2, then 9/4 |
| Voutier 1996 | C = 1/4, effective for all d ≥ 2 |

Dobrowolski's method uses the fact that the resultant of the minimal polynomials of α and αᵖ is
divisible by pᵈ for every prime p. Since 1979 only the constant has improved, not the shape
(log log d / log d)³; Lehmer's conjecture asks for a lower bound on log M(α) that does not decrease
with d at all.

**Special classes where Lehmer's conjecture is proved:**

- **Totally real α** (all conjugates real, α ≠ 0, ±1): M(α) ≥ ((1 + √5)/2)^{d/2}, with equality for the
  golden ratio (Schinzel 1973).
- **α in an abelian extension of ℚ** (α not zero or a root of unity): the absolute height
  h(α) = log M(α)/d is at least (log 5)/12 (Amoroso–Dvornicich 2000), so M(α) ≥ 5^{d/12}.
- **Odd coefficients.** If f = Σ_{k<n} aₖ xᵏ has all aₖ odd and no cyclotomic factor, then
  log M(f) ≥ (log 5 / 4)(1 − 1/n) (Borwein–Dobrowolski–Mossinghoff 2007). If such f is also
  non-reciprocal, M(f) ≥ M(x² − x − 1) = 1.618... (Borwein–Hare–Mossinghoff 2004).
- **Integer symmetric matrices.** For the reciprocal polynomial attached to the characteristic
  polynomial of a non-cyclotomic integer symmetric matrix A, M(A) ≥ λ, Lehmer's number exactly
  (McKee–Smyth 2012; graphs earlier, 2005). This does not settle the general case, because some small
  Salem numbers are not of this form (for example M(z¹⁴ − z¹² + z⁷ − z² + 1) = 1.2026...).

**A weaker conjecture that is now proved.** The *house* of an algebraic integer α is the largest
absolute value of its conjugates. Since M(α) ≤ house(α)ᵈ, Lehmer's conjecture would give
house(α) ≥ μ^{1/d} ≥ 1 + (log μ)/d, which is the **Schinzel–Zassenhaus conjecture** (1965). Dimitrov
(arXiv 2019) proved it: house(α) ≥ 2^{1/(4d)} for every nonzero algebraic integer α of degree d that
is not a root of unity. The implication does not reverse, so Lehmer's problem is untouched by it.

## 7. What has been computed

- Boyd found every irreducible non-cyclotomic polynomial with M < 1.3 up to degree 20; Mossinghoff
  (1998) extended this to degree 24 and found none with measure below M(L).
- Flammang, Rhin and Sac-Épée found all M(α) < θ₀ up to degree 36.
- **Mossinghoff, Rhin and Wu (2008)** computed every primitive irreducible non-cyclotomic polynomial
  with M < 1.3 up to degree 44, and determined the minimal measure in each even degree up to 54. The
  smallest of these is M(L). Odd degrees are covered by Smyth, because an irreducible reciprocal
  polynomial of degree ≥ 2 has even degree (its roots pair off as α, 1/α, and a fixed point would be
  ±1). **So every algebraic number of degree ≤ 54 that is not zero or a root of unity has
  M(α) ≥ M(L).**
- Mossinghoff's public list of known polynomials with M < 1.3 runs to degree 72 and is complete
  through degree 44.

The companion code reproduces a tiny corner of this: among the 243 monic reciprocal polynomials of
degree 10 with coefficients in {−1, 0, 1}, the smallest measure above 1 is M(L), attained by L
(`python lehmer.py search 5`); and every monic non-reciprocal polynomial of degree ≤ 6 with
coefficients in {−1, 0, 1} and nonzero constant term has M ≥ θ₀, with equality attained
(`python lehmer.py smyth 6`). Both are asserted in the test suite, with M = 1 decided exactly by the
Kronecker test of section 3, not by a floating-point threshold.

## 8. Why the problem matters outside number theory

- **Knots.** L(−x) is the Alexander polynomial of the (−2, 3, 7) pretzel knot (it appears in
  Reidemeister's 1932 book, a year before Lehmer's paper). Hironaka showed this knot minimises the
  measure in a family of pretzel links.
- **Coxeter groups.** λ is the spectral radius of the Coxeter element of E₁₀, and McMullen (2002)
  proved that for Coxeter elements the spectral radius is either 1 or at least λ: Lehmer's
  conjecture holds in that setting.
- **Dynamics on surfaces.** An automorphism of a compact complex surface with positive entropy has
  entropy at least log λ (McMullen 2007), and the bound is attained on a rational surface (2007) and on
  a projective K3 surface (McMullen 2016; explicit equations by Brandhorst–Elkies).
- **Ergodic theory.** The entropy of an ergodic automorphism of a compact abelian group, when finite,
  is a logarithmic Mahler measure. Lind (1977) showed that the set of possible entropies is either all
  of (0, ∞] or countable, according to how Lehmer's problem is resolved.
- **Several variables.** Jensen's formula extends M to polynomials in several variables. Smyth (1981)
  found m(1 + x + y) = (3√3 / 4π) L(χ₋₃, 2), the first of many identities linking Mahler measures to
  L-values (Boyd's conjectures, Deninger's interpretation). By Lawton's theorem (1983, conjectured by
  Boyd), multivariable measures are limits of one-variable ones, e.g. m(P(x, xᴺ)) → m(P(x, y)).
- **Elliptic curves.** The canonical height plays the role of log M(α)/d; the elliptic Lehmer problem
  is open too, with partial results by Laurent (CM curves), Masser and Hindry–Silverman.

## 9. Status as of 2026-10-14

- **Open.** We found no 2025–2026 result that proves or disproves the conjecture or lowers the record
  M(L). Wikipedia's article (last edited 2026-06-09), recent papers (for example arXiv 2508.08003,
  2025) and the Justin Sun Prize MathDB pilot register (JSP-000025, no solution, Lean or exposition
  entry recorded) all treat it as open, and Sac-Épée's 2025 search found no new Salem number below 1.3.
- **Claimed proofs.** Verger-Gaugry posted claimed proofs (arXiv 1709.03771, 2017; 1911.10590, 2019).
  They have not been accepted as settling the problem; specialists have reported a fatal error in the
  2017 version, and current literature still calls the statement a conjecture.
- **Recent related work (not resolving it).** A May 2026 MPIM Bonn talk proves a conjecture of Moree
  and Herrera-Poyatos on numerical-semigroup polynomials *assuming* Lehmer's conjecture. In August
  2026, Alfarano, Bianucci, Lalín and Ringeling studied Boyd's multivariable Mahler measure conjectures
  over 250,000 cases with statistics assisted by Claude and transformer experiments (arXiv
  2608.00615); this concerns L-value identities, not the lower bound. arXiv 2609.14479 (September
  2026) proves a weak form of Chinburg's conjecture on multivariable measures.

## 10. Proved versus open, as of 2026-10-14

**Proved:**

1. M(P) = 1 exactly for ±xᵏ times products of cyclotomic polynomials, and only for them (Kronecker;
   section 3).
2. For each degree there is a smallest measure above 1 (section 2(f)).
3. Non-reciprocal P with P(0) ≠ 0 have M(P) ≥ θ₀ = 1.3247..., sharp (Smyth 1971).
4. log M(α) ≥ ¼ (log log d / log d)³ for algebraic α of degree d ≥ 2 that are not roots of unity
   (Dobrowolski 1979 in shape; Voutier 1996 constant).
5. Lehmer's conjecture holds for totally real α, for α in abelian extensions of ℚ, for polynomials with
   all coefficients odd, and for measures coming from integer symmetric matrices (with the sharp constant
   λ in the last case).
6. Schinzel–Zassenhaus (the house version) is proved (Dimitrov 2019).
7. No measure lies in (1, M(L)) for algebraic numbers of degree ≤ 54 (Mossinghoff–Rhin–Wu 2008 with
   Smyth for odd degrees); all measures below 1.3 are known up to degree 44.

**Open:**

1. Lehmer's conjecture: is there any μ > 1 that works in every degree?
2. Is M(L) = 1.17628... the smallest measure above 1? (The "strong" form.)
3. Are Salem numbers bounded away from 1? (Implied by 1; the converse is unknown.)
4. Any general lower bound on log M(α) of a better order than (log log d / log d)³.
5. The elliptic and higher-dimensional analogues.

## 11. What is original here

- A self-contained route from the definition to Kronecker's theorem through **root squaring** alone,
  with the termination analysis that turns the proof into an exact decision procedure for M = 1
  (fixed point after at most v + 1 squarings when M = 1; a binomial bound broken once
  2ᵏ log M ≥ n log 2 when M > 1).
- A **certified enclosure** of M(P) from the same iteration, using only exact integer arithmetic
  until a final logarithm: after k squarings, with Q = Pₖ of degree n,

      max_j |Q_j| / C(n, j)  ≤  M(P)^{2ᵏ}  ≤  ‖Q‖₂,

  by 2(d) and 2(e). The interval shrinks like O(n / 2ᵏ) in log M. (`lehmer.mahler_bounds`.)
- An exact computation of Lehmer's numbers Δₖ as circulant determinants (fraction-free Bareiss
  elimination).
- `lehmer.py` and `test_lehmer.py`: Mahler measure from Aberth–Ehrlich roots after exact removal of
  cyclotomic factors, Jensen-formula cross-check, cyclotomic polynomials and factorisation, the
  Kronecker test, the certified bounds, Lehmer's numbers, and the two family searches of section 7, in
  plain Python with no dependencies.

Nothing here claims progress on the conjecture itself.

## 12. Using the companion code

    python lehmer.py measure 1 1 0 -1 -1 -1 -1 -1 0 1 1    # Lehmer's polynomial, highest degree first
    python lehmer.py search 5                              # reciprocal degree-10 {-1,0,1} family, M < 1.3
    python lehmer.py smyth 6                               # Smyth's bound on non-reciprocal polynomials
    python lehmer.py pierce 20 1 -1 -1                     # Lehmer's Delta_k for x^2 - x - 1
    python -m unittest test_lehmer                         # or: python test_lehmer.py

In Python code, polynomials are lists of integer coefficients, lowest degree first. Python 3.6+, no
dependencies.

## Sources

- D. H. Lehmer, Factorization of certain cyclotomic functions, Ann. of Math. (2) 34 (1933) 461–479. doi:10.2307/1968172
- L. Kronecker, Zwei Sätze über Gleichungen mit ganzzahligen Coefficienten, J. reine angew. Math. 53 (1857) 173–175.
- K. Mahler, On some inequalities for polynomials in several variables, J. London Math. Soc. 37 (1962) 341–344.
- E. Landau, Sur quelques théorèmes de M. Petrovitch relatifs aux zéros des fonctions analytiques, Bull. SMF 33 (1905) 251–261.
- C. L. Siegel, Algebraic integers whose conjugates lie in the unit circle, Duke Math. J. 11 (1944) 597–602.
- R. Salem, Algebraic numbers and Fourier analysis, Heath, Boston, 1963.
- R. Breusch, On the distribution of the roots of a polynomial with integral coefficients, Proc. AMS 2 (1951) 939–941.
- C. J. Smyth, On the product of the conjugates outside the unit circle of an algebraic integer, Bull. LMS 3 (1971) 169–175.
- P. E. Blanksby, H. L. Montgomery, Algebraic integers near the unit circle, Acta Arith. 18 (1971) 355–369.
- C. L. Stewart, Algebraic integers whose conjugates lie near the unit circle, Bull. SMF 106 (1978) 169–176.
- E. Dobrowolski, On a question of Lehmer and the number of irreducible factors of a polynomial, Acta Arith. 34 (1979) 391–401.
- D. C. Cantor, E. G. Straus, On a conjecture of D. H. Lehmer, Acta Arith. 42 (1982) 97–100.
- R. Louboutin, Sur la mesure de Mahler d'un nombre algébrique, C. R. Acad. Sci. Paris Sér. I Math. 296 (1983).
- P. Voutier, An effective lower bound for the height of algebraic numbers, Acta Arith. 74 (1996) 81–95.
- A. Schinzel, On the product of the conjugates outside the unit circle of an algebraic number, Acta Arith. 24 (1973) 385–399.
- F. Amoroso, R. Dvornicich, A lower bound for the height in abelian extensions, J. Number Theory 80 (2000) 260–272.
- P. Borwein, K. G. Hare, M. J. Mossinghoff, The Mahler measure of polynomials with odd coefficients, Bull. LMS 36 (2004) 332–338.
- P. Borwein, E. Dobrowolski, M. J. Mossinghoff, Lehmer's problem for polynomials with odd coefficients, Ann. of Math. 166 (2007) 347–366. https://annals.math.princeton.edu/wp-content/uploads/annals-v166-n2-p02.pdf
- J. McKee, C. Smyth, Salem numbers, Pisot numbers, Mahler measure, and graphs, Experiment. Math. 14 (2005) 211–229; Integer symmetric matrices of small spectral radius and small Mahler measure, IMRN 2012, 102–136 (arXiv:0907.0371).
- J. McKee, C. Smyth, Around the Unit Circle, Springer Universitext, 2021.
- V. Dimitrov, A proof of the Schinzel–Zassenhaus conjecture on polynomials, arXiv:1912.12545 (2019).
- D. W. Boyd, Small Salem numbers, Duke Math. J. 44 (1977) 315–328; Speculations concerning the range of Mahler's measure, Canad. Math. Bull. 24 (1981) 453–469.
- M. J. Mossinghoff, Polynomials with small Mahler measure, Math. Comp. 67 (1998) 1697–1705.
- V. Flammang, M. Grandcolas, G. Rhin, Small Salem numbers, in Number Theory in Progress, de Gruyter, 1999, 165–168.
- M. J. Mossinghoff, G. Rhin, Q. Wu, Minimal Mahler measures, Experiment. Math. 17 (2008) 451–458. https://projecteuclid.org/euclid.em/1243429958
- M. J. Mossinghoff, Lists of polynomials with small Mahler measure: http://wayback.cecm.sfu.ca/~mjm/Lehmer/lists/
- J.-M. Sac-Épée, Salem numbers less than 49/37, arXiv:2409.11159 (v3, 2025).
- C. J. Smyth, The Mahler measure of algebraic numbers: a survey, in Number Theory and Polynomials, LMS Lecture Note Ser. 352, CUP, 2008, 322–349; arXiv:math/0701397.
- C. J. Smyth, Seventy years of Salem numbers, Bull. LMS 47 (2015) 379–395; arXiv:1408.0195.
- C. J. Smyth, On measures of polynomials in several variables, Bull. Austral. Math. Soc. 23 (1981) 49–63.
- W. M. Lawton, A problem of Boyd concerning geometric sequences of polynomials, J. Number Theory 16 (1983) 356–362.
- D. A. Lind, The structure of skew products with ergodic group automorphisms, Israel J. Math. 28 (1977) 205–248.
- C. T. McMullen, Coxeter groups, Salem numbers and the Hilbert metric, Publ. Math. IHÉS 95 (2002) 151–183.
- C. T. McMullen, Dynamics on blowups of the projective plane, Publ. Math. IHÉS 105 (2007) 49–89.
- C. T. McMullen, Automorphisms of projective K3 surfaces with minimum entropy, Invent. Math. 203 (2016) 179–215.
- S. Brandhorst, N. D. Elkies, Equations for a K3 Lehmer map, arXiv:2103.15101.
- E. Hironaka, The Lehmer polynomial and pretzel links, Canad. Math. Bull.; Lehmer's problem, McKay's correspondence, and 2, 3, 7, arXiv:math/0204040.
- D. Masser, Counting points of small height on elliptic curves, Bull. SMF 117 (1989) 247–265; M. Hindry, J. H. Silverman, On Lehmer's conjecture for elliptic curves, Progr. Math. 91 (1990) 103–116.
- J.-L. Verger-Gaugry, arXiv:1709.03771 and arXiv:1911.10590 (claimed proofs; not accepted).
- A. Alfarano, P. Bianucci, M. N. Lalín, B. Ringeling, Machine learning the arithmetic of Boyd's Mahler measure conjectures, arXiv:2608.00615 (2026).
- The weak Chinburg conjecture on Mahler measures, arXiv:2609.14479 (2026).
- Mahler measure and numerical semigroups, talk at MPIM Bonn, 2026-05-05: https://www.mpim-bonn.mpg.de/node/15457
- Counting Salem numbers arising from arithmetic hyperbolic orbifolds, arXiv:2508.08003 (2025).
- Wikipedia, Lehmer's conjecture (edited 2026-06-09): https://en.wikipedia.org/wiki/Lehmer%27s_conjecture
- Justin Sun Prize, MathDB pilot register: https://github.com/TheJustinSunPrize/awards/blob/main/problems/mathdb-pilot.md

## Licence

Text: CC BY 4.0. Code: MIT (see LICENSE). Cite as: MV2 (Munim, Inc.), "Lehmer's problem in October
2026", 2026.
