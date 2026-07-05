import Mathlib

/-!
# NitroSAT proof skeletons in Lean 4 / Mathlib

This file formalizes the parts of the NitroSAT note that can be stated as
ordinary finite-dimensional algebra/calculus-free claims, and it exposes the
analytic-number-theory and spectral-genericity inputs as hypotheses.

The point is deliberate: Lean should not be asked to magically prove BV/BDH,
Riemann-zero estimates, Lambert-W asymptotics, or the Spectral Genericity
assumption unless those are separately imported as theorems. Once those inputs
are given, the downstream solver-stability claims are clean algebra.
-/

noncomputable section

open scoped BigOperators
open Real

namespace NitroSAT

/-! ## Common `O(1)` wrapper

Using a tiny boundedness predicate avoids relying on mathlib's asymptotics API.
For the condition-number theorem, `IsO1 f` means exactly that `f` is bounded by
some absolute constant independent of the problem size.
-/

def IsO1 (f : ℕ → ℝ) : Prop :=
  ∃ C : ℝ, 0 ≤ C ∧ ∀ K : ℕ, |f K| ≤ C

lemma IsO1_const (C : ℝ) : IsO1 (fun _ : ℕ => C) := by
  refine ⟨|C|, abs_nonneg C, ?_⟩
  intro K
  simp

lemma IsO1_mul {f g : ℕ → ℝ} (hf : IsO1 f) (hg : IsO1 g) :
    IsO1 (fun K => f K * g K) := by
  rcases hf with ⟨Cf, hCf_nonneg, hf_bound⟩
  rcases hg with ⟨Cg, hCg_nonneg, hg_bound⟩
  refine ⟨Cf * Cg, mul_nonneg hCf_nonneg hCg_nonneg, ?_⟩
  intro K
  calc
    |f K * g K| = |f K| * |g K| := by rw [abs_mul]
    _ ≤ Cf * Cg := by
      exact mul_le_mul (hf_bound K) (hg_bound K) (abs_nonneg (g K)) hCf_nonneg

/-! ## Prime weights -/

def primeWeight (p : ℝ) : ℝ :=
  (1 : ℝ) / (1 + Real.log p)

lemma primeWeight_pos {p : ℝ} (hp : 1 ≤ p) :
    0 < primeWeight p := by
  unfold primeWeight
  have hlog : 0 ≤ Real.log p := Real.log_nonneg hp
  have hden : 0 < 1 + Real.log p := by linarith
  positivity

lemma primeWeight_le_one {p : ℝ} (hp : 1 ≤ p) :
    primeWeight p ≤ 1 := by
  unfold primeWeight
  have hlog : 0 ≤ Real.log p := Real.log_nonneg hp
  have hden_pos : 0 < 1 + Real.log p := by linarith
  have hden_ge_one : (1 : ℝ) ≤ 1 + Real.log p := by linarith
  rw [div_le_iff₀ hden_pos]
  linarith

/-! ## Theorem 1: quadratic gradient vanishing

The inverse conformal factor in the inverted Poincare disk is

`r^2 * (1-r^2)^2 / 4`.

The Lean theorem below proves the useful quantitative form:
for `0 ≤ r ≤ 1`, this factor is bounded above by `r^2 / 4`. Hence any
coordinate gradient whose norm is bounded by `C` has Riemannian gradient norm
bounded by `(r^2 / 4) C`, i.e. it vanishes quadratically as `r → 0`.
-/

def invMetricFactor (r : ℝ) : ℝ :=
  (r ^ 2 * (1 - r ^ 2) ^ 2) / 4

lemma invMetricFactor_nonneg (r : ℝ) :
    0 ≤ invMetricFactor r := by
  unfold invMetricFactor
  positivity

lemma invMetricFactor_le_quadratic {r : ℝ} (h0 : 0 ≤ r) (h1 : r ≤ 1) :
    invMetricFactor r ≤ r ^ 2 / 4 := by
  unfold invMetricFactor
  have hr2_nonneg : 0 ≤ r ^ 2 := sq_nonneg r
  have h_one_minus_r_nonneg : 0 ≤ 1 - r := sub_nonneg.mpr h1
  have hmul : 0 ≤ r * (1 - r) := mul_nonneg h0 h_one_minus_r_nonneg
  have hr2_le_r : r ^ 2 ≤ r := by nlinarith
  have hr2_le_one : r ^ 2 ≤ 1 := by linarith
  have ha0 : 0 ≤ 1 - r ^ 2 := by linarith
  have ha1 : 1 - r ^ 2 ≤ 1 := by nlinarith [hr2_nonneg]
  have hamul : 0 ≤ (1 - r ^ 2) * (1 - (1 - r ^ 2)) := by
    exact mul_nonneg ha0 (sub_nonneg.mpr ha1)
  have ha2_le_a : (1 - r ^ 2) ^ 2 ≤ 1 - r ^ 2 := by nlinarith
  have ha2_le_one : (1 - r ^ 2) ^ 2 ≤ 1 := by linarith
  have hprod : r ^ 2 * (1 - r ^ 2) ^ 2 ≤ r ^ 2 * 1 := by
    exact mul_le_mul_of_nonneg_left ha2_le_one hr2_nonneg
  nlinarith

theorem quadratic_gradient_vanishing_bound
    {r C gradNorm : ℝ}
    (h0 : 0 ≤ r) (h1 : r ≤ 1)
    (hgrad_nonneg : 0 ≤ gradNorm)
    (hgrad_bound : gradNorm ≤ C)
    (hC : 0 ≤ C) :
    invMetricFactor r * gradNorm ≤ (r ^ 2 / 4) * C := by
  have hf_nonneg : 0 ≤ invMetricFactor r := invMetricFactor_nonneg r
  calc
    invMetricFactor r * gradNorm ≤ invMetricFactor r * C := by
      exact mul_le_mul_of_nonneg_left hgrad_bound hf_nonneg
    _ ≤ (r ^ 2 / 4) * C := by
      exact mul_le_mul_of_nonneg_right (invMetricFactor_le_quadratic h0 h1) hC

/-! ## Interior convexity algebra

The analytic proof that the Hessian has the stated lower bound is not encoded
here; the file proves the downstream algebra: if the stated curvature inequality
holds and diffusion contributes a nonnegative spectral term, then the convergence
rate `mu` is positive.
-/

def mu (β Wmax kmax dClause δ λ λ₂ : ℝ) : ℝ :=
  4 / β - (Wmax * kmax ^ 2 * dClause) / δ ^ 2 + λ * λ₂

theorem mu_pos_of_convexity_bound
    {β Wmax kmax dClause δ λ λ₂ : ℝ}
    (hcurv : (Wmax * kmax ^ 2 * dClause) / δ ^ 2 < 4 / β)
    (hdiff : 0 ≤ λ * λ₂) :
    0 < mu β Wmax kmax dClause δ λ λ₂ := by
  unfold mu
  linarith

lemma exp_decay_factor_le_one {μ t : ℝ} (hμ : 0 ≤ μ) (ht : 0 ≤ t) :
    Real.exp ((-μ) * t) ≤ 1 := by
  have hmul : 0 ≤ μ * t := mul_nonneg hμ ht
  have hle : (-μ) * t ≤ 0 := by nlinarith
  rw [← Real.exp_zero]
  exact Real.exp_le_exp.mpr hle

theorem contraction_factor_not_expansive
    {μ t d0 : ℝ} (hμ : 0 ≤ μ) (ht : 0 ≤ t) (hd0 : 0 ≤ d0) :
    Real.exp ((-μ) * t) * d0 ≤ d0 := by
  exact mul_le_of_le_one_left hd0 (exp_decay_factor_le_one hμ ht)

/-! ## Lyapunov structure

This is the finite-dimensional algebraic core of
`dF/dt = - ∇Fᵀ g⁻¹ ∇F ≤ 0`.
-/

def quadForm {ι : Type} [Fintype ι] (Ginv : ι → ι → ℝ) (v : ι → ℝ) : ℝ :=
  ∑ i : ι, ∑ j : ι, v i * Ginv i j * v j

def PositiveSemidefiniteKernel {ι : Type} [Fintype ι]
    (Ginv : ι → ι → ℝ) : Prop :=
  ∀ v : ι → ℝ, 0 ≤ quadForm Ginv v

theorem lyapunov_decrease
    {ι : Type} [Fintype ι]
    {Ginv : ι → ι → ℝ}
    (hpsd : PositiveSemidefiniteKernel Ginv)
    (gradF : ι → ℝ) :
    -quadForm Ginv gradF ≤ 0 := by
  have hq : 0 ≤ quadForm Ginv gradF := hpsd gradF
  linarith

theorem lyapunov_equality_only_at_stationary
    {ι : Type} [Fintype ι]
    {Ginv : ι → ι → ℝ}
    (hzero : ∀ v : ι → ℝ, quadForm Ginv v = 0 → v = 0)
    {gradF : ι → ℝ}
    (hdF_zero : -quadForm Ginv gradF = 0) :
    gradF = 0 := by
  apply hzero
  linarith

/-! ## O(M) scaling from clause locality

A CNF incidence relation is represented as `Inc v c`, meaning variable `v`
occurs in clause `c`. The proof shows exactly the double-counting identity

`sum_v degree(v) = sum_c width(c)`

and then the bounded-width linear bound.
-/

section Incidence

variable {Var Clause : Type}
variable [Fintype Var] [Fintype Clause]
variable (Inc : Var → Clause → Prop)
variable [DecidableRel Inc]

def degree (v : Var) : ℕ :=
  (Finset.univ.filter (fun c : Clause => Inc v c)).card

def width (c : Clause) : ℕ :=
  (Finset.univ.filter (fun v : Var => Inc v c)).card

theorem sum_degrees_eq_sum_widths :
    (∑ v : Var, degree Inc v) = ∑ c : Clause, width Inc c := by
  classical
  simp [degree, width, Finset.card_filter, Finset.sum_comm]

theorem gradient_work_linear_in_clause_count
    {k : ℕ} (hwidth : ∀ c : Clause, width Inc c ≤ k) :
    (∑ v : Var, degree Inc v) ≤ Fintype.card Clause * k := by
  classical
  rw [sum_degrees_eq_sum_widths (Inc := Inc)]
  calc
    (∑ c : Clause, width Inc c) ≤ ∑ c : Clause, k := by
      exact Finset.sum_le_sum (fun c _ => hwidth c)
    _ = Fintype.card Clause * k := by
      simp [Finset.card_univ, Nat.mul_comm]

end Incidence

/-! ## Heat damping and implicit spectral preconditioning

The full theorem needs spectral genericity / BDH-type estimates as inputs. In
Lean, those estimates are hypotheses; the proof below shows how the downstream
`O(1)` condition-number conclusion follows once the three factors in the paper's
formula are uniformly bounded:

`κ_eff = spectralRatio * varianceRatio * heatRatio`.
-/

def heatRatio (t λmax λ₂ : ℝ) : ℝ :=
  Real.exp ((-t) * λmax) / Real.exp ((-t) * λ₂)

lemma heatRatio_le_one {t λmax λ₂ : ℝ}
    (ht : 0 ≤ t) (hλ : λ₂ ≤ λmax) :
    heatRatio t λmax λ₂ ≤ 1 := by
  unfold heatRatio
  have hden_pos : 0 < Real.exp ((-t) * λ₂) := Real.exp_pos _
  have hmul : t * λ₂ ≤ t * λmax := mul_le_mul_of_nonneg_left hλ ht
  have hexp : Real.exp ((-t) * λmax) ≤ Real.exp ((-t) * λ₂) := by
    apply Real.exp_le_exp.mpr
    nlinarith
  rw [div_le_iff₀ hden_pos]
  simpa using hexp

def kappaEffSeq
    (spectralRatio varianceRatio heatRatioSeq : ℕ → ℝ) : ℕ → ℝ :=
  fun K => spectralRatio K * varianceRatio K * heatRatioSeq K

theorem implicit_spectral_preconditioning_O1
    {spectralRatio varianceRatio heatRatioSeq : ℕ → ℝ}
    (hspec : IsO1 spectralRatio)
    (hvar : IsO1 varianceRatio)
    (hheat : IsO1 heatRatioSeq) :
    IsO1 (kappaEffSeq spectralRatio varianceRatio heatRatioSeq) := by
  unfold kappaEffSeq
  exact IsO1_mul (IsO1_mul hspec hvar) hheat

/-! ## Lambert-W / BAHA interface

Mathlib does not ship a high-level Lambert-W asymptotic theorem. So the correct
formal style is: assume a branch `LambertW` solves the fold equation and assume
injectivity on the branch being used. Then the BAHA jump is the Lambert-W value.
The large-window approximation is also encoded as an explicit asymptotic input.
-/

def foldEquation (ξ Δ : ℝ) : Prop :=
  Δ * Real.exp Δ = ξ

theorem fold_solution_is_lambertW
    (LambertW : ℝ → ℝ)
    {ξ Δ : ℝ}
    (hΔ : foldEquation ξ Δ)
    (hW : foldEquation ξ (LambertW ξ))
    (hbranch_inj : ∀ {a b : ℝ},
      a * Real.exp a = b * Real.exp b → a = b) :
    Δ = LambertW ξ := by
  apply hbranch_inj
  exact hΔ.trans hW.symm

def fisherHorizon (T : ℝ) : ℝ :=
  T / Real.exp 1

theorem baha_large_window_formula
    (LambertW : ℝ → ℝ)
    {T ξ Δ : ℝ}
    (hΔ : Δ = LambertW ξ)
    (hAsymp : LambertW ξ = fisherHorizon T - Real.log (fisherHorizon T)) :
    Δ = T / Real.exp 1 - Real.log (T / Real.exp 1) := by
  rw [hΔ, hAsymp, fisherHorizon]

/-! ## Exponent competition

The paper's stability condition is the scalar inequality `1 - σ > γ`. These
lemmas formalize the RH-side and off-critical-line-side algebra.
-/

def StableExponent (σ γ : ℝ) : Prop :=
  1 - σ > γ

theorem stable_exponent_under_RH {γ : ℝ} (hγ : γ < (1 : ℝ) / 2) :
    StableExponent ((1 : ℝ) / 2) γ := by
  unfold StableExponent
  nlinarith

theorem not_stable_if_prime_fluctuation_too_slow
    {σ γ : ℝ} (hγ : 1 - σ < γ) :
    ¬ StableExponent σ γ := by
  intro hs
  unfold StableExponent at hs
  linarith

/-! ## Generic finite-size phase-transition interface

The Lambert-W phase transition formula is formalized as a conditional theorem:
if the critical size equation is one of the model assumptions, Lean can transport
it into whatever downstream notation the solver proof uses.
-/

def criticalC (δ kmax dClause β : ℝ) : ℝ :=
  4 * δ ^ 2 / (kmax ^ 2 * dClause * β)

def criticalLogK (LambertW : ℝ → ℝ) (C : ℝ) : ℝ :=
  -C * LambertW (-(1 / C))

theorem lambertW_phase_transition_identity
    (LambertW : ℝ → ℝ)
    {δ kmax dClause β logKstar : ℝ}
    (hcrit : logKstar = criticalLogK LambertW (criticalC δ kmax dClause β)) :
    logKstar =
      -(criticalC δ kmax dClause β) *
        LambertW (-(1 / criticalC δ kmax dClause β)) := by
  simpa [criticalLogK] using hcrit

end NitroSAT
