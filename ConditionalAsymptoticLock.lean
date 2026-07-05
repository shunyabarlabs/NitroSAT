import Mathlib

/-!
# Conditional Asymptotic Lock

This file proves the formal core of the conditional Asymptotic Lock claim.
It does not assume or prove the hard analytic/spectral ingredients. Instead,
it treats them as named hypotheses:

1. the solver's asymptotic stability at spectral exponent `gamma` is exactly
   the scalar exponent threshold `1 - sigma > gamma`;
2. the critical graph family realizes every subcritical exponent `gamma < 1/2`;
3. optionally, RH is equivalent to the exponent bound `sigma <= 1/2`.

Under those hypotheses, Lean proves the downstream lock:

  asymptotic stability on the critical family <-> sigma <= 1/2.

No `sorry`, `axiom`, or `admit` is used.
-/

noncomputable section

namespace NitroSAT
namespace AsymptoticLock

/--
The scalar exponent condition: prime fluctuation decay beats spectral-gap decay.
Here `sigma` is the prime-fluctuation/zeta-zero exponent and `gamma` is the
spectral decay exponent.
-/
def ExponentStable (sigma gamma : Real) : Prop :=
  1 - sigma > gamma

/--
Stability for every graph family whose spectral decay exponent is subcritical,
i.e. every `gamma < 1/2`.
-/
def SubcriticalStable (sigma : Real) : Prop :=
  forall gamma : Real, gamma < (1 : Real) / 2 -> ExponentStable sigma gamma

/--
Pure exponent algebra: stability for all subcritical spectral exponents is
equivalent to `sigma <= 1/2`.
-/
theorem subcritical_stable_iff_sigma_le_half (sigma : Real) :
    SubcriticalStable sigma <-> sigma <= (1 : Real) / 2 := by
  constructor
  · intro h
    by_contra hle
    have hsigma : (1 : Real) / 2 < sigma := not_le.mp hle
    let gamma : Real := ((1 - sigma) + ((1 : Real) / 2)) / 2
    have hgamma_lt : gamma < (1 : Real) / 2 := by
      dsimp [gamma]
      linarith
    have hstable : ExponentStable sigma gamma := h gamma hgamma_lt
    unfold ExponentStable at hstable
    dsimp [gamma] at hstable
    linarith
  · intro hsigma gamma hgamma
    unfold ExponentStable
    linarith

/--
Conditional Asymptotic Lock.

Hypothesis `threshold` packages the linear-coupling/no-threshold-shift claim:
solver stability at exponent `gamma` is equivalent to `1 - sigma > gamma`.

Hypothesis `critical_family` packages the graph-family/genericity claim:
asymptotic stability on the critical family is equivalent to stability for all
subcritical exponents `gamma < 1/2`.

Then the conditional lock follows by real arithmetic.
-/
theorem conditional_asymptotic_lock
    {SolverStableAtExponent : Real -> Real -> Prop}
    {AsymptoticallyStableOnCriticalFamily : Real -> Prop}
    (threshold :
      forall sigma gamma : Real,
        SolverStableAtExponent sigma gamma <-> ExponentStable sigma gamma)
    (critical_family :
      forall sigma : Real,
        AsymptoticallyStableOnCriticalFamily sigma <->
          forall gamma : Real,
            gamma < (1 : Real) / 2 -> SolverStableAtExponent sigma gamma)
    (sigma : Real) :
    AsymptoticallyStableOnCriticalFamily sigma <-> sigma <= (1 : Real) / 2 := by
  rw [critical_family sigma]
  constructor
  · intro h
    have hsub : SubcriticalStable sigma := by
      intro gamma hgamma
      exact (threshold sigma gamma).mp (h gamma hgamma)
    exact (subcritical_stable_iff_sigma_le_half sigma).mp hsub
  · intro hsigma gamma hgamma
    have hsub : SubcriticalStable sigma :=
      (subcritical_stable_iff_sigma_le_half sigma).mpr hsigma
    exact (threshold sigma gamma).mpr (hsub gamma hgamma)

/--
If one separately identifies RH with the exponent bound `sigma <= 1/2`, then
the conditional lock can be restated as stability iff RH.
-/
theorem conditional_asymptotic_lock_iff_RH
    {SolverStableAtExponent : Real -> Real -> Prop}
    {AsymptoticallyStableOnCriticalFamily : Real -> Prop}
    (RH : Prop)
    (threshold :
      forall sigma gamma : Real,
        SolverStableAtExponent sigma gamma <-> ExponentStable sigma gamma)
    (critical_family :
      forall sigma : Real,
        AsymptoticallyStableOnCriticalFamily sigma <->
          forall gamma : Real,
            gamma < (1 : Real) / 2 -> SolverStableAtExponent sigma gamma)
    (sigma : Real)
    (rh_exponent_equivalence : RH <-> sigma <= (1 : Real) / 2) :
    AsymptoticallyStableOnCriticalFamily sigma <-> RH := by
  have hlock :
      AsymptoticallyStableOnCriticalFamily sigma <-> sigma <= (1 : Real) / 2 :=
    conditional_asymptotic_lock
      (SolverStableAtExponent := SolverStableAtExponent)
      (AsymptoticallyStableOnCriticalFamily := AsymptoticallyStableOnCriticalFamily)
      threshold critical_family sigma
  exact hlock.trans rh_exponent_equivalence.symm

end AsymptoticLock
end NitroSAT
