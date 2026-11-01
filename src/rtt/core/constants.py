"""Certified constants with provenance.

Exact values are stored as (numerator, denominator) integer pairs so that rigorous code can
build exact rationals (``flint.fmpq``) or exact Arb balls from them; the float aliases are for
display and for the float proposal stages only.
"""
from __future__ import annotations

from fractions import Fraction

# ---------------------------------------------------------------------------------------------
# Universal ceiling  k_B Tc <= A / A0   (main theorem; numerical hypothesis Proposition 6 = H(a0, eta) certified by wE1 (126 cells,
# audited SOUND), Lean: TcR1.stable_of_ge / TcR1.unstable_imp_lt / TcR1.main in
# lean/TcR1/Basic.lean with H as hypothesis).
# ---------------------------------------------------------------------------------------------
A0_FRACTION = Fraction(433148145, 10**8)          # a0 = 4.33148145 exactly
A0 = float(A0_FRACTION)                           # 4.33148145 (display only)
A0_ETA_FRACTION = Fraction(1, 10**9)              # eta in H(a0, eta) (Prop. 6): 1e-9
C_UNIV = 1.0 / A0                                 # 0.2308678929...  (display only)

# Sharpness (Proposition 7): Tc >= A / 4.33148149 for the extremal Einstein mode.
A0_SHARP_LOWER_FRACTION = Fraction(433148149, 10**8)
EINSTEIN_OPT_OMEGA_OVER_T = Fraction(75748661, 10**7)   # certified near-extremal Omega/T

# Corollary (screening filter): Tc >= 300 K requires A >= 300 * a0 = 1299.444435 K.
DEFAULT_T_STAR_K = 300.0

# Defaults for the certify engine (configurable; recorded in each Verdict).
DEFAULT_MU_HI = 0.15            # mu* range [0, mu_hi] treated for CERTIFIED_CANDIDATE
DEFAULT_OMEGA_C_FACTOR = 10.0   # omega_c = factor * max omega  (ARCH section 4, N_c rule)
ARB_PREC_BITS = 128             # working precision of the Arb checks

PROVENANCE = {
    "A0": {
        "value": "433148145/10^8",
        "statement": "In isotropic Migdal-Eliashberg theory every alpha^2F, mu*>=0 and cutoff give "
                     "k_B Tc <= A/a0 (strictly: every T >= A/a0 is stable).",
        "proof": "main theorem (Lemmas 1-5, Prop. 6 by the wE1 interval certificate, audited); Lean TcR1",
        "lean": "TcR1.stable_of_ge, TcR1.main (hypothesis H(a0, 1e-9) = Prop. 6)",
    },
    "A0_SHARP_LOWER": {
        "value": "433148149/10^8",
        "statement": "Einstein mode with Omega/T = 7.5748661 is unstable at A/T = 4.33148149",
        "proof": "Prop. 7 (interval certificate of the Einstein-mode instability)",
    },
}
