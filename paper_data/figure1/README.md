# figure1 — escape-probability figure

**Author:** Yimin Huang
**Affiliation:** Fudan University
**Email:** huangym23@m.fudan.edu.cn (hyimin0924@gmail.com)

`plot_escape_probability.py` produces `escape_probability.pdf`.

The left panel is a schematic of the cell-by-cell solver. The right panel
plots three escape-probability quantities against line-center optical
depth τ_ℓ. Two of the three are exact line-by-line ports of Cloudy's
C++ implementations in `cloudy/source/rt_escprob.cpp`;
the third is a pedagogical simplification of the CRD+wing case.

This README documents the correspondence so the figure can be cited as
"Cloudy's escape probability" without ambiguity.

## 1. `beta_K2(tau)` ↔ Cloudy `esca0k2(taume)`  (rt_escprob.cpp, lines 425–481)

**Status: exact port for τ ≥ 0.**

| Item | Cloudy (C++) | Python | Match |
|---|---|---|---|
| τ rescale | `tau = taume*SQRTPI` | `tau = taume*SQRTPI` | ✓ |
| Small-τ branch (`tau < 0.01`) | `1 - 2*tau` | `1.0 - 2.0*tau` | ✓ |
| Intermediate-τ (`tau ≤ 11`) | `tau/2.5066283 * log(tau/SQRTPI) + suma/sumb` | identical | ✓ |
| Large-τ (`tau > 11`) | `(sumc/sumd)/(2*tau*sqrt(log(tau/SQRTPI)))` | identical | ✓ |
| Horner polynomials `suma, sumb, sumc, sumd` | as in Cloudy | identical | ✓ |
| Coefficients `a[5], b[6], c[5], d[6]` | (see Cloudy source) | `_A, _B, _C, _D` | ✓ element-for-element |

Coefficients are from Hummer 1981, *JQSRT* 26, 187 (rational approximation
to Hummer's K₂ one-sided escape probability for a pure Doppler profile).
The Hummer & Rybicki 1982, *ApJ* 254, 767 paper plots this same K₂
function as the unmarked lower-envelope curve in their Fig. 2
(their Eqs. 2.11–2.12).

The only deliberate divergence is the negative-τ (maser) branch:
Cloudy calls `escmase(taume)`, the Python returns the small-τ formula.
This is irrelevant for the plotted positive-τ range.

## 2. `beta_PRD(tau, a)` ↔ Cloudy `esc_PRD_1side(tau, a)`  (rt_escprob.cpp, lines 116–162)

**Status: exact port for τ ≥ 0.**

| Item | Cloudy (C++) | Python | Match |
|---|---|---|---|
| Product | `atau = a*tau` | identical | ✓ |
| Prefactor | `3 * pow(2*a, -0.12)` | `3.0 * (2.0*a)**(-0.12)` | ✓ |
| Branch `atau > 1` | `b = 1.6 + pref/(1+atau)` | identical | ✓ |
| Branch `atau ≤ 1` | `b = 1.6 + pref*sqrt(atau)/(1+sqrt(atau))` | identical | ✓ |
| Cap | `b = MIN2(6, b)` | `b = min(6.0, b)` | ✓ |
| Return | `1/(1 + b*tau)` | `1.0/(1.0 + b*tau)` | ✓ |

Negative-τ behavior again differs (Cloudy: `escmase`; Python: NaN),
which does not affect the plotted curve at `a = 1e-3`.

## 3. `β_ℓ = (1 − p_w) β_{K_2}(τ) + p_w`  (right-panel dashed curves)

**Status: pedagogical simplification, NOT a port of `esc_CRDwing_1side`.**

The script blends `beta_K2` with a *constant* wing-escape floor `p_w`,
following the schematic Hummer 1982 / Ferland et al. 2017 picture used
in lectures and textbook discussions.

Cloudy's actual subordinate-line CRD-with-damping-wings routine is
`esc_CRDwing_1side(tau, a)` (rt_escprob.cpp, lines 164–203), which uses
a **τ- and a-dependent** wing fraction:

```
scal  = a*(1+a+τ) / ((1+a)^2 + a*τ)
pwing = scal * sqrt(a) / sqrt(a + 2.25*sqrt(π)*τ)        (τ > 0)
β     = esca0k2(τ) * (1 − pwing) + pwing
```

The three dashed curves with `p_w ∈ {1e-4, 1e-3, 1e-2}` are illustrative
of how a non-zero wing escape probability flattens the large-τ tail of
β_ℓ — they do **not** reproduce Cloudy's `esc_CRDwing_1side` output.
If a faithful reproduction is needed, port the `pwing(τ, a)` expression
above and replace the constant-`p_w` blend.

## Summary

- `beta_K2(tau)` reproduces Cloudy `esca0k2` exactly (τ ≥ 0). ✓
- `beta_PRD(tau, a)` reproduces Cloudy `esc_PRD_1side` exactly (τ ≥ 0). ✓
- The constant-`p_w` curves are schematic, not a port of
  `esc_CRDwing_1side`. ✗ (intentional simplification)
