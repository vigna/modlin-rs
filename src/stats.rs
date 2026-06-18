/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

//! The **F**_ₚ_ null distributions of the matrix rank and of the linear
//! complexity, the one-sided per-sample *p*-values derived from them, and *p*-value
//! formatting.

use std::borrow::Cow;

/// Σ_{*j*=*a*}^{*b*} ln(1 − *p*⁻*ʲ*); terms past *p*⁻*ʲ* ≈ 0 are dropped.
fn sum_ln_1m(p: f64, a: usize, b: usize) -> f64 {
    if a > b {
        return 0.0;
    }
    let inv = 1.0 / p;
    let mut pj = inv.powi(a as i32); // p⁻ᵃ, underflowing to 0 contributes nothing
    let mut s = 0.0;
    for _ in a..=b {
        if pj < 1e-18 {
            break; // ln(1 − p⁻ʲ) ≈ 0 from here on
        }
        s += (-pj).ln_1p(); // ln(1 − p⁻ʲ)
        pj *= inv;
    }
    s
}

/// Pr[corank = *d*] for a uniform random *n* × *n* matrix over **F**_ₚ_, exact for
/// finite *n*.
///
/// In logarithmic terms:
///
/// ln Pr[corank = *d*] = −*d*²·ln *p* + 2·Σ_{*j*=*d*+1}^{*n*} ln(1−*p*⁻*ʲ*) − Σ_{*k*=1}^{*n*−*d*} ln(1−*p*⁻*ᵏ*).
///
/// In particular, Pr[corank = 0] = ∏_{*i*=1}^{*n*}(1−*p*⁻*ⁱ*) and Pr[corank =
/// *d*] ≈ *p*⁻*ᵈ*² for large *p*. Returns 0 on underflow.
pub fn corank_prob(p: f64, n: usize, d: usize) -> f64 {
    if d > n {
        return 0.0;
    }
    let lp = -(d as f64).powi(2) * p.ln() + 2.0 * sum_ln_1m(p, d + 1, n) - sum_ln_1m(p, 1, n - d);
    if lp < -745.0 { 0.0 } else { lp.exp() }
}

/// One-sided *p*-value of a single matrix's corank for the deficiency alternative:
/// the upper tail Pr[corank ≥ *c*] for a uniform random *n* × *n* matrix over **F**_ₚ_.
///
/// Equals 1 for a full-rank matrix (*c* = 0) and ≈ *p*⁻*ᶜ*² for corank *c*.
/// When the tail is smaller than the smallest positive `f64` it is floored to
/// [`f64::MIN_POSITIVE`] rather than 0, so a printed *p*-value reads as
/// astronomically small instead of looking like a broken test.
pub fn corank_tail_pvalue(p: f64, n: usize, c: usize) -> f64 {
    if c == 0 {
        return 1.0;
    }
    let mut s = 0.0;
    let mut d = c;
    while d <= n {
        let pr = corank_prob(p, n, d);
        if pr == 0.0 {
            break; // the tail has underflowed
        }
        s += pr;
        d += 1;
    }
    s.clamp(f64::MIN_POSITIVE, 1.0)
}

/// One-sided *p*-value of a single sequence's linear complexity for the
/// low-complexity alternative: the lower tail Pr[*Lₙ* ≤ *ℓ*] for a uniform random
/// length-*n* sequence over **F**_ₚ_.
///
/// Equals 1 at or above the mode ⌈*n*/2⌉ (not anomalously low) and ≈
/// *p*²*ˡ*⁻*ⁿ*⁺¹/(*p*+1) below it. When that is smaller than the smallest
/// positive `f64` it is floored to [`f64::MIN_POSITIVE`] rather than 0, so a
/// printed *p*-value reads as astronomically small instead of looking like a
/// broken test.
pub fn lc_left_tail_pvalue(p: f64, n: usize, ell: usize) -> f64 {
    if 2 * ell >= n {
        return 1.0; // at or above the mode
    }
    let exponent = 2.0 * ell as f64 - n as f64 + 1.0;
    let log10p = exponent * p.log10() - (p + 1.0).log10();
    if log10p <= f64::MIN_10_EXP as f64 {
        f64::MIN_POSITIVE
    } else {
        10f64.powf(log10p).min(1.0)
    }
}

/// Pretty-prints a *p*-value, writing `0` for 0.0 and `1` for 1.0.
pub fn pretty_p_value(p: f64) -> Cow<'static, str> {
    if p == 0.0 {
        "0".into()
    } else if p == 1.0 {
        "1".into()
    } else {
        format!("{:?}", p).into()
    }
}
