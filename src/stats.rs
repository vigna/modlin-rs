/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! The 𝐅*ₚ* null distributions of the matrix rank and of the linear
//! complexity, the one-sided per-sample *p*-values derived from them, and *p*-value
//! formatting.

use std::borrow::Cow;
use std::f64::consts::LN_2;

/// Readable floor reported for an underflowed *p*-value instead of 0: the
/// smallest power of ten that is still a normal `f64`.
pub const P_FLOOR: f64 = 1e-307;

/// Natural logarithm of a positive integer *p*.
fn ln_u64(p: u64) -> f64 {
    let e = p.ilog2();
    let mantissa = p as f64 / (1u64 << e) as f64; // p / 2^e ∈ [1..2]
    (e as f64).mul_add(LN_2, mantissa.ln())
}

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

/// Pr[corank = *d*] for a uniform random *n* × *n* matrix over 𝐅*ₚ*, exact for
/// finite *n*.
///
/// In logarithmic terms:
///
/// ln Pr[corank = *d*] = −*d*²·ln *p* + 2·Σ_{*j*=*d*+1}^{*n*} ln(1−*p*⁻*ʲ*) − Σ_{*k*=1}^{*n*−*d*} ln(1−*p*⁻*ᵏ*).
///
/// In particular, Pr[corank = 0] = ∏_{*i*=1}^{*n*}(1−*p*⁻*ⁱ*) and Pr[corank =
/// *d*] ≈ *p*⁻*ᵈ*² for large *p*. Returns 0 on underflow.
pub fn corank_prob(p: u64, n: usize, d: usize) -> f64 {
    if d > n {
        return 0.0;
    }
    // The dominant −d²·ln p term uses the integer-exact logarithm; the
    // corrections sum tiny p⁻ʲ terms, where f64's 1/p is far more than enough.
    let pf = p as f64;
    let lp =
        -(d as f64).powi(2) * ln_u64(p) + 2.0 * sum_ln_1m(pf, d + 1, n) - sum_ln_1m(pf, 1, n - d);
    if lp < -745.0 { 0.0 } else { lp.exp() }
}

/// One-sided *p*-value of a single matrix's corank for the deficiency alternative:
/// the upper tail Pr[corank ≥ *c*] for a uniform random *n* × *n* matrix over 𝐅*ₚ*.
///
/// Equals 1 for a full-rank matrix (*c* = 0) and ≈ *p*⁻*ᶜ*² for corank *c*.
/// When the tail underflows it is floored to [`P_FLOOR`] rather than 0.
pub fn corank_tail_pvalue(p: u64, n: usize, c: usize) -> f64 {
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
    s.clamp(P_FLOOR, 1.0)
}

/// One-sided *p*-value of a single sequence's linear complexity for the
/// low-complexity alternative: the lower tail Pr[*Lₙ* ≤ *ℓ*] for a uniform
/// random length-*n* sequence over 𝐅*ₚ*.
///
/// This is the closed-form CDF of the linear-complexity distribution, whose two
/// branches meet at the mode floor ⌊*n*/2⌋:
///
/// Pr[*Lₙ* ≤ *ℓ*] = (1 + *p*²*ˡ*⁺¹) / (*p*ⁿ(*p* + 1))      for 2*ℓ* ≤ *n*,
///
/// Pr[*Lₙ* ≤ *ℓ*] = 1 − (*p*ⁿ⁻²*ˡ* − *p*⁻ⁿ) / (*p* + 1)    for 2*ℓ* > *n*.
///
/// Below the mode the lower branch is ≈ *p*²*ˡ*⁻*ⁿ*, dropping by a factor *p*⁻²
/// for every step away; at or above the mode it is ≈ 1 (so for a large field
/// the mode rounds to exactly 1.0 in `f64`). When the tail underflows it is
/// floored to [`P_FLOOR`] rather than 0.
pub fn lc_left_tail_pvalue(p: u64, n: usize, ell: usize) -> f64 {
    if ell >= n {
        return 1.0; // complexity cannot exceed the sequence length
    }
    let lnp = ln_u64(p);
    let pr = if 2 * ell <= n {
        // Lower branch (at or below the mode floor ⌊n/2⌋), in logs because the
        // p^{2ℓ+1} term can overflow f64 and the whole tail can underflow far
        // below the f64 floor.
        let a = (2 * ell + 1) as f64 * lnp; // ln p^{2ℓ+1}
        let ln_num = if a > 0.0 {
            a + (-a).exp().ln_1p() // ln(1 + p^{2ℓ+1})
        } else {
            a.exp().ln_1p()
        };
        let ln_den = n as f64 * lnp + (p as f64 + 1.0).ln(); // ln(p^n (p+1))
        (ln_num - ln_den).exp()
    } else {
        // Upper branch (above the mode floor): 1 − (p^{n−2ℓ} − p^{−n})/(p+1).
        // Both exponents are negative, so the subtracted term is a small positive.
        let hi = ((n as f64 - 2.0 * ell as f64) * lnp).exp(); // p^{n−2ℓ} < 1
        let lo = (-(n as f64) * lnp).exp(); // p^{−n}
        1.0 - (hi - lo) / (p as f64 + 1.0)
    };
    pr.clamp(P_FLOOR, 1.0)
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
