/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

use modlin::stats::{corank_prob, corank_tail_pvalue, lc_left_tail_pvalue};

const M61: u64 = 2305843009213693951; // 2⁶¹ − 1

#[test]
fn corank_distribution_sane() {
    // Full-rank probability and the 1/p deficiency leading term.
    let p0 = corank_prob(M61, 64, 0);
    assert!((p0 - (1.0 - 1.0 / M61 as f64)).abs() / p0 < 1e-12);
    let p1 = corank_prob(M61, 64, 1);
    assert!((p1 - 1.0 / M61 as f64).abs() / p1 < 1e-6);
    // Over F₂ the full-rank probability tends to ≈ 0.2888.
    assert!((corank_prob(2, 256, 0) - 0.2887880951).abs() < 1e-6);
    // The distribution sums to ≈ 1.
    let s: f64 = (0..=64).map(|d| corank_prob(M61, 64, d)).sum();
    assert!((s - 1.0).abs() < 1e-12);
}

#[test]
fn single_observation_tails() {
    let n = 305;
    assert_eq!(corank_tail_pvalue(M61, n, 0), 1.0); // full rank: unremarkable
    let p1 = corank_tail_pvalue(M61, n, 1);
    assert!((p1 - 1.0 / M61 as f64).abs() / p1 < 1e-6); // one deficiency ≈ 1/p
    // Deep deficiency / far-below-mode underflow to the smallest positive float,
    // never to 0 (so a printed p-value does not look like a broken test).
    assert_eq!(corank_tail_pvalue(M61, n, 33), f64::MIN_POSITIVE);
    // At the mode of a large field the exact tail (≈ 1 − 1/(p+1)) rounds to 1.0.
    assert_eq!(lc_left_tail_pvalue(M61, 10_000, 5000), 1.0);
    assert_eq!(lc_left_tail_pvalue(M61, 10_000, 272), f64::MIN_POSITIVE);
    let p = lc_left_tail_pvalue(2, 100, 40);
    assert!(p > 0.0 && p < 1.0); // a moderate left-tail value over F₂
}

#[test]
fn lc_left_tail_is_exact_finite_n() {
    // Exact finite-n values where the old asymptotic/“1.0 at the mode” shortcut
    // was wrong: over F₂, Pr[L₂ ≤ 1] = 3/4, Pr[L₅ ≤ 2] = 11/32, Pr[L₅ ≤ 3] = 27/32.
    assert!((lc_left_tail_pvalue(2, 2, 1) - 0.75).abs() < 1e-12);
    assert!((lc_left_tail_pvalue(2, 5, 2) - 0.343_75).abs() < 1e-12);
    assert!((lc_left_tail_pvalue(2, 5, 3) - 0.843_75).abs() < 1e-12);

    // Agree with a direct sum of the linear-complexity pmf for every ℓ, across a
    // few small fields and lengths spanning both branches of the CDF.
    fn pmf(q: f64, n: usize, l: usize) -> f64 {
        if l == 0 {
            q.powi(-(n as i32))
        } else if 2 * l <= n {
            (q - 1.0) * q.powi(2 * l as i32 - n as i32 - 1)
        } else {
            (q - 1.0) * q.powi(n as i32 - 2 * l as i32)
        }
    }
    for &q in &[2u64, 3, 5, 7] {
        for n in 1..=24usize {
            let mut cdf = 0.0;
            for ell in 0..=n {
                cdf += pmf(q as f64, n, ell);
                let got = lc_left_tail_pvalue(q, n, ell);
                assert!(
                    (got - cdf.min(1.0)).abs() <= 1e-9 * cdf + 1e-15,
                    "q={q} n={n} ell={ell}: got {got}, pmf sum {cdf}"
                );
            }
        }
    }
}
