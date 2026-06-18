/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

use modlin::stats::{corank_prob, corank_tail_pvalue, lc_left_tail_pvalue};

const M61: f64 = 2305843009213693951.0; // 2⁶¹ − 1

#[test]
fn corank_distribution_sane() {
    // Full-rank probability and the 1/p deficiency leading term.
    let p0 = corank_prob(M61, 64, 0);
    assert!((p0 - (1.0 - 1.0 / M61)).abs() / p0 < 1e-12);
    let p1 = corank_prob(M61, 64, 1);
    assert!((p1 - 1.0 / M61).abs() / p1 < 1e-6);
    // Over F₂ the full-rank probability tends to ≈ 0.2888.
    assert!((corank_prob(2.0, 256, 0) - 0.2887880951).abs() < 1e-6);
    // The distribution sums to ≈ 1.
    let s: f64 = (0..=64).map(|d| corank_prob(M61, 64, d)).sum();
    assert!((s - 1.0).abs() < 1e-12);
}

#[test]
fn single_observation_tails() {
    let n = 305;
    assert_eq!(corank_tail_pvalue(M61, n, 0), 1.0); // full rank: unremarkable
    let p1 = corank_tail_pvalue(M61, n, 1);
    assert!((p1 - 1.0 / M61).abs() / p1 < 1e-6); // one deficiency ≈ 1/p
    // Deep deficiency / far-below-mode underflow to the smallest positive float,
    // never to 0 (so a printed p-value does not look like a broken test).
    assert_eq!(corank_tail_pvalue(M61, n, 33), f64::MIN_POSITIVE);
    assert_eq!(lc_left_tail_pvalue(M61, 10_000, 5000), 1.0); // at the mode
    assert_eq!(lc_left_tail_pvalue(M61, 10_000, 272), f64::MIN_POSITIVE);
    let p = lc_left_tail_pvalue(2.0, 100, 40);
    assert!(p > 0.0 && p < 1.0); // a moderate left-tail value over F₂
}
