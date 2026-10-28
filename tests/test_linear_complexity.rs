/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

use dsi_progress_logger::prelude::*;

use modlin::fp::{Field, HGCD_THRESHOLD, berlekamp_massey, linear_complexity};
use modlin::poly;

mod common;
use common::{addmod, mulmod, powmod, recurrence, stream, submod};

fn lc_ref(p: u64, s: &[u64]) -> usize {
    let n = s.len();
    let mut c = vec![0u64; n];
    let mut b = vec![0u64; n];
    c[0] = 1;
    b[0] = 1;
    let mut l = 0;
    let mut m = 1;
    let mut bb = 1;
    for i in 0..n {
        let mut d = s[i];
        for j in 1..=l {
            d = addmod(d, mulmod(c[j], s[i - j], p), p);
        }
        if d == 0 {
            m += 1;
            continue;
        }
        let coef = mulmod(d, powmod(bb, p - 2, p), p);
        let t = c.clone();
        for j in 0..(n - m) {
            c[j + m] = submod(c[j + m], mulmod(coef, b[j], p), p);
        }
        if 2 * l <= i {
            l = i + 1 - l;
            b = t;
            bb = d;
            m = 1;
        } else {
            m += 1;
        }
    }
    l
}

#[test]
fn test_linear_complexity_of_lfsr_sequence() {
    let p = 2305843009213693951; // 2⁶¹ − 1
    let f = Field::new(p);
    // An order-3 recurrence s[i] = 3·s[i-1] + 5·s[i-2] + 7·s[i-3] has LC ≤ 3.
    let s = recurrence(&[3, 5, 7], &[1, 2, 3], 40, p);
    assert_eq!(linear_complexity(&f, &s, no_logging![]), 3);
    // The all-equal sequence has linear complexity 1.
    assert_eq!(linear_complexity(&f, &[4; 20], no_logging![]), 1);
    // The zero sequence has linear complexity 0.
    assert_eq!(linear_complexity(&f, &[0; 20], no_logging![]), 0);
    // So does the empty sequence.
    assert_eq!(linear_complexity(&f, &[], no_logging![]), 0);
    assert_eq!(poly::half_gcd(&f, &[]), 0);
}

#[test]
fn test_matches_reference() {
    let mut next = stream(0);

    for &p in &[2u64, 3, 97, 65537, 2305843009213693951] {
        let f = Field::new(p);
        // Fully random streams of various lengths.
        for &n in &[1usize, 2, 3, 5, 16, 50, 100, 200] {
            let s: Vec<u64> = (0..n).map(|_| next() % p).collect();
            assert_eq!(
                linear_complexity(&f, &s, no_logging![]),
                lc_ref(p, &s),
                "random p={p} n={n}"
            );
        }
        // Order-k LFSR streams with random coefficients: low complexity (≤ k).
        for &k in &[1usize, 4, 10, 25] {
            let coeffs: Vec<u64> = (0..k).map(|_| next() % p).collect();
            let init: Vec<u64> = (0..k).map(|_| next() % p).collect();
            let s = recurrence(&coeffs, &init, 150, p);
            let got = linear_complexity(&f, &s, no_logging![]);
            assert_eq!(got, lc_ref(p, &s), "lfsr p={p} k={k}");
            assert!(got <= k, "lfsr complexity {got} exceeds order {k}");
        }
    }
}

/// Sequences on which the half-GCD algorithm is compared with the quadratic
/// Berlekamp–Massey algorithm: random streams, low-order recurrences, sparse
/// streams, and streams with long runs of leading or trailing zeros.
fn sample_sequences(p: u64, n: usize, seed: u64) -> Vec<(String, Vec<u64>)> {
    let mut next = stream(seed);
    let mut out = Vec::new();
    out.push(("random".to_string(), (0..n).map(|_| next() % p).collect()));
    for k in [1, 7, n / 5, n / 2 - 1, n / 2 + 3] {
        let coeffs: Vec<u64> = (0..k).map(|_| next() % p).collect();
        let init: Vec<u64> = (0..k.max(1)).map(|_| next() % p).collect();
        out.push((
            format!("recurrence of order {k}"),
            recurrence(&coeffs, &init, n, p),
        ));
    }
    out.push((
        "sparse".to_string(),
        (0..n)
            .map(|_| if next() % 8 == 0 { next() % p } else { 0 })
            .collect(),
    ));
    for z in [n / 3, n - 1] {
        let mut s = vec![0; z];
        s.extend((z..n).map(|_| next() % p));
        out.push((format!("{z} leading zeros"), s));
    }
    let mut s: Vec<u64> = (0..n / 2).map(|_| next() % p).collect();
    s.resize(n, 0);
    out.push(("trailing zeros".to_string(), s));
    out
}

#[test]
fn test_hgcd_matches_berlekamp_massey() {
    let lengths: &[usize] = if cfg!(feature = "slow_tests") {
        &[2, 65, 130, 257, 1000, 3001]
    } else {
        &[2, 65, 130, 257, 1000]
    };
    for p in [
        2,
        3,
        65537,
        (1 << 31) - 1,
        2305843009213693951,
        9223372036854775783,
    ] {
        let f = Field::new(p);
        for (i, &n) in lengths.iter().enumerate() {
            for (name, s) in sample_sequences(p, n, p.wrapping_mul(1000).wrapping_add(i as u64)) {
                assert_eq!(
                    poly::half_gcd(&f, &s),
                    berlekamp_massey(&f, &s, no_logging![]),
                    "p={p} n={n} {name}"
                );
            }
        }
    }
}

#[test]
fn test_dispatch_to_hgcd() {
    // From HGCD_THRESHOLD on, linear_complexity uses the half-GCD algorithm.
    for p in [3, 2305843009213693951] {
        let f = Field::new(p);
        for (name, s) in sample_sequences(p, HGCD_THRESHOLD + 501, p) {
            assert_eq!(
                linear_complexity(&f, &s, no_logging![]),
                berlekamp_massey(&f, &s, no_logging![]),
                "p={p} {name}"
            );
        }
    }
}
