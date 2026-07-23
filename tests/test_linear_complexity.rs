/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

use dsi_progress_logger::prelude::*;

use modlin::fp::{Field, linear_complexity};

fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    (a as u128 * b as u128 % p as u128) as u64
}

fn addmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 + b as u128) % p as u128) as u64
}

fn submod(a: u64, b: u64, p: u64) -> u64 {
    if a >= b { a - b } else { a + p - b }
}

fn powmod(a: u64, mut e: u64, p: u64) -> u64 {
    let pm = p as u128;
    let mut r: u128 = 1 % pm;
    let mut base = a as u128 % pm;
    while e > 0 {
        if e & 1 == 1 {
            r = r * base % pm;
        }
        base = base * base % pm;
        e >>= 1;
    }
    r as u64
}

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
fn linear_complexity_of_lfsr_sequence() {
    let p = 2305843009213693951; // 2⁶¹ − 1
    let f = Field::new(p);
    // An order-3 recurrence s[i] = 3·s[i-1] + 5·s[i-2] + 7·s[i-3] has LC ≤ 3.
    let mut s = vec![1, 2, 3];
    for i in 3..40 {
        let v = addmod(
            addmod(mulmod(3, s[i - 1], p), mulmod(5, s[i - 2], p), p),
            mulmod(7, s[i - 3], p),
            p,
        );
        s.push(v);
    }
    assert_eq!(linear_complexity(&f, &s, no_logging![]), 3);
    // The all-equal sequence has linear complexity 1.
    assert_eq!(linear_complexity(&f, &[4; 20], no_logging![]), 1);
    // The zero sequence has linear complexity 0.
    assert_eq!(linear_complexity(&f, &[0; 20], no_logging![]), 0);
}

#[test]
fn matches_reference() {
    fn mix(x: u64) -> u64 {
        let mut z = x.wrapping_mul(0x9E3779B97F4A7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z ^ (z >> 31)
    }
    let mut seed = 0u64;
    let mut next = || {
        seed = seed.wrapping_add(1);
        mix(seed)
    };

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
            let mut s: Vec<u64> = (0..k).map(|_| next() % p).collect();
            for i in k..150 {
                let mut v = 0;
                for (j, &cj) in coeffs.iter().enumerate() {
                    v = addmod(v, mulmod(cj, s[i - 1 - j], p), p);
                }
                s.push(v);
            }
            let got = linear_complexity(&f, &s, no_logging![]);
            assert_eq!(got, lc_ref(p, &s), "lfsr p={p} k={k}");
            assert!(got <= k, "lfsr complexity {got} exceeds order {k}");
        }
    }
}
