/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Helpers shared by the integration tests.

// Each test file uses a different subset of the helpers.
#![allow(dead_code)]

/// A bijective mixing function on 64-bit integers (the SplitMix finalizer
/// applied to a multiple of the golden ratio), used to generate test data.
pub fn mix(x: u64) -> u64 {
    let mut z = x.wrapping_mul(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// A stream of test data: the values of [`mix`] on `seed` + 1, `seed` + 2, ….
pub fn stream(mut seed: u64) -> impl FnMut() -> u64 {
    move || {
        seed = seed.wrapping_add(1);
        mix(seed)
    }
}

pub fn addmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 + b as u128) % p as u128) as u64
}

pub fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    (a as u128 * b as u128 % p as u128) as u64
}

pub fn submod(a: u64, b: u64, p: u64) -> u64 {
    if a >= b { a - b } else { a + p - b }
}

pub fn powmod(a: u64, mut e: u64, p: u64) -> u64 {
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

/// The first `n` values of the sequence over 𝐅*ₚ* that starts with `init`
/// and then satisfies *xᵢ* = Σⱼ `coeffs`ⱼ *xᵢ*₋₁₋ⱼ, where `init` has at least
/// as many values as `coeffs`.
pub fn recurrence(coeffs: &[u64], init: &[u64], n: usize, p: u64) -> Vec<u64> {
    let mut s = init.to_vec();
    while s.len() < n {
        let i = s.len();
        let v = coeffs
            .iter()
            .enumerate()
            .fold(0, |v, (j, &c)| addmod(v, mulmod(c, s[i - 1 - j], p), p));
        s.push(v);
    }
    s.truncate(n);
    s
}
