/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Helpers shared by unit tests.

/// A bijective mixing function on 64-bit integers (the SplitMix finalizer
/// applied to a multiple of the golden ratio), used to generate test data.
pub fn mix(x: u64) -> u64 {
    let mut z = x.wrapping_mul(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// A pseudorandom polynomial over 𝐅*ₚ* with `len` coefficients (from the
/// lowest degree) and a nonzero leading coefficient.
pub fn random_poly(len: usize, p: u64, seed: u64) -> Vec<u64> {
    let mut v: Vec<u64> = (0..len as u64)
        .map(|i| mix(seed ^ i.wrapping_mul(0x9E3779B97F4A7C15)) % p)
        .collect();
    if let Some(last) = v.last_mut() {
        if *last == 0 {
            *last = 1;
        }
    }
    v
}
