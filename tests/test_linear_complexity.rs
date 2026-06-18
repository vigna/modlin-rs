/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

//! Black-box tests of the Berlekamp–Massey linear-complexity routine
//! ([`modlin::gf::linear_complexity`]).

use dsi_progress_logger::prelude::*;
use modlin::gf::{Field, linear_complexity};

fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    (a as u128 * b as u128 % p as u128) as u64
}

fn addmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 + b as u128) % p as u128) as u64
}

#[test]
fn linear_complexity_of_lfsr_sequence() {
    let p = 2305843009213693951u64; // 2⁶¹ − 1
    let f = Field::new(p);
    // An order-3 recurrence s[i] = 3·s[i-1] + 5·s[i-2] + 7·s[i-3] has LC ≤ 3.
    let mut s = vec![1u64, 2, 3];
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
    assert_eq!(linear_complexity(&f, &[4u64; 20], no_logging![]), 1);
    // The zero sequence has linear complexity 0.
    assert_eq!(linear_complexity(&f, &[0u64; 20], no_logging![]), 0);
}
