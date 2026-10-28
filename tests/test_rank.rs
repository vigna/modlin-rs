/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

use dsi_progress_logger::prelude::*;

use modlin::fp::{Field, rank};

mod common;
use common::{addmod, mulmod, powmod, stream, submod};

fn rank_ref(p: u64, m: &mut [u64], n: usize) -> usize {
    assert_eq!(m.len(), n * n);
    let mut rank = 0;
    for col in 0..n {
        if rank == n {
            break;
        }
        let Some(pr) = (rank..n).find(|&r| m[r * n + col] != 0) else {
            continue;
        };
        if pr != rank {
            for c in col..n {
                m.swap(pr * n + c, rank * n + c);
            }
        }
        let inv = powmod(m[rank * n + col], p - 2, p);
        for r in (rank + 1)..n {
            let head = m[r * n + col];
            if head != 0 {
                let f = mulmod(head, inv, p);
                for c in col..n {
                    m[r * n + c] = submod(m[r * n + c], mulmod(f, m[rank * n + c], p), p);
                }
            }
        }
        rank += 1;
    }
    rank
}

#[test]
fn test_rank_of_known_matrices() {
    let f = Field::new(2305843009213693951); // 2⁶¹ − 1

    let mut id = vec![0; 9];
    for i in 0..3 {
        id[i * 3 + i] = 1;
    }
    assert_eq!(rank(&f, &mut id, 3, no_logging![]), 3);

    // r2 = r0 + r1 → rank 2.
    let mut m = vec![1, 2, 3, 4, 5, 6, 5, 7, 9];
    assert_eq!(rank(&f, &mut m, 3, no_logging![]), 2);

    let mut z = vec![0; 16];
    assert_eq!(rank(&f, &mut z, 4, no_logging![]), 0);
}

#[test]
fn test_rank_over_gf2() {
    let f = Field::new(2);
    let mut m = vec![1, 1, 0, 0, 1, 1, 1, 0, 1];
    assert_eq!(rank(&f, &mut m, 3, no_logging![]), 2);
}

#[test]
fn test_blocked_matches_reference() {
    let mut next = stream(0);

    for &p in &[2, 97, 2305843009213693951, 9223372036854775783] {
        let f = Field::new(p);
        let sizes: &[usize] = if cfg!(feature = "slow_tests") {
            &[1, 2, 5, 63, 64, 65, 129, 200, 257, 300]
        } else {
            &[1, 2, 5, 63, 64, 65, 129, 200]
        };
        for &n in sizes {
            for &target in &[n, n / 2, n / 3, 1, 0] {
                let target = target.min(n);
                // Build an n×n matrix of exact rank ≤ target: target random
                // independent rows, the rest random linear combinations of them.
                let mut a = vec![0; n * n];
                for r in 0..target {
                    for c in 0..n {
                        a[r * n + c] = next() % p;
                    }
                }
                for r in target..n {
                    for k in 0..target {
                        let coeff = next() % p;
                        if coeff == 0 {
                            continue;
                        }
                        for c in 0..n {
                            a[r * n + c] = addmod(a[r * n + c], mulmod(coeff, a[k * n + c], p), p);
                        }
                    }
                }
                // Shuffle rows so pivots are not pre-sorted.
                for r in (1..n).rev() {
                    let s = (next() as usize) % (r + 1);
                    if s != r {
                        for c in 0..n {
                            a.swap(r * n + c, s * n + c);
                        }
                    }
                }

                let mut b = a.clone();
                let expected = rank_ref(p, &mut a, n);
                let got = rank(&f, &mut b, n, no_logging![]);
                assert_eq!(got, expected, "p={p} n={n} target={target}");
                assert!(expected <= target);
            }
        }
    }
}
