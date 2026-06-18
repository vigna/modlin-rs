/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

//! Modular matrix rank and linear complexity on **F**_ₚ_.

use dsi_progress_logger::prelude::*;
use rayon::prelude::*;

/// High 128 bits of the 256-bit product *a* · *b*.
#[inline(always)]
fn mul_hi_128(a: u128, b: u128) -> u128 {
    let (a1, a0) = ((a >> 64) as u64 as u128, a as u64 as u128);
    let (b1, b0) = ((b >> 64) as u64 as u128, b as u64 as u128);
    let lo = a0 * b0;
    let mid1 = a1 * b0;
    let mid2 = a0 * b1;
    let hi = a1 * b1;
    // Carry out of bits [64..128) into the high half.
    let carry = (lo >> 64) + (mid1 & 0xFFFF_FFFF_FFFF_FFFF) + (mid2 & 0xFFFF_FFFF_FFFF_FFFF);
    hi + (mid1 >> 64) + (mid2 >> 64) + (carry >> 64)
}

/// Column-panel width for the blocked elimination: the trailing update is
/// deferred and applied once per panel of this many columns, cutting passes
/// over the (out-of-cache) trailing submatrix by roughly this factor.
const PANEL: usize = 64;
/// Trailing column-tile width: the panel's PANEL reduced pivot rows, restricted
/// to a tile, stay resident in cache while every far row's tile is updated.
const TILE: usize = 256;

/// A prime field **F**_ₚ_ with *p* < 2⁶³.
#[derive(Clone, Copy)]
pub struct Field {
    p: u64,
    /// Barrett reciprocal μ = ⌊(2¹²⁸ − 1)/*p*⌋ ≤ ⌊2¹²⁸/*p*⌋, used for
    /// division-free reduction of any `u128` modulo *p* (see
    /// [Self::reduce_u128]).
    mu: u128,
    /// Largest number of products that can be summed in a `u128` without
    /// overflow; the blocked trailing update accumulates this many before a
    /// single reduction (see [Self::reduce_u128]).
    safe_batch: usize,
}

impl Field {
    pub fn new(p: u64) -> Self {
        assert!(p >= 2, "modulus must be at least 2");
        assert!(p < (1 << 63), "modulus must be below 2⁶³");
        let mu = u128::MAX / p as u128;
        let max_prod = (p as u128 - 1) * (p as u128 - 1);
        // How many such products fit under 2¹²⁸.
        let safe_batch = (u128::MAX / max_prod).min(usize::MAX as u128) as usize;
        Self {
            p,
            mu,
            safe_batch: safe_batch.max(1),
        }
    }

    /// Reduces an arbitrary 64-bit value into [0 . . *p*).
    #[inline(always)]
    pub fn reduce(&self, x: u64) -> u64 {
        self.reduce_u128(x as u128)
    }

    /// Reduces a 128-bit value into [0 . . *p*) by Barrett reduction.
    #[inline(always)]
    pub fn reduce_u128(&self, t: u128) -> u64 {
        let q = mul_hi_128(t, self.mu); // ⌊t·μ / 2¹²⁸⌋ ≤ ⌊t/p⌋, so q·p ≤ t
        let p = self.p as u128;
        let mut r = t - q * p; // r = (t mod p) + (⌊t/p⌋ − q)·p < 3p
        while r >= p {
            r -= p;
        }
        r as u64
    }

    #[inline(always)]
    fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p { s - self.p } else { s }
    }

    #[inline(always)]
    fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b { a - b } else { a + self.p - b }
    }

    #[inline(always)]
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.reduce_u128(a as u128 * b as u128)
    }

    fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1;
        while e != 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }

    /// Multiplicative inverse via Fermat: *aᵖ*⁻² mod *p* (*a* ≠ 0).
    #[inline(always)]
    fn inv(&self, a: u64) -> u64 {
        self.pow(a, self.p - 2)
    }
}

/// Rank over **F**_ₚ_ of an *n* × *n* matrix stored row-major in `m`,
/// by blocked Gaussian elimination.
///
/// This function uses the default Rayon thread pool for parallelization. To
/// customize the thread pool, you can a new thread pool and then call `install`
/// to run the rank computation within that pool, or set the `RAYON_NUM_THREADS`
/// environment variable to control the number of threads.
///
/// `pl` reports progress; pass `no_logging![]` to disable.
pub fn rank(field: &Field, m: &mut [u64], n: usize, pl: &mut impl ProgressLog) -> usize {
    debug_assert_eq!(m.len(), n * n);
    // Reduced pivot rows of the current panel (full width; only the trailing part
    // is used). Reused across panels.
    let mut u = vec![0; PANEL * n].into_boxed_slice();
    // Column-major repacking of the pivot rows for one trailing tile, so the
    // accumulation loop reads the pb factors contiguously.
    let mut utile = vec![0; PANEL * TILE].into_boxed_slice();
    // Pivot columns found in the current panel.
    let mut pcols: Vec<usize> = Vec::with_capacity(PANEL);

    let mut rank = 0;
    let mut col = 0;
    while col < n && rank < n {
        let jend = (col + PANEL).min(n);
        let p0 = rank; // first pivot row of this panel

        // Panel factorization over columns [col..jend) for rows >= rank.
        // Eliminate within the panel only; store each multiplier in place at its
        // pivot column (LAPACK-style L), deferring the trailing columns.
        pcols.clear();
        for j in col..jend {
            if rank == n {
                break;
            }
            let mut pr = None;
            for r in rank..n {
                if m[r * n + j] != 0 {
                    pr = Some(r);
                    break;
                }
            }
            let Some(pr) = pr else { continue };
            if pr != rank {
                // Swap whole rows (the stored multipliers travel with the row).
                let (lo, hi) = m.split_at_mut(pr * n); // rank < pr
                lo[rank * n..rank * n + n].swap_with_slice(&mut hi[0..n]);
            }
            let inv = field.inv(m[rank * n + j]);
            let (head, tail) = m.split_at_mut((rank + 1) * n);
            let prow = &head[rank * n..rank * n + n];
            tail.par_chunks_mut(n).for_each(|rr| {
                let e = rr[j];
                if e != 0 {
                    let f = field.mul(e, inv);
                    rr[j] = f; // store the multiplier instead of zero
                    for c in (j + 1)..jend {
                        rr[c] = field.sub(rr[c], field.mul(f, prow[c]));
                    }
                }
            });
            pcols.push(j);
            rank += 1;
        }

        let pb = pcols.len();
        if pb > 0 && jend < n {
            // Build U: the pivot rows' trailing parts, reduced among
            // themselves (the small triangular solve that couples them).
            for k in 0..pb {
                let src = (p0 + k) * n;
                u[k * n + jend..k * n + n].copy_from_slice(&m[src + jend..src + n]);
            }
            for k in 1..pb {
                for k2 in 0..k {
                    let f = m[(p0 + k) * n + pcols[k2]];
                    if f != 0 {
                        let (lo, hi) = u.split_at_mut(k * n);
                        let uk2 = &lo[k2 * n..k2 * n + n];
                        let uk = &mut hi[0..n];
                        for c in jend..n {
                            uk[c] = field.sub(uk[c], field.mul(f, uk2[c]));
                        }
                    }
                }
            }

            // Trailing update (rank-pb update) of the far rows [rank..n).
            // For each column-tile we pack the pb pivot rows column-major and
            // keep them cache-resident while every far row's tile is updated.
            // Each far-row element sums its pb products in a u128, reducing only
            // once per safe_batch products (a single reduction when safe_batch ≥
            // pb) — the blocked structure's real payoff.
            let sb = field.safe_batch;
            let mut c0 = jend;
            while c0 < n {
                let c1 = (c0 + TILE).min(n);
                let w = c1 - c0;
                // Pack U[*][c0..c1] column-major: utile[ci * pb + k] = U[k][c0+ci].
                for k in 0..pb {
                    let urow = &u[k * n + c0..k * n + c1];
                    for (ci, &val) in urow.iter().enumerate() {
                        utile[ci * pb + k] = val;
                    }
                }
                let utile = &utile[..w * pb];
                let pcols = &pcols[..];
                m[rank * n..].par_chunks_mut(n).for_each(|rr| {
                    let mut fs = [0; PANEL];
                    for k in 0..pb {
                        fs[k] = rr[pcols[k]];
                    }
                    for ci in 0..w {
                        let base = ci * pb;
                        let mut total = 0;
                        let mut acc = 0u128;
                        let mut cnt = 0;
                        for k in 0..pb {
                            acc += fs[k] as u128 * utile[base + k] as u128;
                            cnt += 1;
                            if cnt == sb {
                                total = field.add(total, field.reduce_u128(acc));
                                acc = 0;
                                cnt = 0;
                            }
                        }
                        if cnt > 0 {
                            total = field.add(total, field.reduce_u128(acc));
                        }
                        rr[c0 + ci] = field.sub(rr[c0 + ci], total);
                    }
                });
                c0 = c1;
            }
        }

        pl.update_with_count(jend - col);
        col = jend;
    }
    rank
}

/// Linear complexity of a sequence over **F**_ₚ_ (the length of the shortest
/// linear-feedback shift register that generates it) by the Berlekamp–Massey
/// algorithm. The connection-polynomial update touches only the deg(*b*) + 1
/// nonzero coefficients of the stored polynomial *b*, so the cost is O(*n* · *L*)
/// in the recovered complexity *L* rather than O(*n*²).
pub fn linear_complexity(field: &Field, s: &[u64], pl: &mut impl ProgressLog) -> usize {
    let n = s.len();
    let mut c = vec![0; n].into_boxed_slice(); // current connection polynomial, c[0] = 1
    let mut b = vec![0; n].into_boxed_slice(); // last connection polynomial before a length change
    let mut t = vec![0; n].into_boxed_slice(); // reusable scratch (avoids per-step allocation)
    c[0] = 1;
    b[0] = 1;
    let mut l = 0; // current linear complexity
    let mut bl = 1; // number of meaningful coefficients of b (deg b + 1)
    let mut m = 1; // steps since the last length change
    let mut bb_inv = 1; // inverse of the discrepancy at that last change (bb starts at 1)

    for i in 0..n {
        pl.update(); // once per outer step (n total), cheap relative to the inner loop
        // Discrepancy d = s[i] + Σ_{j=1}^{l} c[j]·s[i−j].
        let mut d = s[i];
        for j in 1..=l {
            d = field.add(d, field.mul(c[j], s[i - j]));
        }
        if d == 0 {
            m += 1;
            continue;
        }
        // c(x) ← c(x) − (d/bb)·xᵐ·b(x). Only b's deg(b) + 1 nonzero coefficients
        // contribute, so we update c[m..m + bl) instead of all of c.
        let coef = field.mul(d, bb_inv);
        let jmax = bl.min(n - m);
        if 2 * l <= i {
            // Length change: the pre-update c (degree ≤ l) becomes the new b, so
            // snapshot it before the update.
            let old_l = l;
            t[..old_l + 1].copy_from_slice(&c[..old_l + 1]);
            for j in 0..jmax {
                c[j + m] = field.sub(c[j + m], field.mul(coef, b[j]));
            }
            std::mem::swap(&mut b, &mut t); // b ← old c; t reused next time
            bl = old_l + 1;
            bb_inv = field.inv(d);
            l = i + 1 - old_l;
            m = 1;
        } else {
            for j in 0..jmax {
                c[j + m] = field.sub(c[j + m], field.mul(coef, b[j]));
            }
            m += 1;
        }
    }
    l
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn barrett_matches_division() {
        // The Barrett reduction must equal the % p for every modulus, over the
        // full u128 input range and the worst-case corners.
        fn mix(x: u64) -> u64 {
            let mut z = x.wrapping_mul(0x9E3779B97F4A7C15);
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z ^ (z >> 31)
        }
        let primes = [
            2,
            3,
            97,
            65537,
            (1 << 31) - 1,
            (1 << 61) - 1,
            2305843009213693951, // 2⁶¹ − 1 (MIXMAX)
            (1 << 62) + 135,
            (1 << 63) - 25,
        ];
        for &p in &primes {
            let f = Field::new(p);
            let pm = p as u128;
            // reduce_u128 over corners and a pseudo-random sweep up to ~2¹²⁸.
            let mut corners = vec![
                0u128,
                1,
                pm - 1,
                pm,
                pm + 1,
                u128::MAX,
                u128::MAX - 1,
                1u128 << 127,
                (1u128 << 127) | 1,
            ];
            for i in 0..256 {
                let hi = mix(i) as u128;
                let lo = mix(i ^ 0xABCD) as u128;
                corners.push((hi << 64) | lo);
            }
            for &t in &corners {
                assert_eq!(f.reduce_u128(t), (t % pm) as u64, "p={p} t={t}");
            }
            // reduce(x) for 64-bit inputs, and mul for products.
            for &x in &[0, 1, p - 1, p.wrapping_sub(0), u64::MAX, u64::MAX - 1] {
                assert_eq!(f.reduce(x), (x as u128 % pm) as u64, "p={p} x={x}");
            }
            for i in 0..256 {
                let a = mix(i) % p;
                let b = mix(i ^ 0x1234) % p;
                assert_eq!(f.mul(a, b), ((a as u128 * b as u128) % pm) as u64, "p={p}");
            }
        }
    }
}
