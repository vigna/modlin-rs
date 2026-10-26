/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Fast linear complexity over 𝐅*ₚ* by the half-GCD algorithm.
//!
//! The linear complexity of a sequence *s*₀, …, *s*ₙ₋₁ can be read off the
//! remainder sequence of the Euclidean algorithm on *x*ⁿ and *S*(*x*) =
//! Σ *sᵢ* *x*ⁿ⁻¹⁻ⁱ: if *Rⱼ*₋₁ and *Rⱼ* are the two consecutive remainders such
//! that 2 deg *Rⱼ* < *n* ≤ 2 deg *Rⱼ*₋₁, the linear complexity is *n* − deg
//! *Rⱼ*₋₁ if deg *Rⱼ* + deg *Rⱼ*₋₁ < *n*, and *n* − deg *Rⱼ* otherwise (the
//! zero polynomial has degree −∞). The half-GCD algorithm computes these two
//! remainders with O(*M*(*n*) log *n*) field operations, where *M*(*n*) is the
//! cost of multiplying two polynomials of degree *n*; with the number-theoretic
//! transforms of this crate, the total cost is O(*n* log² *n*).
//!
//! # Implementation Notes
//!
//! Polynomials are vectors of coefficients from the lowest degree, with no
//! trailing zeros (so the zero polynomial is the empty vector). The half-GCD
//! follows the formulation of Thull and Yap, “A unified approach to HGCD
//! algorithms for polynomials and integers” (1990), with a schoolbook
//! Euclidean base case.

use rayon::prelude::*;

use crate::fp::Field;
use crate::ntt::{Multiplier, Spectrum, trim};

/// The maximum number of threads used by [`linear_complexity`].
pub const MAX_THREADS: usize = 16;

/// The maximum length of a sequence accepted by [`linear_complexity`]: the
/// products of the half-GCD algorithm must fit in a transform of at most 2³²
/// points.
pub const MAX_LENGTH: usize = ((1u64 << 32) - 2) as usize;

/// Below this degree, the half-GCD runs the Euclidean algorithm directly.
const HGCD_BASE: usize = 64;

/// Below this length of the quotient (or of the divisor), division with
/// remainder uses the schoolbook algorithm.
const DIV_CLASSICAL: usize = 64;

/// From this length of all nonzero operands on, [`Ring::apply`] and
/// [`Ring::mat_mul`] transform each operand once and add the products in the
/// transform domain.
const SHARE_CUTOFF: usize = 256;

/// The length of the shortest nonzero polynomial among `polys`.
fn shortest(polys: &[&[u64]]) -> usize {
    polys
        .iter()
        .map(|p| p.len())
        .filter(|&l| l > 0)
        .min()
        .unwrap_or(0)
}

/// The maximum number of coefficients of a product computed by
/// [`linear_complexity`] on a sequence of length `n`.
///
/// # Implementation Notes
///
/// The longest products are those of the division with remainder at the top
/// level of the recursion: the product of the divisor and the quotient, and
/// the truncated products of the Newton iteration, have at most *n* + 2
/// coefficients. All other products are shorter.
const fn max_product_len(n: usize) -> usize {
    n + 2
}

/// Degree of a polynomial, with −1 for the zero polynomial.
#[inline(always)]
const fn deg(a: &[u64]) -> isize {
    a.len() as isize - 1
}

/// The quotient of `a` by *x*ᵏ.
#[inline(always)]
fn shift_right(a: &[u64], k: usize) -> &[u64] {
    if k >= a.len() { &[] } else { &a[k..] }
}

/// The remainder of `a` by *x*ᵏ, with no trailing zeros.
#[inline(always)]
fn low(a: &[u64], k: usize) -> &[u64] {
    let mut l = a.len().min(k);
    while l > 0 && a[l - 1] == 0 {
        l -= 1;
    }
    &a[..l]
}

/// `a` · *x*ᵏ.
fn shift_left(a: &[u64], k: usize) -> Vec<u64> {
    if a.is_empty() {
        return Vec::new();
    }
    let mut c = vec![0; k];
    c.extend_from_slice(a);
    c
}

/// A 2 × 2 polynomial matrix.
struct Mat([Vec<u64>; 4]);

impl Mat {
    fn identity() -> Self {
        Self([vec![1], Vec::new(), Vec::new(), vec![1]])
    }
}

/// Polynomial arithmetic over 𝐅*ₚ*.
struct Ring {
    mult: Multiplier,
}

impl Ring {
    const fn field(&self) -> &Field {
        self.mult.field()
    }

    fn add(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
        let mut c = long.to_vec();
        for (x, &y) in c.iter_mut().zip(short) {
            *x = self.field().add(*x, y);
        }
        trim(&mut c);
        c
    }

    fn sub(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        let mut c = a.to_vec();
        if c.len() < b.len() {
            c.resize(b.len(), 0);
        }
        for (x, &y) in c.iter_mut().zip(b) {
            *x = self.field().sub(*x, y);
        }
        trim(&mut c);
        c
    }

    fn mul(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        self.mult.mul(a, b)
    }

    /// `a` · `b` modulo *x*ᵏ, padded with zeros to length *k*.
    fn mul_low(&self, a: &[u64], b: &[u64], k: usize) -> Vec<u64> {
        let a = &a[..a.len().min(k)];
        let b = &b[..b.len().min(k)];
        let mut c = self.mul(a, b);
        c.resize(k, 0);
        c
    }

    /// The inverse of `f` modulo *x*ᵏ by Newton iteration; `f`(0) must be
    /// nonzero. The result is padded with zeros to length *k*.
    fn inv_series(&self, f: &[u64], k: usize) -> Vec<u64> {
        let mut g = vec![self.field().inv(f[0])];
        let mut l = 1;
        while l < k {
            let l2 = (2 * l).min(k);
            // e = f·g ≡ 1 mod xˡ, so 1 − e = −xˡ·h; then g ← g + xˡ·g·h.
            let e = self.mul_low(f, &g, l2);
            let h: Vec<u64> = e[l..].iter().map(|&x| self.field().sub(0, x)).collect();
            let mut gh = self.mul(&g, &h);
            gh.resize(l2 - l, 0);
            g.extend_from_slice(&gh);
            l = l2;
        }
        g
    }

    /// Schoolbook division with remainder of `a` by the nonzero `b`.
    fn divrem_classical(&self, a: &[u64], b: &[u64]) -> (Vec<u64>, Vec<u64>) {
        debug_assert!(!b.is_empty(), "division by zero");
        let f = self.field();
        let db = b.len() - 1;
        if a.len() < b.len() {
            return (Vec::new(), a.to_vec());
        }
        let mut r = a.to_vec();
        let mut q = vec![0; a.len() - db];
        let lead_inv = f.inv(b[db]);
        for i in (0..q.len()).rev() {
            let c = f.mul(r[i + db], lead_inv);
            q[i] = c;
            if c != 0 {
                for (x, &y) in r[i..i + db].iter_mut().zip(&b[..db]) {
                    *x = f.sub(*x, f.mul(c, y));
                }
            }
        }
        r.truncate(db);
        trim(&mut r);
        trim(&mut q);
        (q, r)
    }

    /// Division with remainder of `a` by the nonzero `b`.
    fn divrem(&self, a: &[u64], b: &[u64]) -> (Vec<u64>, Vec<u64>) {
        debug_assert!(!b.is_empty(), "division by zero");
        if a.len() < b.len() {
            return (Vec::new(), a.to_vec());
        }
        let qlen = a.len() - b.len() + 1;
        if qlen < DIV_CLASSICAL || b.len() < DIV_CLASSICAL {
            return self.divrem_classical(a, b);
        }
        // The reversal of the quotient is the reversal of a divided by the
        // reversal of b, modulo x^qlen.
        let rev_a: Vec<u64> = a.iter().rev().take(qlen).copied().collect();
        let rev_b: Vec<u64> = b.iter().rev().take(qlen).copied().collect();
        let rev_q = self.mul_low(&rev_a, &self.inv_series(&rev_b, qlen), qlen);
        let mut q: Vec<u64> = rev_q.into_iter().rev().collect();
        trim(&mut q);
        let mut r = self.sub(a, &self.mul(b, &q));
        r.truncate(b.len() - 1);
        trim(&mut r);
        (q, r)
    }

    /// The transforms of `polys`, all of length `size`.
    fn transforms(&self, polys: &[&[u64]], size: usize) -> Vec<Spectrum> {
        polys
            .par_iter()
            .map(|p| self.mult.transform(p, size))
            .collect()
    }

    /// (`m`₀`a` + `m`₁`b`, `m`₂`a` + `m`₃`b`).
    fn apply(&self, m: &Mat, a: &[u64], b: &[u64]) -> (Vec<u64>, Vec<u64>) {
        let [m0, m1, m2, m3] = &m.0;
        if shortest(&[m0, m1, m2, m3, a, b]) >= SHARE_CUTOFF {
            let len = (m0.len().max(m2.len()) + a.len()).max(m1.len().max(m3.len()) + b.len()) - 1;
            let t = self.transforms(&[a, b, m0, m1, m2, m3], self.mult.size(len));
            return rayon::join(
                || self.mult.dot(&[(&t[2], &t[0]), (&t[3], &t[1])]),
                || self.mult.dot(&[(&t[4], &t[0]), (&t[5], &t[1])]),
            );
        }
        let ((x0, x1), (x2, x3)) = rayon::join(
            || rayon::join(|| self.mul(m0, a), || self.mul(m1, b)),
            || rayon::join(|| self.mul(m2, a), || self.mul(m3, b)),
        );
        (self.add(&x0, &x1), self.add(&x2, &x3))
    }

    /// The product `s` · `r`.
    fn mat_mul(&self, s: &Mat, r: &Mat) -> Mat {
        let [s0, s1, s2, s3] = &s.0;
        let [r0, r1, r2, r3] = &r.0;
        if shortest(&[s0, s1, s2, s3, r0, r1, r2, r3]) >= SHARE_CUTOFF {
            let ls = s0.len().max(s1.len()).max(s2.len()).max(s3.len());
            let lr = r0.len().max(r1.len()).max(r2.len()).max(r3.len());
            let t = self.transforms(
                &[s0, s1, s2, s3, r0, r1, r2, r3],
                self.mult.size(ls + lr - 1),
            );
            let ((a, b), (c, d)) = rayon::join(
                || {
                    rayon::join(
                        || self.mult.dot(&[(&t[0], &t[4]), (&t[1], &t[6])]),
                        || self.mult.dot(&[(&t[0], &t[5]), (&t[1], &t[7])]),
                    )
                },
                || {
                    rayon::join(
                        || self.mult.dot(&[(&t[2], &t[4]), (&t[3], &t[6])]),
                        || self.mult.dot(&[(&t[2], &t[5]), (&t[3], &t[7])]),
                    )
                },
            );
            return Mat([a, b, c, d]);
        }
        let (((a0, a1), (b0, b1)), ((c0, c1), (d0, d1))) = rayon::join(
            || {
                rayon::join(
                    || rayon::join(|| self.mul(s0, r0), || self.mul(s1, r2)),
                    || rayon::join(|| self.mul(s0, r1), || self.mul(s1, r3)),
                )
            },
            || {
                rayon::join(
                    || rayon::join(|| self.mul(s2, r0), || self.mul(s3, r2)),
                    || rayon::join(|| self.mul(s2, r1), || self.mul(s3, r3)),
                )
            },
        );
        Mat([
            self.add(&a0, &a1),
            self.add(&b0, &b1),
            self.add(&c0, &c1),
            self.add(&d0, &d1),
        ])
    }

    /// Left-multiplies `r` by the matrix [[0, 1], [1, −`q`]] of one Euclidean
    /// step.
    fn step(&self, r: Mat, q: &[u64]) -> Mat {
        let [r0, r1, r2, r3] = r.0;
        // Parallelism pays only for long products (e.g., not in the base case).
        let (t0, t1) = if r2.len().max(r3.len()) < SHARE_CUTOFF {
            (self.mul(q, &r2), self.mul(q, &r3))
        } else {
            rayon::join(|| self.mul(q, &r2), || self.mul(q, &r3))
        };
        let n2 = self.sub(&r0, &t0);
        let n3 = self.sub(&r1, &t1);
        Mat([r2, r3, n2, n3])
    }

    /// The half-GCD of `a` and `b`, where deg `a` > deg `b`: a matrix *M*,
    /// product of Euclidean steps, and the pair (*c*, *d*) = *M* (`a`, `b`),
    /// which satisfies deg *c* ≥ ⌈deg `a` / 2⌉ > deg *d*.
    ///
    /// # Implementation Notes
    ///
    /// Writing `a` = *a*ₕ *x*ᵐ + *a*ₗ and `b` = *b*ₕ *x*ᵐ + *b*ₗ, a recursive
    /// call on (*a*ₕ, *b*ₕ) returns *R* and *R* (*a*ₕ, *b*ₕ), so *R* (`a`, `b`)
    /// = *x*ᵐ *R* (*a*ₕ, *b*ₕ) + *R* (*a*ₗ, *b*ₗ) requires only products by
    /// the low parts; the same holds for the second recursive call.
    fn hgcd(&self, a: &[u64], b: &[u64]) -> (Mat, Vec<u64>, Vec<u64>) {
        let n = a.len() - 1;
        let m = n.div_ceil(2);
        if deg(b) < m as isize {
            return (Mat::identity(), a.to_vec(), b.to_vec());
        }
        if n < HGCD_BASE {
            return self.hgcd_base(a, b, m);
        }
        let (r, a1, b1) = self.hgcd_shifted(a, b, m);
        if deg(&b1) < m as isize {
            return (r, a1, b1);
        }
        let (q, d) = self.divrem(&a1, &b1);
        let r = self.step(r, &q);
        if deg(&d) < m as isize {
            return (r, b1, d);
        }
        let k = 2 * m - (b1.len() - 1);
        let (s, c, d) = self.hgcd_shifted(&b1, &d, k);
        (self.mat_mul(&s, &r), c, d)
    }

    /// The half-GCD *M* of the quotients of `a` and `b` by *x*ᵏ, with *M*
    /// (`a`, `b`).
    fn hgcd_shifted(&self, a: &[u64], b: &[u64], k: usize) -> (Mat, Vec<u64>, Vec<u64>) {
        let (r, ah, bh) = self.hgcd(shift_right(a, k), shift_right(b, k));
        let (al, bl) = self.apply(&r, low(a, k), low(b, k));
        (
            r,
            self.add(&shift_left(&ah, k), &al),
            self.add(&shift_left(&bh, k), &bl),
        )
    }

    /// The half-GCD by the Euclidean algorithm, stopping at the first
    /// remainder of degree less than `m`.
    fn hgcd_base(&self, a: &[u64], b: &[u64], m: usize) -> (Mat, Vec<u64>, Vec<u64>) {
        let mut r = Mat::identity();
        let (mut a, mut b) = (a.to_vec(), b.to_vec());
        while deg(&b) >= m as isize {
            let (q, rem) = self.divrem_classical(&a, &b);
            r = self.step(r, &q);
            a = b;
            b = rem;
        }
        (r, a, b)
    }
}

/// An estimate of the peak memory, in bytes, used by [`linear_complexity`] on
/// a sequence of length `n`: the transform tables (two words per point for
/// each of the three primes), the six transforms of the largest matrix-vector
/// product, and the polynomials of the recursion.
pub const fn memory_bytes(n: usize) -> usize {
    let tables = 6 * max_product_len(n).next_power_of_two();
    let spectra = 18 * (3 * n / 4 + 1).next_power_of_two();
    (tables + spectra + 24 * n) * size_of::<u64>()
}

/// Linear complexity of a sequence over 𝐅*ₚ* by the half-GCD algorithm,
/// which requires O(*n* log² *n*) field operations.
///
/// The result is the same as that of [`berlekamp_massey`], which is faster
/// for short sequences; [`linear_complexity`] chooses between the two.
///
/// This function uses the current Rayon thread pool for parallelization, but
/// at most [`MAX_THREADS`] threads, as beyond that the overhead of the
/// parallel tasks outweighs their benefits.
///
/// # Panics
///
/// If the sequence is longer than [`MAX_LENGTH`].
///
/// [`berlekamp_massey`]: crate::fp::berlekamp_massey
/// [`linear_complexity`]: crate::fp::linear_complexity
pub fn linear_complexity(field: &Field, s: &[u64]) -> usize {
    assert!(s.len() <= MAX_LENGTH, "sequence longer than {MAX_LENGTH}");
    if rayon::current_num_threads() > MAX_THREADS {
        if let Ok(pool) = rayon::ThreadPoolBuilder::new()
            .num_threads(MAX_THREADS)
            .build()
        {
            return pool.install(|| linear_complexity_in_pool(field, s));
        }
    }
    linear_complexity_in_pool(field, s)
}

/// [`linear_complexity`] in the current thread pool.
fn linear_complexity_in_pool(field: &Field, s: &[u64]) -> usize {
    let n = s.len();
    // S(x) = Σ sᵢ xⁿ⁻¹⁻ⁱ.
    let mut b: Vec<u64> = s.iter().rev().copied().collect();
    trim(&mut b);
    if b.is_empty() {
        return 0;
    }
    let mut a = vec![0; n + 1];
    a[n] = 1;
    let ring = Ring {
        mult: Multiplier::new(*field, max_product_len(n)),
    };
    let (_, c, d) = ring.hgcd(&a, &b);
    let (dc, dd) = (deg(&c), deg(&d));
    debug_assert!(2 * dd < n as isize && n as isize <= 2 * dc);
    if dd + dc < n as isize {
        n - dc as usize
    } else {
        n - dd as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::{mix, random_poly};

    fn ring(p: u64, max_len: usize) -> Ring {
        Ring {
            mult: Multiplier::new(Field::new(p), max_len),
        }
    }

    #[test]
    fn test_divrem_matches_schoolbook() {
        for p in [3, (1 << 61) - 1, 9223372036854775783] {
            let r = ring(p, 1 << 14);
            for (la, lb) in [(10, 3), (200, 100), (1000, 300), (3000, 2900), (4000, 70)] {
                let a = random_poly(la, p, la as u64);
                let b = random_poly(lb, p, 7 * lb as u64);
                assert_eq!(
                    r.divrem(&a, &b),
                    r.divrem_classical(&a, &b),
                    "p={p} {la}/{lb}"
                );
            }
        }
    }

    /// The remainders of the Euclidean algorithm straddling ⌈deg a / 2⌉.
    fn euclid_straddle(r: &Ring, a: &[u64], b: &[u64]) -> (Vec<u64>, Vec<u64>) {
        let m = (a.len() - 1).div_ceil(2) as isize;
        let (mut a, mut b) = (a.to_vec(), b.to_vec());
        while deg(&b) >= m {
            let (_, rem) = r.divrem_classical(&a, &b);
            a = b;
            b = rem;
        }
        (a, b)
    }

    #[test]
    fn test_hgcd_matches_euclid() {
        for p in [2, 3, 7, (1 << 61) - 1] {
            let r = ring(p, 1 << 14);
            for (la, lb) in [(65, 64), (130, 129), (300, 250), (1000, 999), (2049, 1500)] {
                let a = random_poly(la, p, 3 * la as u64);
                let b = random_poly(lb, p, 5 * lb as u64);
                let (m, c, d) = r.hgcd(&a, &b);
                // The returned pair must be the image of the input.
                assert_eq!(r.apply(&m, &a, &b), (c.clone(), d.clone()));
                let (ec, ed) = euclid_straddle(&r, &a, &b);
                // Remainders are determined up to a constant factor.
                assert_eq!(deg(&c), deg(&ec), "p={p} {la},{lb}");
                assert_eq!(deg(&d), deg(&ed), "p={p} {la},{lb}");
                if !ec.is_empty() {
                    let k = r
                        .field()
                        .mul(c[c.len() - 1], r.field().inv(ec[ec.len() - 1]));
                    let scaled: Vec<u64> = ec.iter().map(|&x| r.field().mul(x, k)).collect();
                    assert_eq!(c, scaled, "p={p} {la},{lb}");
                }
            }
        }
    }

    #[test]
    fn test_thread_cap() -> anyhow::Result<()> {
        // Run in a pool with more than MAX_THREADS threads, so that
        // linear_complexity builds its own pool.
        let p = (1 << 61) - 1;
        let field = Field::new(p);
        let s: Vec<u64> = (0..5000).map(|i| mix(i) % p).collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(MAX_THREADS + 1)
            .build()?;
        assert_eq!(pool.install(|| linear_complexity(&field, &s)), 2500);
        Ok(())
    }
}
