/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Multiplication of polynomials over 𝐅*ₚ*, for any prime *p* < 2⁶³, by
//! number-theoretic transforms modulo three primes below 2⁶² and the Chinese
//! remainder theorem.
//!
//! # Implementation Notes
//!
//! Transforms have at most 2³² points, so a polynomial with coefficients in
//! [0 . . *p*) appearing in a product has fewer than 2³² coefficients, and a
//! sum of at most three products (see [`Multiplier::dot`]) has integer
//! coefficients smaller than 3 · 2³² · 2¹²⁶ < 2¹⁶⁰, which is less than the
//! product of the three primes, so the coefficients modulo the three primes
//! determine the integer coefficients, and thus the coefficients modulo *p*.
//!
//! The forward transform is a decimation-in-frequency transform that maps
//! natural order to bit-reversed order, and the inverse transform a
//! decimation-in-time transform that maps bit-reversed order back to natural
//! order, so no permutation is ever needed. Butterflies use Shoup's
//! precomputed-quotient multiplication and Harvey's lazy reduction, keeping
//! values in [0 . . 2*m*), which is possible because 4*m* < 2⁶⁴.

use rayon::prelude::*;

use crate::fp::Field;

/// Primes *m* = *c* · 2³² + 1 < 2⁶², with primitive 2³²-th roots of unity
/// modulo them.
const PRIMES: [(u64, u64); 3] = [
    (4611685941117976577, 69433692538710738),
    (4611685692009873409, 458164920477615602),
    (4611685606110527489, 3318345213167893729),
];

/// The base-2 logarithm of the order of the roots of unity in [`PRIMES`].
const ROOT_LOG: u32 = 32;

/// Below this length of the shorter operand, [`Multiplier::mul`] uses
/// schoolbook multiplication.
const CLASSICAL_CUTOFF: usize = 48;

/// Up to this length, transforms perform their stages iteratively, as the data
/// fit in cache; longer transforms recurse on the two halves.
const ITERATIVE_CUTOFF: usize = 1 << 12;

/// From this length on, the halves of a recursive transform, and the
/// butterflies of its outer stage, are processed in parallel.
const PARALLEL_CUTOFF: usize = 1 << 15;

/// The number of butterflies of a stage processed by a parallel task.
const STAGE_CHUNK: usize = 1 << 12;

/// The number of coefficients reconstructed by a parallel task of the Chinese
/// remainder theorem.
const CRT_CHUNK: usize = 1 << 12;

/// `x` − `m` if `x` ≥ `m`, and `x` otherwise.
#[inline(always)]
const fn csub(x: u64, m: u64) -> u64 {
    if x >= m { x - m } else { x }
}

/// Shoup's precomputed quotient ⌊*w* · 2⁶⁴ / *m*⌋.
#[inline(always)]
const fn shoup(w: u64, m: u64) -> u64 {
    (((w as u128) << 64) / m as u128) as u64
}

/// *x* · *w* modulo *m* in [0 . . 2*m*), for any *x*, given *w* < *m* and its
/// Shoup quotient `wp`.
#[inline(always)]
const fn mul_shoup(x: u64, w: u64, wp: u64, m: u64) -> u64 {
    let q = ((x as u128 * wp as u128) >> 64) as u64;
    x.wrapping_mul(w).wrapping_sub(q.wrapping_mul(m))
}

/// Applies `butterflies` to the halves `x` and `y` of a block, passing the
/// offset of the first butterfly, in parallel chunks for long blocks.
#[inline(always)]
fn for_each_chunk(
    x: &mut [u64],
    y: &mut [u64],
    butterflies: impl Fn(&mut [u64], &mut [u64], usize) + Send + Sync,
) {
    if 2 * x.len() >= PARALLEL_CUTOFF {
        x.par_chunks_mut(STAGE_CHUNK)
            .zip(y.par_chunks_mut(STAGE_CHUNK))
            .enumerate()
            .for_each(|(c, (x, y))| butterflies(x, y, c * STAGE_CHUNK));
    } else {
        butterflies(x, y, 0);
    }
}

/// Transform tables for one prime.
struct Tables {
    /// Arithmetic modulo the prime *m*.
    field: Field,
    /// For each power of two *h* less than the length of this vector, the
    /// twiddles *W*₂*ₕʲ*, *j* < *h*, at positions [*h* . . 2*h*), where
    /// *W*₂*ₕ* is a primitive 2*h*-th root of unity, so that each stage of a
    /// transform reads its twiddles sequentially.
    tw: Vec<u64>,
    /// The Shoup quotients of the entries of `tw`.
    tw_p: Vec<u64>,
}

impl Tables {
    fn new(m: u64, root: u64, log_max: u32) -> Self {
        let field = Field::new(m);
        let size = 1usize << log_max;
        let half = size / 2;
        let mut tw = vec![0; size];
        let mut tw_p = vec![0; size];
        let w = field.pow(root, 1 << (ROOT_LOG - log_max));
        let mut x = 1;
        for j in 0..half {
            tw[half + j] = x;
            tw_p[half + j] = shoup(x, m);
            x = field.mul(x, w);
        }
        // W₂ₕʲ = W₄ₕ²ʲ, so each stage is the even half of the next one.
        let mut h = half / 2;
        while h >= 1 {
            for j in 0..h {
                tw[h + j] = tw[2 * h + 2 * j];
                tw_p[h + j] = tw_p[2 * h + 2 * j];
            }
            h /= 2;
        }
        Self { field, tw, tw_p }
    }

    /// Forward transform of `a` (whose length is a power of two not larger
    /// than the length of the twiddle tables), from natural to bit-reversed
    /// order; entries must be in [0 . . 2*m*), and stay so.
    ///
    /// After the outer stage, the two halves are independent transforms, so
    /// the transform recurses on them: this way, the stages of a half that
    /// fits in cache are performed in cache.
    fn forward(&self, a: &mut [u64]) {
        let n = a.len();
        if n <= ITERATIVE_CUTOFF {
            return self.forward_iterative(a);
        }
        let (x, y) = a.split_at_mut(n / 2);
        self.forward_stage(x, y);
        if n >= PARALLEL_CUTOFF {
            rayon::join(|| self.forward(x), || self.forward(y));
        } else {
            self.forward(x);
            self.forward(y);
        }
    }

    /// The butterflies of a forward stage between the halves `x` and `y` of
    /// a block.
    fn forward_stage(&self, x: &mut [u64], y: &mut [u64]) {
        let m = self.field.p();
        let m2 = 2 * m;
        let h = x.len();
        let butterflies = |x: &mut [u64], y: &mut [u64], j0: usize| {
            let w = &self.tw[h + j0..h + j0 + x.len()];
            let wp = &self.tw_p[h + j0..h + j0 + x.len()];
            for ((u, v), (&w, &wp)) in x.iter_mut().zip(y.iter_mut()).zip(w.iter().zip(wp)) {
                let s = *u + *v;
                let d = *u + m2 - *v;
                *u = csub(s, m2);
                *v = mul_shoup(d, w, wp, m);
            }
        };
        for_each_chunk(x, y, butterflies);
    }

    /// Iterative forward transform, for short inputs.
    fn forward_iterative(&self, a: &mut [u64]) {
        let mut h = a.len() / 2;
        while h >= 1 {
            for block in a.chunks_exact_mut(2 * h) {
                let (x, y) = block.split_at_mut(h);
                self.forward_stage(x, y);
            }
            h /= 2;
        }
    }

    /// Inverse transform of `a`, from bit-reversed to natural order, without
    /// the division by the length; entries must be in [0 . . 2*m*), and stay
    /// so. As for [`Self::forward`], the transform recurses on the two
    /// halves, and then performs the outer stage.
    fn inverse(&self, a: &mut [u64]) {
        let n = a.len();
        if n <= ITERATIVE_CUTOFF {
            return self.inverse_iterative(a);
        }
        let (x, y) = a.split_at_mut(n / 2);
        if n >= PARALLEL_CUTOFF {
            rayon::join(|| self.inverse(x), || self.inverse(y));
        } else {
            self.inverse(x);
            self.inverse(y);
        }
        self.inverse_stage(x, y);
    }

    /// The butterflies of an inverse stage between the halves `x` and `y` of
    /// a block.
    ///
    /// The inverse twiddles are derived from the forward ones, as *W*₂*ₕ*⁻*ʲ* =
    /// −*W*₂*ₕʰ*⁻*ʲ*, and the Shoup quotient of *m* − *w* is the complement of
    /// that of *w*.
    fn inverse_stage(&self, x: &mut [u64], y: &mut [u64]) {
        let m = self.field.p();
        let m2 = 2 * m;
        let h = x.len();
        let butterfly = |u: &mut u64, v: &mut u64, t: u64| {
            let s = *u + t;
            let d = *u + m2 - t;
            *u = csub(s, m2);
            *v = csub(d, m2);
        };
        let butterflies = |x: &mut [u64], y: &mut [u64], j0: usize| {
            let mut start = 0;
            if j0 == 0 {
                // The twiddle of index 0 is 1.
                let t = y[0];
                butterfly(&mut x[0], &mut y[0], t);
                start = 1;
            }
            // Indices 2h − j for j ∈ [j0 + start..j0 + len).
            let hi = 2 * h - j0 - start;
            let lo = 2 * h + 1 - j0 - x.len();
            let w = self.tw[lo..=hi].iter().rev();
            let wp = self.tw_p[lo..=hi].iter().rev();
            for ((u, v), (&w, &wp)) in x[start..]
                .iter_mut()
                .zip(y[start..].iter_mut())
                .zip(w.zip(wp))
            {
                let t = mul_shoup(*v, m - w, !wp, m);
                butterfly(u, v, t);
            }
        };
        for_each_chunk(x, y, butterflies);
    }

    /// Iterative inverse transform, for short inputs.
    fn inverse_iterative(&self, a: &mut [u64]) {
        let mut h = 1;
        while h < a.len() {
            for block in a.chunks_exact_mut(2 * h) {
                let (x, y) = block.split_at_mut(h);
                self.inverse_stage(x, y);
            }
            h *= 2;
        }
    }

    /// The forward transform of `a`, padded with zeros to the given
    /// power-of-two length, with entries in [0 . . 2*m*).
    fn transform(&self, a: &[u64], len: usize) -> Vec<u64> {
        let m2 = 2 * self.field.p();
        // Field elements are below 2⁶³, so one conditional subtraction of 2m
        // brings them into [0..2m).
        let mut v = vec![0; len];
        for (d, &s) in v.iter_mut().zip(a) {
            *d = csub(s, m2);
        }
        self.forward(&mut v);
        v
    }

    /// Σ *xᵢ* · *yᵢ*, pointwise, for at most three pairs of transforms with
    /// entries in [0 . . 2*m*); the result has entries in [0 . . *m*).
    fn dot(&self, terms: &[(&[u64], &[u64])]) -> Vec<u64> {
        // Each product is below 4m² < 2¹²⁶, so three of them fit in a u128.
        debug_assert!((1..=3).contains(&terms.len()));
        let len = terms[0].0.len();
        (0..len)
            .map(|i| {
                let acc: u128 = terms.iter().map(|(x, y)| x[i] as u128 * y[i] as u128).sum();
                self.field.reduce_u128(acc)
            })
            .collect()
    }

    /// The inverse transform of `a`, divided by the length, with entries in
    /// [0 . . *m*).
    fn untransform(&self, mut a: Vec<u64>) -> Vec<u64> {
        let m = self.field.p();
        self.inverse(&mut a);
        let len_inv = self.field.inv(a.len() as u64 % m);
        let len_inv_p = shoup(len_inv, m);
        for x in a.iter_mut() {
            *x = csub(mul_shoup(*x, len_inv, len_inv_p, m), m);
        }
        a
    }
}

/// The transforms of a polynomial modulo the three primes.
pub(crate) struct Spectrum([Vec<u64>; 3]);

/// Multiplier of polynomials over 𝐅*ₚ*, for products with at most the number
/// of coefficients given at construction.
pub(crate) struct Multiplier {
    field: Field,
    /// The maximum transform length.
    max_len: usize,
    tables: [Tables; 3],
    /// *m*₁⁻¹ mod *m*₂.
    c12: u64,
    /// *m*₁ mod *m*₃.
    m1_mod_m3: u64,
    /// (*m*₁*m*₂)⁻¹ mod *m*₃.
    c123: u64,
    /// *m*₁ mod *p*.
    m1_mod_p: u64,
    /// *m*₁*m*₂ mod *p*.
    m12_mod_p: u64,
}

impl Multiplier {
    /// Creates a multiplier over `field` for products with at most `max_len`
    /// coefficients.
    pub(crate) fn new(field: Field, max_len: usize) -> Self {
        let log_max = max_len.next_power_of_two().trailing_zeros().max(1);
        assert!(
            log_max <= ROOT_LOG,
            "the products must have at most 2³² coefficients"
        );
        let tables = PRIMES.map(|(m, root)| Tables::new(m, root, log_max));
        let [(m1, _), (m2, _), (m3, _)] = PRIMES;
        let (f2, f3) = (&tables[1].field, &tables[2].field);
        let p = field.p();
        let m1_mod_p = m1 % p;
        Self {
            field,
            max_len: 1 << log_max,
            c12: f2.inv(m1 % m2),
            m1_mod_m3: m1 % m3,
            c123: f3.inv(f3.mul(m1 % m3, m2 % m3)),
            m1_mod_p,
            m12_mod_p: field.mul(m1_mod_p, m2 % p),
            tables,
        }
    }

    /// The transform length for products with at most `len` coefficients.
    pub(crate) fn size(&self, len: usize) -> usize {
        let size = len.next_power_of_two();
        assert!(
            size <= self.max_len,
            "the product must fit the transforms of this multiplier"
        );
        size
    }

    /// The transforms of `a` of length `size`, as returned by [`Self::size`].
    pub(crate) fn transform(&self, a: &[u64], size: usize) -> Spectrum {
        let t: Vec<Vec<u64>> = self
            .tables
            .par_iter()
            .map(|t| t.transform(a, size))
            .collect();
        Spectrum(t.try_into().unwrap())
    }

    /// The polynomial Σ *xᵢ* · *yᵢ*, for at most three pairs of transforms of
    /// the same length, with no trailing zeros. The length of the transforms
    /// must be at least the number of coefficients of the result.
    pub(crate) fn dot(&self, terms: &[(&Spectrum, &Spectrum)]) -> Vec<u64> {
        let residues: Vec<Vec<u64>> = (0..3usize)
            .into_par_iter()
            .map(|k| {
                let pairs: Vec<(&[u64], &[u64])> = terms
                    .iter()
                    .map(|(x, y)| (x.0[k].as_slice(), y.0[k].as_slice()))
                    .collect();
                let t = &self.tables[k];
                t.untransform(t.dot(&pairs))
            })
            .collect();
        self.crt(&residues[0], &residues[1], &residues[2])
    }

    /// The field of the coefficients.
    pub(crate) const fn field(&self) -> &Field {
        &self.field
    }

    /// The product of `a` and `b` (coefficients from the lowest degree,
    /// possibly with trailing zeros), with no trailing zeros.
    pub(crate) fn mul(&self, a: &[u64], b: &[u64]) -> Vec<u64> {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        if a.len().min(b.len()) < CLASSICAL_CUTOFF {
            return mul_classical(&self.field, a, b);
        }
        let size = self.size(a.len() + b.len() - 1);
        let (ta, tb) = rayon::join(|| self.transform(a, size), || self.transform(b, size));
        self.dot(&[(&ta, &tb)])
    }

    /// The polynomial with the given residues modulo the three primes, with
    /// no trailing zeros.
    fn crt(&self, r1: &[u64], r2: &[u64], r3: &[u64]) -> Vec<u64> {
        let [_, (m2, _), (m3, _)] = PRIMES;
        let (f2, f3) = (&self.tables[1].field, &self.tables[2].field);
        let mut c = vec![0; r1.len()];
        // Garner's algorithm: the integer coefficient is t₁ + m₁t₂ + m₁m₂t₃,
        // with t₁ < m₁, t₂ < m₂, and t₃ < m₃.
        c.par_chunks_mut(CRT_CHUNK)
            .enumerate()
            .for_each(|(chunk, out)| {
                let base = chunk * CRT_CHUNK;
                for (i, x) in out.iter_mut().enumerate() {
                    let k = base + i;
                    let t1 = r1[k];
                    let t1_m2 = csub(t1, m2);
                    let t2 = f2.mul(f2.sub(r2[k], t1_m2), self.c12);
                    let t1_m3 = csub(t1, m3);
                    let t2_m3 = csub(t2, m3);
                    let u = f3.add(t1_m3, f3.mul(self.m1_mod_m3, t2_m3));
                    let t3 = f3.mul(f3.sub(r3[k], u), self.c123);
                    *x = self.field.reduce_u128(
                        t1 as u128
                            + self.m1_mod_p as u128 * t2 as u128
                            + self.m12_mod_p as u128 * t3 as u128,
                    );
                }
            });
        trim(&mut c);
        c
    }
}

/// Schoolbook product of the nonempty `a` and `b` (possibly with trailing
/// zeros), with no trailing zeros.
fn mul_classical(field: &Field, a: &[u64], b: &[u64]) -> Vec<u64> {
    debug_assert!(
        !a.is_empty() && !b.is_empty(),
        "the factors must be nonempty"
    );
    // With b reversed, each coefficient is a dot product of contiguous slices,
    // as b[k − i] = rb[lb − 1 − k + i].
    let lb = b.len();
    let rb: Vec<u64> = b.iter().rev().copied().collect();
    let mut c: Vec<u64> = (0..a.len() + lb - 1)
        .map(|k| {
            let lo = k.saturating_sub(lb - 1);
            let hi = k.min(a.len() - 1);
            let start = lb - 1 + lo - k;
            field.dot(&a[lo..=hi], &rb[start..=start + hi - lo])
        })
        .collect();
    trim(&mut c);
    c
}

/// Removes trailing zeros.
#[inline(always)]
pub(crate) fn trim(a: &mut Vec<u64>) {
    while a.last() == Some(&0) {
        a.pop();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::{mix, random_poly};

    #[test]
    fn test_primes_and_roots() {
        for (m, root) in PRIMES {
            let f = Field::new(m);
            assert!(m < 1 << 62);
            assert_eq!(f.pow(root, 1 << (ROOT_LOG - 1)), m - 1);
        }
    }

    #[test]
    fn test_mul_matches_schoolbook() {
        for p in [
            2,
            3,
            65537,
            (1 << 31) - 1,
            (1 << 61) - 1,
            9223372036854775783,
        ] {
            let field = Field::new(p);
            let mult = Multiplier::new(field, 1 << 12);
            for (la, lb) in [
                (1, 1),
                (47, 48),
                (48, 48),
                (100, 300),
                (1000, 1000),
                (1500, 2500),
            ] {
                let a = random_poly(la, p, 1 + la as u64);
                let b = random_poly(lb, p, 2 + lb as u64);
                assert_eq!(
                    mult.mul(&a, &b),
                    mul_classical(&field, &a, &b),
                    "p={p} {la}x{lb}"
                );
            }
        }
    }

    #[test]
    fn test_mul_extreme_coefficients() {
        // All coefficients equal to p − 1 maximize the integer coefficients
        // that the Chinese remainder theorem must recover.
        let p = 9223372036854775783; // largest prime below 2⁶³
        let field = Field::new(p);
        let mult = Multiplier::new(field, 1 << 13);
        let a = vec![p - 1; 4000];
        assert_eq!(mult.mul(&a, &a), mul_classical(&field, &a, &a));
    }

    #[test]
    fn test_mul_parallel_by_evaluation() {
        // A product long enough to use the parallel branches of the transforms
        // (at least PARALLEL_CUTOFF points, several stage chunks) and of the
        // Chinese remainder theorem, checked by evaluation at random points.
        let p = 9223372036854775783;
        let field = Field::new(p);
        let (la, lb) = (10_000, 9_000);
        assert!(Multiplier::new(field, la + lb).size(la + lb - 1) >= PARALLEL_CUTOFF);
        let mult = Multiplier::new(field, la + lb);
        let a = random_poly(la, p, 1);
        let b = random_poly(lb, p, 2);
        let c = mult.mul(&a, &b);
        assert_eq!(c.len(), la + lb - 1);
        let eval = |poly: &[u64], x: u64| {
            poly.iter()
                .rev()
                .fold(0, |acc, &coef| field.add(field.mul(acc, x), coef))
        };
        for t in 0..10 {
            let x = mix(t) % p;
            assert_eq!(eval(&c, x), field.mul(eval(&a, x), eval(&b, x)));
        }
    }
}
