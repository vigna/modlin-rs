/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

//! Pseudorandom number generators selected at build time via Cargo features.
//!
//! Exactly one feature must be enabled when building the crate; each variant
//! exposes a single `Prng` type with a `new(seed: u64) -> Self` constructor and
//! a `next_u64(&mut self) -> u64` step function.
//!
//! The test reduces each output modulo the field prime, so a generator must
//! emits a 64-bit integer output: for example, SplitMix emits its full 64-bit
//! word, and MIXMAX (61 bits) is left-justified into 64 bits.

// When no PRNG feature is selected, the _prng marker (enabled by every PRNG
// feature in Cargo.toml) is off. Emit one clear error AND expose a placeholder Prng
// so the rest of the crate still type-checks.
#[cfg(not(feature = "_prng"))]
compile_error!(
    "no PRNG selected: enable exactly one PRNG feature, as in --features splitmix \
     (see the [features] table in Cargo.toml)"
);

#[cfg(not(feature = "_prng"))]
mod placeholder {
    #[derive(Clone, Copy)]
    pub struct Prng;
    impl Prng {
        pub const NAME: &str = "(no generator selected)";
        pub fn new(_seed: u64) -> Self {
            Self
        }
        #[inline(always)]
        pub fn next_u64(&mut self) -> u64 {
            0
        }
    }
}

#[cfg(not(feature = "_prng"))]
pub use placeholder::Prng;

// ----- SplitMix -------------------------------------------------------------------

#[cfg(feature = "splitmix")]
#[derive(Clone, Copy)]
pub struct Prng {
    x: u64,
}

#[cfg(feature = "splitmix")]
impl Prng {
    pub const NAME: &str = "SplitMix";
    pub fn new(seed: u64) -> Self {
        Self { x: seed }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        const PHI: u64 = 0x9e3779b97f4a7c15;
        let mut z = self.x;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        self.x = self.x.wrapping_add(PHI);
        z ^ (z >> 31)
    }
}

// ----- MIXMAX (matrix generator over the field of size 2⁶¹ − 1) ----------------------
//
// Bit-for-bit port of CERN ROOT's MIXMAX. One generic engine carries the
// four template parameters of ROOT's MixMaxEngine<N, SkipNumber> as const generics
// (SPECIALMUL/SPECIAL are derived from N in the C #if ladder; here they are
// explicit const parameters), and the three ROOT typedefs are feature-gated Prng
// wrappers. Seeding uses ROOT's seed_spbox (its SetSeedFast routine).

#[cfg(any(feature = "mixmax", feature = "mixmax17", feature = "mixmax256"))]
mod mixmax {
    // mod_mulspec/fmodmul_m61 are exercised only by the SPECIAL≠0 variants
    // (N=240, N=256), so they look dead in an N=17 build.
    #![allow(dead_code)]

    /// The prime modulus *p* = 2⁶¹ − 1.
    const MM_P: u64 = 2305843009213693951;
    const MM_BITS: u32 = 61;

    /// Payne reduction modulo 2⁶¹ − 1 (non-canonical: the result may slightly
    /// exceed *p*, exactly as ROOT's MOD_MERSENNE/MOD_PAYNE does).
    #[inline(always)]
    fn mod_mersenne(k: u64) -> u64 {
        (k & MM_P) + (k >> MM_BITS)
    }

    /// (*a* + *b*) mod *p*, non-canonical (ROOT modadd).
    #[inline(always)]
    fn modadd(a: u64, b: u64) -> u64 {
        mod_mersenne(a.wrapping_add(b))
    }

    /// (*a*·*b* + cum) mod *p* via 128-bit arithmetic (ROOT fmodmulM61/mod128).
    #[inline(always)]
    fn fmodmul_m61(cum: u64, a: u64, b: u64) -> u64 {
        let s = (a as u128)
            .wrapping_mul(b as u128)
            .wrapping_add(cum as u128);
        let lo = s as u64;
        let hi = (s >> 64) as u64;
        let s1 = (lo & MM_P)
            .wrapping_add(hi.wrapping_mul(8))
            .wrapping_add(lo >> MM_BITS);
        mod_mersenne(s1)
    }

    /// MixMax<N, SPECIALMUL, SPECIAL, SKIP> mirrors ROOT's MixMaxEngine<N, SkipNumber>:
    /// SPECIALMUL/SPECIAL encode the per-N matrix tweak, SKIP is the SkipNumber
    /// thinning. SPECIAL == u64::MAX represents the C value -1 (the N=256 variant),
    /// whose special-entry map is *p* − *k*; any other nonzero SPECIAL uses SPECIAL·*k*.
    #[derive(Clone, Copy)]
    pub struct MixMax<const N: usize, const SPECIALMUL: u32, const SPECIAL: u64, const SKIP: usize> {
        v: [u64; N],
        sumtot: u64,
        counter: usize,
    }

    impl<const N: usize, const SPECIALMUL: u32, const SPECIAL: u64, const SKIP: usize>
        MixMax<N, SPECIALMUL, SPECIAL, SKIP>
    {
        /// Knuth's 64-bit LCG multiplier, used by seed_spbox.
        const MULT64: u64 = 6364136223846793005;

        /// Seeds via ROOT's seed_spbox: a 64-bit LCG plus a half-word swap fills the
        /// state, then counter = N so the first draw re-iterates the whole vector.
        /// A zero seed (which the C code rejects with exit) is remapped to a fixed
        /// nonzero constant.
        pub fn new(seed: u64) -> Self {
            let mut l = if seed == 0 { 0x9e3779b97f4a7c15 } else { seed };
            let mut v = [0; N];
            let mut sumtot: u64 = 0;
            let mut ovflow: u64 = 0;
            for slot in v.iter_mut() {
                l = l.wrapping_mul(Self::MULT64);
                l = (l << 32) ^ (l >> 32);
                *slot = l & MM_P;
                let (s, c) = sumtot.overflowing_add(*slot);
                sumtot = s;
                ovflow += c as u64;
            }
            let sumtot = mod_mersenne(mod_mersenne(sumtot).wrapping_add(ovflow << 3));
            Self {
                v,
                sumtot,
                counter: N,
            }
        }

        /// The special-entry map MOD_MULSPEC for the SPECIAL≠0 variants.
        #[inline(always)]
        fn mod_mulspec(k: u64) -> u64 {
            if SPECIAL == u64::MAX {
                MM_P.wrapping_sub(k) // N=256: MERSBASE − k (C SPECIAL == -1)
            } else {
                fmodmul_m61(0, SPECIAL, k) // N=240: SPECIAL · k mod p
            }
        }

        /// One full vector application — ROOT's iterate_raw_vec. sumtot is a raw
        /// u64 allowed to wrap; each 2⁶⁴ wrap is reinjected as 2⁶⁴ mod *p* = 8
        /// through ovflow << 3.
        fn iterate(&mut self) {
            let temp2 = self.v[1]; // used only when SPECIAL != 0
            let mut temp_v = self.sumtot;
            self.v[0] = temp_v;
            let mut sumtot = temp_v;
            let mut ovflow: u64 = 0;
            let mut temp_p: u64 = 0;
            for i in 1..N {
                if SPECIALMUL != 0 {
                    // MULWU(temp_p): 61-bit rotate left by SPECIALMUL = ·2^SPECIALMUL mod p.
                    let temp_po =
                        ((temp_p << SPECIALMUL) & MM_P) | (temp_p >> (MM_BITS - SPECIALMUL));
                    temp_p = modadd(temp_p, self.v[i]);
                    temp_v = mod_mersenne(temp_v.wrapping_add(temp_p).wrapping_add(temp_po));
                } else {
                    temp_p = modadd(temp_p, self.v[i]);
                    temp_v = modadd(temp_v, temp_p);
                }
                self.v[i] = temp_v;
                let (s, c) = sumtot.overflowing_add(temp_v);
                sumtot = s;
                ovflow += c as u64;
            }
            if SPECIAL != 0 {
                let t2 = Self::mod_mulspec(temp2);
                self.v[2] = modadd(self.v[2], t2);
                let (s, c) = sumtot.overflowing_add(t2);
                sumtot = s;
                ovflow += c as u64;
            }
            self.sumtot = mod_mersenne(mod_mersenne(sumtot).wrapping_add(ovflow << 3));
        }

        /// ROOT's get_next/GET_BY_MACRO: hands out entry counter of V, re-iterating
        /// the whole vector and resuming at entry 1 on exhaustion.
        #[inline(always)]
        fn get_next(&mut self) -> u64 {
            let i = self.counter;
            if i < N {
                self.counter += 1;
                self.v[i]
            } else {
                self.iterate();
                self.counter = 2;
                self.v[1]
            }
        }

        /// One raw 61-bit output in the range 1 to 2⁶¹ − 1, reproducing ROOT's
        /// MixMaxEngine<N, SkipNumber>::IntRndm: when a re-iterate is imminent
        /// (counter == N), perform SKIP extra iterations first, then draw.
        #[inline(always)]
        pub fn next_raw(&mut self) -> u64 {
            let counter = self.counter;
            if counter >= N {
                for _ in 0..SKIP {
                    self.iterate();
                }
            }
            self.counter = counter;
            self.get_next()
        }
    }
}

#[cfg(feature = "mixmax")]
#[derive(Clone, Copy)]
pub struct Prng(mixmax::MixMax<240, 51, 487013230256099140, 0>);

#[cfg(feature = "mixmax")]
impl Prng {
    pub const NAME: &str = "MIXMAX (TRandomMixMax, N=240, s=487013230256099140)";
    pub fn new(seed: u64) -> Self {
        Self(mixmax::MixMax::new(seed))
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        self.0.next_raw() << 3
    }
}

#[cfg(feature = "mixmax17")]
#[derive(Clone, Copy)]
pub struct Prng(mixmax::MixMax<17, 36, 0, 0>);

#[cfg(feature = "mixmax17")]
impl Prng {
    pub const NAME: &str = "MIXMAX (TRandomMixMax17, N=17)";
    pub fn new(seed: u64) -> Self {
        Self(mixmax::MixMax::new(seed))
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        self.0.next_raw() << 3
    }
}

#[cfg(feature = "mixmax256")]
#[derive(Clone, Copy)]
pub struct Prng(mixmax::MixMax<256, 0, { u64::MAX }, 2>);

#[cfg(feature = "mixmax256")]
impl Prng {
    pub const NAME: &str = "MIXMAX (TRandomMixMax256, N=256, skip=2)";
    pub fn new(seed: u64) -> Self {
        Self(mixmax::MixMax::new(seed))
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        self.0.next_raw() << 3
    }
}
