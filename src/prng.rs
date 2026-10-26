/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

//! Pseudorandom number generators selected at build time via Cargo features.
//!
//! Exactly one generator feature must be enabled when building the crate; each
//! variant exposes a single `Prng` type with a `NAME` constant, a
//! `try_new(seed: u64) -> Result<Self, &'static str>` constructor (which
//! rejects seeds a generator cannot accept), and a `next_u64(&mut self) -> u64`
//! step function. Some variants have additional constructors from an explicit
//! state, used to check them against reference outputs.
//!
//! The test reduces each output modulo the field prime. The 64-bit generators
//! emit their full 64-bit output, and MIXMAX (61 bits) is left-justified into
//! 64 bits; the generators whose output has 32 bits or fewer (MINSTD, the
//! fifth-order MRG, and MRG32k3a) emit it unshifted, so it is uniform over
//! their own field, but not over larger ones.

// When no PRNG feature is selected, the _prng marker (enabled by every PRNG
// feature in Cargo.toml) is off. Emit one clear error AND expose a placeholder Prng
// so the rest of the crate still type-checks.
#[cfg(not(feature = "_prng"))]
compile_error!(
    "no PRNG selected: enable exactly one PRNG feature, as in --features xoroshiro128pp \
     (see the [features] table in Cargo.toml)"
);

#[cfg(all(
    feature = "mixmax-lux",
    not(any(feature = "mixmax", feature = "mixmax17", feature = "mixmax256"))
))]
compile_error!(
    "mixmax-lux is an output modifier of the MIXMAX generators: combine it with \
     mixmax, mixmax17, or mixmax256, as in --features mixmax17,mixmax-lux"
);

#[cfg(not(feature = "_prng"))]
mod placeholder {
    #[derive(Clone, Copy)]
    pub struct Prng;
    impl Prng {
        pub const NAME: &str = "(no generator selected)";
        pub fn try_new(_seed: u64) -> Result<Self, &'static str> {
            Ok(Self)
        }
        #[inline(always)]
        pub fn next_u64(&mut self) -> u64 {
            0
        }
    }
}

#[cfg(not(feature = "_prng"))]
pub use placeholder::Prng;

// ----- SplitMix step ---------------------------------------------------------
//
// One step of SplitMix: used as the SplitMix generator, and to expand the
// 64-bit seed into the state of other generators.

#[cfg(any(
    feature = "xoroshiro128pp",
    feature = "splitmix",
    feature = "minstd",
    feature = "mrg5",
    feature = "mrg32k3a",
    feature = "mrg32k3a-c1"
))]
#[inline(always)]
const fn splitmix(z: &mut u64) -> u64 {
    *z = z.wrapping_add(0x9e3779b97f4a7c15);
    let mut w = *z;
    w = (w ^ (w >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    w = (w ^ (w >> 27)).wrapping_mul(0x94d049bb133111eb);
    w ^ (w >> 31)
}

// ----- xoroshiro128++ --------------------------------------------------------

#[cfg(feature = "xoroshiro128pp")]
#[derive(Clone, Copy)]
pub struct Prng {
    s: [u64; 2],
}

#[cfg(feature = "xoroshiro128pp")]
impl Prng {
    pub const NAME: &str = "xoroshiro128++";

    pub fn try_new(seed: u64) -> Result<Self, &'static str> {
        // SplitMix expansion of the 64-bit seed into the 128-bit state; every
        // seed (including 0) maps to a valid nonzero state, so none is rejected.
        let mut z = seed;
        let s0 = splitmix(&mut z);
        let s1 = splitmix(&mut z);
        Ok(Self { s: [s0, s1] })
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        let s0 = self.s[0];
        let mut s1 = self.s[1];
        let result = s0.wrapping_add(s1).rotate_left(17).wrapping_add(s0);
        s1 ^= s0;
        self.s[0] = s0.rotate_left(49) ^ s1 ^ (s1 << 21);
        self.s[1] = s1.rotate_left(28);
        result
    }
}

// ----- MIXMAX (matrix generator over the field of size 2⁶¹ − 1) --------------
//
// Bit-for-bit port of CERN ROOT's MIXMAX. One generic engine carries the four
// template parameters of ROOT's MixMaxEngine<N, SkipNumber> as const generics.
// Seeding uses ROOT's seed_spbox (its SetSeedFast routine).

// The parameters of the mixmax-lux decimation, as literals, so that they can
// appear both in mixmax::LUX_KEEP/LUX_DROP and in the generator names.
#[cfg(any(feature = "mixmax", feature = "mixmax17", feature = "mixmax256"))]
macro_rules! lux_keep {
    () => {
        5
    };
}
#[cfg(any(feature = "mixmax", feature = "mixmax17", feature = "mixmax256"))]
macro_rules! lux_drop {
    () => {
        9
    };
}

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
        #[cfg(feature = "mixmax-lux")]
        lux_kept: u32,
    }

    impl<const N: usize, const SPECIALMUL: u32, const SPECIAL: u64, const SKIP: usize>
        MixMax<N, SPECIALMUL, SPECIAL, SKIP>
    {
        /// Knuth's 64-bit LCG multiplier, used by seed_spbox.
        const MULT64: u64 = 6364136223846793005;

        /// Seeds via ROOT's seed_spbox: a 64-bit LCG plus a half-word swap fills the
        /// state, then counter = N so the first draw re-iterates the whole vector.
        /// The zero seed is rejected (ROOT's C code likewise rejects it, with exit),
        /// since it would seed the all-zero state.
        pub fn try_new(seed: u64) -> Result<Self, &'static str> {
            if seed == 0 {
                return Err("the zero seed is not allowed for MIXMAX");
            }
            let mut l = seed;
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
            Ok(Self {
                v,
                sumtot,
                counter: N,
                #[cfg(feature = "mixmax-lux")]
                lux_kept: 0,
            })
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

        /// One full vector application (ROOT's iterate_raw_vec). sumtot is a raw
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

        /// Final 64-bit output: the raw 61-bit value left-justified by 3. With
        /// `mixmax-lux`, decimation keeps `LUX_KEEP` outputs then discards `LUX_DROP`.
        #[inline(always)]
        pub fn output(&mut self) -> u64 {
            #[cfg(feature = "mixmax-lux")]
            {
                if self.lux_kept == LUX_KEEP {
                    for _ in 0..LUX_DROP {
                        let _ = self.next_raw();
                    }
                    self.lux_kept = 0;
                }
                self.lux_kept += 1;
            }
            self.next_raw() << 3
        }
    }

    /// `mixmax-lux` decimation pattern: keep `LUX_KEEP` consecutive outputs, then
    /// discard `LUX_DROP` before resuming (RANLUX-style "luxury"). Decimation is a
    /// linear operation, so the kept stream stays linear over the field; it only
    /// raises the linear complexity by a bounded factor (≤ LUX_KEEP).
    #[cfg(feature = "mixmax-lux")]
    pub const LUX_KEEP: u32 = lux_keep!();
    #[cfg(feature = "mixmax-lux")]
    pub const LUX_DROP: u32 = lux_drop!();
}

#[cfg(feature = "mixmax")]
#[derive(Clone, Copy)]
pub struct Prng(mixmax::MixMax<240, 51, 487013230256099140, 0>);

#[cfg(feature = "mixmax")]
impl Prng {
    pub const NAME: &str = if cfg!(feature = "mixmax-lux") {
        concat!(
            "MIXMAX (TRandomMixMax, N=240, s=487013230256099140), lux: keep ",
            lux_keep!(),
            " drop ",
            lux_drop!()
        )
    } else {
        "MIXMAX (TRandomMixMax, N=240, s=487013230256099140)"
    };
    pub fn try_new(seed: u64) -> Result<Self, &'static str> {
        Ok(Self(mixmax::MixMax::try_new(seed)?))
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        self.0.output()
    }
}

#[cfg(feature = "mixmax17")]
#[derive(Clone, Copy)]
pub struct Prng(mixmax::MixMax<17, 36, 0, 0>);

#[cfg(feature = "mixmax17")]
impl Prng {
    pub const NAME: &str = if cfg!(feature = "mixmax-lux") {
        concat!(
            "MIXMAX (TRandomMixMax17, N=17), lux: keep ",
            lux_keep!(),
            " drop ",
            lux_drop!()
        )
    } else {
        "MIXMAX (TRandomMixMax17, N=17)"
    };
    pub fn try_new(seed: u64) -> Result<Self, &'static str> {
        Ok(Self(mixmax::MixMax::try_new(seed)?))
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        self.0.output()
    }
}

#[cfg(feature = "mixmax256")]
#[derive(Clone, Copy)]
pub struct Prng(mixmax::MixMax<256, 0, { u64::MAX }, 2>);

#[cfg(feature = "mixmax256")]
impl Prng {
    pub const NAME: &str = if cfg!(feature = "mixmax-lux") {
        concat!(
            "MIXMAX (TRandomMixMax256, N=256, skip=2), lux: keep ",
            lux_keep!(),
            " drop ",
            lux_drop!()
        )
    } else {
        "MIXMAX (TRandomMixMax256, N=256, skip=2)"
    };
    pub fn try_new(seed: u64) -> Result<Self, &'static str> {
        Ok(Self(mixmax::MixMax::try_new(seed)?))
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        self.0.output()
    }
}

// ----- MINSTD (Park–Miller LCG over the field of size 2³¹ − 1) ---------------
//
// x ← 16807·x mod 2³¹ − 1 (C++'s minstd_rand0). The output is the new state, a
// value in [1..2³¹ − 1), emitted unshifted, so that reduction modulo 2³¹ − 1
// is the identity.

#[cfg(feature = "minstd")]
#[derive(Clone, Copy)]
pub struct Prng {
    x: u64,
}

#[cfg(feature = "minstd")]
impl Prng {
    pub const NAME: &str = "MINSTD (x ← 16807·x mod 2³¹ − 1)";
    const M: u64 = (1 << 31) - 1;
    const A: u64 = 16807;

    pub fn try_new(seed: u64) -> Result<Self, &'static str> {
        let mut z = seed;
        Ok(Self::from_state(1 + splitmix(&mut z) % (Self::M - 1)))
    }

    /// Starts from the given state.
    ///
    /// # Panics
    ///
    /// If the state is not in [1..2³¹ − 1).
    pub const fn from_state(x: u64) -> Self {
        assert!(0 < x && x < Self::M);
        Self { x }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        self.x = self.x * Self::A % Self::M;
        self.x
    }
}

// ----- Fifth-order MRG (L'Ecuyer, Blouin, and Couture, 1993) -----------------
//
// xₙ = (107374182·xₙ₋₁ + 104480·xₙ₋₅) mod 2³¹ − 1 (GSL's gsl_rng_mrg). The
// output is xₙ, unshifted.

#[cfg(feature = "mrg5")]
#[derive(Clone, Copy)]
pub struct Prng {
    /// The last five values, oldest first.
    s: [u64; 5],
}

#[cfg(feature = "mrg5")]
impl Prng {
    pub const NAME: &str = "MRG5 (xₙ = 107374182·xₙ₋₁ + 104480·xₙ₋₅ mod 2³¹ − 1)";
    const M: u64 = (1 << 31) - 1;
    const A1: u64 = 107374182;
    const A5: u64 = 104480;

    pub fn try_new(seed: u64) -> Result<Self, &'static str> {
        let mut z = seed;
        loop {
            let s = std::array::from_fn(|_| splitmix(&mut z) % Self::M);
            if s.iter().any(|&x| x != 0) {
                return Ok(Self::from_state(s));
            }
        }
    }

    /// Starts from the given state, oldest value first.
    ///
    /// # Panics
    ///
    /// If an entry of the state is not below 2³¹ − 1, or all entries are zero.
    pub fn from_state(s: [u64; 5]) -> Self {
        assert!(s.iter().all(|&x| x < Self::M) && s.iter().any(|&x| x != 0));
        Self { s }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        // Both products are below 2⁵⁸, so the sum cannot overflow.
        let x = (Self::A1 * self.s[4] + Self::A5 * self.s[0]) % Self::M;
        self.s.copy_within(1.., 0);
        self.s[4] = x;
        x
    }
}

// ----- MRG32k3a (L'Ecuyer, 1999) ---------------------------------------------
//
// Two third-order MRGs modulo the distinct primes m₁ = 2³² − 209 and m₂ = 2³² −
// 22853, combined by subtraction modulo m₁. The feature mrg32k3a emits the
// combined output in [1..m₁], as in L'Ecuyer's reference implementation
// (before normalization to (0..1)); the feature mrg32k3a-c1 emits instead
// the first component alone, in [0..m₁), which is linear over the field of
// size m₁. Outputs are unshifted.

#[cfg(any(feature = "mrg32k3a", feature = "mrg32k3a-c1"))]
mod mrg32k3a {
    /// The modulus of the first component, 2³² − 209.
    pub const M1: i64 = 4294967087;
    /// The modulus of the second component, 2³² − 22853.
    #[cfg(feature = "mrg32k3a")]
    pub const M2: i64 = 4294944443;
    const A12: i64 = 1403580;
    const A13N: i64 = 810728;
    #[cfg(feature = "mrg32k3a")]
    const A21: i64 = 527612;
    #[cfg(feature = "mrg32k3a")]
    const A23N: i64 = 1370589;

    /// The state of a component, oldest value first.
    ///
    /// # Panics
    ///
    /// If an entry is not below `m`, or all entries are zero.
    pub fn component(s: [u64; 3], m: i64) -> [i64; 3] {
        assert!(s.iter().all(|&x| x < m as u64) && s.iter().any(|&x| x != 0));
        s.map(|x| x as i64)
    }

    /// A component state for modulus `m` expanded from `z` by SplitMix.
    pub fn seed_component(z: &mut u64, m: i64) -> [i64; 3] {
        loop {
            let s: [u64; 3] = std::array::from_fn(|_| super::splitmix(z) % m as u64);
            if s.iter().any(|&x| x != 0) {
                return component(s, m);
            }
        }
    }

    /// Advances the first component and returns its new value.
    #[inline(always)]
    pub fn next1(s: &mut [i64; 3]) -> u64 {
        // |A12·s − A13N·s| < 2⁵³, so i64 arithmetic is exact.
        let x = (A12 * s[1] - A13N * s[0]).rem_euclid(M1);
        *s = [s[1], s[2], x];
        x as u64
    }

    /// Advances the second component and returns its new value.
    #[cfg(feature = "mrg32k3a")]
    #[inline(always)]
    pub fn next2(s: &mut [i64; 3]) -> u64 {
        let x = (A21 * s[2] - A23N * s[0]).rem_euclid(M2);
        *s = [s[1], s[2], x];
        x as u64
    }
}

#[cfg(feature = "mrg32k3a")]
#[derive(Clone, Copy)]
pub struct Prng {
    s1: [i64; 3],
    s2: [i64; 3],
}

#[cfg(feature = "mrg32k3a")]
impl Prng {
    pub const NAME: &str = "MRG32k3a (combined, m₁ = 2³² − 209, m₂ = 2³² − 22853)";

    pub fn try_new(seed: u64) -> Result<Self, &'static str> {
        let mut z = seed;
        let s1 = mrg32k3a::seed_component(&mut z, mrg32k3a::M1);
        let s2 = mrg32k3a::seed_component(&mut z, mrg32k3a::M2);
        Ok(Self { s1, s2 })
    }

    /// Starts from the given component states, oldest value first.
    ///
    /// # Panics
    ///
    /// If an entry of a state is not below the modulus of its component, or all
    /// entries of a state are zero.
    pub fn from_state(s1: [u64; 3], s2: [u64; 3]) -> Self {
        Self {
            s1: mrg32k3a::component(s1, mrg32k3a::M1),
            s2: mrg32k3a::component(s2, mrg32k3a::M2),
        }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        let p1 = mrg32k3a::next1(&mut self.s1);
        let p2 = mrg32k3a::next2(&mut self.s2);
        if p1 > p2 {
            p1 - p2
        } else {
            p1 + mrg32k3a::M1 as u64 - p2
        }
    }
}

#[cfg(feature = "mrg32k3a-c1")]
#[derive(Clone, Copy)]
pub struct Prng {
    s: [i64; 3],
}

#[cfg(feature = "mrg32k3a-c1")]
impl Prng {
    pub const NAME: &str = "MRG32k3a, first component only (m₁ = 2³² − 209)";

    /// The state is expanded from the seed as the first component of
    /// MRG32k3a, so the output is that of the first component of MRG32k3a with
    /// the same seed.
    pub fn try_new(seed: u64) -> Result<Self, &'static str> {
        let mut z = seed;
        Ok(Self {
            s: mrg32k3a::seed_component(&mut z, mrg32k3a::M1),
        })
    }

    /// Starts from the given state, oldest value first.
    ///
    /// # Panics
    ///
    /// If an entry of the state is not below 2³² − 209, or all entries are
    /// zero.
    pub fn from_state(s: [u64; 3]) -> Self {
        Self {
            s: mrg32k3a::component(s, mrg32k3a::M1),
        }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        mrg32k3a::next1(&mut self.s)
    }
}

// ----- MT19937-64 (Nishimura, 2000) ------------------------------------------
//
// Port of Matsumoto and Nishimura's mt19937-64.c, seeded by init_genrand64. The
// generator is linear over the field with two elements; its full 64-bit output
// is emitted.

/// The number of 64-bit words of the state of MT19937-64.
#[cfg(feature = "mt19937-64")]
const MT_NN: usize = 312;

#[cfg(feature = "mt19937-64")]
#[derive(Clone)]
pub struct Prng {
    mt: [u64; MT_NN],
    mti: usize,
}

#[cfg(feature = "mt19937-64")]
impl Prng {
    pub const NAME: &str = "MT19937-64";
    const MM: usize = 156;
    const MATRIX_A: u64 = 0xB5026F5AA96619E9;
    const UM: u64 = 0xFFFFFFFF80000000;
    const LM: u64 = 0x7FFFFFFF;

    pub fn try_new(seed: u64) -> Result<Self, &'static str> {
        let mut mt = [0; MT_NN];
        mt[0] = seed;
        for i in 1..MT_NN {
            mt[i] = 6364136223846793005u64
                .wrapping_mul(mt[i - 1] ^ (mt[i - 1] >> 62))
                .wrapping_add(i as u64);
        }
        Ok(Self { mt, mti: MT_NN })
    }

    #[inline(always)]
    fn twist(&mut self, i: usize, j: usize, k: usize) {
        let x = (self.mt[i] & Self::UM) | (self.mt[j] & Self::LM);
        let mag = if x & 1 == 0 { 0 } else { Self::MATRIX_A };
        self.mt[i] = self.mt[k] ^ (x >> 1) ^ mag;
    }

    fn refill(&mut self) {
        for i in 0..MT_NN - Self::MM {
            self.twist(i, i + 1, i + Self::MM);
        }
        for i in MT_NN - Self::MM..MT_NN - 1 {
            self.twist(i, i + 1, i + Self::MM - MT_NN);
        }
        self.twist(MT_NN - 1, 0, Self::MM - 1);
        self.mti = 0;
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        if self.mti >= MT_NN {
            self.refill();
        }
        let mut x = self.mt[self.mti];
        self.mti += 1;
        x ^= (x >> 29) & 0x5555555555555555;
        x ^= (x << 17) & 0x71D67FFFEDA60000;
        x ^= (x << 37) & 0xFFF7EEE000000000;
        x ^ (x >> 43)
    }
}

// ----- SplitMix (Steele, Lea, and Flood, 2014) -------------------------------
//
// A Weyl sequence followed by a nonlinear mixing function; the seed is the
// initial state.

#[cfg(feature = "splitmix")]
#[derive(Clone, Copy)]
pub struct Prng {
    z: u64,
}

#[cfg(feature = "splitmix")]
impl Prng {
    pub const NAME: &str = "SplitMix";

    pub const fn try_new(seed: u64) -> Result<Self, &'static str> {
        Ok(Self { z: seed })
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        splitmix(&mut self.z)
    }
}

// ----- Philox4x64-10 (Salmon, Moraes, Dror, and Shaw, 2011) ------------------
//
// Counter-based generator of Random123 (also in cuRAND and NumPy): ten rounds
// of a multiply-and-xor bijection applied to a 256-bit counter under a 128-bit
// key. The key is (seed, 0), the counter starts at zero and is incremented
// after each block; the four 64-bit words of each block are emitted in order.

#[cfg(feature = "philox4x64")]
#[derive(Clone, Copy)]
pub struct Prng {
    key: [u64; 2],
    ctr: [u64; 4],
    buf: [u64; 4],
    pos: usize,
}

#[cfg(feature = "philox4x64")]
impl Prng {
    pub const NAME: &str = "Philox4x64-10";
    const M0: u64 = 0xD2E7470EE14C6C93;
    const M1: u64 = 0xCA5A826395121157;
    const W0: u64 = 0x9E3779B97F4A7C15;
    const W1: u64 = 0xBB67AE8584CAA73B;

    pub const fn try_new(seed: u64) -> Result<Self, &'static str> {
        Ok(Self::from_state([0; 4], [seed, 0]))
    }

    /// Starts from the given counter and key.
    pub const fn from_state(ctr: [u64; 4], key: [u64; 2]) -> Self {
        Self {
            key,
            ctr,
            buf: [0; 4],
            pos: 4,
        }
    }

    /// The ten-round Philox4x64 bijection of `ctr` under `key`.
    pub fn block(ctr: [u64; 4], key: [u64; 2]) -> [u64; 4] {
        let mulhilo = |a: u64, b: u64| {
            let p = a as u128 * b as u128;
            ((p >> 64) as u64, p as u64)
        };
        let mut c = ctr;
        let mut k = key;
        for round in 0..10 {
            if round > 0 {
                k = [k[0].wrapping_add(Self::W0), k[1].wrapping_add(Self::W1)];
            }
            let (hi0, lo0) = mulhilo(Self::M0, c[0]);
            let (hi1, lo1) = mulhilo(Self::M1, c[2]);
            c = [hi1 ^ c[1] ^ k[0], lo1, hi0 ^ c[3] ^ k[1], lo0];
        }
        c
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        if self.pos == 4 {
            self.buf = Self::block(self.ctr, self.key);
            // 256-bit increment, least significant word first.
            for w in self.ctr.iter_mut() {
                *w = w.wrapping_add(1);
                if *w != 0 {
                    break;
                }
            }
            self.pos = 0;
        }
        let x = self.buf[self.pos];
        self.pos += 1;
        x
    }
}

#[cfg(all(test, feature = "mixmax-lux"))]
mod tests {
    use super::*;

    /// The decimated stream must be the raw stream with LUX_DROP values
    /// discarded after every LUX_KEEP values.
    #[test]
    fn test_lux_decimates_raw_stream() -> anyhow::Result<()> {
        let mut lux = Prng::try_new(1234567890123456789).map_err(anyhow::Error::msg)?;
        let mut raw = lux;
        for _ in 0..20 {
            for _ in 0..mixmax::LUX_KEEP {
                assert_eq!(lux.next_u64(), raw.0.next_raw() << 3);
            }
            for _ in 0..mixmax::LUX_DROP {
                raw.0.next_raw();
            }
        }
        Ok(())
    }
}
