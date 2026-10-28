/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

// Bit-exactness check against the official MIXMAX C source. The expected values
// are the first 10 raw 61-bit outputs printed by the unmodified ROOT
// mixmax.h/.icc (spbox-seeded) for the seed below. next_u64() left-shifts the
// raw output by 3, so we compare next_u64() >> 3.
// The mixmax-lux decimation discards part of the raw output stream, so these
// bit-exactness vectors (the first ten consecutive ROOT outputs) do not apply to it.
#[cfg(all(
    any(feature = "mixmax", feature = "mixmax17", feature = "mixmax256"),
    not(feature = "mixmax-lux")
))]
mod mixmax_ref {
    use modlin::prng::Prng;

    const SEED: u64 = 1234567890123456789;

    #[cfg(feature = "mixmax17")]
    const EXPECTED: [u64; 10] = [
        725153038902651861,
        1378221679953707950,
        1286954255613532202,
        2183378240693018162,
        480095223529858829,
        1028486193651132422,
        1683776858464191402,
        752698145398996242,
        1196878046144236647,
        1221112387665998636,
    ];

    #[cfg(feature = "mixmax")]
    const EXPECTED: [u64; 10] = [
        1894724080392768937,
        2178403735632896525,
        82974682673747095,
        1362029167703491378,
        175589691633294316,
        1981764336260928646,
        1644676922156291251,
        463159379824582110,
        73102843550620337,
        783499878333038877,
    ];

    #[cfg(feature = "mixmax256")]
    const EXPECTED: [u64; 10] = [
        379984832372095102,
        2262036764019001978,
        952890261108261993,
        1848869284066170746,
        412703182074144958,
        1796627384723021327,
        1959941208714421405,
        1149684300799040022,
        525039624349349232,
        110320404022327900,
    ];

    #[test]
    fn test_matches_official_mixmax() -> anyhow::Result<()> {
        let mut p = Prng::try_new(SEED)?;
        for (i, &expected) in EXPECTED.iter().enumerate() {
            let raw = p.next_u64() >> 3;
            assert_eq!(raw, expected, "MIXMAX output {i} disagrees with ROOT");
        }
        Ok(())
    }
}

// MIXMAX rejects the zero seed, which would yield the all-zero state.
#[cfg(any(feature = "mixmax", feature = "mixmax17", feature = "mixmax256"))]
mod mixmax_seed {
    use modlin::prng::Prng;

    #[test]
    fn test_mixmax_rejects_zero_seed() {
        assert_eq!(
            Prng::try_new(0).unwrap_err().to_string(),
            "the zero seed is not allowed for MIXMAX"
        );
        assert!(Prng::try_new(1).is_ok());
    }
}

// Bit-exactness check against the reference xoroshiro128++ algorithm
// (<https://prng.di.unimi.it/xoroshiro128plusplus.c>). The expected values are
// the first 10 outputs of the algorithm whose 128-bit state is SplitMix-seeded
// from the seed below, computed from the published recurrence independently of
// this crate. (The internal recurrence is anchored separately: from the raw
// state [1, 0] the first output is 131073 and the next state [562949955518465,
// 268435456].)
#[cfg(feature = "xoroshiro128pp")]
mod xoroshiro_ref {
    use modlin::prng::Prng;

    const SEED: u64 = 0;

    const EXPECTED: [u64; 10] = [
        8027914721839836897,
        13805533416164201645,
        5256508173613850168,
        7973558954284022901,
        8526501294691771125,
        6116102375994396471,
        16028966417245382669,
        12808598746819302742,
        15824426267781808726,
        5829521525559713354,
    ];

    #[test]
    fn test_matches_reference_xoroshiro() -> anyhow::Result<()> {
        let mut p = Prng::try_new(SEED)?;
        for (i, &expected) in EXPECTED.iter().enumerate() {
            assert_eq!(
                p.next_u64(),
                expected,
                "xoroshiro128++ output {i} disagrees with the reference"
            );
        }
        Ok(())
    }
}

// C++11 fixes the 10000th output of a default-constructed std::minstd_rand0
// (state 1) at 1043618065; the first outputs below were printed by libc++.
#[cfg(feature = "minstd")]
mod minstd_ref {
    use modlin::prng::Prng;

    #[test]
    fn test_minstd_matches_cpp() {
        let mut p = Prng::from_state(1);
        let expected = [16807, 282475249, 1622650073, 984943658, 1144108930];
        for (i, &e) in expected.iter().enumerate() {
            assert_eq!(p.next_u64(), e, "MINSTD output {i}");
        }
        let mut p = Prng::from_state(1);
        let mut x = 0;
        for _ in 0..10000 {
            x = p.next_u64();
        }
        assert_eq!(x, 1043618065);
    }
}

// The expected values were computed by an independent transcription of the
// recurrence xₙ = (107374182·xₙ₋₁ + 104480·xₙ₋₅) mod 2³¹ − 1 from the state
// (1, 2, 3, 4, 5), oldest value first.
#[cfg(feature = "mrg5")]
mod mrg5_ref {
    use modlin::prng::Prng;

    #[test]
    fn test_mrg5_matches_recurrence() {
        let mut p = Prng::from_state([1, 2, 3, 4, 5]);
        let expected = [
            536975390, 886009397, 1515571251, 651083988, 631636463, 279519110, 1807621113,
            793306269, 104290668, 2027673655,
        ];
        for (i, &e) in expected.iter().enumerate() {
            assert_eq!(p.next_u64(), e, "MRG5 output {i}");
        }
    }
}

// The expected values are those of L'Ecuyer's reference implementation of
// MRG32k3a with all six seeds equal to 12345, before normalization (the first
// normalized output is 545508589 · 2.328306549295728e-10 = 0.127011122...).
#[cfg(feature = "mrg32k3a")]
mod mrg32k3a_ref {
    use modlin::prng::Prng;

    #[test]
    fn test_mrg32k3a_matches_reference() {
        let mut p = Prng::from_state([12345; 3], [12345; 3]);
        let expected = [
            545508589, 1368065410, 1327943761, 3546985096, 951893194, 2290915636, 2064909380,
            1527117980, 584065747, 3246360482,
        ];
        for (i, &e) in expected.iter().enumerate() {
            assert_eq!(p.next_u64(), e, "MRG32k3a output {i}");
        }
    }
}

#[cfg(feature = "mrg32k3a-c1")]
mod mrg32k3a_c1_ref {
    use modlin::prng::Prng;

    #[test]
    fn test_mrg32k3a_c1_matches_reference() {
        let mut p = Prng::from_state([12345; 3]);
        let expected = [
            3023790853, 3023790853, 3385359573, 1322208174, 2930192941, 2462079208, 2386811717,
            2989318136, 3378525425, 1773647758,
        ];
        for (i, &e) in expected.iter().enumerate() {
            assert_eq!(p.next_u64(), e, "MRG32k3a first-component output {i}");
        }
    }
}

// C++11 fixes the 10000th output of a default-constructed std::mt19937_64 (seed
// 5489) at 9981545732273789042; the first outputs below were printed by libc++.
#[cfg(feature = "mt19937-64")]
mod mt19937_64_ref {
    use modlin::prng::Prng;

    #[test]
    fn test_mt19937_64_matches_cpp() -> anyhow::Result<()> {
        let mut p = Prng::try_new(5489)?;
        let expected = [
            14514284786278117030,
            4620546740167642908,
            13109570281517897720,
            17462938647148434322,
            355488278567739596,
        ];
        for (i, &e) in expected.iter().enumerate() {
            assert_eq!(p.next_u64(), e, "MT19937-64 output {i}");
        }
        let mut p = Prng::try_new(5489)?;
        let mut x = 0;
        for _ in 0..10000 {
            x = p.next_u64();
        }
        assert_eq!(x, 9981545732273789042);
        Ok(())
    }
}

// The first outputs of SplitMix from the zero state.
#[cfg(feature = "splitmix")]
mod splitmix_ref {
    use modlin::prng::Prng;

    #[test]
    fn test_splitmix_matches_reference() -> anyhow::Result<()> {
        let mut p = Prng::try_new(0)?;
        for (i, &e) in [0xe220a8397b1dcdaf, 0x6e789e6aa1b965f4, 0x06c45d188009454f]
            .iter()
            .enumerate()
        {
            assert_eq!(p.next_u64(), e, "SplitMix output {i}");
        }
        Ok(())
    }
}

// The expected values were printed by NumPy's Philox (Philox4x64-10) with key
// (1, 0): NumPy increments the counter before each block, so its first block is
// that of counter 1, which is the second block of Prng::try_new(1).
#[cfg(feature = "philox4x64")]
mod philox4x64_ref {
    use modlin::prng::Prng;

    const EXPECTED: [u64; 8] = [
        5599841837815857887,
        15655913098571550255,
        2880178291573394738,
        573812481542357666,
        16607021752199172844,
        960460491066566553,
        13744443775079107998,
        4542806826184744765,
    ];

    #[test]
    fn test_philox4x64_matches_numpy() -> anyhow::Result<()> {
        assert_eq!(Prng::block([1, 0, 0, 0], [1, 0]), EXPECTED[..4]);
        assert_eq!(Prng::block([2, 0, 0, 0], [1, 0]), EXPECTED[4..]);
        let mut p = Prng::from_state([1, 0, 0, 0], [1, 0]);
        for (i, &e) in EXPECTED.iter().enumerate() {
            assert_eq!(p.next_u64(), e, "Philox4x64-10 output {i}");
        }
        let mut p = Prng::try_new(1)?;
        for _ in 0..4 {
            p.next_u64();
        }
        for (i, &e) in EXPECTED.iter().enumerate() {
            assert_eq!(p.next_u64(), e, "Philox4x64-10 output {i} after seeding");
        }
        Ok(())
    }
}

// Seeding through SplitMix accepts every seed, and distinct seeds yield
// distinct streams.
#[cfg(any(
    feature = "minstd",
    feature = "mrg5",
    feature = "mrg32k3a",
    feature = "mrg32k3a-c1"
))]
mod seeding {
    use modlin::prng::Prng;

    #[test]
    fn test_try_new_accepts_all_seeds() -> anyhow::Result<()> {
        let mut streams = Vec::new();
        for seed in [0, 1, 2, u64::MAX] {
            let mut p = Prng::try_new(seed)?;
            streams.push((0..8).map(|_| p.next_u64()).collect::<Vec<_>>());
        }
        for (i, a) in streams.iter().enumerate() {
            for b in &streams[i + 1..] {
                assert_ne!(a, b);
            }
        }
        Ok(())
    }
}
