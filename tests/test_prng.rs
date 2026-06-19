/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

// Bit-exactness check against the official MIXMAX C source. The expected values
// are the first 10 raw 61-bit outputs printed by the unmodified ROOT
// mixmax.h/.icc (spbox-seeded) for the seed below. next_u64() left-shifts the
// raw output by 3, so we compare next_u64() >> 3.
// The mixmax-star output scrambler deliberately alters the raw output, so these
// bit-exactness vectors (which check the unscrambled ROOT output) do not apply to it.
#[cfg(all(
    any(feature = "mixmax", feature = "mixmax17", feature = "mixmax256"),
    not(feature = "mixmax-star")
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
    fn matches_official_mixmax() {
        let mut p = Prng::try_new(SEED).unwrap();
        for (i, &expected) in EXPECTED.iter().enumerate() {
            let raw = p.next_u64() >> 3;
            assert_eq!(raw, expected, "MIXMAX output {i} disagrees with ROOT");
        }
    }
}

// Bit-exactness check against the reference xoroshiro128++ algorithm
// (<https://prng.di.unimi.it/xoroshiro128plusplus.c>). The expected values are
// the first 10 outputs of the algorithm whose 128-bit state is SplitMix64-seeded
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
    fn matches_reference_xoroshiro() {
        let mut p = Prng::try_new(SEED).unwrap();
        for (i, &expected) in EXPECTED.iter().enumerate() {
            assert_eq!(
                p.next_u64(),
                expected,
                "xoroshiro128++ output {i} disagrees with the reference"
            );
        }
    }
}
