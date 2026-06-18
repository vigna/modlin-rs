/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

use modlin::prng::Prng;

#[test]
fn skip_matches_repeated_next() {
    let seed = 0x0123_4567_89ab_cdef;
    for &n in &[0u64, 1, 2, 7, 1000, 100_000] {
        let mut a = Prng::new(seed);
        if a.try_skip(n).is_err() {
            return; // generator has no jump-ahead; nothing to check
        }
        let mut b = Prng::new(seed);
        for _ in 0..n {
            b.next_u64();
        }
        for k in 0..64 {
            assert_eq!(
                a.next_u64(),
                b.next_u64(),
                "skip({n}) disagreed with stepping at output {k}"
            );
        }
    }
    // Composition: two successive skips must equal one combined skip.
    let (x, y) = (1u64 << 40, (1u64 << 41) + 12_345);
    let mut p = Prng::new(seed);
    p.try_skip(x).unwrap();
    p.try_skip(y).unwrap();
    let mut q = Prng::new(seed);
    q.try_skip(x + y).unwrap();
    for k in 0..64 {
        assert_eq!(
            p.next_u64(),
            q.next_u64(),
            "skip composition failed at output {k}"
        );
    }
}

// Bit-exactness check against the official MIXMAX C source. The expected values
// are the first 10 raw 61-bit outputs printed by the unmodified ROOT
// mixmax.h/.icc (spbox-seeded) for the seed below. next_u64() left-shifts the
// raw output by 3, so we compare next_u64() >> 3.
#[cfg(any(feature = "mixmax", feature = "mixmax17", feature = "mixmax256"))]
mod mixmax_ref {
    use super::*;

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
        let mut p = Prng::new(SEED);
        for (i, &expected) in EXPECTED.iter().enumerate() {
            let raw = p.next_u64() >> 3;
            assert_eq!(raw, expected, "MIXMAX output {i} disagrees with ROOT");
        }
    }
}
