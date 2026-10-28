# Change Log

## [0.2.0] - 2026-10-28

### New

- New generators: MINSTD (`minstd`), the fifth-order MRG of L'Ecuyer, Blouin,
  and Couture (`mrg5`), MRG32k3a (`mrg32k3a`) and its first component alone
  (`mrg32k3a-c1`), MT19937-64 (`mt19937-64`), SplitMix (`splitmix`), and
  Philox4x64-10 (`philox4x64`). MINSTD, MRG5, both MRG32k3a variants, and
  Philox also have a `from_state` constructor, and Philox exposes its `block`
  function.

- New module `poly`, with the half-GCD linear-complexity algorithm
  `half_gcd`, `half_gcd_memory`, and the constants `MAX_THREADS` and
  `MAX_LENGTH`.

- New items in `fp`: `berlekamp_massey` (the quadratic algorithm, previously
  the body of `linear_complexity`), `HGCD_THRESHOLD`,
  `linear_complexity_memory`, and `Field::p`.

- New error type `prng::SeedError`.

- New feature `slow_tests`, which runs the slowest tests.

### Changed

- `stats::corank_prob`, `stats::corank_tail_pvalue`, and
  `stats::lc_left_tail_pvalue` are now `stats::corank_probability`,
  `stats::corank_p_value`, and `stats::linear_complexity_p_value`.

- `Prng::try_new` returns a `SeedError` instead of a `&'static str`.

- `Field` implements `Debug`, `PartialEq`, `Eq`, and `Hash`, and the
  generators implement `Debug`; `Field::new`, `Field::reduce`,
  `Field::reduce_u128`, and several methods of the generators are `const`.

- The generators no longer implement `Copy`, as an accidental copy repeats
  the outputs.

- `fp::rank` panics if the length of the matrix is not the square of its side
  (the check was a debug assertion).

- `fp::linear_complexity` panics on sequences longer than `poly::MAX_LENGTH`.

- `stats::corank_probability`, `stats::corank_p_value`, and
  `stats::linear_complexity_p_value` panic if the modulus is smaller than 2,
  and the last two also if the corank or the linear complexity exceeds the side
  or the length.

- The linear-complexity test uses at most `poly::MAX_THREADS` (16) threads, as
  more threads slow it down.

- Combining `mixmax-lux` with a generator other than MIXMAX is a compile-time
  error (it was silently ignored).

- The command line rejects lengths larger than `poly::MAX_LENGTH` (2³² − 2),
  reports a seed rejected by the generator (e.g., the default zero seed for
  MIXMAX) as an argument error, and accepts the modulus in the same formats as
  the seed.

### Improved

- The linear-complexity test uses, from length 2048 on, a half-GCD algorithm
  with polynomial multiplication by number-theoretic transforms modulo three
  primes, instead of the quadratic Berlekamp–Massey algorithm: a sequence of
  length 10⁶ now requires seconds instead of about half an hour.

- The trailing update of the rank test removes bounds checks and branches from
  its inner products, breaks their carry chains with independent accumulators,
  and, for primes up to 2⁶¹, updates two rows at a time: depending on size and
  number of threads, the rank test is up to twice as fast for primes up to
  2⁶¹, and 10–20% faster for larger primes.

### Fixed

- `berlekamp_massey` (previously `linear_complexity`) no longer panics on an
  empty sequence.

- `stats::corank_probability` (previously `stats::corank_prob`) no longer
  returns NaN for coranks from about 2³¹ on, which made
  `stats::corank_p_value` run for a very long time.

- The command line rejects matrix sides from 2³⁰ on, whose matrices do not fit
  in the address space (sides below 2³² caused a panic).

## [0.1.2] - 2026-07-28

### Changed

- Changes to documentation.

## [0.1.1] - 2026-07-23

### Changed

- License is now MIT or Apache.

## [0.1.0] - 2026-07-01

### New

- First release.
