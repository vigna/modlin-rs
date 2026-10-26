# Change Log

## [0.2.0] - unreleased

### New

- New generators: MINSTD (`minstd`), the fifth-order MRG of L'Ecuyer, Blouin,
  and Couture (`mrg5`), MRG32k3a (`mrg32k3a`) and its first component alone
  (`mrg32k3a-c1`), MT19937-64 (`mt19937-64`), SplitMix (`splitmix`), and
  Philox4x64-10 (`philox4x64`). MINSTD, MRG5, both MRG32k3a variants, and
  Philox also have a `from_state` constructor, and Philox exposes its `block`
  function.

- New module `poly`, with the half-GCD `linear_complexity`, `memory_bytes`,
  and the constants `MAX_THREADS` and `MAX_LENGTH`.

- New items in `fp`: `berlekamp_massey` (the quadratic algorithm, previously
  the body of `linear_complexity`), `HGCD_THRESHOLD`,
  `linear_complexity_algorithm`, `linear_complexity_memory`, and `Field::p`.

- New feature `slow_tests`, which runs the slowest tests.

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

### Changed

- The linear-complexity test uses at most `poly::MAX_THREADS` (16) threads, as
  more threads slow it down.

- Combining `mixmax-lux` with a generator other than MIXMAX is a compile-time
  error (it was silently ignored).

- The command line rejects lengths larger than `poly::MAX_LENGTH` (2³² − 2).

- `cli::init_env_logger` returns the result of `try_init`.

### Fixed

- `berlekamp_massey` (previously `linear_complexity`) no longer panics on an
  empty sequence.

## [0.1.2] - 2026-07-28

### Changed

- Changes to documentation.

## [0.1.1] - 2026-07-23

### Changed

- License is now MIT or Apache.

## [0.1.0] - 2026-07-01

### New

- First release.
