# Modular rank and linear-complexity tests for pseudorandom number generators

This crate implements two empirical tests for pseudorandom number generators
(PRNGs), the _modular rank test_ and the _modular linear-complexity test_. These
tests are generalizations of the standard binary-rank and (binary)
linear-complexity tests from 𝐅₂ to an arbitrary finite field.

While the _binary_ rank and linear-complexity tests are widely used in the
literature (e.g., [dieharder], the [NIST suite], [TestU01], [PractRand]
implement them), they can only detect linearity in generators that are linear
over 𝐅₂. The modular rank and linear-complexity tests can find bias in
generators that are linear over any finite field. This includes linear
congruential generators, single multiple-recursive recurrences, and matrix
generators such as [MIXMAX]; indeed, implementing tests detecting the
statistical defects of [MIXMAX] was the original motivation for this crate.

As in the binary case, the two tests are different measures of the linear
complexity of the output of the generator. The rank test measure it by the rank
of a matrix of outputs, while the linear-complexity test measures it directly
from the output stream using the Berlekamp–Massey algorithm. The first is more
expensive, but robust, as it can find bias even in presence of moderate
scrambling, whereas the second is cheaper but fragile, as it requires the output
to obey a single linear recurrence of low degree.

## Testing linearity over 𝐅*ₚ*

Both tests are implemented for fields 𝐅*ₚ* of prime order _p_, with _p_ <
2⁶³. A generator that is linear over 𝐅*ₚ* has a finite, usually small linear
complexity _L_; a generator with no such structure does not. For a generator
that emits _b_ values per step from a _k_-dimensional state _L_ ≤ _b_ · _k_; for
example, [MIXMAX]-_N_ emits _N_ − 1 values per step, giving _L_ = _N_(_N_ − 1).

The _modular rank test_ reads *n*² successive outputs into an _n_ × _n_ matrix
over 𝐅*ₚ*. Every length-_n_ window of an _L_-linear stream is fixed by the
recurrence, so all the rows lie in a subspace of dimension at most _L_, and the
matrix has rank at most _L_. A side _n_ > _L_ then forces the matrix to be _rank
deficient_ (in practice, the deficiency already appears at far smaller sides). A
generator with no linear structure instead yields full-rank matrices: a uniform
random _n_ × _n_ matrix over 𝐅*ₚ* has corank (i.e., side minus rank) _d_
with probability ≈ *p*⁻ᵈ², so it is singular only with probability ≈ 1/_p_. The
test reports, for each matrix, its own _p_-value.

The _modular linear-complexity test_ reads _L_ off the stream directly, with the
Berlekamp–Massey algorithm, which returns the order of the shortest linear
recurrence a sequence obeys. A length-_n_ stream from an _L_-linear generator
has complexity _L_ once _n_ passes about 2*L*, far below the expected ⌈_n_/2⌉.
Also in this case the test reports, for each sequence, its own _p_-value.

Note that for large fields deviation from the typical case is astronomically rare,
so a single anomalous matrix or sequence already has a per-sample
_p_-value very close to 0.

## Build and run

The generator is selected at build time via a Cargo feature. Exactly one of
`-R`/`-L` selects the test. For example, we can find bias in the 17-dimensional
[MIXMAX] generator from [CERN's ROOT], for which _p_ = 2⁶¹ − 1, using a 500×500
matrix in milliseconds:

```bash
cargo run -r -F mixmax17 -- -R 500 -p 2305843009213693951 -S 1
Generator: MIXMAX (TRandomMixMax17, N=17)
Seed: 0x0000000000000001
Running a modular rank test using 1.907 MiB of RAM: 1 500×500 matrix over the field of size 2305843009213693951
2026-06-21 21:15:42.580 6ms INFO [ThreadId(1)] modlin - Generating matrix entries...
2026-06-21 21:15:42.581 7ms INFO [ThreadId(1)] modlin - Completed.
2026-06-21 21:15:42.581 8ms INFO [ThreadId(1)] modlin - Elapsed: 1ms [250,000 outputs, 223955739.18 outputs/s, 4.47 ns/output]; res/vir/avail/free/total mem 8.96MB/420.73GB/32.78GB/1.93GB/68.72GB
2026-06-21 21:15:42.581 8ms INFO [ThreadId(1)] modlin - Matrix 1/1: ranking (blocked Gaussian elimination over Fₚ)...
2026-06-21 21:15:42.590 17ms INFO [ThreadId(1)] modlin - Completed.
2026-06-21 21:15:42.590 17ms INFO [ThreadId(1)] modlin - Elapsed: 8ms [500 columns, 55700.74 columns/s, 17.95 μs/column]; res/vir/avail/free/total mem 9.91MB/420.79GB/32.78GB/1.93GB/68.72GB
Matrix 1/1	corank=432	p=1e-307
```

The corank should be zero, and the _p_-value should be 1, for a generator with
no linear dependencies.

Finding bias using the Berlekamp–Massey algorithm to measure the linear complexity of a
sequence of 1000 elements is even faster:

```bash
cargo run -r -F mixmax17 -- -L 1000 -p 2305843009213693951 -S 1
Generator: MIXMAX (TRandomMixMax17, N=17)
Seed: 0x0000000000000001
Running a modular linear-complexity test using 0.031 MiB of RAM: 1 sequence of length 1000 over the field of size 2305843009213693951
2026-06-21 21:15:50.127 7ms INFO [ThreadId(1)] modlin - Sequence 1/1: Berlekamp–Massey over Fₚ...
2026-06-21 21:15:50.128 8ms INFO [ThreadId(1)] modlin - Completed.
2026-06-21 21:15:50.129 8ms INFO [ThreadId(1)] modlin - Elapsed: 1ms [1,000 steps, 854335.75 steps/s, 1.17 μs/step]; res/vir/avail/free/total mem 7.03MB/420.73GB/32.76GB/1.94GB/68.72GB
Sequence 1/1	linear complexity=272	p=1e-307
```

The linear complexity here should be approximately 500, and the _p_-value should
be 1, for a generator with no linear dependencies.

The same test finds bias on the largest provided [MIXMAX] generator in [CERN's
ROOT] (256-dimensional) in less than a minute:

```bash
cargo run -r -F mixmax256 -- -L 200000 -p 2305843009213693951 -S 1
Generator: MIXMAX (TRandomMixMax256, N=256, skip=2)
Seed: 0x0000000000000001
Running a modular linear-complexity test using 6.104 MiB of RAM: 1 sequence of length 200000 over the field of size 2305843009213693951
2026-06-21 21:16:06.814 9ms INFO [ThreadId(1)] modlin - Sequence 1/1: Berlekamp–Massey over Fₚ...
2026-06-21 21:16:16.814 10s9ms INFO [ThreadId(1)] modlin - 89,274 steps, 10s, 8927.06 steps/s, 112.02 μs/step; 44.64% done, 12s to end; res/vir/avail/free/total mem 13.43MB/420.73GB/32.80GB/1.91GB/68.72GB
[...]
2026-06-21 21:16:40.334 33s529ms INFO [ThreadId(1)] modlin - Completed.
2026-06-21 21:16:40.335 33s529ms INFO [ThreadId(1)] modlin - Elapsed: 33s [200,000 steps, 5966.55 steps/s, 167.60 μs/step]; res/vir/avail/free/total mem 13.52MB/420.74GB/32.80GB/1.91GB/68.72GB
Sequence 1/1	linear complexity=65280	p=1e-307
```

Running a modular rank test capable of detecting bias on the same generator
requires instead a couple of hours, albeit the time depends on the amount of
cores, as the Gaussian elimination is parallelized. You can customize the amount
of parallism with the environment variable `RAYON_NUM_THREADS`.

Running the same test on [`xoroshiro128++`], a generator without linear
dependencies will find no bias for any _p_:

```bash
cargo run -r -F xoroshiro128pp -- -R 1000 -p 2
Generator: xoroshiro128++
Seed: 0x0000000000000000
Running a modular rank test using 7.629 MiB of RAM: 1 1000×1000 matrix over the field of size 2
2026-06-21 21:15:50.337 8ms INFO [ThreadId(1)] modlin - Generating matrix entries...
2026-06-21 21:15:50.346 17ms INFO [ThreadId(1)] modlin - Completed.
2026-06-21 21:15:50.346 17ms INFO [ThreadId(1)] modlin - Elapsed: 8ms [1,000,000 outputs, 118677923.19 outputs/s, 8.43 ns/output]; res/vir/avail/free/total mem 14.88MB/420.73GB/32.76GB/1.93GB/68.72GB
2026-06-21 21:15:50.346 17ms INFO [ThreadId(1)] modlin - Matrix 1/1: ranking (blocked Gaussian elimination over Fₚ)...
2026-06-21 21:15:50.454 125ms INFO [ThreadId(1)] modlin - Completed.
2026-06-21 21:15:50.455 125ms INFO [ThreadId(1)] modlin - Elapsed: 108ms [1,000 columns, 9252.37 columns/s, 108.08 μs/column]; res/vir/avail/free/total mem 16.30MB/420.79GB/32.76GB/1.93GB/68.72GB
Matrix 1/1	corank=0	p=1

Generator: xoroshiro128++
Seed: 0x0000000000000000
Running a modular linear-complexity test using 0.031 MiB of RAM: 1 sequence of length 1000 over the field of size 2
2026-06-21 21:15:50.668 6ms INFO [ThreadId(1)] modlin - Sequence 1/1: Berlekamp–Massey over Fₚ...
2026-06-21 21:15:50.669 8ms INFO [ThreadId(1)] modlin - Completed.
2026-06-21 21:15:50.669 8ms INFO [ThreadId(1)] modlin - Elapsed: 1ms [1,000 steps, 798057.85 steps/s, 1.25 μs/step]; res/vir/avail/free/total mem 6.86MB/420.73GB/32.76GB/1.93GB/68.72GB
Sequence 1/1	linear complexity=503	p=0.9947916666666666

cargo run -r -F xoroshiro128pp -- -R 1000 -p 2305843009213693951
Generator: xoroshiro128++
Seed: 0x0000000000000000
Running a modular rank test using 7.629 MiB of RAM: 1 1000×1000 matrix over the field of size 2305843009213693951
2026-06-21 21:15:59.086 9ms INFO [ThreadId(1)] modlin - Generating matrix entries...
2026-06-21 21:15:59.089 12ms INFO [ThreadId(1)] modlin - Completed.
2026-06-21 21:15:59.089 12ms INFO [ThreadId(1)] modlin - Elapsed: 3ms [1,000,000 outputs, 305429031.58 outputs/s, 3.27 ns/output]; res/vir/avail/free/total mem 14.22MB/420.59GB/32.74GB/1.93GB/68.72GB
2026-06-21 21:15:59.089 12ms INFO [ThreadId(1)] modlin - Matrix 1/1: ranking (blocked Gaussian elimination over Fₚ)...
2026-06-21 21:15:59.225 149ms INFO [ThreadId(1)] modlin - Completed.
2026-06-21 21:15:59.226 149ms INFO [ThreadId(1)] modlin - Elapsed: 136ms [1,000 columns, 7331.94 columns/s, 136.39 μs/column]; res/vir/avail/free/total mem 15.58MB/420.63GB/32.74GB/1.93GB/68.72GB
Matrix 1/1	corank=0	p=1

cargo run -r -F xoroshiro128pp -- -L 10000 -p 2305843009213693951
Generator: xoroshiro128++
Seed: 0x0000000000000000
Running a modular linear-complexity test using 0.305 MiB of RAM: 1 sequence of length 10000 over the field of size 2305843009213693951
2026-06-21 21:15:59.449 9ms INFO [ThreadId(1)] modlin - Sequence 1/1: Berlekamp–Massey over Fₚ...
2026-06-21 21:15:59.579 138ms INFO [ThreadId(1)] modlin - Completed.
2026-06-21 21:15:59.579 138ms INFO [ThreadId(1)] modlin - Elapsed: 129ms [10,000 steps, 77110.19 steps/s, 12.97 μs/step]; res/vir/avail/free/total mem 8.01MB/420.88GB/32.78GB/1.97GB/68.72GB
Sequence 1/1	linear complexity=5000	p=1
```

Note that _p_-values are _one-sided_: a _p_-value near 0 thus flags an anomaly,
whereas a _p_-value near 1 is just the generic case—not a failure.

You can also repeat the test multiple times: the test is simply run again on
disjoint, contiguous stretches of the orbit, printing one _p_-value per matrix
(or sequence), and you can decide how to combine them (e.g., a simple Bonferroni
correction, or a full-scale χ² test). For large _p_, however, the probability of
a deficient matrix under the randomness hypothesis is so low that a single test
is normally enough to rule out the randomness hypothesis (as above).

There is some progress logging during the computation, and you can adjust the
frequency with a command-line option or the logging level using the `RUST_LOG`
environment variable.

## Adding your own generator

To add a new generator, add a feature in `Cargo.toml` and a corresponding
implementation in the [`prng`] module.

[dieharder]: https://webhome.phy.duke.edu/~rgb/General/dieharder.php
[NIST suite]: https://csrc.nist.gov/projects/random-bit-generation/documentation-and-software
[TestU01]: https://dl.acm.org/doi/10.1145/1268776.1268777
[PractRand]: https://pracrand.sourceforge.net/
[MIXMAX]: https://doi.org/10.1016/j.cpc.2015.06.003
[`prng`]: https://docs.rs/modlin/latest/modlin/prng/index.html
[CERN's ROOT]: https://root.cern/doc/v628/classROOT_1_1Math_1_1MixMaxEngine.html
[`xoroshiro128++`]: https://prng.di.unimi.it/xoroshiro128plusplus.c
