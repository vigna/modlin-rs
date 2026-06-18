# Modular rank and linear-complexity tests for pseudorandom number generators

This crate implements two empirical tests for pseudorandom number generators
(PRNGs), the _modular rank test_ and the _modular linear-complexity test_. These
tests are generalizations of the standard binary-rank and (binary)
linear-complexity tests from **F**₂ to an arbitrary finite field.

While the _binary_ rank and linear-complexity tests are widely used in the
literature (e.g., [dieharder], the [NIST suite], [TestU01], [PractRand]
implement them), they can only detect linearity in generators that are linear
over **F**₂. The modular rank and linear-complexity tests can find bias in
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

## Testing linearity over **F**_ₚ_

Both tests are implemented for fields **F**_ₚ_ of prime order _p_, with _p_ <
2⁶³. A generator that is linear over **F**_ₚ_ has a finite, usually small linear
complexity _L_; a generator with no such structure does not. For a generator
that emits _b_ values per step from a _k_-dimensional state _L_ ≤ _b_ · _k_; for
example, [MIXMAX]-_N_ emits _N_ − 1 values per step, giving _L_ = _N_(_N_ − 1).

The _modular rank test_ reads *n*² successive outputs into an _n_ × _n_ matrix
over **F**_ₚ_. Every length-_n_ window of an _L_-linear stream is fixed by the
recurrence, so all the rows lie in a subspace of dimension at most _L_, and the
matrix has rank at most _L_. A side _n_ > _L_ then forces the matrix to be _rank
deficient_ (in practice, the deficiency already appears at far smaller sides). A
generator with no linear structure instead yields full-rank matrices: a uniform
random _n_ × _n_ matrix over **F**_ₚ_ has corank (i.e., side minus rank) _d_
with probability ≈ *p*⁻ᵈ², so it is singular only with probability ≈ 1/_p_. The
test reports, for each matrix, its own _p_-value.

The _modular linear-complexity test_ reads _L_ off the stream directly, with the
Berlekamp–Massey algorithm, which returns the order of the shortest linear
recurrence a sequence obeys. A length-_n_ stream from an _L_-linear generator
has complexity _L_ once _n_ passes about 2*L*, far below the expected ⌈_n_/2⌉.
Also in this case the test reports, for each sequence, its own _p_-value.

Note that for large fields deviation from the typical case is astronomically rare,
so a single anomalous matrix or sequence already has a per-sample
_p_-value of essentially 0.

## Build and run

The generator is selected at build time via a Cargo feature. Exactly one of
`-R`/`-L` selects the test. For example, we can find bias in the 17-dimensional
[MIXMAX] generator from [CERN's ROOT], for which _p_ = 2⁶¹ − 1, using a 500×500
matrix in milliseconds:

```bash
cargo run --release --features mixmax17 --  -R 500 -p 2305843009213693951
Generator: MIXMAX (TRandomMixMax17, N=17)
Seed: 0x18ba1ef3747cbab0
Running a modular rank test: 1 500×500 matrix over the field of size 2305843009213693951
2026-06-18 08:17:13.570 12ms INFO [ThreadId(1)] modlin - Generating matrix entries...
2026-06-18 08:17:13.572 13ms INFO [ThreadId(1)] modlin - Completed.
2026-06-18 08:17:13.572 14ms INFO [ThreadId(1)] modlin - Elapsed: 1ms [250,000 outputs, 150492894.33 outputs/s, 6.64 ns/output]; res/vir/avail/free/total mem 10.34MB/420.89GB/35.25GB/17.18GB/68.72GB
2026-06-18 08:17:13.572 14ms INFO [ThreadId(1)] modlin - Matrix 1/1: ranking (blocked Gaussian elimination over Fₚ)...
2026-06-18 08:17:13.584 25ms INFO [ThreadId(1)] modlin - Completed.
2026-06-18 08:17:13.584 25ms INFO [ThreadId(1)] modlin - Elapsed: 11ms [500 columns, 44010.21 columns/s, 22.72 μs/column]; res/vir/avail/free/total mem 11.34MB/421.07GB/35.25GB/17.18GB/68.72GB
Matrix 1/1	corank=432	p=2.2250738585072014e-308
```

The corank should be zero, and the _p_-value should be 1, for a generator with
no linear dependencies.

Finding bias using the Berlekamp–Massey algorithm to measure the linear complexity of a
sequence of 1000 elements is even faster:

```bash
cargo run --release --features mixmax17 -- -L 1000 -p 2305843009213693951
Generator: MIXMAX (TRandomMixMax17, N=17)
Seed: 0x0000000000000000
Running a modular linear-complexity test: 1 sequence of length 1000 over the field of size 2305843009213693951
2026-06-18 15:52:08.823 6ms INFO [ThreadId(1)] modlin - Sequence 1/1: Berlekamp–Massey over Fₚ...
2026-06-18 15:52:08.824 7ms INFO [ThreadId(1)] modlin - Completed.
2026-06-18 15:52:08.825 7ms INFO [ThreadId(1)] modlin - Elapsed: 0ms [1,000 steps, 1073777.08 steps/s, 931.29 ns/step]; res/vir/avail/free/total mem 6.80MB/420.59GB/10.91GB/414.12MB/68.72GB
Sequence 1/1	linear complexity=272	p=2.2250738585072014e-308
```

The linear complexity here should be approximately 500, and the _p_-value should
be 1, for a generator with no linear dependencies.

The same test finds bias on the largest provided [MIXMAX] generator in [CERN's
ROOT] (256-dimensional) in less than a minute:

```bash
cargo run --release --features mixmax256 -- -L 200000 -p 2305843009213693951
Generator: MIXMAX (TRandomMixMax256, N=256, skip=2)
Seed: 0x0000000000000000
Running a modular linear-complexity test: 1 sequence of length 200000 over the field of size 2305843009213693951
2026-06-18 15:52:35.390 9ms INFO [ThreadId(1)] modlin - Sequence 1/1: Berlekamp–Massey over Fₚ...
2026-06-18 15:52:45.391 10s10ms INFO [ThreadId(1)] modlin - 74,268 steps, 10s, 7426.50 steps/s, 134.65 μs/step; 37.13% done, 16s to end; res/vir/avail/free/total mem 15.14MB/421.02GB/10.71GB/508.49MB/68.72GB
[...]
2026-06-18 15:53:21.981 46s599ms INFO [ThreadId(1)] modlin - Completed.
2026-06-18 15:53:21.981 46s600ms INFO [ThreadId(1)] modlin - Elapsed: 46s [200,000 steps, 4292.76 steps/s, 232.95 μs/step]; res/vir/avail/free/total mem 15.29MB/421.04GB/10.71GB/508.49MB/68.72GB
Sequence 1/1	linear complexity=65280	p=2.2250738585072014e-308
```

Running a modular rank test capable of detecting bias on the same generator
requires instead a couple of hours, albeit the time depends on the amount of
cores, as the Gaussian elimination is parallelized. You can customize the amount
of parallism with the environment variable `RAYON_NUM_THREADS`.

Running the same test on a generator without linear dependencies will find no bias
for any _p_:

```bash
cargo run --release --features splitmix -- -R 1000 -p 2
Generator: SplitMix
Seed: 0x0000000000000000
Running a modular rank test: 1 1000×1000 matrix over the field of size 2
2026-06-18 16:09:49.192 7ms INFO [ThreadId(1)] modlin - Generating matrix entries...
2026-06-18 16:09:49.208 23ms INFO [ThreadId(1)] modlin - Completed.
2026-06-18 16:09:49.208 23ms INFO [ThreadId(1)] modlin - Elapsed: 15ms [1,000,000 outputs, 64040305.43 outputs/s, 15.62 ns/output]; res/vir/avail/free/total mem 14.94MB/420.74GB/8.80GB/74.35MB/68.72GB
2026-06-18 16:09:49.208 23ms INFO [ThreadId(1)] modlin - Matrix 1/1: ranking (blocked Gaussian elimination over Fₚ)...
2026-06-18 16:09:49.969 785ms INFO [ThreadId(1)] modlin - Completed.
2026-06-18 16:09:49.970 785ms INFO [ThreadId(1)] modlin - Elapsed: 761ms [1,000 columns, 1313.77 columns/s, 761.17 μs/column]; res/vir/avail/free/total mem 16.42MB/420.92GB/8.80GB/74.66MB/68.72GB
Matrix 1/1	corank=1	p=0.7112119049133976

cargo run --release --features splitmix -- -L 10000 -p 2
Generator: SplitMix
Seed: 0x0000000000000000
Running a modular linear-complexity test: 1 sequence of length 10000 over the field of size 2
2026-06-18 16:03:59.777 11ms INFO [ThreadId(1)] modlin - Sequence 1/1: Berlekamp–Massey over Fₚ...
2026-06-18 16:03:59.878 112ms INFO [ThreadId(1)] modlin - Completed.
2026-06-18 16:03:59.878 112ms INFO [ThreadId(1)] modlin - Elapsed: 100ms [10,000 steps, 99169.54 steps/s, 10.08 μs/step]; res/vir/avail/free/total mem 8.70MB/420.88GB/11.23GB/172.15MB/68.72GB
Sequence 1/1	linear complexity=4999	p=0.16666666666666666

cargo run --release --features splitmix -- -R 1000 -p 2305843009213693951
Generator: SplitMix
Seed: 0x0000000000000000
Running a modular rank test: 1 1000×1000 matrix over the field of size 2305843009213693951
2026-06-18 16:10:52.572 7ms INFO [ThreadId(1)] modlin - Generating matrix entries...
2026-06-18 16:10:52.580 14ms INFO [ThreadId(1)] modlin - Completed.
2026-06-18 16:10:52.580 14ms INFO [ThreadId(1)] modlin - Elapsed: 7ms [1,000,000 outputs, 136794227.28 outputs/s, 7.31 ns/output]; res/vir/avail/free/total mem 15.14MB/420.74GB/8.57GB/636.93MB/68.72GB
2026-06-18 16:10:52.580 14ms INFO [ThreadId(1)] modlin - Matrix 1/1: ranking (blocked Gaussian elimination over Fₚ)...
2026-06-18 16:10:54.547 1s982ms INFO [ThreadId(1)] modlin - Completed.
2026-06-18 16:10:54.547 1s982ms INFO [ThreadId(1)] modlin - Elapsed: 1s [1,000 columns, 508.31 columns/s, 1.97 ms/column]; res/vir/avail/free/total mem 16.42MB/420.79GB/8.57GB/637.09MB/68.72GB
Matrix 1/1	corank=0	p=1

cargo run --release --features splitmix -- -L 10000 -p 2305843009213693951
Generator: SplitMix
Seed: 0x0000000000000000
Running a modular linear-complexity test: 1 sequence of length 10000 over the field of size 2305843009213693951
2026-06-18 16:05:49.498 7ms INFO [ThreadId(1)] modlin - Sequence 1/1: Berlekamp–Massey over Fₚ...
2026-06-18 16:05:49.771 280ms INFO [ThreadId(1)] modlin - Completed.
2026-06-18 16:05:49.772 281ms INFO [ThreadId(1)] modlin - Elapsed: 273ms [10,000 steps, 36572.38 steps/s, 27.34 μs/step]; res/vir/avail/free/total mem 8.42MB/420.89GB/11.66GB/60.96MB/68.72GB
Sequence 1/1	linear complexity=5000	p=1
```

Note that _p_-values are _one-sided_: a _p_-value near 0 thus flags an anomaly,
whereas _p_ = 1 is just the generic case—not a failure.

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
