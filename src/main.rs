/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

use std::time::Duration;

use clap::Parser;
use dsi_progress_logger::prelude::*;
use pluralizer::pluralize;

use modlin::cli::{self, Args};
use modlin::fp::{self, Field};
use modlin::prng::Prng;
use modlin::stats::{corank_p_value, linear_complexity_p_value, pretty_p_value};

fn main() -> Result<(), log::SetLoggerError> {
    let args = Args::parse();
    args.validate();
    // The generator may reject the seed (e.g., MIXMAX forbids the zero seed,
    // which is the default).
    let rng = Prng::try_new(args.seed)
        .unwrap_or_else(|e| Args::die(format!("{e}; choose another seed with -S")));
    cli::init_env_logger()?;

    println!("Generator: {}", Prng::NAME);
    println!("Seed: {:#018x}", args.seed);

    let reps = args.reps;
    let field = Field::new(args.modulus);

    if let Some(length) = args.linear_complexity {
        run_linear_complexity(&field, length, reps, rng, args.log_interval);
    } else if let Some(side) = args.rank {
        run_rank(&field, side, reps, rng, args.log_interval);
    } else {
        unreachable!();
    }
    Ok(())
}

/// Formats a byte count in GiB from 1 GiB on, and in MiB below, so that small
/// amounts don't degenerate into `0.000 GiB`.
fn format_memory(bytes: usize) -> String {
    const GIB: f64 = (1u64 << 30) as f64;
    const MIB: f64 = (1u64 << 20) as f64;
    let bytes = bytes as f64;
    if bytes >= GIB {
        format!("{:.3} GiB", bytes / GIB)
    } else {
        format!("{:.3} MiB", bytes / MIB)
    }
}

/// A progress logger displaying memory, for `expected_updates` updates on
/// items named `item_name`.
fn new_logger(item_name: &str, expected_updates: usize, log_interval: Duration) -> ProgressLogger {
    progress_logger![
        item_name = item_name,
        display_memory = true,
        log_interval = log_interval,
        expected_updates = Some(expected_updates),
    ]
}

/// Fills `buf` with the next outputs of `rng`, reduced modulo the size of
/// `field`.
fn fill(field: &Field, rng: &mut Prng, buf: &mut [u64]) {
    for x in buf {
        *x = field.reduce(rng.next_u64());
    }
}

/// Modular rank test: rank `reps` independent *n* × *n* matrices, one at a
/// time, each under its own progress line, printing for each its one-sided
/// *p*-value Pr[corank ≥ its corank] under the randomness hypothesis. With
/// `reps` > 1 the test is simply repeated on disjoint stretches of the orbit.
fn run_rank(field: &Field, n: usize, reps: usize, mut rng: Prng, log_interval: Duration) {
    println!(
        "Running a modular rank test using {} of RAM: {reps} {n}×{n} {} over the field of size {}",
        format_memory(n * n * size_of::<u64>()),
        pluralize("matrix", reps as isize, false),
        field.p(),
    );

    let mut data = vec![0; n * n].into_boxed_slice();

    for i in 0..reps {
        // Fill this matrix from the orbit (disjoint, contiguous outputs). When the
        // whole job is a single large matrix the fill is itself slow, so give it a
        // progress line.
        if reps == 1 {
            let mut gpl = new_logger("output", n * n, log_interval);
            gpl.start("Generating matrix entries...");
            for row in data.chunks_mut(n) {
                fill(field, &mut rng, row);
                gpl.update_with_count(row.len());
            }
            gpl.done();
        } else {
            fill(field, &mut rng, &mut data);
        }
        let mut pl = new_logger("column", n, log_interval);
        pl.start(format!(
            "Matrix {}/{reps}: ranking (blocked Gaussian elimination over Fₚ)...",
            i + 1
        ));
        let r = fp::rank(field, &mut data, n, &mut pl);
        pl.done();

        let corank = n - r;
        let p = corank_p_value(field.p(), n, corank);
        println!(
            "Matrix {}/{reps}\tcorank={corank}\tp={}",
            i + 1,
            pretty_p_value(p)
        );
    }
}

/// The name of the algorithm used by [`fp::linear_complexity`] for sequences of
/// length `n`.
const fn linear_complexity_algorithm(n: usize) -> &'static str {
    if n < fp::HGCD_THRESHOLD {
        "Berlekamp–Massey"
    } else {
        "half-GCD"
    }
}

/// Modular linear-complexity test (see [`fp::linear_complexity`]) over `reps`
/// independent length-*n* sequences, one at a time, each under its own progress
/// line, printing for each its one-sided *p*-value Pr[*Lₙ* ≤ its complexity]
/// under the randomness hypothesis. With `reps` > 1 the test is simply repeated
/// on disjoint stretches of the orbit.
fn run_linear_complexity(
    field: &Field,
    n: usize,
    reps: usize,
    mut rng: Prng,
    log_interval: Duration,
) {
    println!(
        "Running a modular linear-complexity test using {} of RAM: {reps} {} of length {n} \
         over the field of size {}",
        // The sequence plus the memory of the algorithm.
        format_memory(n * size_of::<u64>() + fp::linear_complexity_memory(n)),
        pluralize("sequence", reps as isize, false),
        field.p(),
    );

    let mut seq = vec![0; n].into_boxed_slice();

    for t in 0..reps {
        fill(field, &mut rng, &mut seq);
        let mut pl = new_logger("step", n, log_interval);
        pl.start(format!(
            "Sequence {}/{reps}: {} over Fₚ...",
            t + 1,
            linear_complexity_algorithm(n)
        ));
        let lc = fp::linear_complexity(field, &seq, &mut pl);
        pl.done();

        let p = linear_complexity_p_value(field.p(), n, lc);
        println!(
            "Sequence {}/{reps}\tlinear complexity={lc}\tp={}",
            t + 1,
            pretty_p_value(p)
        );
    }
}
