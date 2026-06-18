/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

use clap::Parser;
use dsi_progress_logger::prelude::*;
use modlin::cli::{self, Args};
use modlin::fp::{self, Field};
use modlin::prng::Prng;
use modlin::stats::{corank_tail_pvalue, lc_left_tail_pvalue, pretty_p_value};
use pluralizer::pluralize;
use std::time::Duration;

fn main() {
    let args = Args::parse();
    args.validate();
    cli::init_env_logger();

    println!("Generator: {}", Prng::NAME);

    let seed = args.seed;
    println!("Seed: {:#018x}", seed);

    let reps = args.reps;
    let field = Field::new(args.modulus);

    if let Some(length) = args.linear_complexity {
        run_linear_complexity(&field, args.modulus, length, reps, seed, args.log_interval);
    } else if let Some(side) = args.rank {
        run_rank(&field, args.modulus, side, reps, seed, args.log_interval);
    } else {
        unreachable!();
    }
}

/// Modular rank test: rank `reps` independent *n* × *n* matrices, one at a time,
/// each under its own progress line, printing for each its one-sided *p*-value
/// Pr[corank ≥ its corank] under the null. With `reps` > 1 the test is simply
/// repeated on disjoint stretches of the orbit.
fn run_rank(field: &Field, modulus: u64, n: usize, reps: usize, seed: u64, log_interval: Duration) {
    println!(
        "Running a modular rank test: {reps} {n}×{n} {} over the field of size {modulus}",
        pluralize("matrix", reps as isize, false)
    );

    let mut data = vec![0; n * n];
    let mut rng = Prng::new(seed);

    for i in 0..reps {
        // Fill this matrix from the orbit (disjoint, contiguous outputs). When the
        // whole job is a single large matrix the fill is itself slow, so give it a
        // progress line.
        if reps == 1 {
            let mut gpl = progress_logger![
                item_name = "output",
                display_memory = true,
                log_interval = log_interval,
                expected_updates = Some(n * n),
            ];
            gpl.start("Generating matrix entries...");
            for row in data.chunks_mut(n) {
                for x in row.iter_mut() {
                    *x = field.reduce(rng.next_u64());
                }
                gpl.update_with_count(row.len());
            }
            gpl.done();
        } else {
            for x in data.iter_mut() {
                *x = field.reduce(rng.next_u64());
            }
        }
        let mut pl = progress_logger![
            item_name = "column",
            display_memory = true,
            log_interval = log_interval,
            expected_updates = Some(n),
        ];
        pl.start(format!(
            "Matrix {}/{reps}: ranking (blocked Gaussian elimination over Fₚ)...",
            i + 1
        ));
        let r = fp::rank(field, &mut data, n, &mut pl);
        pl.done();

        let corank = n - r;
        let p = corank_tail_pvalue(modulus as f64, n, corank);
        println!(
            "Matrix {}/{reps}\tcorank={corank}\tp={}",
            i + 1,
            pretty_p_value(p)
        );
    }
}

/// Modular linear-complexity test: Berlekamp–Massey over `reps` independent
/// length-*n* sequences, one at a time, each under its own progress line,
/// printing for each its one-sided *p*-value Pr[*Lₙ* ≤ its complexity] under the
/// null. With `reps` > 1 the test is simply repeated on disjoint stretches of
/// the orbit.
fn run_linear_complexity(
    field: &Field,
    modulus: u64,
    n: usize,
    reps: usize,
    seed: u64,
    log_interval: Duration,
) {
    println!(
        "Running a modular linear-complexity test: {reps} {} of length {n} \
         over the field of size {modulus}",
        pluralize("sequence", reps as isize, false),
    );

    let mut rng = Prng::new(seed);
    let mut seq = vec![0; n];

    for t in 0..reps {
        for x in seq.iter_mut() {
            *x = field.reduce(rng.next_u64());
        }
        let mut pl = progress_logger![
            item_name = "step",
            display_memory = true,
            log_interval = log_interval,
            expected_updates = Some(n),
        ];
        pl.start(format!(
            "Sequence {}/{reps}: Berlekamp–Massey over Fₚ...",
            t + 1
        ));
        let lc = fp::linear_complexity(field, &seq, &mut pl);
        pl.done();

        let p = lc_left_tail_pvalue(modulus as f64, n, lc);
        println!(
            "Sequence {}/{reps}\tlinear complexity={lc}\tp={}",
            t + 1,
            pretty_p_value(p)
        );
    }
}
