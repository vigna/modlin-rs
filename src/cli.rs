/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

//! Command-line argument definitions.

use clap::{ArgGroup, Parser};
use std::time::Duration;

/// Parses a duration: the suffixes s, m, h, d mean seconds, minutes, hours,
/// days; a number without a suffix is milliseconds. For example "1d2h3m4s567"
/// is 1 day + 2 hours + 3 minutes + 4 seconds + 567 ms.
fn parse_duration(value: &str) -> Result<Duration, String> {
    if value.is_empty() {
        return Err("empty duration (use \"0\" for a zero interval)".to_owned());
    }
    let mut duration = Duration::from_secs(0);
    let mut acc = String::new();
    for c in value.chars() {
        if c.is_ascii_digit() {
            acc.push(c);
        } else if c.is_whitespace() {
            continue;
        } else {
            let dur: u64 = acc
                .parse()
                .map_err(|_| format!("invalid number before '{c}'"))?;
            match c {
                's' => duration += Duration::from_secs(dur),
                'm' => duration += Duration::from_secs(dur * 60),
                'h' => duration += Duration::from_secs(dur * 60 * 60),
                'd' => duration += Duration::from_secs(dur * 60 * 60 * 24),
                _ => return Err(format!("invalid duration suffix: '{c}'")),
            }
            acc.clear();
        }
    }
    if !acc.is_empty() {
        let dur: u64 = acc
            .parse()
            .map_err(|_| "invalid trailing milliseconds".to_owned())?;
        duration += Duration::from_millis(dur);
    }
    Ok(duration)
}

/// Initializes env_logger with a custom format that prefixes each line with the
/// wall-clock time and the elapsed time since initialization. The default level
/// is info, overridable through RUST_LOG.
pub fn init_env_logger() {
    use jiff::{
        SpanRound,
        fmt::friendly::{Designator, Spacing, SpanPrinter},
    };
    use std::io::Write;

    let mut builder =
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"));

    let start = std::time::Instant::now();
    let printer = SpanPrinter::new()
        .spacing(Spacing::None)
        .designator(Designator::Compact);
    let span_round = SpanRound::new()
        .largest(jiff::Unit::Day)
        .smallest(jiff::Unit::Millisecond)
        .days_are_24_hours();

    builder.format(move |buf, record| {
        let Ok(ts) = jiff::Timestamp::try_from(std::time::SystemTime::now()) else {
            return Err(std::io::Error::other("Failed to get timestamp"));
        };
        let style = buf.default_level_style(record.level());
        let elapsed = start.elapsed();
        let span = jiff::Span::new()
            .seconds(elapsed.as_secs() as i64)
            .milliseconds(elapsed.subsec_millis() as i64);
        let span = span.round(span_round).expect("Failed to round span");
        writeln!(
            buf,
            "{} {} {style}{}{style:#} [{:?}] {} - {}",
            ts.strftime("%F %T%.3f"),
            printer.span_to_string(&span),
            record.level(),
            std::thread::current().id(),
            record.target(),
            record.args()
        )
    });
    builder.init();
}

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Runs a modular rank or linear-complexity test. Please use RAYON_NUM_THREADS to customize the number of threads, and RUST_LOG to customize the logging.",
    next_line_help = true,
    max_term_width = 100
)]
#[command(group(ArgGroup::new("test").required(true).args(["rank", "linear_complexity"])))]
pub struct Args {
    /// Run the modular rank test on SIDE×SIDE matrices; mutually exclusive with -L.​
    #[arg(short = 'R', long, value_name = "SIDE")]
    pub rank: Option<usize>,

    /// Run the modular linear-complexity test on sequences of the given length; mutually exclusive with -R.​
    #[arg(short = 'L', long, value_name = "LENGTH")]
    pub linear_complexity: Option<usize>,

    /// Number of matrices or sequences to test.​
    #[arg(default_value_t = 1)]
    pub reps: usize,

    /// Prime field modulus p: a prime below 2⁶³.​
    #[arg(short = 'p', long)]
    pub modulus: u64,

    /// PRNG seed.​
    #[arg(short = 'S', long, default_value_t = 0)]
    pub seed: u64,

    /// How often to log progress. Suffixes: s (seconds), m
    /// (minutes), h (hours), d (days); without a suffix, the value is
    /// milliseconds. Example: "1d2h3m4s567".​
    #[arg(long, value_parser = parse_duration, default_value = "10s")]
    pub log_interval: Duration,
}

impl Args {
    /// The test size: the matrix side (rank test) or the sequence length
    /// (linear-complexity test). The clap group guarantees exactly one is set.
    pub fn size(&self) -> usize {
        self.rank
            .or(self.linear_complexity)
            .expect("clap requires exactly one of -R/-L")
    }

    pub fn validate(&self) {
        let size = self.size();
        if size < 2 {
            Self::die("SIDE/LENGTH must be at least 2");
        }
        if self.rank.is_some() && size.checked_mul(size).is_none() {
            Self::die("SIDE is too large: SIDE² overflows a 64-bit index");
        }
        if self.reps < 1 {
            Self::die("reps must be at least 1");
        }
        if self.modulus < 2 {
            Self::die("modulus must be at least 2");
        }
        if self.modulus >= (1 << 63) {
            Self::die("modulus must be below 2⁶³");
        }
        if !primal::is_prime(self.modulus) {
            Self::die("modulus must be prime");
        }
    }

    /// Reports an argument error in clap's own style and exits with status 2.
    pub fn die(msg: &str) -> ! {
        use clap::CommandFactory;
        Args::command()
            .error(clap::error::ErrorKind::ValueValidation, msg)
            .exit()
    }
}
