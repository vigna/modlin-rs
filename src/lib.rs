/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR MIT
 */

#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]
#![warn(missing_docs)]

#[doc(hidden)]
pub mod cli;
pub mod fp;
mod ntt;
pub mod poly;
pub mod prng;
pub mod stats;
#[cfg(test)]
mod test_utils;
