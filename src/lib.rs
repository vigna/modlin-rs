/*
 * SPDX-FileCopyrightText: 2026 Sebastiano Vigna
 *
 * SPDX-License-Identifier: Apache-2.0 OR LGPL-2.1-or-later
 */

#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

#[doc(hidden)]
pub mod cli;
pub mod fp;
#[allow(dead_code)]
pub mod prng;
pub mod stats;
