// Copyright 2021 Davide Mancusi <davide.mancusi@cea.fr>
//
// This file is part of fredholm.
//
// fredholm is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// fredholm is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// fredholm.  If not, see <https://www.gnu.org/licenses/>.

use std::sync::Arc;
use structopt::StructOpt;

mod sample;
use sample::run;

mod histo;
mod state;
mod transition;

/// Parameters for the Monte Carlo simulation
#[derive(StructOpt)]
pub struct Params {
    /// the coefficient for the linear term in the kernel
    #[structopt(long, default_value = "1.0")]
    alpha: f64,

    /// the cross section
    #[structopt(long, default_value = "1.0")]
    sigma: f64,

    /// the absorption probability
    #[structopt(long, default_value = "0.9")]
    lambda: f64,
}

/// Solve a Fredholm equation of the second kind with Monte Carlo
#[derive(StructOpt)]
pub struct FredholmConfig {
    /// the number of replicas
    #[structopt(short, long, default_value = "1000")]
    replicas: usize,

    /// the maximum packet size
    #[structopt(short, long, default_value = "1000")]
    packet_size: usize,

    /// the cutoff for the coordinate
    #[structopt(long, default_value = "10.0")]
    cutoff: f64,

    /// the Russian roulette threshold
    #[structopt(long, default_value = "0.5")]
    roulette_threshold: f64,

    /// the number of bins for the output histogram
    #[structopt(long, default_value = "100")]
    bins: usize,

    /// the output file
    #[structopt(short, long, default_value = "fredholm.dat")]
    output: String,

    /// debug flag
    #[structopt(short, long)]
    debug: bool,

    /// number of threads
    #[structopt(short, long, default_value = "1")]
    threads: usize,

    #[structopt(flatten)]
    params: Params,
}

fn main() {
    let config = Arc::new(FredholmConfig::from_args());
    run(config)
}
