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
mod options;
mod state;
mod transition;
use options::FredholmConfig;

fn main() {
    let config = Arc::new(FredholmConfig::from_args());
    run(config)
}
