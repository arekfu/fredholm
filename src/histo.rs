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

use super::sample::Sample;
use std::fmt::Write;
use std::fs;
use std::ops::{Add, AddAssign};

#[derive(Clone)]
pub struct Histo {
    /// the accumulated content of the histogram
    pub contents: Vec<f64>,
    /// the accumulated square of the histogram
    pub contents2: Vec<f64>,
    /// the minimum bound for the histogram
    pub min: f64,
    /// the width of the histogram
    pub width: f64,
    /// the norm of the histogram
    pub norm: usize,
}

impl Histo {
    pub fn new(min: f64, width: f64, n_bins: usize) -> Histo {
        Histo {
            contents: vec![0.0; n_bins],
            contents2: vec![0.0; n_bins],
            min,
            width,
            norm: 0,
        }
    }

    pub fn score<'a, T>(&mut self, iter: T)
    where
        T: Iterator<Item = &'a Sample>,
    {
        iter.for_each(|sample| {
            sample.0.iter().for_each(|state| {
                let bin = (self.contents.len() as f64) * (state.position - self.min) / self.width;
                let ubin = bin as usize;
                self.contents[ubin] += state.weight;
                self.contents2[ubin] += state.weight.powi(2);
            });
            self.norm += 1;
        })
    }

    pub fn write(&self, path: &str) {
        let flen = self.contents.len() as f64;
        let fnorm = self.norm as f64;
        let bin_width = self.width / flen;
        let data = (0..)
            .zip(self.contents.iter())
            .zip(self.contents2.iter())
            .fold(String::new(), |mut s, ((i, &f), &f2)| {
                let bin = self.min + (i as f64) * self.width / flen;
                let mut av = f / fnorm;
                let mut err = (1.0 / (self.norm - 1) as f64) * (f2 / fnorm - av.powi(2));
                err = err.sqrt() / bin_width;
                av /= bin_width;
                writeln!(s, "{} {} {}", bin, av, err).ok();
                s
            });
        fs::write(path, data).expect("Unable to write output file");
    }
}

impl AddAssign<&Histo> for Histo {
    fn add_assign(&mut self, rhs: &Histo) {
        assert_eq!(self.min, rhs.min);
        assert_eq!(self.width, rhs.width);
        self.contents
            .iter_mut()
            .zip(rhs.contents.iter())
            .for_each(|(this, other)| *this += *other);
        self.contents2
            .iter_mut()
            .zip(rhs.contents.iter())
            .for_each(|(this, other)| *this += *other);
        self.norm += rhs.norm;
    }
}

impl Add<&Histo> for &Histo {
    type Output = Histo;
    fn add(self, rhs: &Histo) -> Histo {
        let mut clone = self.clone();
        clone += rhs;
        clone
    }
}

#[cfg(test)]
mod tests {
    use super::super::sample::Samples;
    use super::super::state::State;
    use super::Histo;
    use super::Sample;

    #[test]
    pub fn one_bin() {
        let mut histo = Histo::new(0.0, 10.0, 1);
        assert_eq!(histo.min, 0.0);
        assert_eq!(histo.width, 10.0);
        assert_eq!(histo.contents.len(), 1);
        assert_eq!(histo.contents2.len(), 1);
        assert_eq!(histo.norm, 0);

        let states = vec![State {
            position: 5.0,
            weight: 1.0,
        }];
        let sample = Sample(states);
        let samples = Samples(vec![sample]);
        histo.score(samples.0.iter());
        assert_eq!(histo.norm, 1);
        assert_eq!(histo.contents, vec![1.0]);
        assert_eq!(histo.contents2, vec![1.0]);
    }

    #[test]
    pub fn weight() {
        let mut histo = Histo::new(0.0, 10.0, 1);
        let states = vec![State {
            position: 5.0,
            weight: 0.5,
        }];
        let sample = Sample(states);
        let samples = Samples(vec![sample]);
        histo.score(samples.0.iter());
        assert_eq!(histo.norm, 1);
        assert_eq!(histo.contents, vec![0.5]);
        assert_eq!(histo.contents2, vec![0.25]);
    }

    #[test]
    pub fn locate() {
        let mut histo = Histo::new(0.0, 5.0, 5);
        assert_eq!(histo.min, 0.0);
        assert_eq!(histo.width, 5.0);
        assert_eq!(histo.contents.len(), 5);
        assert_eq!(histo.contents2.len(), 5);
        assert_eq!(histo.norm, 0);

        let sample = Sample(vec![State {
            position: 2.3,
            weight: 1.0,
        }]);
        let samples = Samples(vec![sample]);
        histo.score(samples.0.iter());
        assert_eq!(histo.norm, 1);
        assert_eq!(histo.contents, vec![0.0, 0.0, 1.0, 0.0, 0.0]);
        assert_eq!(histo.contents2, vec![0.0, 0.0, 1.0, 0.0, 0.0]);
    }

    #[test]
    pub fn add() {
        let mut histo = Histo::new(0.0, 5.0, 5);

        let sample = Sample(vec![State {
            position: 2.3,
            weight: 1.0,
        }]);
        let samples = Samples(vec![sample]);
        histo.score(samples.0.iter());
        let histo2 = histo.clone();
        let histo_sum = &histo + &histo2;
        assert_eq!(histo_sum.norm, 2);
        assert_eq!(histo_sum.contents, vec![0.0, 0.0, 2.0, 0.0, 0.0]);
        assert_eq!(histo_sum.contents2, vec![0.0, 0.0, 2.0, 0.0, 0.0]);

        histo += &histo_sum;
        assert_eq!(histo.norm, 3);
        assert_eq!(histo.contents, vec![0.0, 0.0, 3.0, 0.0, 0.0]);
        assert_eq!(histo.contents2, vec![0.0, 0.0, 3.0, 0.0, 0.0]);
    }
}
