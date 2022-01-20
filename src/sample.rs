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

use rand::rngs::ThreadRng;
use rand::{thread_rng, Rng};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time;

use crate::histo::Histo;
use crate::options::SamplingAlgorithm;
use crate::state::State;
use crate::transition::{ELTransition, PNTransition, Transition};
use crate::FredholmConfig;

#[derive(Debug)]
pub struct Sample(pub Vec<State>);

#[derive(Debug)]
pub struct Samples(pub Vec<Sample>);

pub fn run(config: Arc<FredholmConfig>) {
    let now = time::Instant::now();
    let transition: Arc<dyn Transition<ThreadRng>> = match config.sampling_algorithm {
        SamplingAlgorithm::PositiveNegativeDecomposition => {
            Arc::new(PNTransition::new(&config.params))
        }
        SamplingAlgorithm::ExpLinearDecomposition => Arc::new(ELTransition::new(&config.params)),
    };

    let mut histo = Arc::new(Mutex::new(Histo::new(0.0, config.cutoff, config.bins)));

    let mut handles = Vec::new();
    for i in 0..config.threads {
        println!("spawning thread n. {}...", i);
        let config = Arc::clone(&config);
        let transition = Arc::clone(&transition);
        let histo = Arc::clone(&mut histo);
        let handle = thread::spawn(move || {
            let mut rng = thread_rng();
            let n_histories = config.replicas / config.threads;
            let n_packets = n_histories / config.packet_size;
            for _ in 0..n_packets {
                let samples = generate_samples(
                    &transition,
                    config.packet_size,
                    config.cutoff,
                    config.roulette_threshold,
                    &mut rng,
                );

                if config.debug {
                    println!("{:#?}", samples);
                }

                let mut histo = histo.lock().unwrap();
                histo.score(samples.0.iter());

                // let length = average_length(samples.0.iter(), samples.0.len());
                // println!("[{}] average history length: {}", i, length);
            }
        });
        handles.push(handle);
    }

    println!("joining...");

    for handle in handles {
        handle.join().unwrap();
    }

    println!("joined!");

    histo.lock().unwrap().write(&config.output);

    let elapsed = now.elapsed().as_millis();
    let elapsed_s = 1e-3 * elapsed as f64;
    println!(
        "elapsed: {} s ({} us / replica)",
        elapsed_s,
        1000.0 * elapsed as f64 / config.replicas as f64,
    );
}

fn generate_samples<T>(
    transition: &Arc<dyn Transition<T>>,
    replicas: usize,
    cutoff: f64,
    roulette_threshold: f64,
    rng: &mut T,
) -> Samples
where
    T: Rng,
{
    let vec: Vec<Sample> = (0..replicas)
        .map(|_| generate_sample(transition, cutoff, roulette_threshold, rng))
        .collect();
    Samples(vec)
}

fn generate_sample<T>(
    transition: &Arc<dyn Transition<T>>,
    cutoff: f64,
    roulette_threshold: f64,
    rng: &mut T,
) -> Sample
where
    T: Rng,
{
    let mut states = Vec::new();

    let mut state = Some(State {
        position: 0.0,
        weight: 1.0,
    });

    while let Some(some_state) = state {
        state = transition
            .sample(&some_state, cutoff, rng)
            .map(|state| {
                states.push(state);
                state
            })
            .and_then(|state| roulette(state, roulette_threshold, rng));
    }
    Sample(states)
}

fn roulette<T>(mut state: State, threshold: f64, rng: &mut T) -> Option<State>
where
    T: Rng,
{
    if state.weight.abs() > threshold {
        return Some(state);
    }
    let xi: f64 = rng.gen();
    if xi > state.weight.abs() {
        None
    } else {
        state.weight = 1.0_f64.copysign(state.weight);
        Some(state)
    }
}

#[allow(dead_code)]
fn average_length<'a, T>(iter: T, len: usize) -> f64
where
    T: Iterator<Item = &'a Sample>,
{
    let average = iter.fold(0, |acc, sample| acc + sample.0.len());
    average as f64 / len as f64
}
