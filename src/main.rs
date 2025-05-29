use std::{
    path::Path,
    sync::{LazyLock, mpsc},
    thread,
    time::Duration,
};

extern crate nalgebra as na;

mod display;
mod numeric;
mod poission;
mod sim;

use display::DisplayPacket;
use indicatif::ProgressIterator;
use na::{DMatrix, dmatrix};
use num_traits::ToPrimitive;

type ScalarField = DMatrix<f32>;
type VectorField = [ScalarField; 2];

static FRAMES_PATH: LazyLock<&Path> = LazyLock::new(|| Path::new("sim-frames"));

fn main() {
    let (sender, receiver) = mpsc::channel();

    thread::spawn(|| {
        display::image_io_loop(receiver).unwrap();
    });

    let simtime = 15.00;

    let sim = sim::NewtonianSim::new(
        1.,
        0.002,
        (8., 0.),
        &sim::NewtonianSim::sample_shape_mask(150, 150),
        (5., 5.),
        0.0025,
        simtime,
    );

    let iter_count = sim.iter_count();
    for (i, u) in sim.enumerate().progress_count(iter_count) {
        let (ux, uy) = (u[0].to_owned(), u[1].to_owned());

        sender
            .send(DisplayPacket {
                velocity_x: ux,
                velocity_y: uy,
                i,
            })
            .unwrap();
    }

    display::play_video((iter_count as f32 / simtime).floor() as usize, *FRAMES_PATH).unwrap();
}
