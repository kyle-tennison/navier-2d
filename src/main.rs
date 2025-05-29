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
use na::{DMatrix, dmatrix};

type ScalarField = DMatrix<f32>;
type VectorField = [ScalarField; 2];

static FRAMES_PATH: LazyLock<&Path> = LazyLock::new(|| Path::new("sim-frames"));

fn main() {
    println!("Hello, world!");

    let (sender, receiver) = mpsc::channel();

    thread::spawn(|| {
        display::image_io_loop(receiver).unwrap();
    });

    let sim = sim::NewtonianSim::new(
        1.,
        0.002,
        (8., 0.),
        &sim::NewtonianSim::sample_shape_mask(50, 50),
        (5., 5.),
        0.0025,
        2.00,
    );

    for (i, u) in sim.enumerate() {
        let (ux, uy) = (u[0].to_owned(), u[1].to_owned());

        sender
            .send(DisplayPacket {
                velocity_x: ux,
                velocity_y: uy,
                i,
            })
            .unwrap();
    }

    display::play_video(10, *FRAMES_PATH).unwrap();
}
