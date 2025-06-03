use std::{
    path::Path,
    sync::{LazyLock, mpsc},
    thread,
};

extern crate nalgebra as na;

mod display;
mod numeric;
mod poission;
mod sim;

use display::DisplayPacket;
use indicatif::{ProgressBar, ProgressStyle};
use na::DMatrix;
use num_traits::Zero;

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
        &sim::NewtonianSim::sample_shape_mask(200, 200),
        (5., 5.),
        simtime,
        1.,
    );

    // let mut pbar = ProgressBar::new((simtime*100.).floor() as u64);
    let mut iter_count = 0; // note, this will be nonlinear-timing rn
    let mut t_outer = 0.;

    let bar = ProgressBar::new(10_000);
    bar.set_style(
        ProgressStyle::with_template(
            "[Elapsed: {elapsed_precise}] [{bar:40.cyan/blue}] {percent}% (Remaining: {eta_precise})"
        )
        .unwrap()
        .progress_chars("##-"),
    );

    for (i, (u, t)) in sim.enumerate() {
        let (ux, uy) = (u[0].to_owned(), u[1].to_owned());

        sender
            .send(DisplayPacket {
                velocity_x: ux,
                velocity_y: uy,
                i,
            })
            .unwrap();

        iter_count += 1;
        t_outer = t;

        let progress = ((t / simtime) * 10_000.0).round() as u64;
        bar.set_position(progress);
    }

    let mut fps = (iter_count as f32 / simtime).floor() as usize;

    fps.is_zero().then(|| fps += 1);

    println!("fps: {}; t={}", fps, t_outer);

    display::play_video(fps, *FRAMES_PATH).unwrap();
}
