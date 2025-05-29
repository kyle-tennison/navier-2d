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

    let vel_x: DMatrix<f32> = dmatrix![
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
        0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0;
        0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0;
        0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0;
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    ];

    let vel_y: DMatrix<f32> = vel_x.clone();

    let mut display_packet = DisplayPacket {
        velocity_x: vel_x,
        velocity_y: vel_y,
        i: 0,
    };

    sender.send(display_packet.clone()).unwrap();

    for _ in 0..30 {
        let new_packet = DisplayPacket {
            velocity_x: display_packet.velocity_x.add_scalar(0.1),
            velocity_y: display_packet.velocity_y.add_scalar(0.1),
            i: display_packet.i + 1,
        };

        sender.send(new_packet.clone()).unwrap();
        display_packet = new_packet;
    }

    thread::sleep(Duration::from_secs(1));

    display::play_video(10, *FRAMES_PATH).unwrap();
}
