use std::fs;

extern crate nalgebra as na;

mod observers;
mod postprocessing;
mod preprocessing;
mod sim;

use clap::Parser;
use na::DMatrix;
use num_traits::Zero;
use tracing::{Level, warn};
use tracing_subscriber::{self, fmt::format::FmtSpan};

use crate::{
    postprocessing::display,
    preprocessing::{InterfaceMode, cli::CliArgs},
    sim::task::spawn_sim_thread,
};

type ScalarField = DMatrix<f32>;
type VectorField = [ScalarField; 2];

fn main() {
    // setup logging
    tracing_subscriber::fmt()
        .with_span_events(FmtSpan::ENTER | FmtSpan::EXIT)
        .with_ansi(true)
        .with_max_level(Level::DEBUG)
        .init();

    // read args
    let args = CliArgs::parse();
    let sim_input = args.crate_input();
    sim_input.log();

    let sim_thread = spawn_sim_thread(sim_input.clone());

    if let InterfaceMode::ImageStream(settings) = sim_input.mode {
        let output = sim_thread.join().expect("Sim thread panicked");

        if settings.display_video {
            display::play_video(
                sim_input.simulation_time,
                60,
                &output.temporal_map,
                &settings.frames_dir,
            )
            .unwrap();
        }

        if !settings.retain_frames {
            _ = fs::remove_dir_all(settings.frames_dir)
                .inspect_err(|err| warn!("Unable to cleanup frames output: {:?}", err));
        }
    };
}
