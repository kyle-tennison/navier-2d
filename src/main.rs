extern crate nalgebra as na;

mod observers;
mod postprocessing;
mod preprocessing;
mod sim;

use clap::Parser;
use na::DMatrix;
use tracing::Level;
use tracing_subscriber::{self, fmt::format::FmtSpan};

use crate::{
    postprocessing::postprocess, preprocessing::cli::CliArgs, sim::task::spawn_sim_thread,
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

    // run sim
    let sim_thread = spawn_sim_thread(sim_input.clone());
    let sim_output = sim_thread.join().expect("simulation thread panicked");

    // post process
    postprocess(sim_input, sim_output);
}
