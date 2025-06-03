use std::{
    path::{Path, PathBuf}, process::exit, sync::{mpsc, LazyLock}, thread
};

extern crate nalgebra as na;

mod display;
mod numeric;
mod poission;
mod sim;
mod preprocessor;

use clap::Parser;
use display::DisplayPacket;
use indicatif::{ProgressBar, ProgressStyle};
use na::DMatrix;
use num_traits::Zero;
use tracing::{Level, debug, error, info, warn};
use tracing_subscriber::{self, fmt::format::FmtSpan};

type ScalarField = DMatrix<f32>;
type VectorField = [ScalarField; 2];

static DEFAULT_FRAMES_PATH: LazyLock<&Path> = LazyLock::new(|| Path::new("sim-frames2"));

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[arg(short, long)]
    frames_dir: Option<String>,

    #[arg(short, long)]
    show_video: bool,

    #[arg(help = "The path to a PNG image to use as the mask. White pixels will permit fluid flow, black pixels will be treated as solid.")]
    mask_path: Option<String>
}

fn main() {
    // setup logging
    tracing_subscriber::fmt()
        .with_span_events(FmtSpan::ENTER | FmtSpan::EXIT)
        .with_ansi(true)
        .with_max_level(Level::DEBUG)
        .init();

    let args = Args::parse();

    let frames_path = args
        .frames_dir
        .and_then(|p| Some(PathBuf::from(p)))
        .unwrap_or((*DEFAULT_FRAMES_PATH).into());

    debug!("Using frames path {:?}", frames_path);
    debug!("Showing video? {:?}", args.show_video);

    let mask: DMatrix<bool>;

    if let Some(mask_path) = args.mask_path {
        info!("Loading mask path {}...", &mask_path);
        mask = match preprocessor::mask_from_image(PathBuf::from(&mask_path).as_path()){
            Ok(m) => m,
            Err(err) => {
                error!("Failed to load mask from provided image {}, because: {}", mask_path, err);
                exit(1);
            }
        };
        debug!("Sucessfully created mask from image");
    }
    else{
        warn!("Using example shape mask.");
        mask = sim::NewtonianSim::sample_shape_mask(200, 200);
    }



    let (sender, receiver) = mpsc::channel();

    let frames_path_tread = frames_path.clone();
    thread::spawn(move || {
        display::image_io_loop(receiver, &frames_path_tread).unwrap();
    });

    let simtime = 15.00;

    let sim = sim::NewtonianSim::new(
        1.,
        0.002,
        (8., 0.),
        &mask,
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

    display::play_video(fps, frames_path.as_path()).unwrap();
}
