use std::{
    fs, path::{Path, PathBuf}, process::exit, sync::{mpsc, LazyLock}, thread
};

extern crate nalgebra as na;

mod observers;
mod postprocessing;
mod preprocessing;
mod sim;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use na::DMatrix;
use num_traits::Zero;
use tracing::{Level, debug, error, info, warn};
use tracing_subscriber::{self, fmt::format::FmtSpan};

use crate::{
    postprocessing::display,
    preprocessing::{InterfaceMode, cli::CliArgs, preprocessor},
    sim::{sim::NewtonianSim, task::spawn_sim_thread},
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

    let sim_thread = spawn_sim_thread(sim_input.clone());

    match sim_input.mode {
        InterfaceMode::ImageStream(settings) => {
            let output = sim_thread.join().expect("Sim thread panicked");

            if settings.display_video {
                let mut fps =
                    (output.total_iter as f32 / sim_input.simulation_time).floor() as usize;
                fps.is_zero().then(|| fps += 1);
                display::play_video(fps, &settings.frames_dir).unwrap();
            }

            if !settings.retain_frames {
                _ = fs::remove_dir_all(settings.frames_dir).inspect_err(|err| warn!("Unable to cleanup frames output: {:?}", err));
            }
        }
        _ => (),
    };

    todo!();

    // let args = Args::parse();

    // let frames_path = args
    //     .frames_dir
    //     .map(PathBuf::from)
    //     .unwrap_or((*DEFAULT_FRAMES_PATH).into());

    // debug!("Using frames path {:?}", frames_path);
    // debug!("Showing video? {:?}", args.show_video);

    // let mask: DMatrix<bool>;

    // if let Some(mask_path) = args.mask_path {
    //     info!("Loading mask path {}...", &mask_path);
    //     mask = match preprocessor::mask_from_image(PathBuf::from(&mask_path).as_path()) {
    //         Ok(m) => m,
    //         Err(err) => {
    //             error!(
    //                 "Failed to load mask from provided image {}, because: {}",
    //                 mask_path, err
    //             );
    //             exit(1);
    //         }
    //     };
    //     debug!("Sucessfully created mask from image");
    // } else {
    //     warn!("Using example shape mask.");
    //     mask = NewtonianSim::sample_shape_mask(200, 200);
    // }

    // let (sender, receiver) = mpsc::channel();

    // let frames_path_tread = frames_path.clone();
    // thread::spawn(move || {
    //     imgstream::image_io_loop(receiver, &frames_path_tread).unwrap();
    // });

    // let simtime = args.simtime;

    // let sim = NewtonianSim::new(1., 0.002, (8., 0.), &mask, (5., 5.), simtime, 1.);

    // // let mut pbar = ProgressBar::new((simtime*100.).floor() as u64);
    // let mut iter_count = 0; // note, this will be nonlinear-timing rn
    // let mut t_outer = 0.;

    // let bar = ProgressBar::new(10_000);
    // bar.set_style(
    //     ProgressStyle::with_template(
    //         "[Elapsed: {elapsed_precise}] [{bar:40.cyan/blue}] {percent}% (Remaining: {eta_precise})"
    //     )
    //     .unwrap()
    //     .progress_chars("##-"),
    // );

    // for (i, (u, t)) in sim.enumerate() {
    //     let (ux, uy) = (u[0].to_owned(), u[1].to_owned());

    //     sender
    //         .send(DisplayPacket {
    //             velocity_x: ux,
    //             velocity_y: uy,
    //             i,
    //         })
    //         .unwrap();

    //     iter_count += 1;
    //     t_outer = t;

    //     let progress = ((t / simtime) * 10_000.0).round() as u64;
    //     bar.set_position(progress);
    // }
}
