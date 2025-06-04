use std::{
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
    process::exit,
    sync::LazyLock,
};

use clap::{Parser, command};
use na::DMatrix;
use tracing::{error, info};

use crate::preprocessing::{
    ImageStreamSettings, InterfaceMode, SimulationInput, WebStreamSettings,
    preprocessor::mask_from_image, serial_mask::SerialMask,
};

static DEFAULT_FRAMES_PATH: LazyLock<&Path> = LazyLock::new(|| Path::new("sim-frames2"));

// Raw, CLI input
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct CliArgs {

    #[arg(help = "The path to a PNG image to use as the solid object.")]
    mask_path: Option<PathBuf>,

    #[arg(long, help = "An input file with pre-loaded parameters.")]
    input_json: Option<PathBuf>,

    #[arg(
        long,
        help = "The mode to run the simulation in: `stream` or `video`",
        default_value = "video"
    )]
    mode: String,

    #[arg(long, help = "Optional directiory to save input file to.")]
    input_json_savepath: Option<PathBuf>,

    #[arg(
        long,
        help = "An optional directory pointing to where frames should be saved."
    )]
    frames_dir: Option<PathBuf>,

    #[arg(
        short,
        long,
        help = "Whether or not frames should be retained after saving. Always false if streaming.",
        default_value = "false"
    )]
    retain_frames: bool,

    #[arg(
        short,
        long,
        help = "Whether the frame animation should play after solving. Disabled if streaming."
    )]
    display_video: bool,

    #[arg(long, help = "Inflow x velocity.", default_value = "8.0")]
    inflow_x: f32,

    #[arg(long, help = "Inflow y velocity.", default_value = "0.0")]
    inflow_y: f32,

    #[arg(long, help = "Domain length in x axis.", default_value = "5.0")]
    length_x: f32,

    #[arg(long, help = "Domain length in y axis.", default_value = "5.0")]
    length_y: f32,

    #[arg(
        short,
        long,
        default_value = "10",
        help = "Simulation time in seconds."
    )]
    simtime: f32,

    #[arg(long, default_value = "1.0", help = "Fluid density in g/ml")]
    density: f32,

    #[arg(long, default_value = "0.002", help = "Shear viscosity")]
    viscosity: f32,

    #[arg(long, default_value = "1.0", help = "Maximum CFL")]
    cfl: f32,
}

impl CliArgs {
    pub fn crate_input(&self) -> SimulationInput {
        // if the input file is supplied, just use that
        if let Some(input_filepath) = &self.input_json {
            if !input_filepath.exists() {
                error!("Input file {:?} does not exist.", input_filepath)
            }
            if input_filepath.is_dir() {
                error!("Input file {:?} is a directory.", input_filepath)
            }

            info!(
                "Using input file {}",
                input_filepath.to_str().unwrap_or("<unknown>")
            );

            let input_file = File::open(input_filepath)
                .inspect_err(|err| {
                    error!("Failed to open input file: {:?}", err);
                    exit(1);
                })
                .unwrap();

            let reader = BufReader::new(input_file);
            let mut loaded_input: SimulationInput = serde_json::from_reader(reader)
                .inspect_err(|err| error!("Failed to deserialize input file: {:?}", err))
                .unwrap();

            if loaded_input.mask.is_none() {
                if let Some(mask_path) = &self.mask_path {
                    let mask: DMatrix<bool> = mask_from_image(&mask_path)
                        .inspect_err(|f| {
                            error!("{:?}", f);
                            exit(1);
                        })
                        .unwrap();
                    loaded_input.mask = Some(SerialMask::from_mask(&mask));
                } else {
                    error!(
                        "The provided input file does not define a mask; this must be provided with the `--mask-path <png>` argument."
                    );
                    exit(1);
                }
            }

            return loaded_input;
        }

        // otherwise, build the input from the other arguments
        let mode = match self.mode.as_str() {
            "video" => {
                let frames_dir = self
                    .frames_dir
                    .as_ref()
                    .map(PathBuf::from)
                    .unwrap_or((*DEFAULT_FRAMES_PATH).into());

                InterfaceMode::ImageStream(ImageStreamSettings {
                    frames_dir,
                    retain_frames: self.retain_frames,
                    display_video: self.display_video,
                })
            }
            "stream" => InterfaceMode::WebStream(WebStreamSettings {}),
            _ => {
                error!(
                    "'{}' is not a valid interface mode. Use --help for info.",
                    self.mode
                );
                exit(1);
            }
        };

        let mask: SerialMask;
        if let Some(mask_path) = &self.mask_path {
            let na_mask: DMatrix<bool> = mask_from_image(&mask_path)
                .inspect_err(|f| {
                    error!("{:?}", f);
                    exit(1);
                })
                .unwrap();
            mask = SerialMask::from_mask(&na_mask);
        } else {
            error!("A png file depicting the solid object needs to be provided.");
            exit(1);
        }

        SimulationInput {
            mask: Some(mask),
            mode,
            simulation_time: self.simtime,
            inflow: (self.inflow_x, self.inflow_y),
            density: self.density,
            viscosity: self.viscosity,
            length: (self.length_x, self.length_y),
            cfl: self.cfl,
        }
    }
}
