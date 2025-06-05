// Contains post-processers for analyzing simulation results

pub mod display;

use crate::{
    preprocessing::{InterfaceMode, SimulationInput},
    sim::task::SimulationOutput,
};
use std::fs;
use tracing::warn;

pub fn postprocess(sim_input: SimulationInput, sim_output: SimulationOutput) {
    if let InterfaceMode::ImageStream(settings) = sim_input.mode {
        if settings.display_video {
            display::play_video(
                sim_input.simulation_time,
                60,
                &sim_output.temporal_map,
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
