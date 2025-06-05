// Handles preprocessing -- i.e. everything that happens prior to the simulation

use std::{any, path::PathBuf};

use na::DMatrix;
use serde::{Deserialize, Serialize};
use tracing::info;

use crate::preprocessing::serial_mask::SerialMask;

pub mod cli;
pub mod image_input;
pub mod serial_mask;

/// Available modes for the simulator to run in
#[derive(Serialize, Deserialize, Clone)]
pub enum InterfaceMode {
    ImageStream(ImageStreamSettings),
    WebStream(WebStreamSettings),
}

/// Settings when running in ImageStream mode
#[derive(Serialize, Deserialize, Clone)]
pub struct ImageStreamSettings {
    pub frames_dir: PathBuf,
    pub retain_frames: bool,
    pub display_video: bool,
}

/// Settings for when running in WebStream mode
#[derive(Serialize, Deserialize, Clone)]
pub struct WebStreamSettings {}

/// All the simulation parameters passed to the solver thread
#[derive(Serialize, Deserialize, Clone)]
pub struct SimulationInput {
    pub mode: InterfaceMode,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub mask: Option<SerialMask>, // may be None when saved, but needs to be loaded to solve
    pub simulation_time: f32,
    pub inflow: (f32, f32),
    pub density: f32,
    pub viscosity: f32,
    pub length: (f32, f32),
    pub cfl: f32,
}

impl SimulationInput {
    /// Convert the SerialMask attribute into a nalgebra DMatrix
    pub fn get_mask(&self) -> DMatrix<bool> {
        self.mask.as_ref().unwrap().to_mask()
    }
}

impl SimulationInput {
    /// Log the simulation input configuration over tracing.
    pub fn log(&self) {
        info!(
            "Simulation is shown below:\n\n\
        \t mode:       {}\n\
        \t time range: {} s\n\
        \t inflow:     < {} m/s, {} m/s >\n\
        \t density:    {} g/ml\n\
        \t viscosity:  {} mPa s \n\
        \t length:     <{} m, {} m > \n\
        \t CFL:        {} (max)\n\n\
        ",
            any::type_name_of_val(&self.mode),
            self.simulation_time,
            self.inflow.0,
            self.inflow.1,
            self.density,
            self.viscosity,
            self.length.0,
            self.length.1,
            self.cfl,
        );

        let mode_str = serde_json::to_string_pretty(&self.mode).unwrap();

        info!("Mode parameters are:\n\n{}", mode_str);
    }
}
