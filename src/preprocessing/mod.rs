use std::{any, path::PathBuf};

use na::DMatrix;
use serde::{Deserialize, Serialize};
use tracing::info;

use crate::preprocessing::serial_mask::SerialMask;

pub mod cli;
pub mod preprocessor;
pub mod serial_mask;

#[derive(Serialize, Deserialize, Clone)]
pub struct ImageStreamSettings {
    pub frames_dir: PathBuf,
    pub retain_frames: bool,
    pub display_video: bool,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct WebStreamSettings {}

#[derive(Serialize, Deserialize, Clone)]
pub enum InterfaceMode {
    ImageStream(ImageStreamSettings),
    WebStream(WebStreamSettings),
}

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
    pub fn get_mask(&self) -> DMatrix<bool> {
        self.mask.as_ref().unwrap().to_mask()
    }
}

impl SimulationInput {
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
