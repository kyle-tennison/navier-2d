use std::path::PathBuf;

use na::DMatrix;
use serde::{Deserialize, Serialize};

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
    mask: Option<SerialMask>, // may be None when saved, but needs to be loaded to solve
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
