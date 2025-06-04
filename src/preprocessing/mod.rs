use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::preprocessing::serial_mask::SerialMask;

pub mod cli;
pub mod preprocessor;
pub mod serial_mask;

#[derive(Serialize, Deserialize)]
pub struct ImageStreamSettings {
    pub frames_dir: PathBuf,
    pub retain_frames: bool,
    pub display_video: bool,
}

#[derive(Serialize, Deserialize)]
pub struct WebStreamSettings {}

#[derive(Serialize, Deserialize)]
pub enum InterfaceMode {
    ImageStream(ImageStreamSettings),
    WebStream(WebStreamSettings),
}

#[derive(Serialize, Deserialize)]
pub struct SimulationInput {
    pub mode: InterfaceMode,
    pub mask: Option<SerialMask>, // may be None when saved, but needs to be loaded to solve
    pub simulation_time: f32,
    pub inflow: (f32, f32),
    pub density: f32,
    pub viscosity: f32,
}
