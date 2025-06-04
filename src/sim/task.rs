use std::{
    sync::mpsc,
    thread::{self, JoinHandle},
};

use indicatif::{ProgressBar, ProgressStyle};

use crate::{
    observers::imgstream::{self, DisplayPacket},
    preprocessing::{ImageStreamSettings, InterfaceMode, SimulationInput},
    sim::sim::NewtonianSim,
};

pub struct TaskOutput {
    pub total_iter: usize,
}

pub fn imgstream_task(
    settings: &ImageStreamSettings,
    sim: NewtonianSim,
    simulation_input: &SimulationInput,
) -> TaskOutput {
    let bar = ProgressBar::new(10_000);
    bar.set_style(
        ProgressStyle::with_template(
            "[Elapsed: {elapsed_precise}] [{bar:40.cyan/blue}] {percent}% (Remaining: {eta_precise})"
        )
        .unwrap()
        .progress_chars("##-"),
    );

    let (sender, receiver) = mpsc::channel();

    // spawn image io thread
    let settings_clone = (*settings).clone();
    thread::spawn(move || {
        imgstream::image_io_loop(receiver, &settings_clone.frames_dir).unwrap();
    });

    let mut total_iter = 0;
    for (i, (u, t)) in sim.enumerate() {
        let (ux, uy) = (u[0].to_owned(), u[1].to_owned());

        sender
            .send(DisplayPacket {
                velocity_x: ux,
                velocity_y: uy,
                i,
            })
            .unwrap();

        let progress = ((t / simulation_input.simulation_time) * 10_000.0).round() as u64;
        bar.set_position(progress);
        total_iter += 1;
    }

    return TaskOutput { total_iter };
}

pub fn spawn_sim_thread(simulation_input: SimulationInput) -> JoinHandle<TaskOutput> {
    thread::spawn(move || {
        let sim = NewtonianSim::new(
            simulation_input.density,
            simulation_input.viscosity,
            simulation_input.inflow,
            &simulation_input.get_mask(),
            simulation_input.length,
            simulation_input.simulation_time,
            simulation_input.cfl,
        );

        match &simulation_input.mode {
            InterfaceMode::ImageStream(settings) => {
                imgstream_task(settings, sim, &simulation_input)
            }
            InterfaceMode::WebStream(_) => {
                panic!("not implemented");
            }
        }
    })
}
