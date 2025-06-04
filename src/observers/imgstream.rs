use na::DMatrix;
use plotters::prelude::*;
use std::{error::Error, fs, path::Path, sync::mpsc};

#[derive(Clone)]
pub struct DisplayPacket {
    pub velocity_x: DMatrix<f32>,
    pub velocity_y: DMatrix<f32>,
    pub i: usize,
}

pub fn image_save(
    bitmap: &DMatrix<f32>,
    filename: &str,
    frames_dir: &Path,
) -> Result<(), Box<dyn Error>> {
    let (rows, cols) = bitmap.shape();

    let filename = frames_dir.join(filename);

    let root = BitMapBackend::new(&filename, (cols as u32, rows as u32)).into_drawing_area();
    root.fill(&WHITE)?;

    let bitmap = bitmap / bitmap.max();

    for i in 0..rows {
        for j in 0..cols {
            let pixel_mag = bitmap.get((i, j)).ok_or("Pixel not on velocity field")?;
            let pixel_intensity = (254.0 * pixel_mag).floor() as u8;
            let pixel_color = &RGBColor(pixel_intensity, pixel_intensity, pixel_intensity);

            root.draw_pixel((j as i32, i as i32), pixel_color)?;
        }
    }
    root.present().expect("Failed to present bitmap");

    Ok(())
}

pub fn image_io_loop(
    inbound_bitmaps: mpsc::Receiver<DisplayPacket>,
    frames_dir: &Path,
) -> Result<(), Box<dyn Error>> {
    if (frames_dir).exists() {
        fs::remove_dir_all(frames_dir)?;
    }
    fs::create_dir(frames_dir)?;

    loop {
        // read inbound packet
        let inbound = inbound_bitmaps
            .recv()
            .map_err(|err| panic!("Couldn't read inbound bitmap: {:?}", err))
            .unwrap();

        let velocity_magnitude = (inbound.velocity_x.map(|x| x.powi(2))
            + inbound.velocity_y.map(|y| y.powi(2)))
        .map(|k| k.sqrt());

        image_save(
            &velocity_magnitude,
            format!("{}.png", inbound.i).as_str(),
            frames_dir,
        )
        .expect("Image save failed");
    }
}
