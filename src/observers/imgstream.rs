/// Handles image rendering and IO
use na::DMatrix;
use plotters::prelude::*;
use tracing::debug;
use std::{error::Error, fs, path::Path, sync::mpsc};

/// A display packet that streams from the solver thread to this IO thread
#[derive(Clone)]
pub struct DisplayPacket {
    pub velocity_x: DMatrix<f32>,
    pub velocity_y: DMatrix<f32>,
    pub i: usize,
}

/// Saves an bitmap (velocity) with a blacked out mask as a png image
///
/// Parameters
/// - `bitmap` - The bitmap to be shown as the color gradient
/// - `solid_bitmap` - The bitmap to be shown as a black solid
/// - `filename` - The name of the file to save
pub fn image_save(
    bitmap: &DMatrix<f32>,
    solid_bitmap: &DMatrix<bool>,
    filename: &Path,
) -> Result<(), Box<dyn Error>> {
    let bitmap = bitmap.normalize();
    let (rows, cols) = bitmap.shape();

    let root = BitMapBackend::new(&filename, (cols as u32, rows as u32)).into_drawing_area();
    root.fill(&WHITE)?;

    for i in 0..rows {
        for j in 0..cols {
            let pixel_mag = bitmap.get((i, j)).ok_or("Pixel not on velocity field")?;
            let pixel_intensity = (254.0 * pixel_mag).floor() as u8;

            let pixel_color: RGBColor = if *solid_bitmap.index((i, j)) {
                RGBColor(0, 0, 0)
            } else {
                let t = pixel_intensity as f32 / 255.0;
                RGBColor(
                    (if t < 0.5 { (0.5 + t) * 255.0 } else { 255.0 }) as u8,
                    (if t < 0.5 {
                        (0.5 + t) * 255.0
                    } else {
                        (1.5 - t) * 255.0
                    }) as u8,
                    (if t < 0.5 { 255.0 } else { (1.5 - t) * 255.0 }) as u8,
                )
            };

            root.draw_pixel((j as i32, i as i32), &pixel_color)?;
        }
    }
    root.present().expect("Failed to present bitmap");

    Ok(())
}

/// Runs the image IO task
///
/// Parameters
/// - `inbound_bitmaps` - mpsc receiver for inbound DisplayPackets
/// - `solid_mask` - The solid mask that is assumed to be constant between all inbount bitmaps
/// - `frames_dir` - The directory that frames should be saved to
pub fn image_io_loop(
    inbound_bitmaps: mpsc::Receiver<DisplayPacket>,
    solid_mask: DMatrix<bool>,
    frames_dir: &Path,
) -> Result<(), Box<dyn Error>> {
    if (frames_dir).exists() {
        fs::remove_dir_all(frames_dir)?;
    }
    fs::create_dir(frames_dir)?;

    loop {
        // read inbound packet
        let inbound = match inbound_bitmaps
            .recv() {
                Ok(i) => i,
                Err(_) => {
                    debug!("Closing image IO loop");
                    break;
                }
            };

        

        let velocity_magnitude = (inbound.velocity_x.map(|x| x.powi(2))
            + inbound.velocity_y.map(|y| y.powi(2)))
        .map(|k| k.sqrt());

        image_save(
            &velocity_magnitude,
            &solid_mask,
            frames_dir.join(format!("{}.png", inbound.i)).as_path(),
        )
        .expect("Image save failed");
    }

    Ok(())
}
