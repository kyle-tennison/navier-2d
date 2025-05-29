use std::{
    error::Error,
    fs,
    path::Path,
    sync::mpsc,
    time::{Duration, Instant},
};

use na::DMatrix;
use plotters::prelude::*;

use image::{DynamicImage, GenericImageView, imageops::FilterType};
use minifb::{Key, Window, WindowOptions};
use screen_size::get_primary_screen_size as get_screen_size;

use crate::FRAMES_PATH;

#[derive(Clone)]
pub struct DisplayPacket {
    pub velocity_x: DMatrix<f32>,
    pub velocity_y: DMatrix<f32>,
    pub i: usize,
}

pub fn image_save(bitmap: &DMatrix<f32>, filename: &str) -> Result<(), Box<dyn Error>> {
    let (rows, cols) = bitmap.shape();

    let filename = FRAMES_PATH.join(filename);

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

pub fn image_io_loop(inbound_bitmaps: mpsc::Receiver<DisplayPacket>) -> Result<(), Box<dyn Error>> {
    if (*FRAMES_PATH).exists() {
        fs::remove_dir_all(*FRAMES_PATH)?;
    }
    fs::create_dir(*FRAMES_PATH)?;

    loop {
        // read inbound packet
        let inbound = inbound_bitmaps
            .recv()
            .map_err(|err| panic!("Couldn't read inbound bitmap: {:?}", err))
            .unwrap();

        let velocity_magnitude = (inbound.velocity_x.map(|x| x.powi(2))
            + inbound.velocity_y.map(|y| y.powi(2)))
        .map(|k| k.sqrt());

        image_save(&velocity_magnitude, format!("{}.png", inbound.i).as_str())
            .expect("Image save failed");
    }
}

pub fn play_video(fps: usize, frames_dir: &Path) -> Result<(), Box<dyn Error>> {
    // collect & sort frame paths
    let mut paths: Vec<_> = fs::read_dir(frames_dir)?
        .filter_map(Result::ok)
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|e| e.to_str()) == Some("png"))
        .collect();
    paths.sort_by_key(|p| {
        p.file_stem()
            .and_then(|s| s.to_str())
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0)
    });
    if paths.is_empty() {
        return Err("no PNG frames found".into());
    }

    // load all frames as DynamicImage
    let originals: Vec<DynamicImage> = paths
        .iter()
        .map(|p| image::open(p))
        .collect::<Result<_, _>>()?;

    // determine base dimensions
    let (w, h) = originals[0].dimensions();
    let (screen_w, _) = get_screen_size().expect("failed to get screen size");
    let init_w = screen_w / 2;
    let init_h = (init_w as f32 * (h as f32 / w as f32)) as u32;

    // create window
    let mut window = Window::new(
        "Navier 2D",
        init_w as usize,
        init_h as usize,
        WindowOptions {
            resize: true,
            ..WindowOptions::default()
        },
    )
    .expect("Could not create window for video");

    let frame_time = Duration::from_secs_f64(1.0 / fps as f64);
    let start = Instant::now();

    while window.is_open() && !window.is_key_down(Key::Escape) {
        // get current window size
        let (win_w, win_h) = window.get_size();
        // current frame index
        let elapsed = Instant::now().duration_since(start);
        let idx = ((elapsed.as_secs_f64() * fps as f64) as usize) % originals.len();

        // resize & convert to RGBA buffer
        let img = originals[idx]
            .resize_exact(win_w as u32, win_h as u32, FilterType::Nearest)
            .to_rgba8();

        let buffer: Vec<u32> = img
            .pixels()
            .map(|px| {
                ((px[3] as u32) << 24)
                    | ((px[0] as u32) << 16)
                    | ((px[1] as u32) << 8)
                    | (px[2] as u32)
            })
            .collect();

        window
            .update_with_buffer(&buffer, win_w, win_h)
            .expect("Could not update window with buffer");
        // throttle to fps
        let next = start + frame_time * (idx + 1) as u32;
        if let Some(d) = next.checked_duration_since(Instant::now()) {
            std::thread::sleep(d);
        }
    }

    Ok(())
}
