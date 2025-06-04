use std::{
    error::Error,
    fs,
    path::Path,
    time::{Duration, Instant},
};

use plotters::prelude::*;

use image::{DynamicImage, GenericImageView, imageops::FilterType};
use minifb::{Key, Window, WindowOptions};
use screen_size::get_primary_screen_size as get_screen_size;

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
    let originals: Vec<DynamicImage> = paths.iter().map(image::open).collect::<Result<_, _>>()?;

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
