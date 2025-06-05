/// Handles PNG simulation input

use image::{GenericImageView, ImageReader, Pixel};
use na::DMatrix;
use std::{error::Error, path::Path};

const THRESHOLD_LUMA: u8 = 127;

/// Load a DMatrix boolean mask from a PNG image by looking at pixel luminosity.
/// 
/// Parameters
/// - `image` - The path to the image to process
/// 
/// Returns
/// - The boolean mask as a Result
pub fn mask_from_image(image: &Path) -> Result<DMatrix<bool>, Box<dyn Error>> {
    let image = ImageReader::open(image)?.decode()?;

    let (nrows, ncols) = (image.height(), image.width());

    let mut mask: DMatrix<bool> = DMatrix::from_element(nrows as usize, ncols as usize, false);

    // load mask
    image.pixels().for_each(|(x, y, color)| {
        *(mask.index_mut((y as usize, x as usize))) = color.to_luma().0[0] < THRESHOLD_LUMA
    });

    Ok(mask)
}
