// Numeric approximations

use std::mem::MaybeUninit;

use na::{DMatrix, Dim, Dyn, Matrix, RawStorage};

type ScalarField = DMatrix<f32>;
type VectorField = [ScalarField; 2];

/// Compute the divergence of some vector field F=<u,v>.
///
/// Mathematically, this is du/dx + dv/dy
///
/// Parameters:
/// - `field` - The `VectorField` to take the divergence of
/// - `dy` - The y-axis step size
/// - `dx` - THe x-axis step size
///
/// Returns:
///     A `ScalarField` of the divergence.
pub fn divergence(field: VectorField, dy: f32, dx: f32) -> ScalarField {
    let (rows, cols) = field[0].shape();

    let mut du_dx: DMatrix<f32> = DMatrix::zeros(rows, cols);
    let mut dv_dy: DMatrix<f32> = DMatrix::zeros(rows, cols);

    // -- x-axis --

    // set interior nodes
    for r in 0..rows {
        for c in 1..(cols - 1) {
            let dui_dx = (field[0].index((r, c + 1)) - field[0].index((r, c - 1))) / (2.0 * dx);
            *(unsafe { du_dx.get_unchecked_mut((r, c)) }) = dui_dx
        }
    }

    // set edge nodes
    for r in 0..rows {
        let dui_dx_left = (field[0].index((r, 1)) - field[0].index((r, 0))) / dx;
        let dui_dx_right = (field[0].index((r, cols - 1)) - field[0].index((r, cols - 2))) / dx;

        *(unsafe { du_dx.get_unchecked_mut((r, 0)) }) = dui_dx_left;
        *(unsafe { du_dx.get_unchecked_mut((r, cols - 1)) }) = dui_dx_right;
    }

    // -- y-axis --

    // set interior nodes
    for r in 1..(rows - 1) {
        for c in 0..cols {
            let dui_dy = (field[1].index((r + 1, c)) - field[1].index((r - 1, c))) / (2.0 * dy);
            *(unsafe { dv_dy.get_unchecked_mut((r, c)) }) = dui_dy
        }
    }

    // set edge nodes
    for c in 0..cols {
        let dui_dy_left = (field[1].index((1, c)) - field[1].index((0, c))) / dy;
        let dui_dy_right = (field[1].index((rows - 1, c)) - field[1].index((rows - 2, c))) / dy;

        *(unsafe { dv_dy.get_unchecked_mut((0, c)) }) = dui_dy_left;
        *(unsafe { dv_dy.get_unchecked_mut((rows - 1, c)) }) = dui_dy_right;
    }

    // -- sum matrices --

    let div = du_dx + dv_dy;
    return div;
}

#[cfg(test)]
mod tests {
    use na::dmatrix;

    use super::*;

    #[test]
    fn test_divergence() {
        let field_x: DMatrix<f32> = dmatrix![
            1., 5., 2., 1., 6.;
            5., 4., 3., 0., 2.;
            2., 8., 2., 9., 8.;
            9., 2., 5., 2., 3.;
            5., 2., 2., 1., 7.;
        ];

        let field_y: DMatrix<f32> = dmatrix![
            7., 0., 3., 8., 1.;
            4., 2., 9., 6., 0.;
            5., 1., 4., 7., 3.;
            2., 8., 0., 5., 9.;
            6., 3., 7., 2., 1.;
        ];

        let field = [field_x, field_y];

        let expected_div: DMatrix<f32> = dmatrix![
              2.,   5.,   8.,   0.,   8.;
             -4.,  -1.,  -3.,  -2.,   6.;
             10.,   6.,  -8.,   5.,   7.;
            -13.,  -2.,   3.,  -7.,   0.;
              2., -13.,  13.,  -1.,  -4.;
        ];

        let (dy, dx) = (0.5, 0.5);

        let actual_div = divergence(field, dy, dx);

        println!(
            "Expected div:\n{}\n\nActual div:\n{}\n\n",
            expected_div, actual_div
        );

        assert_eq!(expected_div, actual_div);
    }
}
