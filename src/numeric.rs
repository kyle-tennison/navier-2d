// Numeric approximations

use std::mem::MaybeUninit;

use na::{DMatrix, Dim, Dyn, Matrix, RawStorage};

type ScalarField = DMatrix<f32>;
type VectorField = [ScalarField; 2];



fn gradient_x(field: &ScalarField, dx: f32) -> ScalarField {

    let (rows, cols) = field.shape();

    let mut df_dx: DMatrix<f32> = DMatrix::zeros(rows, cols);

    // set interior nodes
    for r in 0..rows {
        for c in 1..(cols - 1) {
            let dfi_dx = (field.index((r, c + 1)) - field.index((r, c - 1))) / (2.0 * dx);
            *(unsafe { df_dx.get_unchecked_mut((r, c)) }) = dfi_dx
        }
    };

    // set edge nodes
    for r in 0..rows {
        let dfi_dx_left = (field.index((r, 1)) - field.index((r, 0))) / dx;
        let dfi_dx_right = (field.index((r, cols - 1)) - field.index((r, cols - 2))) / dx;

        *(unsafe { df_dx.get_unchecked_mut((r, 0)) }) = dfi_dx_left;
        *(unsafe { df_dx.get_unchecked_mut((r, cols - 1)) }) = dfi_dx_right;
    };

    return df_dx;

}

fn gradient_y(field: &ScalarField, dy: f32) -> ScalarField {

    let (rows, cols) = field.shape();

    let mut df_dy: DMatrix<f32> = DMatrix::zeros(rows, cols);

    // set interior nodes
    for r in 1..(rows - 1) {
        for c in 0..cols {
            let dui_dy = (field.index((r + 1, c)) - field.index((r - 1, c))) / (2.0 * dy);
            *(unsafe { df_dy.get_unchecked_mut((r, c)) }) = dui_dy
        }
    };

    // set edge nodes
    for c in 0..cols {
        let dui_dy_left = (field.index((1, c)) - field.index((0, c))) / dy;
        let dui_dy_right = (field.index((rows - 1, c)) - field.index((rows - 2, c))) / dy;

        *(unsafe { df_dy.get_unchecked_mut((0, c)) }) = dui_dy_left;
        *(unsafe { df_dy.get_unchecked_mut((rows - 1, c)) }) = dui_dy_right;
    };

    return df_dy;

}



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
pub fn divergence(field: &VectorField, dy: f32, dx: f32) -> ScalarField {

    let du_dx: DMatrix<f32> = gradient_x(&field[0], dx);
    let dv_dy: DMatrix<f32> = gradient_y(&field[1], dy);

    let div = du_dx + dv_dy;
    return div;
}


pub fn laplacian(field: &VectorField, dy: f32, dx: f32) -> VectorField {

    let (u,v) = (&field[0], &field[1]);

    // first-order derivatives
    let du_dy = gradient_y(u, dy);
    let du_dx = gradient_x(u, dx);

    let dv_dy = gradient_y(v, dy);
    let dv_dx = gradient_x(v, dx);

    // second-order derivatives
    let d2u_dy2 = gradient_y(&du_dy, dy);
    let d2u_dx2 = gradient_x(&du_dx, dx);

    let d2v_dy2 = gradient_y(&dv_dy, dy);
    let d2v_dx2 = gradient_x(&dv_dx, dx);

    // form field from components
    let lap_u = d2u_dx2 + d2u_dy2;
    let lap_v = d2v_dx2 + d2v_dy2;
    let laplacian = [lap_u, lap_v];

    return laplacian;

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

        let actual_div = divergence(&field, dy, dx);

        println!(
            "Expected div:\n{}\n\nActual div:\n{}\n\n",
            expected_div, actual_div
        );

        assert_eq!(expected_div, actual_div);
    }


    #[test]
    fn test_laplacian() {
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

        let expected_lap_x: DMatrix<f32> = dmatrix![
            -28.,  -2.,  -1.,  34.,  32.;
             -4.,  -2.,   1.,  12.,  19.;
            -22., -20.,   6., -19., -19.;
              8.,  16.,  -6.,  -2.,  15.;
            -16.,  17.,  -4.,  25.,  32.;
        ];

        let expected_lap_y: DMatrix<f32> = dmatrix![
             28.,  16., -20., -16., -16.;
             22.,  10., -35., -13.,   5.;
             17.,  15.,   2., -18., -18.;
            -18., -31.,  34.,   6., -27.;
             28., -19.,  15.,  -3., -20.;
        ];

        let expected_lap = [expected_lap_x, expected_lap_y];

        let actual_lap = laplacian(&field, 0.5, 0.5);

        assert_eq!(expected_lap, actual_lap);

    }
}
