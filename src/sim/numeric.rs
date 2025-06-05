// Numeric approximations

use na::DMatrix;

use crate::{ScalarField, VectorField};

/// Compute the finite-difference gradient of a scalar field **in the x axis**.
/// Uses a central finite difference for interior nodes and a first-order
/// forward/backward (depends on side) finite difference for edge nodes.
///
/// Parameters:
/// - `field` - A reference to the scalar field to take the gradient of
/// - `dx` - The size of the elements in the x-axis
///
/// Returns:
/// - A `ScalarField` with the finite difference of each element
pub fn gradient_x(field: &ScalarField, dx: f32) -> ScalarField {
    let (rows, cols) = field.shape();

    let mut df_dx: DMatrix<f32> = DMatrix::zeros(rows, cols);

    // set interior nodes
    for r in 0..rows {
        for c in 1..(cols - 1) {
            let dfi_dx = (field.index((r, c + 1)) - field.index((r, c - 1))) / (2.0 * dx);
            *(df_dx.index_mut((r, c))) = dfi_dx
        }
    }

    // set edge nodes
    for r in 0..rows {
        let dfi_dx_left = (field.index((r, 1)) - field.index((r, 0))) / dx;
        let dfi_dx_right = (field.index((r, cols - 1)) - field.index((r, cols - 2))) / dx;

        *(df_dx.index_mut((r, 0))) = dfi_dx_left;
        *(df_dx.index_mut((r, cols - 1))) = dfi_dx_right;
    }

    df_dx
}

/// Compute the finite-difference gradient of a scalar field **in the y axis**.
/// Uses a central finite difference for interior nodes and a first-order
/// forward/backward (depends on side) finite difference for edge nodes.
///
/// Parameters:
/// - `field` - A reference to the scalar field to take the gradient of
/// - `dy` - The size of the elements in the y-axis
///
/// Returns:
/// - A `ScalarField` with the finite difference of each element
pub fn gradient_y(field: &ScalarField, dy: f32) -> ScalarField {
    let (rows, cols) = field.shape();

    let mut df_dy: DMatrix<f32> = DMatrix::zeros(rows, cols);

    // set interior nodes
    for r in 1..(rows - 1) {
        for c in 0..cols {
            let dui_dy = (field.index((r + 1, c)) - field.index((r - 1, c))) / (2.0 * dy);
            *(df_dy.index_mut((r, c))) = dui_dy
        }
    }

    // set edge nodes
    for c in 0..cols {
        let dui_dy_left = (field.index((1, c)) - field.index((0, c))) / dy;
        let dui_dy_right = (field.index((rows - 1, c)) - field.index((rows - 2, c))) / dy;

        *(df_dy.index_mut((0, c))) = dui_dy_left;
        *(df_dy.index_mut((rows - 1, c))) = dui_dy_right;
    }

    df_dy
}

/// Compute the gradient of a scalar field **in the x axis**. Interior
/// nodes are computed using an upwind scheme, while edge nodes use
/// finite differences.
///
/// Parameters:
/// - `field` - A reference to the scalar field to take the gradient of
/// - `sign_field` - The field to reference for upwind/downwind direction
/// - `dx` - The size of the elements in the x-axis
///
/// Returns:
/// - A `ScalarField` with the x-gradient of each element
fn gradeint_x_upwind(field: &ScalarField, sign_field: &ScalarField, dx: f32) -> ScalarField {
    let (rows, cols) = field.shape();

    let mut df_dx: DMatrix<f32> = DMatrix::zeros(rows, cols);

    // set interior nodes
    let mut dui_dx: f32;
    for r in 0..rows {
        for c in 1..(cols - 1) {
            if *sign_field.index((r, c)) > 0. {
                dui_dx = (field.index((r, c)) - field.index((r, c - 1))) / dx;
            } else {
                dui_dx = (field.index((r, c + 1)) - field.index((r, c))) / dx;
            }

            *(df_dx.index_mut((r, c))) = dui_dx
        }
    }

    // set edge nodes (using 1st order finite-diff)
    for r in 0..rows {
        let dfi_dx_left = (field.index((r, 1)) - field.index((r, 0))) / dx;
        let dfi_dx_right = (field.index((r, cols - 1)) - field.index((r, cols - 2))) / dx;

        *(df_dx.index_mut((r, 0))) = dfi_dx_left;
        *(df_dx.index_mut((r, cols - 1))) = dfi_dx_right;
    }

    df_dx
}

/// Compute the gradient of a scalar field **in the y axis**. Interior
/// nodes are computed using an upwind scheme, while edge nodes use
/// finite differences.
///
/// Parameters:
/// - `field` - A reference to the scalar field to take the gradient of
/// - `sign_field` - The field to reference for upwind/downwind direction
/// - `dy` - The size of the elements in the y-axis
///
/// Returns:
/// - A `ScalarField` with the y-gradient of each element
fn gradeint_y_upwind(field: &ScalarField, sign_field: &ScalarField, dy: f32) -> ScalarField {
    let (rows, cols) = field.shape();

    let mut df_dy: DMatrix<f32> = DMatrix::zeros(rows, cols);

    // set interior nodes
    let mut dui_dy: f32;
    for r in 1..(rows - 1) {
        for c in 0..cols {
            if *sign_field.index((r, c)) > 0. {
                dui_dy = (field.index((r, c)) - field.index((r - 1, c))) / dy;
            } else {
                dui_dy = (field.index((r + 1, c)) - field.index((r, c))) / dy;
            }

            *(df_dy.index_mut((r, c))) = dui_dy
        }
    }

    // set edge nodes
    for c in 0..cols {
        let dui_dy_left = (field.index((1, c)) - field.index((0, c))) / dy;
        let dui_dy_right = (field.index((rows - 1, c)) - field.index((rows - 2, c))) / dy;

        *(df_dy.index_mut((0, c))) = dui_dy_left;
        *(df_dy.index_mut((rows - 1, c))) = dui_dy_right;
    }

    df_dy
}

/// Compute the divergence of some vector field F=<u,v>. That is, ∇⋅F
///
/// Mathematically, this is du/dx + dv/dy
///
/// Parameters:
/// - `field` - The `VectorField` to take the divergence of
/// - `dy` - The y-axis step size
/// - `dx` - The x-axis step size
///
/// Returns:
///     A `ScalarField` of the divergence.
pub fn divergence(field: &VectorField, dy: f32, dx: f32) -> ScalarField {
    let du_dx: DMatrix<f32> = gradient_x(&field[0], dx);
    let dv_dy: DMatrix<f32> = gradient_y(&field[1], dy);

    du_dx + dv_dy
}

/// Compute the laplacian of a scalar field f. That is ∇²f
///
/// Mathematically, this is ∇²f = ∂²f/∂x² + ∂²f/∂y²
///
/// Parameters
/// - `field` - The vector field to take the laplaican of
/// - `dy` - The y-axis step size
/// - `dx` - THe x-axis step size
///
/// Returns:
///     A `VectorField` of the laplacian.
pub fn laplacian(field: &ScalarField, dy: f32, dx: f32) -> ScalarField {
    // first-order derivatives
    let df_dx = gradient_x(field, dx);
    let df_dy = gradient_y(field, dy);

    // second-order derivatives
    let d2f_dx2 = gradient_x(&df_dx, dx);
    let d2f_dy2 = gradient_y(&df_dy, dy);

    d2f_dx2 + d2f_dy2
}

/// Compute the laplacian of some (cartesian) vector field F=<u,v>. That is ∇²F
///
/// Mathematically, this is ∇²F = <∇²u, ∇²v>
///
/// Parameters
/// - `field` - The vector field to take the laplaican of
/// - `dy` - The y-axis step size
/// - `dx` - THe x-axis step size
///
/// Returns:
///     A `VectorField` of the laplacian.
pub fn laplacian_vf(field: &VectorField, dy: f32, dx: f32) -> VectorField {
    let (u, v) = (&field[0], &field[1]);

    [laplacian(u, dy, dx), laplacian(v, dy, dx)]
}

pub fn zero_where_mask(field: &mut ScalarField, mask: &DMatrix<bool>) {
    // zero-out velocity in object shape
    let mask_flat = mask.as_slice();

    for (i, msk) in mask_flat.iter().enumerate() {
        if *msk {
            field.as_mut_slice()[i] = 0.;
        }
    }
}

/// Computes the upwind advection (u⋅∇)u, where u is a vector field.
///
/// Parameters:
/// - `field` - The field to take the advection of (i.e. u shown above)
/// - `dy` - The y-axis step size
/// - `dx` - THe x-axis step size
///
/// Returns:
///     A `VectorField` of the advection.
pub fn advection_upwind(
    field: &VectorField,
    mask: &DMatrix<bool>,
    dy: f32,
    dx: f32,
) -> VectorField {
    let (u, v) = (&field[0], &field[1]);

    let mut du_dx = gradeint_x_upwind(u, u, dx);
    let mut du_dy = gradeint_y_upwind(u, v, dy);

    let mut dv_dx = gradeint_x_upwind(v, u, dx);
    let mut dv_dy = gradeint_y_upwind(v, v, dy);

    zero_where_mask(&mut du_dx, mask);
    zero_where_mask(&mut du_dy, mask);
    zero_where_mask(&mut dv_dx, mask);
    zero_where_mask(&mut dv_dy, mask);

    let adv_u = u.component_mul(&du_dx) + v.component_mul(&du_dy);
    let adv_v = u.component_mul(&dv_dx) + v.component_mul(&dv_dy);

    [adv_u, adv_v]
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

        dbg!(format!(
            "Expected div:\n{}\n\nActual div:\n{}\n\n",
            expected_div, actual_div
        ));

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

        let actual_lap = laplacian_vf(&field, 0.5, 0.5);

        assert_eq!(expected_lap, actual_lap);
    }

    #[test]
    fn test_advection() {
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

        let expected_adv_x = dmatrix![
         -66.,   40.,  -12.,   -2.,   58.;
          62.,  -12.,   12.,  -12.,    8.;
         -54.,  104.,  -32.,  252.,   20.;
         136., -124.,   30.,  -82.,  -84.;
         -68.,  -12.,  -42.,   -6.,   92.;
        ];

        let expected_adv_y = dmatrix![
              26.,  -70.,  -12.,  106.,  -84.;
              16.,   -8.,  150.,  -24.,  -24.;
              18.,  -66.,  -28.,   68.,  -46.;
            -138.,  136.,  -80.,    0.,  132.;
              98.,  -42.,  114.,  -22.,  -30.;
        ];

        let (nrows, ncols) = &expected_adv_x.shape();
        let expected_adv = [expected_adv_x, expected_adv_y];

        let mask = DMatrix::from_element(*nrows, *ncols, false);
        let actual_adv = advection_upwind(&field, &mask, 0.5, 0.5);

        println!(
            "Expected adv_x:\n{}\n\nActual adv_x:\n{}\n\n",
            &expected_adv[0], &actual_adv[0]
        );

        println!(
            "Expected adv_y:\n{}\n\nActual adv_y:\n{}\n\n",
            &expected_adv[1], &actual_adv[1]
        );

        println!("{}", (&expected_adv[0] - &actual_adv[0]));
        assert!((&expected_adv[0] - &actual_adv[0]).norm() < 1.);
        assert!((&expected_adv[1] - &actual_adv[1]).norm() < 1.);
    }
}
