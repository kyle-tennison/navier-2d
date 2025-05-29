use na::DMatrix;

use crate::{numeric, poission, ScalarField, VectorField };

struct NewtonianSim {
    density: f32,
    shear_viscosity: f32,
    solid_mask: DMatrix<bool>,
    rows: usize,
    cols: usize,
    dx: f32,
    dy: f32,
    dt: f32,
    simtime: f32,
    t: f32,
    u: VectorField,
    p: ScalarField,
    f: VectorField,
}

impl NewtonianSim {
    /// Create a new Newtonian fluid simulation instance
    ///
    /// Parameters
    /// - `domain` - The domain (x, y) of the 2D simulation
    /// - `step` - The step size per cell
    pub fn new(density: f32, shear_viscosity: f32, solid_mask: &DMatrix<bool>, length: (f32, f32), timestep: f32, simtime: f32) -> Self {

        let (rows, cols) = solid_mask.shape();

        let (lx, ly) = length;
        let dy = ly / (rows as f32);
        let dx = lx / (cols as f32);

        let ux: DMatrix<f32> = DMatrix::zeros(rows, cols);
        let uy: DMatrix<f32> = DMatrix::zeros(rows, cols);
        let p: DMatrix<f32> = DMatrix::zeros(rows, cols);
        let fx: DMatrix<f32> = DMatrix::zeros(rows, cols);
        let fy: DMatrix<f32> = DMatrix::zeros(rows, cols);

        NewtonianSim {
            density,
            shear_viscosity,
            solid_mask: solid_mask.to_owned(),
            rows,
            cols,
            dx,
            dy,
            dt: timestep,
            simtime,
            t: 0.,
            u: [ux, uy],
            p,
            f: [fx, fy],
        }
    }

    /// Set the boundary values on the velocity field
    fn set_bv_u(&mut self) {
        // set inflow on left
        self.u[0].rows_mut(0, 3).fill(0.);
        self.u[1].rows_mut(0, 3).fill(0.);

        // set zero-derivative at outflow
        let xcol = self.u[0].column(self.cols - 2).into_owned();
        self.u[0].set_column(self.cols - 1, &xcol);

        let ycol = self.u[1].column(self.cols - 2).into_owned();
        self.u[1].set_column(self.cols - 1, &ycol);

        // set no-slip on top & bottom
        self.u[0].row_mut(0).fill(0.);
        self.u[1].row_mut(0).fill(0.);
        self.u[0].row_mut(self.rows - 1).fill(0.);
        self.u[1].row_mut(self.rows - 1).fill(0.);

        ()
    }

    /// Predict u*
    fn predict_u_star(&self) -> VectorField {

        let laplacian_u: VectorField = numeric::laplacian_vf(&self.u, self.dy, self.dx);
        let advection_u: VectorField = numeric::advection_upwind(&self.u, self.dy, self.dx);

        let ustar_x: ScalarField = {
            &self.u[0] + (self.dt / self.density) * (
                self.shear_viscosity * &laplacian_u[0] + &self.f[0] - self.density * &advection_u[0]
            )
        };
        let ustar_y: ScalarField = {
            &self.u[1] + (self.dt / self.density) * (
                self.shear_viscosity * &laplacian_u[1] + &self.f[1] - self.density * &advection_u[1]
            )
        };

        [ustar_x, ustar_y]

    }


}

impl Iterator for NewtonianSim {
    type Item = VectorField;

    fn next(&mut self) -> Option<Self::Item> {
        if self.t > self.simtime {
            return None;
        }

        self.set_bv_u();

        // zero-out velocity in object shape
        let mask_flat = self.solid_mask.as_slice();

        for i in 0..mask_flat.len() {
            if mask_flat[i] {
                self.u[0].as_mut_slice()[i] = 0.;
                self.u[1].as_mut_slice()[i] = 0.;
            }
        }

        // predict ustar
        let u_star: VectorField = self.predict_u_star();

        // solve pressure gradient
        let poission_rhs: ScalarField = (self.density / self.dt) * numeric::divergence(&u_star, self.dy, self.dx);
        let p: ScalarField = poission::poission_solve(&poission_rhs, &self.solid_mask, self.dx);

        let dp_dx: ScalarField = numeric::gradient_x(&p, self.dx);
        let dp_dy: ScalarField = numeric::gradient_y(&p, self.dy);

        // compute next velocity field
        let next_ux: ScalarField = &u_star[0] - (self.dt / self.density) * dp_dx;
        let next_uy: ScalarField = &u_star[1] - (self.dt / self.density) * dp_dy;

        let next_u: VectorField = [next_ux, next_uy];

        // normalize velocity to prevent explosion
        let next_u: VectorField = [next_u[0].normalize(), next_u[1].normalize()];

        return Some(next_u);

    }
}
