use na::DMatrix;

use crate::{display::image_save, numeric, poission, ScalarField, VectorField};

pub struct NewtonianSim {
    density: f32,
    shear_viscosity: f32,
    inflow: (f32, f32), // (x,y)
    solid_mask: DMatrix<bool>,
    dt: f32,
    simtime: f32,
    rows: usize,
    cols: usize,
    dx: f32,
    dy: f32,
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
    pub fn new(
        density: f32,
        shear_viscosity: f32,
        inflow: (f32, f32),
        solid_mask: &DMatrix<bool>,
        length: (f32, f32),
        timestep: f32,
        simtime: f32,
    ) -> Self {
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
            inflow,
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
        self.u[0].columns_mut(0, 3).fill(self.inflow.0);
        self.u[1].columns_mut(0, 3).fill(self.inflow.1);

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
            &self.u[0]
                + (self.dt / self.density)
                    * (self.shear_viscosity * &laplacian_u[0] + &self.f[0]
                        - self.density * &advection_u[0])
        };
        let ustar_y: ScalarField = {
            &self.u[1]
                + (self.dt / self.density)
                    * (self.shear_viscosity * &laplacian_u[1] + &self.f[1]
                        - self.density * &advection_u[1])
        };

        [ustar_x, ustar_y]
    }

    pub fn sample_shape_mask(ny: usize, nx: usize) -> DMatrix<bool> {
        let mut mask = DMatrix::from_element(ny, nx, false);

        let height = (2 * ny) / 2 - ny / 2;
        let width = (2 * nx) / 3 - nx / 3;

        let new_height = height / 2;
        let new_width = width / 2;

        let y_start = ny / 2 - new_height / 2;
        let y_end = y_start + new_height;

        let x_start = nx / 2 - 2 * new_width;
        let x_end = x_start + new_width;

        for y in y_start..y_end {
            for x in x_start..x_end {
                mask[(y, x)] = true;
            }
        }

        mask
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
        let poission_rhs: ScalarField =
        (self.density / self.dt) * numeric::divergence(&u_star, self.dy, self.dx);

        let u_mag = &self.u[0].pow(2) + &self.u[1].pow(2);
        let u_star_mag = u_star[0].pow(2) +u_star[1].pow(2);
        image_save(
            &u_mag, format!("u_mag-{}.png", self.t).as_str()
        ).unwrap();
        image_save(
            &u_star_mag, format!("u_star_mag-{}.png", self.t).as_str()
        ).unwrap();


        let ustar_div = &numeric::divergence(&u_star, self.dy, self.dx);
        image_save(
            ustar_div, format!("u_star_div-{}.png", self.t).as_str()
        ).unwrap();

        assert_ne!(ustar_div.norm(), 0.);

        println!("ustar norm: {}", ustar_div.norm());

        println!("scale: {}", self.density / self.dt);
        println!("expected norm: {}", (self.density / self.dt) * ustar_div.norm());

        let this = (self.density / self.dt) * ustar_div;
        println!("actual norm: {}", this.norm());
        println!("Passed norm: {}", poission_rhs.norm());

        
        
        #[cfg(debug_assertions)]
        {
            let u_star_norm = u_star[0].norm() + u_star[1].norm();
            debug_assert_ne!(u_star_norm, 0.);
            
            let poission_rhs_norm = poission_rhs.norm();
            debug_assert_ne!(poission_rhs_norm, 0.);
        }
        
        // panic!();

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
