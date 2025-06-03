use na::DMatrix;
use rand::{rngs::ThreadRng, Rng};

use crate::{display::image_save, numeric, poission, ScalarField, VectorField};

const MAX_VELOCITY: f32 = 1000.;

pub struct NewtonianSim {
    pub density: f32,
    pub shear_viscosity: f32,
    pub inflow: (f32, f32), // (x,y)
    pub solid_mask: DMatrix<bool>,
    pub simtime: f32,
    pub cfl: f32,
    rows: usize,
    cols: usize,
    dx: f32,
    dy: f32,
    pub t: f32,
    u: VectorField,
    f: VectorField,
    i: usize,
    rng: ThreadRng
}

impl<'a> NewtonianSim {
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
        simtime: f32,
        cfl: f32,
    ) -> Self {
        let (rows, cols) = solid_mask.shape();

        let (lx, ly) = length;
        let dy = ly / (rows as f32);
        let dx = lx / (cols as f32);

        let ux: DMatrix<f32> = DMatrix::zeros(rows, cols);
        let uy: DMatrix<f32> = DMatrix::zeros(rows, cols);
        let fx: DMatrix<f32> = DMatrix::zeros(rows, cols);
        let fy: DMatrix<f32> = DMatrix::zeros(rows, cols);

        let rng = rand::rng();

        NewtonianSim {
            density,
            shear_viscosity,
            inflow,
            solid_mask: solid_mask.to_owned(),
            simtime,
            cfl,
            rows,
            cols,
            dx,
            dy,
            t: 0.,
            u: [ux, uy],
            f: [fx, fy],
            i: 0,
            rng
        }
    }

    /// Set the boundary values on the velocity field
    fn set_bv_u(&mut self) {

        let velocity_noise: f32 = self.rng.random_range(0.0..0.1);

        // set inflow on left
        self.u[0].columns_mut(0, 3).fill(self.inflow.0 * (1.+velocity_noise));
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
    }

    /// Predict u*
    fn predict_u_star(&self, dt: f32) -> VectorField {
        let laplacian_u: VectorField = numeric::laplacian_vf(&self.u, self.dy, self.dx);
        let advection_u: VectorField = numeric::advection_upwind(&self.u, &self.solid_mask, self.dy, self.dx);

        let ustar_x: ScalarField = {
            &self.u[0]
                + (dt / self.density)
                    * (self.shear_viscosity * &laplacian_u[0] + &self.f[0]
                        - self.density * &advection_u[0])
        };
        let ustar_y: ScalarField = {
            &self.u[1]
                + (dt / self.density)
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
    type Item = (VectorField, f32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.t > self.simtime {
            return None;
        }
        self.set_bv_u();

        let max_ux = *(&self.u[0].iter().fold(0.0f32, |m,&x| m.max(x.abs())));
        let max_uy = *(&self.u[1].iter().fold(0.0f32, |m,&x| m.max(x.abs())));
        let dt = self.cfl / (max_ux/self.dx + max_uy/self.dy);

        // zero-out velocity in object shape
        let mask_flat = self.solid_mask.as_slice();

        for (i, msk) in mask_flat.iter().enumerate() {
            if *msk {
                self.u[0].as_mut_slice()[i] = 0.;
                self.u[1].as_mut_slice()[i] = 0.;
            }
        }

        // predict ustar
        let u_star: VectorField = self.predict_u_star(dt);

        // solve pressure gradient
        let poission_rhs: ScalarField =
            (self.density / dt) * numeric::divergence(&u_star, self.dy, self.dx);

        #[cfg(debug_assertions)]
        {
            let u_star_norm = u_star[0].norm() + u_star[1].norm();
            debug_assert_ne!(u_star_norm, 0.);

            let poission_rhs_norm = poission_rhs.norm();
            debug_assert_ne!(poission_rhs_norm, 0.);
        }

        let p: ScalarField = poission::poission_solve(&poission_rhs, &self.solid_mask, self.dx);

        let dp_dx: ScalarField = numeric::gradient_x(&p, self.dx);
        let dp_dy: ScalarField = numeric::gradient_y(&p, self.dy);

        // compute next velocity field
        let mut next_ux: ScalarField = &u_star[0] - (dt / self.density) * dp_dx;
        let mut next_uy: ScalarField = &u_star[1] - (dt / self.density) * dp_dy;

        // cap velocity

        if next_ux.iter().any(|ux| *ux > MAX_VELOCITY) || next_uy.iter().any(|uy| *uy > MAX_VELOCITY) {
            eprintln!("Velocity exceeded maximum; simulation exploded :(");
            return None;
        }

        next_ux.iter_mut().for_each(|f| {
            if *f > MAX_VELOCITY {
                *f = MAX_VELOCITY
            }
        });
        next_uy.iter_mut().for_each(|f| {
            if *f > MAX_VELOCITY {
                *f = MAX_VELOCITY
            }
        });

        let next_u: VectorField = [next_ux, next_uy];

        // normalize velocity to prevent explosion
        // let next_u: VectorField = [next_u[0].normalize(), next_u[1].normalize()];

        self.i += 1;
        self.t += dt;
        self.u = next_u.clone();
        

        // let velocity_magnitude = (self.u[0].map(|x| x.powi(2))
        //     + self.u[1].map(|y| y.powi(2)))
        // .map(|k| k.sqrt());


        // let advection = (numeric::advection_upwind(&self.u, &self.solid_mask, self.dy, self.dx));

        // let advection_magnitude = (advection[0].map(|x| x.powi(2))
        //     + advection[1].map(|y| y.powi(2)))
        // .map(|k| k.sqrt());
    
        // let pressure_gradient = [numeric::gradient_x(&p, self.dx), numeric::gradient_y(&p, self.dy)];
        // let pressure_gradient_magnitude = (pressure_gradient[0].map(|x| x.powi(2))
        //     + pressure_gradient[1].map(|y| y.powi(2)))
        // .map(|k| k.sqrt());
    
        // let laplacian = numeric::laplacian_vf(&self.u, self.dy, self.dx);
        // let laplacian_magnitude = (laplacian[0].map(|x| x.powi(2))
        //     + laplacian[1].map(|y| y.powi(2)))
        // .map(|k| k.sqrt());
    
        
        // image_save(&p, format!("pressure-{}.png", self.i).as_str()).unwrap();
        // image_save(&p, format!("velocity-{}.png", self.i).as_str()).unwrap();
        // image_save(&velocity_magnitude, format!("velocity-{}.png", self.i).as_str()).unwrap();
        // image_save(&advection_magnitude, format!("advection-{}.png", self.i).as_str()).unwrap();
        // image_save(&laplacian_magnitude, format!("laplacian-{}.png", self.i).as_str()).unwrap();
        // image_save(&pressure_gradien`t_magnitude, format!("pressure-gradient-{}.png", self.i).as_str()).unwrap();


        Some((next_u, self.t))
    }
}
