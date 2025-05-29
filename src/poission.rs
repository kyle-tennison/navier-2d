use std::usize;

use argmin::{
    core::{
        Executor, KV, Operator, State,
        observers::{Observe, ObserverMode},
    },
    solver::conjugategradient::ConjugateGradient,
};
use indicatif::ProgressBar;
use na::{DMatrix, DVector};
use nalgebra_sparse::{CooMatrix, CsrMatrix};
use num_traits::cast::ToPrimitive;

use crate::ScalarField;

const SOLVE_BAR_TOTAL: usize = 1000; // steps in the progress bar
pub const TARGET_CG_COST: f32 = 1e-4; // maximum error in conjugate-gradient solve
pub const MAX_CG_ITER: u64 = 10_000_000; // maximum number of iterations for conjugate-gradient

/// Runs multiplication for Conjugate Gradient Solver
struct ConjugateGradientOperator<'a> {
    a: &'a CsrMatrix<f32>,
}

impl<'a> Operator for ConjugateGradientOperator<'a> {
    type Param = Vec<f32>;
    type Output = Vec<f32>;

    fn apply(&self, x: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        Ok((self.a * DVector::from_vec(x.to_vec()))
            .data
            .as_vec()
            .clone())
    }
}

pub fn poission_solve(field: &ScalarField, mask: &DMatrix<bool>, step_size: f32) -> ScalarField {
    let (rows, cols) = field.shape();

    // let ij_to_k = |i: usize, j: usize| {i + (j - 1) * cols};

    let mut a_coo: CooMatrix<f32> = CooMatrix::new(rows * cols, rows * cols);
    let mut field_flat: DVector<f32> = DVector::from_row_slice(field.as_slice());

    let ij_to_k = { |(i, j): (i32, i32)| (i + j * (cols as i32)) as usize };

    for i in 0..rows {
        for j in 0..cols {
            let k = ij_to_k((i as i32, j as i32));

            // Dirichlet condition that p=0 inside the mask
            if *mask.index((i, j)) {
                a_coo.push(k, k, 1.);
                *(field_flat.index_mut(k)) = 0.;
                continue;
            }

            // Neumann condition that boundary dp/dn = 0
            let mut diag: usize = 0;
            for (di, dj) in [(-1, 0), (1, 0), (0, -1), (0, 1)] {
                let (ni, nj) = ((i as i32) + di, (j as i32) + dj);

                if (0 <= ni && ni < rows as i32)
                    && (0 <= nj && nj < cols as i32)
                    && !(*mask.index((ni as usize, nj as usize)))
                {
                    let nk = ij_to_k((ni, nj));

                    a_coo.push(k, nk, -1.);
                    diag += 1;
                } else {
                    diag += 1;
                }
            }

            // Add center coefficient at the end
            a_coo.push(k, k, diag as f32);
            *(field_flat.index_mut(k)) = -field.index((i, j)) * (step_size.powi(2));
        }
    }

    *(field_flat.index_mut(0)) = 0.; // pin pressure

    let a_csr = CsrMatrix::from(&a_coo); // convert to csr sparse

    // solve system
    let b: Vec<f32> = field_flat.iter().map(|f| *f).collect();

    #[cfg(debug_assertions)]
    {
        let norm_b = b.iter().map(|v| v * v).sum::<f32>().sqrt();

        assert_ne!(norm_b, 0.);

        debug_assert_eq!(b.iter().any(|i| i.is_nan() | i.is_infinite()), false);
        debug_assert_eq!(
            a_coo
                .triplet_iter()
                .any(|i| i.2.is_nan() | i.2.is_infinite()),
            false
        );
    }

    let solver: ConjugateGradient<Vec<f32>, f32> = ConjugateGradient::new(b);
    let initial_guess: Vec<f32> = vec![0.0; field_flat.nrows()];

    let operator = ConjugateGradientOperator { a: &a_csr };

    // Now run the solver as before, but with DVector<f32> everywhere
    let res = match Executor::new(operator, solver)
        .configure(|state| {
            state
                .param(initial_guess)
                .max_iters(MAX_CG_ITER)
                .target_cost(TARGET_CG_COST)
        })
        .run()
    {
        Ok(r) => r,
        Err(err) => {
            panic!("Conjugate Gradient error: {err}");
        }
    };

    let best_param = res
        .state()
        .best_param
        .as_ref()
        .expect("Conjugate Gradient failed.")
        .to_owned();

    let p = DMatrix::from_vec(rows, cols, best_param);

    return p;
}
