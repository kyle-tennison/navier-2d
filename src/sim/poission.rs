// Iterative solver for a generic poission equation.

use argmin::{
    core::{Executor, Operator},
    solver::conjugategradient::ConjugateGradient,
};
use na::{DMatrix, DVector};
use nalgebra_sparse::{CooMatrix, CsrMatrix};

use crate::ScalarField;

pub const TARGET_CG_COST: f32 = 1e-4; // maximum error in conjugate-gradient solve
pub const MAX_CG_ITER: u64 = 100_000; // maximum number of iterations for conjugate-gradient

/// Runs multiplication for Conjugate Gradient Solver
struct ConjugateGradientOperator<'a> {
    a: &'a CsrMatrix<f32>,
}

impl Operator for ConjugateGradientOperator<'_> {
    type Param = Vec<f32>;
    type Output = Vec<f32>;

    fn apply(&self, x: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        Ok((self.a * DVector::from_vec(x.to_vec()))
            .data
            .as_vec()
            .clone())
    }
}

/// Iteratively solve the poission equation using the Conjugate Gradient method.
/// Mathematically, this is: ∇²p = g
///
/// Parameters
/// - `field` - The field (g); i.e. the RHS of the poission equation
/// - `mask` - A boolean mask representing regions to exclude from the solution (e.g., walls, solids, etc)
/// - `step_size` - The space step size (dx or dy) that is assumed to be uniform in both axes
pub fn poission_solve(field: &ScalarField, mask: &DMatrix<bool>, step_size: f32) -> ScalarField {
    // create a mapping between the two coordinate systems
    let (rows, cols) = field.shape();
    let ij_to_k = { |(i, j): (i32, i32)| (i + j * (rows as i32)) as usize };

    // start with a coordinate sparse rep for easy loading; convert to compressed sparse-row after
    let mut a_coo: CooMatrix<f32> = CooMatrix::new(rows * cols, rows * cols);

    // flattern RHS into a vector
    let mut field_flat: DVector<f32> = DVector::from_row_slice(field.as_slice());

    // iterate over field, load coordinate matrix
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
            for (di, dj) in [(-1, 0), (1, 0), (0, -1), (0, 1)] {
                let (ni, nj) = ((i as i32) + di, (j as i32) + dj);

                if (0 <= ni && ni < rows as i32)
                    && (0 <= nj && nj < cols as i32)
                    && !(*mask.index((ni as usize, nj as usize)))
                {
                    let nk = ij_to_k((ni, nj));
                    a_coo.push(k, nk, -1.);
                } 
            }

            // add center coefficient at the end
            a_coo.push(k, k, 4 as f32);
            *(field_flat.index_mut(k)) = -field.index((i, j)) * (step_size.powi(2));
        }
    }

    *(field_flat.index_mut(0)) = 0.; // pin pressure

    let a_csr = CsrMatrix::from(&a_coo); // convert to csr sparse

    // solve system
    let b: Vec<f32> = field_flat.iter().copied().collect();
    let solver: ConjugateGradient<Vec<f32>, f32> = ConjugateGradient::new(b);
    let initial_guess: Vec<f32> = vec![0.0; field_flat.nrows()];
    let operator = ConjugateGradientOperator { a: &a_csr };

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

    DMatrix::from_vec(rows, cols, best_param)
}
