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

/// Observer bar for argmin solver
struct ConjugateGradientObserverBar {
    bar: ProgressBar,
    final_mag: f32,
    done_solve: bool,
}

impl ConjugateGradientObserverBar {
    fn new() -> ConjugateGradientObserverBar {
        ConjugateGradientObserverBar {
            bar: ProgressBar::new(SOLVE_BAR_TOTAL as u64),
            final_mag: TARGET_CG_COST.log10().floor(),
            done_solve: false,
        }
    }
}

impl<I> Observe<I> for ConjugateGradientObserverBar
where
    I: State,
{
    fn observe_init(
        &mut self,
        _name: &str,
        _state: &I,
        _kv: &KV,
    ) -> Result<(), argmin::core::Error> {
        Ok(())
    }

    fn observe_iter(&mut self, state: &I, _kv: &KV) -> Result<(), argmin::core::Error> {
        let value = state.get_cost().to_f32();

        let cost_mag = state.get_cost().to_f32().unwrap().ln().floor();

        let mut progress = ((SOLVE_BAR_TOTAL as f32) / (cost_mag - self.final_mag).sqrt())
            .ceil()
            .to_usize()
            .unwrap();

        if progress > SOLVE_BAR_TOTAL {
            progress = SOLVE_BAR_TOTAL
        }

        if progress == SOLVE_BAR_TOTAL {
            self.done_solve = true;
            self.bar.set_position(progress as u64);
        }

        if !self.done_solve {
            self.bar.set_position(progress as u64);
        }

        Ok(())
    }

    fn observe_final(&mut self, state: &I) -> Result<(), argmin::core::Error> {
        self.bar.finish();
        let iterations = state.get_iter();

        println!(
            "info: finished conjugate gradient approximation in {} iterations",
            iterations
        );
        Ok(())
    }
}

pub fn poission_solve(field: &ScalarField, mask: DMatrix<bool>, step_size: f32) {
    let (rows, cols) = field.shape();

    // let ij_to_k = |i: usize, j: usize| {i + (j - 1) * cols};

    let mut a_coo: CooMatrix<f32> = CooMatrix::new(rows, cols);
    let mut field_flat: DVector<f32> = DVector::from_row_slice(field.as_slice());

    for i in 0..rows {
        for j in 0..cols {
            let k = i + (j - 1) * cols;

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
                    let nk = (ni + (nj - 1) * (cols as i32)) as usize;

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
    let solver: ConjugateGradient<Vec<f32>, f32> = ConjugateGradient::new(b);
    let initial_guess: Vec<f32> = vec![0.0; field_flat.nrows()];

    let operator = ConjugateGradientOperator { a: &a_csr };
    let observer = ConjugateGradientObserverBar::new();

    // let executor = Executor::new(operator, solver);

    // Now run the solver as before, but with DVector<f32> everywhere
    let res = match Executor::new(operator, solver)
        .configure(|state| {
            state
                .param(initial_guess)
                .max_iters(MAX_CG_ITER)
                .target_cost(TARGET_CG_COST)
        })
        .add_observer(observer, ObserverMode::NewBest)
        .run()
    {
        Ok(r) => r,
        Err(err) => {
            panic!("Conjugate Gradient error: {err}");
        }
    };
    //         state
    //             .param(initial_guess)
    //             .max_iters(MAX_CG_ITER)
    //             .target_cost(TARGET_CG_COST)
    //     })
    //     .add_observer(observer, ObserverMode::NewBest)
    //     .run()
    // {
    //     Ok(r) => r,
    //     Err(err) => {
    //         panic!("Conjugate Gradient error: {err}");
    //     }
    // };
}
