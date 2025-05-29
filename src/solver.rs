use na::DMatrix;

struct FluidSimulation {
    object_mask: DMatrix<f32>,
    u: [DMatrix<f32>; 2],
    p: DMatrix<f32>,
}
