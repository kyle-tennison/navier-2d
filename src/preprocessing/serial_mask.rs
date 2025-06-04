use na::DMatrix;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone)]
pub struct SerialMask {
    data: Vec<u8>,
    nrows: usize,
    ncols: usize,
}

impl SerialMask {
    pub fn from_mask(mask: &DMatrix<bool>) -> Self {
        let (nrows, ncols) = mask.shape();
        let bool_data = mask.data.as_slice().to_vec();

        Self {
            data: bool_data.iter().map(|b| (*b) as u8).collect(),
            nrows,
            ncols,
        }
    }

    pub fn to_mask(&self) -> DMatrix<bool> {
        DMatrix::from_iterator(
            self.nrows,
            self.ncols,
            self.data.clone().into_iter().map(|b| b != 0),
        )
    }
}

#[cfg(test)]
mod tests {
    use na::DMatrix;
    use rand::Rng;
    

    use crate::preprocessing::serial_mask::SerialMask;

    #[test]
    pub fn test_save_load() {
        let mut rng = rand::rng();

        let origional_mask: DMatrix<bool> = DMatrix::from_fn(5, 5, |_, _| rng.random_bool(0.5));

        let serial_mask = SerialMask::from_mask(&origional_mask);

        let serialized = serde_json::to_string(&serial_mask).unwrap();

        println!("Serialized mask:\n\n{serialized}");

        let deserialized: SerialMask = serde_json::from_str(&serialized).unwrap();

        let final_mask = deserialized.to_mask();

        assert_eq!(origional_mask, final_mask);
    }
}
