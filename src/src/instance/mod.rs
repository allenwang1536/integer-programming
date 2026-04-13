pub mod parse;

pub struct IPInstance {
    pub num_tests: usize,    // N
    pub num_diseases: usize, // M
    pub test_costs: Vec<u8>,
    pub test_matrix: Vec<u8>, // N x M
}

impl IPInstance {
    pub fn at(&self, t: usize, d: usize) -> u8 {
        self.test_matrix[t * self.num_diseases + d]
    }
}
