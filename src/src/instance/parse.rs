use std::{fmt::Display, fs};

use anyhow::{bail, Result};

use crate::instance::IPInstance;

impl Display for IPInstance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Instance with {} tests, {} diseases",
            self.num_tests, self.num_diseases
        )
    }
}

impl IPInstance {
    pub fn from_path(path: &String) -> Result<Self> {
        let bytes = fs::read(path)?;
        let s = std::str::from_utf8(&bytes)?;
        let mut it = s.split_ascii_whitespace();

        let num_tests: usize = next(&mut it)?;
        let num_diseases: usize = next(&mut it)?;

        let mut test_costs = Vec::with_capacity(num_tests);
        for _ in 0..num_tests {
            test_costs.push(next(&mut it)?);
        }

        let mut test_matrix = Vec::with_capacity(num_tests * num_diseases);
        for _ in 0..num_tests * num_diseases {
            let x: u8 = next(&mut it)?;
            if x > 1 {
                bail!("matrix entry must be 0 or 1, got {x}");
            }
            test_matrix.push(x);
        }

        if it.next().is_some() {
            bail!("extra trailing data");
        }

        Ok(Self {
            num_tests,
            num_diseases,
            test_costs,
            test_matrix,
        })
    }
}

fn next<T: std::str::FromStr>(it: &mut std::str::SplitAsciiWhitespace<'_>) -> Result<T>
where
    T::Err: Send + Sync + std::error::Error + 'static,
{
    Ok(it.next().unwrap().parse()?)
}
