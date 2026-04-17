use std::{path::Path, time::Instant};

use anyhow::Result;
use clap::Parser;
use integer_programming::{instance::IPInstance, solver::BNCSolver};

#[derive(Parser)]
struct Args {
    instance_path: String,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let start = Instant::now();
    let instance = IPInstance::from_path(&args.instance_path)?;
    let mut solver = BNCSolver::new(instance);
    let ip_opt = solver.solve();
    let elapsed = start.elapsed();

    let fname = Path::new(&args.instance_path)
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();
    let res_str = format!(
        "{{\"Instance\": \"{}\", \"Time\": {:.2}, \"Result\": {}, \"Solution\": \"OPT\"}}",
        fname,
        elapsed.as_secs_f64(),
        ip_opt
    );
    println!("{}", res_str);

    Ok(())
}
