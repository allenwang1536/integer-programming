use std::{env, path::PathBuf};

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());

    let project_root = manifest_dir
        .parent()
        .expect("CARGO_MANIFEST_DIR has no parent")
        .to_path_buf();
    let default_ortools = project_root.join(".vendor/or-tools/9.15.6755");
    let ortools_root = env::var("ORTOOLS_ROOT")
        .map(PathBuf::from)
        .unwrap_or(default_ortools);

    println!("cargo:rerun-if-changed=cpp/CMakeLists.txt");
    println!("cargo:rerun-if-changed=cpp/ortools_shim.cc");
    println!("cargo:rerun-if-env-changed=ORTOOLS_ROOT");

    let dst = cmake::Config::new("cpp")
        .define("CMAKE_BUILD_TYPE", "Release")
        .define("CMAKE_PREFIX_PATH", &ortools_root)
        .build();

    let built_lib_dir = dst.join("lib");
    let ortools_lib_dir = ortools_root.join("lib");

    println!("cargo:rustc-link-search=native={}", built_lib_dir.display());
    println!("cargo:rustc-link-lib=dylib=ortools_shim");

    println!(
        "cargo:rustc-link-arg=-Wl,-rpath,{}",
        built_lib_dir.display()
    );
    println!(
        "cargo:rustc-link-arg=-Wl,-rpath,{}",
        ortools_lib_dir.display()
    );
}

