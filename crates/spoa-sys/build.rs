use std::fs;
use walkdir::WalkDir;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");

    let use_simde = if std::env::var("CARGO_FEATURE_SIMDE").is_ok() {
        "ON"
    } else {
        "OFF"
    };

    let dst = cmake::Config::new("spoa")
        .define("spoa_install", "OFF")
        .define("spoa_build_exe", "OFF")
        .define("spoa_build_tests", "OFF")
        .define("spoa_optimize_for_native", "OFF")
        .define("spoa_use_simde", use_simde)
        .build_target("spoa")
        .build();
    let build_dir = dst.join("build");
    let mut copied = false;
    for dir in ["lib", "lib64"] {
        let src = build_dir.join(dir).join("libspoa.a");
        if src.exists() {
            let out_dir = dst.join(dir);
            let _ = fs::create_dir_all(&out_dir);
            fs::copy(&src, out_dir.join("libspoa.a"))
                .expect("failed to copy libspoa.a into build output");
            copied = true;
        }
    }
    if !copied {
        panic!("libspoa.a not found under {}", build_dir.display());
    }
    let lib64 = dst.join("lib64");
    let lib = dst.join("lib");
    println!("cargo:rustc-link-search=native={}", lib64.display());
    println!("cargo:rustc-link-search=native={}", lib.display());
    println!("cargo:rustc-link-lib=static=spoa");

    println!("cargo:rerun-if-changed=spoa/CMakeLists.txt");
    for entry in WalkDir::new("spoa").into_iter().filter_map(Result::ok) {
        if let Some(ext) = entry.path().extension() {
            if ext == "hpp" || ext == "cpp" {
                println!("cargo:rerun-if-changed={}", entry.path().display());
            }
        }
    }

    println!("cargo:rerun-if-changed=include/bindings.hpp");
    println!("cargo:rerun-if-changed=src/bindings.cpp");
    println!("cargo:rerun-if-changed=src/lib.rs");

    cxx_build::bridge("src/lib.rs")
        .flag("-std=c++11")
        .file("src/bindings.cpp")
        .compile("spoa_sys");
}
