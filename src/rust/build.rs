use std::env;
use std::fs;
use std::path::Path;

fn main() {
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let out_dir = env::var("OUT_DIR").expect("OUT_DIR not set");
    let path = Path::new(&out_dir).join("threads.rs");
    let contents = format!("pub const DEFAULT_THREADS: usize = {};\n", threads);
    fs::write(&path, contents).expect("Failed to write threads.rs");
    println!("cargo:rerun-if-changed=build.rs");
}
