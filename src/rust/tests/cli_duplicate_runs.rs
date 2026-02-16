use std::process::Command;

fn cli_path() -> std::path::PathBuf {
    if let Ok(bin) = std::env::var("CARGO_BIN_EXE_convergence_cli") {
        return std::path::PathBuf::from(bin);
    }
    let mut path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("debug");
    path.push("convergence_cli");
    if cfg!(windows) {
        path.set_extension("exe");
    }
    path
}

#[test]
fn cli_fails_on_identical_runs() {
    let dir = std::env::temp_dir();
    let file_name = format!(
        "mcmc_cli_identical_reps_{}_{}.log",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    );
    let log_path = dir.join(file_name);
    let data = "Iteration\tReplicate_ID\tx\n1\t0\t0.1\n2\t0\t0.2\n1\t1\t0.1\n2\t1\t0.2\n";
    std::fs::write(&log_path, data).unwrap();

    let output = Command::new(cli_path())
        .arg("--files")
        .arg(log_path.to_string_lossy().to_string())
        .arg("--format")
        .arg("revbayes")
        .output()
        .expect("failed to run convergence_cli");

    let _ = std::fs::remove_file(log_path);

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Detected identical MCMC runs"));
}
