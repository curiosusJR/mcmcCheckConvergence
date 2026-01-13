use std::env;
use std::fs;
use std::path::PathBuf;

use mcmcCheckConvergence::core::{check_convergence, Control};

fn get_arg(args: &[String], flag: &str) -> Option<String> {
    args.iter()
        .position(|a| a == flag)
        .and_then(|idx| args.get(idx + 1))
        .cloned()
}

fn has_flag(args: &[String], flag: &str) -> bool {
    args.iter().any(|a| a == flag)
}

fn json_escape(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    for ch in input.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            _ => out.push(ch),
        }
    }
    out
}

fn tsv_escape(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    for ch in input.chars() {
        match ch {
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            _ => out.push(ch),
        }
    }
    out
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || has_flag(&args, "--help") || has_flag(&args, "-h") {
        println!("Usage:");
        println!("  convergence_cli --files f1,f2[,f3,f4] [options]");
        println!("  convergence_cli --path /path/to/output [options]");
        println!();
        println!("Inputs:");
        println!("  --files <f1,f2,...>   Comma-separated list of input files.");
        println!("                        RevBayes expects paired _run_1/_run_2");
        println!("                        .log and .trees files.");
        println!("  --path <dir>          Directory scan for .log/.trees files.");
        println!("                        RevBayes expects _run_1/_run_2 stems.");
        println!();
        println!("Options:");
        println!("  --format <revbayes>   Input format (currently only revbayes).");
        println!("  --burnin <value>      Fraction (0-1) or percent (> 1).");
        println!("  --precision <float>   ESS precision threshold.");
        println!("  --tracer <true|false> Enable ESS tracer output.");
        println!("  --namesToExclude <re> Regex for column names to ignore.");
        println!("  -j, --threads <n>     Number of threads for Rust parallelism.");
        println!("  --fast-splits         Use a faster bitset-based split backend.");
        println!("  --continuous-only     Only analyze .log files.");
        println!("  --json                Emit a single JSON object.");
        println!("  --tsv                 Emit key/value rows.");
        println!("  --message-only        Print only the full message.");
        println!("  --quiet               Print only converged and burnin.");
        println!("  --help, -h            Show this help text.");
        println!();
        println!("Examples:");
        println!("  convergence_cli --path tests/test_1");
        println!("  convergence_cli --path tests/test_1 --continuous-only");
        println!("  convergence_cli --files a.log,b.log,a.trees,b.trees --json");
        return;
    }

    let path = get_arg(&args, "--path");
    let files = get_arg(&args, "--files");
    let format = get_arg(&args, "--format").unwrap_or_else(|| "revbayes".to_string());
    let burnin = get_arg(&args, "--burnin").and_then(|v| v.parse::<f64>().ok());
    let precision = get_arg(&args, "--precision").and_then(|v| v.parse::<f64>().ok());
    let tracer = get_arg(&args, "--tracer").map(|v| matches!(v.as_str(), "true" | "t" | "1" | "yes"));
    let names_to_exclude = get_arg(&args, "--namesToExclude");
    let threads = get_arg(&args, "--threads")
        .or_else(|| get_arg(&args, "-j"))
        .and_then(|v| v.parse::<usize>().ok());
    let fast_splits = has_flag(&args, "--fast-splits");
    let quiet = has_flag(&args, "--quiet");
    let message_only = has_flag(&args, "--message-only");
    let json_out = has_flag(&args, "--json");
    let tsv_out = has_flag(&args, "--tsv");
    let continuous_only = has_flag(&args, "--continuous-only");

    let list_files: Vec<String> = if let Some(path) = path {
        let mut files = Vec::new();
        let dir = PathBuf::from(path);
        if let Ok(entries) = fs::read_dir(&dir) {
            for entry in entries.flatten() {
                let p = entry.path();
                if let Some(ext) = p.extension().and_then(|e| e.to_str()) {
                    if ext == "log" || ext == "trees" {
                        files.push(p.to_string_lossy().to_string());
                    }
                }
            }
        }
        files.sort();
        if format == "revbayes" {
            let mut run_files = Vec::new();
            for f in &files {
                let path = PathBuf::from(f);
                let stem = path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("");
                if stem.ends_with("_run_1") || stem.ends_with("_run_2") {
                    run_files.push(f.clone());
                }
            }
            let mut log_count = run_files.iter().filter(|f| f.ends_with(".log")).count();
            let mut tree_count = run_files.iter().filter(|f| f.ends_with(".trees")).count();
            if continuous_only {
                run_files.retain(|f| f.ends_with(".log"));
                tree_count = 0;
                log_count = run_files.len();
            }
            if continuous_only {
                if log_count != 2 {
                    eprintln!(
                        "convergence_cli: expected _run_1/_run_2 .log files in --path"
                    );
                    std::process::exit(1);
                }
            } else if log_count != 2 || tree_count != 2 {
                eprintln!(
                    "convergence_cli: expected _run_1/_run_2 .log and .trees files in --path"
                );
                std::process::exit(1);
            }
            run_files
        } else {
            files
        }
    } else if let Some(files) = files {
        let mut out: Vec<String> = files.split(',').map(|s| s.trim().to_string()).collect();
        if continuous_only {
            out.retain(|f| f.ends_with(".log"));
            if out.is_empty() {
                eprintln!("convergence_cli: no .log files provided for --continuous-only");
                std::process::exit(1);
            }
        }
        out.sort();
        out
    } else {
        eprintln!("Provide --path or --files.");
        std::process::exit(1);
    };

    let mut control = Control::default();
    if let Some(val) = burnin {
        control.burnin = val;
    }
    if let Some(val) = precision {
        control.precision = val;
    }
    if let Some(val) = tracer {
        control.tracer = val;
    }
    if let Some(val) = names_to_exclude {
        control.names_to_exclude = val;
    }
    if let Some(val) = threads {
        if val > 0 {
            control.threads = Some(val);
        }
    }
    if fast_splits {
        control.fast_splits = true;
    }
    control.emit_logs = false;

    let result = match check_convergence(&list_files, &format, &control) {
        Ok(r) => r,
        Err(err) => {
            eprintln!("convergence_cli: {}", err);
            std::process::exit(1);
        }
    };

    if json_out {
        println!(
            "{{\"converged\":{},\"burnin\":{},\"message\":\"{}\",\"message_complete\":\"{}\"}}",
            result.converged,
            result.burnin,
            json_escape(&result.message),
            json_escape(&result.message_complete)
        );
        return;
    }

    if tsv_out {
        println!("converged\t{}", result.converged);
        println!("burnin\t{}", result.burnin);
        println!("message\t{}", tsv_escape(&result.message));
        println!(
            "message_complete\t{}",
            tsv_escape(&result.message_complete)
        );
        return;
    }

    if message_only {
        println!("{}", result.message_complete);
        return;
    }

    if quiet {
        println!("converged: {}", result.converged);
        println!("burnin: {}", result.burnin);
        return;
    }

    println!("converged: {}", result.converged);
    println!("burnin: {}", result.burnin);
    println!("message:");
    println!("{}", result.message_complete);
}
