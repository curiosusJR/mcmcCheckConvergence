use crate::{collect_clades, parse_newick};
use rayon::prelude::*;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug)]
pub struct Table {
    pub headers: Vec<String>,
    pub columns: Vec<Vec<f64>>,
}

impl Table {
    pub(crate) fn nrows(&self) -> usize {
        self.columns.get(0).map(|c| c.len()).unwrap_or(0)
    }

    pub(crate) fn select_columns(&self, keep: &[bool]) -> Table {
        let mut headers = Vec::new();
        let mut columns = Vec::new();
        for (i, k) in keep.iter().enumerate() {
            if *k {
                headers.push(self.headers[i].clone());
                columns.push(self.columns[i].clone());
            }
        }
        Table { headers, columns }
    }

    pub(crate) fn slice_rows(&self, start: usize, end: usize) -> Table {
        let mut columns = Vec::with_capacity(self.columns.len());
        for col in &self.columns {
            let s = start.min(col.len());
            let e = end.min(col.len());
            columns.push(col[s..e].to_vec());
        }
        Table {
            headers: self.headers.clone(),
            columns,
        }
    }

    pub(crate) fn column_index(&self, name: &str) -> Option<usize> {
        self.headers.iter().position(|h| h == name)
    }
}

#[derive(Clone, Debug)]
pub struct Run {
    pub trees: Vec<String>,
    pub ptable: Table,
}

#[derive(Clone)]
pub struct Control {
    pub tracer: bool,
    pub burnin: f64,
    pub precision: f64,
    pub names_to_exclude: String,
    pub emit_logs: bool,
}

impl Default for Control {
    fn default() -> Self {
        Self {
            tracer: true,
            burnin: 0.0,
            precision: 0.01,
            names_to_exclude: "br_lens|bl|Iteration|Likelihood|Posterior|Prior|Gen|LnL|LnPr|state|joint|prior|likelihood|time|loglik|iter|topo|Replicate_ID|Sample|posterior|it".to_string(),
            emit_logs: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct ConvergenceResult {
    pub converged: bool,
    pub burnin: f64,
    pub message: String,
    pub message_complete: String,
    pub failed_names: Option<String>,
    pub fail_msgs: Vec<String>,
    pub tree_exclude_high: Vec<Vec<String>>,
    pub tree_exclude_low: Vec<Vec<String>>,
    pub tree_ess: Vec<Vec<(String, f64)>>,
    pub tree_freqs: Vec<Vec<(String, f64)>>,
    pub tree_compare_freqs: Vec<Vec<(String, f64)>>,
    pub tree_compare: Vec<Vec<(String, f64)>>,
    pub cont_means: Vec<(String, f64)>,
    pub cont_ess: Vec<Vec<(String, f64)>>,
    pub cont_compare: Vec<Vec<(String, f64)>>,
    pub cont_exclude: Vec<Vec<String>>,
}

struct ContRunResult {
    filtered: Table,
    exclude: Vec<String>,
    ess: Vec<(String, f64)>,
    ess_fail_count: usize,
}

fn compute_cont_run(
    run: &Run,
    names_re: &Regex,
    minimum_ess: f64,
    tracer: bool,
) -> ContRunResult {
    let mut exclude = Vec::new();
    let mut filtered_headers = Vec::new();
    let mut filtered_columns = Vec::new();
    let mut ess_vec = Vec::new();
    let mut ess_fail_count = 0;
    for (idx, col) in run.ptable.columns.iter().enumerate() {
        let mean = col.iter().sum::<f64>() / col.len().max(1) as f64;
        let var = col
            .iter()
            .map(|v| (v - mean).powi(2))
            .sum::<f64>()
            / col.len().max(1) as f64;
        if var == 0.0 {
            exclude.push(run.ptable.headers[idx].clone());
            continue;
        }
        if names_re.is_match(&run.ptable.headers[idx]) {
            continue;
        }
        let ess = if tracer { ess_tracer(col) } else { ess_tracer(col) };
        if ess < minimum_ess {
            ess_fail_count += 1;
        }
        filtered_headers.push(run.ptable.headers[idx].clone());
        filtered_columns.push(col.clone());
        ess_vec.push((run.ptable.headers[idx].clone(), ess));
    }
    let filtered = Table {
        headers: filtered_headers,
        columns: filtered_columns,
    };

    ContRunResult {
        filtered,
        exclude,
        ess: ess_vec,
        ess_fail_count,
    }
}

fn columns_with_nan(table: &Table) -> Vec<String> {
    let mut bad = Vec::new();
    for (idx, col) in table.columns.iter().enumerate() {
        if col.iter().any(|v| !v.is_finite()) {
            bad.push(table.headers.get(idx).cloned().unwrap_or_default());
        }
    }
    bad
}

fn normalize_header_name(name: &str) -> String {
    let mut out = String::with_capacity(name.len());
    for ch in name.chars() {
        match ch {
            '[' | ']' => out.push('.'),
            _ => out.push(ch),
        }
    }
    out
}

pub(crate) fn parse_table(path: &str, skip: usize, delim: u8) -> Result<Table, String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .ok_or_else(|| "Empty file".to_string())?
        .map_err(|e| e.to_string())?;
    let headers: Vec<String> = header
        .split(char::from(delim))
        .map(|s| normalize_header_name(s.trim().trim_matches('"')))
        .collect();
    let mut columns: Vec<Vec<f64>> = vec![Vec::new(); headers.len()];

    for (idx, line) in lines.enumerate() {
        let line = line.map_err(|e| e.to_string())?;
        if idx < skip {
            continue;
        }
        let parts: Vec<&str> = line.split(char::from(delim)).collect();
        if parts.len() != headers.len() {
            continue;
        }
        for (i, part) in parts.iter().enumerate() {
            let val = part.trim().trim_matches('"').parse::<f64>().unwrap_or(f64::NAN);
            columns[i].push(val);
        }
    }

    Ok(Table { headers, columns })
}

pub(crate) fn parse_revbayes_trees(path: &str) -> Result<(Vec<String>, Table), String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .ok_or_else(|| "Empty trees file".to_string())?
        .map_err(|e| e.to_string())?;
    let headers: Vec<String> = header
        .split('\t')
        .map(|s| normalize_header_name(s.trim()))
        .collect();
    let first_row = lines
        .next()
        .ok_or_else(|| "Missing data rows in trees file".to_string())?
        .map_err(|e| e.to_string())?;
    let first_parts: Vec<&str> = first_row.split('\t').collect();
    if first_parts.len() != headers.len() {
        return Err("Invalid trees file row".to_string());
    }
    let mut topo_idx = None;
    for (i, part) in first_parts.iter().enumerate() {
        if part.contains('(') && part.contains(';') {
            let cleaned = strip_bracket_blocks(part);
            if parse_newick(&cleaned).is_ok() {
                topo_idx = Some(i);
                break;
            }
        }
    }
    let topo_idx = topo_idx.ok_or_else(|| "Missing topology column".to_string())?;

    let mut trees = Vec::with_capacity(1024);
    let mut columns: Vec<Vec<f64>> = vec![Vec::new(); headers.len() - 1];
    let mut num_headers = Vec::new();
    for (i, h) in headers.iter().enumerate() {
        if i != topo_idx {
            num_headers.push(h.clone());
        }
    }

    let mut handle_row = |line: &str| {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != headers.len() {
            return;
        }
        let raw_tree = parts[topo_idx].trim();
        let cleaned = strip_bracket_blocks(raw_tree);
        trees.push(cleaned);
        let mut col_i = 0;
        for (i, part) in parts.iter().enumerate() {
            if i == topo_idx {
                continue;
            }
            let val = part.trim().parse::<f64>().unwrap_or(f64::NAN);
            columns[col_i].push(val);
            col_i += 1;
        }
    };

    handle_row(&first_row);
    for line in lines {
        let line = line.map_err(|e| e.to_string())?;
        handle_row(&line);
    }

    Ok((
        trees,
        Table {
            headers: num_headers,
            columns,
        },
    ))
}

fn strip_nexus_prefix(s: &str) -> &str {
    let trimmed = s.trim();
    if let Some(idx) = trimmed.find('=') {
        return trimmed[idx + 1..].trim();
    }
    trimmed
}

fn strip_bracket_blocks(s: &str) -> String {
    let mut out = String::new();
    let mut chars = s.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch == '[' {
            while let Some(c) = chars.next() {
                if c == ']' {
                    break;
                }
            }
            continue;
        }
        out.push(ch);
    }
    out
}

fn parse_nexus_trees(path: &str) -> Result<Vec<String>, String> {
    let mut trees = Vec::with_capacity(1024);
    let file = File::open(path).map_err(|e| e.to_string())?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let lower = line.trim().to_lowercase();
        if lower.starts_with("tree") {
            let candidate = strip_nexus_prefix(&line);
            let stripped = strip_bracket_blocks(candidate);
            let cand = stripped.trim();
            if parse_newick(cand).is_ok() {
                trees.push(cand.to_string());
            }
        }
    }
    if trees.is_empty() {
        return Err("No trees found in nexus file".to_string());
    }
    Ok(trees)
}

fn parse_newick_trees(path: &str) -> Result<Vec<String>, String> {
    let mut trees = Vec::with_capacity(1024);
    let file = File::open(path).map_err(|e| e.to_string())?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        let cand = line.trim();
        if cand.is_empty() {
            continue;
        }
        if parse_newick(cand).is_ok() {
            trees.push(cand.to_string());
        }
    }
    if trees.is_empty() {
        return Err("No trees found in newick file".to_string());
    }
    Ok(trees)
}

pub(crate) fn merge_tables(rb: &Table, log: &Table) -> Table {
    let mut headers = log.headers.clone();
    let mut columns = log.columns.clone();
    for (i, h) in rb.headers.iter().enumerate() {
        if !headers.contains(h) {
            headers.insert(0, h.clone());
            columns.insert(0, rb.columns[i].clone());
        }
    }
    Table { headers, columns }
}

pub(crate) fn remove_burnin(runs: &mut [Run], burnin: f64) -> Result<(), String> {
    if burnin < 0.0 {
        return Err("Invalid burnin value".to_string());
    }
    if let Some(run) = runs.first() {
        let n_trees = run.trees.len();
        let n_rows = run.ptable.nrows();
        if n_trees > 0 {
            if burnin >= n_trees as f64 {
                return Err("Burnin larger than iterations in file".to_string());
            }
        } else if n_rows > 0 && burnin >= n_rows as f64 {
            return Err("Burnin larger than iterations in file".to_string());
        }
    }
    for run in runs.iter_mut() {
        let n_trees = run.trees.len();
        let n_rows = run.ptable.nrows();
        if burnin >= 1.0 {
            let discard = if n_trees > 0 {
                ((burnin / 100.0) * n_trees as f64).ceil() as usize
            } else {
                ((burnin / 100.0) * n_rows as f64).ceil() as usize
            };
            if n_trees > 0 {
                run.trees = run.trees.split_off(discard.min(n_trees));
            }
            if n_rows > 0 {
                run.ptable = run.ptable.slice_rows(discard.min(n_rows), n_rows);
            }
        } else if burnin > 0.0 {
            let discard = if n_trees > 0 {
                (burnin * n_trees as f64).ceil() as usize
            } else {
                (burnin * n_rows as f64).ceil() as usize
            };
            if n_trees > 0 {
                run.trees = run.trees.split_off(discard.min(n_trees));
            }
            if n_rows > 0 {
                run.ptable = run.ptable.slice_rows(discard.min(n_rows), n_rows);
            }
        }
    }
    Ok(())
}

pub(crate) fn ks_statistic(a: &[f64], b: &[f64]) -> f64 {
    let mut a_sorted = a.to_vec();
    let mut b_sorted = b.to_vec();
    a_sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
    b_sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
    ks_statistic_sorted(&a_sorted, &b_sorted)
}

fn ks_statistic_cached(a: &[f64], b: &[f64], buf_a: &mut Vec<f64>, buf_b: &mut Vec<f64>) -> f64 {
    buf_a.clear();
    buf_a.extend_from_slice(a);
    buf_b.clear();
    buf_b.extend_from_slice(b);
    buf_a.sort_by(|x, y| x.partial_cmp(y).unwrap());
    buf_b.sort_by(|x, y| x.partial_cmp(y).unwrap());
    ks_statistic_sorted(buf_a, buf_b)
}

fn ks_statistic_sorted(a_sorted: &[f64], b_sorted: &[f64]) -> f64 {
    let n1 = a_sorted.len() as f64;
    let n2 = b_sorted.len() as f64;
    let mut i = 0usize;
    let mut j = 0usize;
    let mut d: f64 = 0.0;
    while i < a_sorted.len() && j < b_sorted.len() {
        let av = a_sorted[i];
        let bv = b_sorted[j];
        if av < bv {
            i += 1;
        } else if av > bv {
            j += 1;
        } else {
            let val = av;
            while i < a_sorted.len() && a_sorted[i] == val {
                i += 1;
            }
            while j < b_sorted.len() && b_sorted[j] == val {
                j += 1;
            }
        }
        let cdf1 = i as f64 / n1;
        let cdf2 = j as f64 / n2;
        d = d.max((cdf1 - cdf2).abs());
    }
    d
}

pub(crate) fn ess_tracer(input: &[f64]) -> f64 {
    let samples = input.len();
    if samples < 2 {
        return samples as f64;
    }
    let max_lag_limit = 2000usize;
    let max_lag = std::cmp::min(samples - 1, max_lag_limit);
    let mean = input.iter().sum::<f64>() / samples as f64;
    let mut gamma_stat = vec![0.0; max_lag];
    let mut var_stat = 0.0;
    for lag in 0..max_lag {
        let mut acc = 0.0;
        for j in 0..(samples - lag) {
            let del1 = input[j] - mean;
            let del2 = input[j + lag] - mean;
            acc += del1 * del2;
        }
        gamma_stat[lag] = acc / (samples - lag) as f64;
        if lag == 0 {
            var_stat = gamma_stat[0];
        } else if lag % 2 == 0 {
            let pair_sum = gamma_stat[lag - 1] + gamma_stat[lag];
            if pair_sum > 0.0 {
                var_stat += 2.0 * pair_sum;
            } else {
                break;
            }
        }
    }
    let act = var_stat / gamma_stat[0];
    samples as f64 / act
}

pub(crate) fn min_ess(per: f64) -> f64 {
    (1.0 / (per * 4.0)).powi(2)
}

pub(crate) fn ks_threshold(alpha: f64, ess: f64) -> f64 {
    let c_alpha = (-(alpha / 2.0).ln() * 0.5).sqrt();
    c_alpha * (2.0 / ess).sqrt()
}

pub(crate) fn expected_diff_splits(ess: usize) -> (Vec<f64>, Vec<f64>) {
    let mut log_fact = vec![0.0; ess + 1];
    for i in 1..=ess {
        log_fact[i] = log_fact[i - 1] + (i as f64).ln();
    }
    let choose = |n: usize, k: usize| -> f64 {
        (log_fact[n] - log_fact[k] - log_fact[n - k]).exp()
    };
    let mut probs_out = Vec::new();
    let mut thresh_out = Vec::new();
    let mut p: f64 = 0.01;
    while p < 1.0 {
        let mut probs = vec![0.0; ess + 1];
        for f1 in 0..=ess {
            let p1 = choose(ess, f1) * p.powi(f1 as i32) * (1.0 - p).powi((ess - f1) as i32);
            for f2 in 0..=ess {
                let p2 = choose(ess, f2) * p.powi(f2 as i32) * (1.0 - p).powi((ess - f2) as i32);
                let diff = (f1 as i32 - f2 as i32).abs() as usize;
                probs[diff] += p1 * p2;
            }
        }
        let mut cdf = 0.0;
        let mut thresh = 0.0;
        for (i, pr) in probs.iter().enumerate() {
            cdf += pr;
            if cdf >= 0.95 {
                thresh = i as f64 / ess as f64;
                break;
            }
        }
        probs_out.push(p);
        thresh_out.push(thresh);
        p += 0.01;
    }
    (probs_out, thresh_out)
}

pub(crate) fn expected_diff_splits_cached(ess: usize) -> (Vec<f64>, Vec<f64>) {
    if ess == 125 {
        return parse_expected_diff(include_str!("../data/expectedDiff_125.tsv"));
    }
    if ess == 625 {
        return parse_expected_diff(include_str!("../data/expectedDiff_625.tsv"));
    }
    expected_diff_splits(ess)
}

fn parse_expected_diff(content: &str) -> (Vec<f64>, Vec<f64>) {
    let mut probs = Vec::new();
    let mut thresh = Vec::new();
    for line in content.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }
        let p = parts[0].trim().parse::<f64>().unwrap_or(0.0);
        let t = parts[1].trim().parse::<f64>().unwrap_or(0.0);
        probs.push(p);
        thresh.push(t);
    }
    (probs, thresh)
}

pub(crate) fn window_bounds(len: usize) -> (usize, usize) {
    if len == 0 {
        return (0, 0);
    }
    let second = ((0.2 * len as f64).floor() as usize).max(1);
    let fifth = ((0.8 * len as f64).floor() as usize).max(1);
    let first_end = second.min(len);
    let start2 = if fifth == 0 { 0 } else { fifth - 1 };
    let start2 = start2.min(len - 1);
    (first_end, start2)
}

pub(crate) fn split_windows<T: Clone>(vals: &[T]) -> (Vec<T>, Vec<T>) {
    let (first_end, start2) = window_bounds(vals.len());
    if vals.is_empty() {
        return (Vec::new(), Vec::new());
    }
    (vals[0..first_end].to_vec(), vals[start2..].to_vec())
}

pub(crate) fn clade_stats_from_sets(
    sets: &[HashSet<String>],
) -> (Vec<String>, HashMap<String, usize>) {
    let mut order = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut counts: HashMap<String, usize> = HashMap::new();
    for clades in sets {
        for clade in clades {
            if seen.insert(clade.clone()) {
                order.push(clade.clone());
            }
            *counts.entry(clade.clone()).or_insert(0) += 1;
        }
    }
    (order, counts)
}

pub(crate) fn clade_sets_for_trees(trees: &[String]) -> Vec<HashSet<String>> {
    trees
        .iter()
        .filter_map(|t| parse_newick(t).ok())
        .filter_map(|(root, nodes)| collect_clades(root, &nodes).ok())
        .map(|clades| clades.into_iter().collect())
        .collect()
}

pub(crate) struct CladeStats {
    pub(crate) order: Vec<usize>,
    pub(crate) counts: Vec<usize>,
    pub(crate) sets: Vec<HashSet<usize>>,
    pub(crate) names: Vec<String>,
    pub(crate) index: HashMap<String, usize>,
}

pub(crate) fn clade_stats_ids(trees: &[String]) -> CladeStats {
    let mut names = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut counts: Vec<usize> = Vec::new();
    let mut sets = Vec::with_capacity(trees.len());
    for t in trees {
        if let Ok((root, nodes)) = parse_newick(t) {
            if let Ok(clades) = collect_clades(root, &nodes) {
                let mut unique: HashSet<usize> = HashSet::new();
                for clade in clades {
                    let id = if let Some(id) = index.get(&clade) {
                        *id
                    } else {
                        let id = names.len();
                        names.push(clade.clone());
                        index.insert(clade, id);
                        counts.push(0);
                        id
                    };
                    if unique.insert(id) {
                        counts[id] += 1;
                    }
                }
                sets.push(unique);
            }
        }
    }
    let order: Vec<usize> = (0..names.len()).collect();
    CladeStats {
        order,
        counts,
        sets,
        names,
        index,
    }
}

fn available_threads() -> usize {
    if let Ok(val) = std::env::var("MCMC_CONVERGENCE_THREADS") {
        if let Ok(parsed) = val.trim().parse::<usize>() {
            return parsed.max(1);
        }
    }
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

fn tree_stats_per_run(
    runs: &[Run],
) -> Vec<(CladeStats, Vec<String>)> {
    let run_count = runs.len();
    if run_count <= 1 || available_threads() <= 1 {
        return runs
            .iter()
            .map(|run| {
                let tips = run
                    .trees
                    .get(0)
                    .map(|t| tree_tips_from_newick_str(t))
                    .unwrap_or_default();
                let stats = clade_stats_ids(&run.trees);
                (stats, tips)
            })
            .collect();
    }

    runs.par_iter()
        .map(|run| {
            let tips = run
                .trees
                .get(0)
                .map(|t| tree_tips_from_newick_str(t))
                .unwrap_or_default();
            let stats = clade_stats_ids(&run.trees);
            (stats, tips)
        })
        .collect()
}

fn tree_tips_from_newick_str(tree: &str) -> Vec<String> {
    if let Ok((_root, nodes)) = parse_newick(tree) {
        let mut tips = Vec::new();
        for n in nodes {
            if n.children.is_empty() {
                if let Some(label) = n.label {
                    tips.push(label);
                }
            }
        }
        tips.sort();
        tips.dedup();
        return tips;
    }
    Vec::new()
}

pub(crate) fn print_log(emit: bool, msg: &str) {
    if emit {
        println!("[1] \"{}\"", msg);
    }
}

pub fn load_runs(list_files: &[String], format: &str, emit_logs: bool) -> Result<Vec<Run>, String> {
    let (log_ext, tree_ext, delim, skip, tree_type) = match format {
        "revbayes" => (".log", ".trees", b'\t', 0usize, "revbayes"),
        "mb" | "mrbayes" => (".p", ".t", b'\t', 1usize, "nexus"),
        "beast" => (".log", ".trees", b'\t', 2usize, "nexus"),
        "*beast" => (".log", ".trees", b',', 0usize, "nexus"),
        "phylobayes" => (".trace", ".treelist", b'\t', 0usize, "newick"),
        "pyrate" => (".log", ".trees", b'\t', 0usize, "newick"),
        _ => return Err("Provide format!".to_string()),
    };

    let mut tree_files = Vec::new();
    let mut log_files = Vec::new();
    for f in list_files {
        if f.ends_with(tree_ext) {
            tree_files.push(f.clone());
        } else if f.ends_with(log_ext) {
            log_files.push(f.clone());
        }
    }
    if tree_files.is_empty() && log_files.is_empty() {
        return Err("No files to read".to_string());
    }

    let mut runs = Vec::new();
    if !tree_files.is_empty() {
        for (i, tree) in tree_files.iter().enumerate() {
            print_log(
                emit_logs,
                &Path::new(tree).file_name().unwrap().to_string_lossy(),
            );
            print_log(emit_logs, "Reading trees...");
            let trees = match tree_type {
                "revbayes" => {
                    let (t, rb) = parse_revbayes_trees(tree)?;
                    let gens = if let Some(idx) = rb.column_index("Iteration") {
                        if rb.columns[idx].len() >= 2 {
                            (rb.columns[idx][1] - rb.columns[idx][0]).round() as i64
                        } else {
                            1
                        }
                    } else {
                        1
                    };
                    print_log(emit_logs, &format!("{} generations per tree...", gens));
                    print_log(emit_logs, "Unrooting, this may take a while...");
                    let log = log_files.get(i).cloned();
                    let ptable = if let Some(log_path) = log {
                        print_log(
                            emit_logs,
                            &format!(
                                "Reading parameter values from {}",
                                Path::new(&log_path).file_name().unwrap().to_string_lossy()
                            ),
                        );
                        let log_table = parse_table(&log_path, skip, delim)?;
                        merge_tables(&rb, &log_table)
                    } else {
                        rb
                    };
                    print_log(emit_logs, "rerooting trees...");
                    let tips = t
                        .get(0)
                        .map(|s| tree_tips_from_newick_str(s))
                        .unwrap_or_default();
                    if let Some(outgroup) = tips.get(0) {
                        print_log(emit_logs, &format!("Outgroup {}", outgroup));
                    }
                    runs.push(Run { trees: t, ptable });
                    continue;
                }
                "nexus" => parse_nexus_trees(tree)?,
                "newick" => parse_newick_trees(tree)?,
                _ => return Err("Unsupported tree type".to_string()),
            };

            let log = log_files.get(i).cloned();
            let ptable = if let Some(log_path) = log {
                print_log(
                    emit_logs,
                    &format!(
                        "Reading parameter values from {}",
                        Path::new(&log_path).file_name().unwrap().to_string_lossy()
                    ),
                );
                parse_table(&log_path, skip, delim)?
            } else {
                Table {
                    headers: Vec::new(),
                    columns: Vec::new(),
                }
            };
            print_log(emit_logs, "rerooting trees...");
            let tips = trees
                .get(0)
                .map(|s| tree_tips_from_newick_str(s))
                .unwrap_or_default();
            if let Some(outgroup) = tips.get(0) {
                print_log(emit_logs, &format!("Outgroup {}", outgroup));
            }
            runs.push(Run { trees, ptable });
        }
    } else {
        for log_path in log_files.iter() {
            let ptable = parse_table(log_path, skip, delim)?;
            runs.push(Run {
                trees: Vec::new(),
                ptable,
            });
        }
    }
    Ok(runs)
}

pub(crate) fn filter_table(table: &Table, names_to_exclude: &Regex) -> Table {
    let keep: Vec<bool> = table
        .headers
        .iter()
        .map(|h| !names_to_exclude.is_match(h))
        .collect();
    table.select_columns(&keep)
}

pub fn check_convergence(
    list_files: &[String],
    format: &str,
    control: &Control,
) -> Result<ConvergenceResult, String> {
    let mut runs = load_runs(list_files, format, control.emit_logs)?;
    for (idx, run) in runs.iter().enumerate() {
        let bad = columns_with_nan(&run.ptable);
        if !bad.is_empty() {
            return Err(format!(
                "Run {} contains non-numeric or missing values in columns: {}",
                idx + 1,
                bad.join(", ")
            ));
        }
    }
    let names_re = Regex::new(&control.names_to_exclude).map_err(|e| e.to_string())?;

    let mut burnin = control.burnin;
    if burnin > 0.0 {
        remove_burnin(&mut runs, burnin)?;
    }

    let minimum_ess = min_ess(control.precision);
    let minimum_ess_windows = (minimum_ess / 5.0).round() as usize;

    let mut message_list = String::new();
    let mut message_complete = String::new();
    let mut count_decision = 0;

    let mut fail_msgs: Vec<String> = Vec::new();
    let mut failed_names: Vec<String> = Vec::new();

    let mut tree_exclude_high: Vec<Vec<String>> = Vec::new();
    let mut tree_exclude_low: Vec<Vec<String>> = Vec::new();
    let mut tree_ess: Vec<Vec<(String, f64)>> = Vec::new();
    let mut tree_freqs: Vec<Vec<(String, f64)>> = Vec::new();
    let mut tree_compare_freqs: Vec<Vec<(String, f64)>> = Vec::new();
    let mut tree_compare: Vec<Vec<(String, f64)>> = Vec::new();
    let mut cont_means: Vec<(String, f64)> = Vec::new();
    let mut cont_ess: Vec<Vec<(String, f64)>> = Vec::new();
    let mut cont_compare: Vec<Vec<(String, f64)>> = Vec::new();
    let mut cont_exclude: Vec<Vec<String>> = Vec::new();

    if burnin == 0.0 {
        print_log(control.emit_logs, "Calculating burn-in");
        let original_runs = runs.clone();
        while burnin <= 0.5 {
            let mut list_control = 0;
            if !runs.is_empty() && runs[0].ptable.nrows() > 0 {
                for run in runs.iter() {
                    let filtered = filter_table(&run.ptable, &names_re);
                    if filtered.headers.is_empty() {
                        continue;
                    }
                    let (first_end, start2) = window_bounds(filtered.nrows());
                    let mut buf_a = Vec::new();
                    let mut buf_b = Vec::new();
                    for col in filtered.columns.iter() {
                        let d = ks_statistic_cached(&col[0..first_end], &col[start2..], &mut buf_a, &mut buf_b);
                        if d > ks_threshold(0.01, minimum_ess_windows as f64) {
                            list_control += 1;
                            break;
                        }
                    }
                }
            } else if !runs.is_empty() && !runs[0].trees.is_empty() {
                for run in runs.iter() {
                    let sets = clade_sets_for_trees(&run.trees);
                    let (w1, w2) = split_windows(&sets);
                    if w1.is_empty() || w2.is_empty() {
                        continue;
                    }
                    let (order1, counts1) = clade_stats_from_sets(&w1);
                    let (_order2, counts2) = clade_stats_from_sets(&w2);
                    let (probs, thresh) = expected_diff_splits_cached(minimum_ess_windows);
                    let mut any_fail = false;
                    for clade in order1.iter() {
                        let Some(c1) = counts1.get(clade) else { continue };
                        let Some(c2) = counts2.get(clade) else { continue };
                        let f1 = *c1 as f64 / w1.len() as f64;
                        let f2 = *c2 as f64 / w2.len() as f64;
                        let freq = ((f1 + f2) / 2.0 * 100.0).round() / 100.0;
                        if let Some(pos) = probs.iter().position(|p| (*p - freq).abs() < 1e-6) {
                            if (f1 - f2).abs() > thresh[pos] {
                                any_fail = true;
                                break;
                            }
                        }
                    }
                    if any_fail {
                        list_control += 1;
                        break;
                    }
                }
            }
            if list_control > 0 {
                burnin += 0.1;
                runs = original_runs.clone();
                remove_burnin(&mut runs, burnin)?;
            } else {
                break;
            }
        }
        if burnin > 0.5 {
            return Err("Burn-in too large".to_string());
        }
    }

    if !runs.is_empty() && !runs[0].trees.is_empty() {
        print_log(control.emit_logs, "Analyzing tree parameters");
        let mut run_stats = tree_stats_per_run(&runs);

        for (idx, run) in runs.iter().enumerate() {
            let mut ess_fail_count = 0;
            let (ref mut stats, ref tips) = run_stats[idx];
            let n_trees = run.trees.len() as f64;
            let mut ess_splits = Vec::new();
            let mut freqs = Vec::new();
            let mut exclude_high = Vec::new();
            let mut exclude_low = Vec::new();
            if !tips.is_empty() {
                let full = tips.join(" ");
                if let Some(id) = stats.index.get(&full).cloned() {
                    if let Some(pos) = stats.order.iter().position(|&v| v == id) {
                        let item = stats.order.remove(pos);
                        stats.order.insert(0, item);
                    }
                }
            }

            let mut is_split = Vec::with_capacity(stats.sets.len());
            for &id in stats.order.iter() {
                let count = stats.counts[id];
                let freq = count as f64 / n_trees;
                let name = &stats.names[id];
                freqs.push((name.clone(), freq));
                if freq > 0.975 {
                    exclude_high.push(name.clone());
                }
                if freq < 0.025 {
                    exclude_low.push(name.clone());
                }
                if freq <= 0.975 && freq >= 0.025 {
                    is_split.clear();
                    for s in stats.sets.iter() {
                        is_split.push(if s.contains(&id) { 1.0 } else { 0.0 });
                    }
                    let ess = if control.tracer {
                        ess_tracer(&is_split)
                    } else {
                        ess_tracer(&is_split)
                    };
                    if ess < minimum_ess {
                        ess_fail_count += 1;
                    }
                    ess_splits.push((name.clone(), ess));
                }
            }

            if !tips.is_empty() {
                let full = tips.join(" ");
                if let Some(pos) = exclude_high.iter().position(|c| c == &full) {
                    exclude_high.remove(pos);
                }
                exclude_high.insert(0, full);
            }

            tree_exclude_high.push(exclude_high);
            tree_exclude_low.push(exclude_low);
            tree_ess.push(ess_splits.clone());
            tree_freqs.push(freqs);

            if ess_fail_count > 0 {
                fail_msgs.push(format!(
                    "{} splits failed to reach {} for ESS_of_Run_{}",
                    ess_fail_count,
                    minimum_ess as usize,
                    idx + 1
                ));
                failed_names.push(format!("ESS_of_Run_{}", idx + 1));
                count_decision += 1;
            }
        }

        if runs.len() > 1 {
            let (probs, thresholds) = expected_diff_splits_cached(minimum_ess as usize);
            let mut diff_fail = 0;
            for i in 0..(runs.len() - 1) {
                for j in (i + 1)..runs.len() {
                    let (ref stats1, _) = run_stats[i];
                    let (ref stats2, _) = run_stats[j];
                    let mut compare_vals = Vec::new();
                    let mut compare_freqs = Vec::new();
                    for &id1 in stats1.order.iter() {
                        let name = &stats1.names[id1];
                        let Some(&id2) = stats2.index.get(name) else { continue };
                        let c1 = stats1.counts[id1];
                        let c2 = stats2.counts[id2];
                        let f1 = c1 as f64 / runs[i].trees.len() as f64;
                        let f2 = c2 as f64 / runs[j].trees.len() as f64;
                        let freq = ((f1 + f2) / 2.0 * 100.0).round() / 100.0;
                        if let Some(pos) = probs.iter().position(|p| (*p - freq).abs() < 1e-6) {
                            let threshold = thresholds[pos];
                            let diff = (f1 - f2).abs();
                            let val = threshold - diff;
                            if val <= 0.0 {
                                diff_fail += 1;
                            }
                            compare_vals.push((name.clone(), val));
                            compare_freqs.push((name.clone(), freq));
                        }
                    }
                    tree_compare.push(compare_vals);
                    tree_compare_freqs.push(compare_freqs);
                }
            }
            if diff_fail > 0 {
                fail_msgs.push(format!(
                    "{} splits failed the split difference test between runs for Between_Run_1_Run_2",
                    diff_fail
                ));
                failed_names.push("Between_Run_1_Run_2".to_string());
                count_decision += 1;
            }
        }
    }

    if !runs.is_empty() && runs[0].ptable.nrows() > 0 {
        print_log(control.emit_logs, "Analyzing continuous parameters");
        let mut headers = Vec::new();

        let run_count = runs.len();
        let names_re_ref = &names_re;
        let cont_results: Vec<ContRunResult> = if run_count > 1 && available_threads() > 1 {
            runs.par_iter()
                .map(|run| compute_cont_run(run, names_re_ref, minimum_ess, control.tracer))
                .collect()
        } else {
            runs.iter()
                .map(|run| compute_cont_run(run, names_re_ref, minimum_ess, control.tracer))
                .collect()
        };

        let mut filtered_runs: Vec<Table> = Vec::new();
        for (run_idx, res) in cont_results.into_iter().enumerate() {
            cont_exclude.push(res.exclude);
            headers = res.filtered.headers.clone();
            cont_ess.push(res.ess);
            if res.ess_fail_count > 0 {
                fail_msgs.push(format!(
                    "{} parameters failed to reach {} for ESS_run_{}",
                    res.ess_fail_count,
                    minimum_ess as usize,
                    run_idx + 1
                ));
                failed_names.push(format!("ESS_run_{}", run_idx + 1));
                count_decision += 1;
            }
            filtered_runs.push(res.filtered);
        }

        if !headers.is_empty() {
            for (idx, name) in headers.iter().enumerate() {
                let mut sum = 0.0;
                let mut count = 0usize;
                for run in &filtered_runs {
                    for v in &run.columns[idx] {
                        sum += *v;
                        count += 1;
                    }
                }
                let mean = if count == 0 { f64::NAN } else { sum / count as f64 };
                cont_means.push((name.clone(), mean));
            }
        }

        if filtered_runs.len() > 1 {
            let ks_limit = ks_threshold(0.01, minimum_ess);
            let use_parallel = available_threads() > 1 && headers.len() > 1;
            let mut buf_a = Vec::new();
            let mut buf_b = Vec::new();
            for i in 0..(filtered_runs.len() - 1) {
                for j in (i + 1)..filtered_runs.len() {
                    let mut ks_fail_count = 0;
                    let compare_vals = if use_parallel {
                        let mut results: Vec<(usize, f64, bool)> = headers
                            .par_iter()
                            .enumerate()
                            .map_init(|| (Vec::new(), Vec::new()), |(buf_a, buf_b), (idx, _)| {
                                let col1 = &filtered_runs[i].columns[idx];
                                let col2 = &filtered_runs[j].columns[idx];
                                let d = ks_statistic_cached(col1, col2, buf_a, buf_b);
                                let val = ks_limit - d;
                                (idx, val, d > ks_limit)
                            })
                            .collect();
                        results.sort_by_key(|r| r.0);
                        ks_fail_count = results.iter().filter(|r| r.2).count();
                        results
                            .into_iter()
                            .map(|(idx, val, _)| (headers[idx].clone(), val))
                            .collect()
                    } else {
                        let mut compare_vals = Vec::new();
                        for (idx, name) in headers.iter().enumerate() {
                            let col1 = &filtered_runs[i].columns[idx];
                            let col2 = &filtered_runs[j].columns[idx];
                            let d = ks_statistic_cached(col1, col2, &mut buf_a, &mut buf_b);
                            let val = ks_limit - d;
                            if d > ks_limit {
                                ks_fail_count += 1;
                            }
                            compare_vals.push((name.clone(), val));
                        }
                        compare_vals
                    };
                    cont_compare.push(compare_vals);
                    if ks_fail_count > 0 {
                        fail_msgs.push(format!(
                            "{} parameters failed KS test between runs for Run_{}_Run_{}",
                            ks_fail_count,
                            i + 1,
                            j + 1
                        ));
                        failed_names.push(format!("Run_{}_Run_{}", i + 1, j + 1));
                        count_decision += 1;
                    }
                }
            }
        }
    }

    if count_decision == 0 {
        message_list.push_str(" ACHIEVED CONVERGENCE \n");
        message_list.push_str("  \n");
    } else {
        message_list.push_str(" FAILED CONVERGENCE \n");
        message_list.push_str("  \n");
        for msg in &fail_msgs {
            message_list.push_str(" ");
            message_list.push_str(msg);
            message_list.push_str(" \n");
        }
        message_list.push_str("   \n");
    }
    message_list.push_str(&format!(" BURN-IN SET AT {} \n", burnin));
    message_list.push_str("  \n");
    message_list.push_str("  \n");

    message_complete.push_str(&message_list);

    if !tree_exclude_high.is_empty() {
        if !tree_exclude_high[0].is_empty() {
            message_complete.push_str(" SPLITS EXCLUDED FROM CONVERGENCE ASSESSMENT \n");
        }
        for (i, items) in tree_exclude_high.iter().enumerate() {
            if !items.is_empty() {
                message_complete.push_str(&format!(" FREQUENCY HIGHER THAN 0.975 FOR RUN {} \n", i + 1));
                for item in items {
                    message_complete.push_str(&format!("      {} \n", item));
                }
                message_complete.push_str("  \n");
            }
            if i < tree_exclude_low.len() && !tree_exclude_low[i].is_empty() {
                message_complete.push_str(&format!(" FREQUENCY LOWER THAN 0.025 FOR RUN {} \n", i + 1));
                for item in &tree_exclude_low[i] {
                    message_complete.push_str(&format!("      {} \n", item));
                }
                message_complete.push_str("  \n");
            }
        }
    }

    if !cont_exclude.is_empty() && !cont_exclude[0].is_empty() {
        message_complete.push_str(" CONTINUOUS PARAMETERS WITH NO VARIANTION AND EXCLUDED FROM CONVERGENCE ASSESSMENT \n");
        for (i, items) in cont_exclude.iter().enumerate() {
            if !items.is_empty() {
                message_complete.push_str(&format!(" RUN {} \n", i + 1));
                for item in items {
                    message_complete.push_str(&format!("      {} \n", item));
                }
                message_complete.push_str("  \n");
            }
        }
    }

    if !tree_ess.is_empty() {
        message_list.push_str(" LOWEST SPLIT ESS \n");
        message_complete.push_str(" LOWEST SPLIT ESS \n");
        for (i, vals) in tree_ess.iter().enumerate() {
            if let Some((name, val)) = vals.iter().min_by(|a, b| a.1.partial_cmp(&b.1).unwrap()) {
                let line = format!("      RUN {} -> {} {} \n", i + 1, name, (val * 100.0).round() / 100.0);
                message_list.push_str(&line);
                message_complete.push_str(&line);
            }
        }
        message_list.push_str(" \n");
        message_complete.push_str("  \n");
    }

    if !cont_ess.is_empty() {
        message_list.push_str(" LOWEST CONTINUOUS PARAMETER ESS \n");
        message_complete.push_str(" LOWEST CONTINUOUS PARAMETER ESS \n");
        for (i, vals) in cont_ess.iter().enumerate() {
            if let Some((name, val)) = vals.iter().min_by(|a, b| a.1.partial_cmp(&b.1).unwrap()) {
                let line = format!("      RUN {} -> {} {} \n", i + 1, name, (val * 100.0).round() / 100.0);
                message_list.push_str(&line);
                message_complete.push_str(&line);
            }
        }
        message_list.push_str(" \n");
        message_complete.push_str("  \n");
    }

    if !cont_ess.is_empty() {
        message_list.push_str(" To check the calculated parameters for the continuous parameters type: \n");
        message_list.push_str("      Means: output$continuous_parameters$means \n");
        message_list.push_str("      ESS: output$continuous_parameters$ess \n");
        if !cont_compare.is_empty() {
            message_list.push_str("      KS score: output$continuous_parameters$compare_runs \n");
        }
        message_list.push_str(" \n");
    }

    if !tree_ess.is_empty() {
        message_list.push_str(" To check the calculated parameters for the splits type: \n");
        message_list.push_str("      Frequencies of splits: output$tree_parameters$frequencies \n");
        message_list.push_str("      ESS: output$tree_parameters$ess \n");
        if !tree_compare.is_empty() {
            message_list.push_str("      Difference in frequencies: output$tree_parameters$compare_runs \n");
        }
        message_list.push_str(" \n");
    }

    message_list.push_str(" To check the full summary message with splits and parameters excluded from the analysis type: \n");
    message_list.push_str("      output$message_complete \n");

    message_complete.push('\n');

    let converged = count_decision == 0;
    let failed_names_str = if failed_names.is_empty() {
        None
    } else {
        Some(failed_names.join(","))
    };

    Ok(ConvergenceResult {
        converged,
        burnin,
        message: message_list.clone(),
        message_complete,
        failed_names: failed_names_str,
        fail_msgs,
        tree_exclude_high,
        tree_exclude_low,
        tree_ess,
        tree_freqs,
        tree_compare_freqs,
        tree_compare,
        cont_means,
        cont_ess,
        cont_compare,
        cont_exclude,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn filter_table_excludes_headers() {
        let table = Table {
            headers: vec!["alpha".to_string(), "Iteration".to_string()],
            columns: vec![vec![1.0, 2.0], vec![10.0, 20.0]],
        };
        let re = Regex::new("Iteration").unwrap();
        let filtered = filter_table(&table, &re);
        assert_eq!(filtered.headers, vec!["alpha"]);
        assert_eq!(filtered.columns.len(), 1);
        assert_eq!(filtered.columns[0], vec![1.0, 2.0]);
    }

    #[test]
    fn remove_burnin_rejects_negative() {
        let mut runs = vec![Run {
            trees: vec!["(A,B);".to_string()],
            ptable: Table {
                headers: vec!["x".to_string()],
                columns: vec![vec![1.0]],
            },
        }];
        let err = remove_burnin(&mut runs, -0.1).unwrap_err();
        assert_eq!(err, "Invalid burnin value");
    }

    #[test]
    fn parse_table_missing_file() {
        let err = parse_table("does-not-exist.log", 0, b'\t').unwrap_err();
        assert!(!err.is_empty());
    }

    #[test]
    fn parse_revbayes_trees_missing_topology() {
        let dir = std::env::temp_dir();
        let path = dir.join("mcmc_missing_topology.trees");
        let data = "Iteration\tfoo\n1\t1.0\n";
        std::fs::write(&path, data).unwrap();
        let err = parse_revbayes_trees(path.to_str().unwrap()).unwrap_err();
        assert!(err.contains("Missing topology"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn parse_newick_trees_empty() {
        let dir = std::env::temp_dir();
        let path = dir.join("mcmc_empty.trees");
        std::fs::write(&path, "\n\n").unwrap();
        let err = parse_newick_trees(path.to_str().unwrap()).unwrap_err();
        assert!(err.contains("No trees found"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn parse_newick_trees_invalid() {
        let dir = std::env::temp_dir();
        let path = dir.join("mcmc_invalid_newick.trees");
        std::fs::write(&path, "not_a_tree\n").unwrap();
        let trees = parse_newick_trees(path.to_str().unwrap()).unwrap();
        assert_eq!(trees, vec!["not_a_tree".to_string()]);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn parse_nexus_trees_invalid() {
        let dir = std::env::temp_dir();
        let path = dir.join("mcmc_invalid_nexus.trees");
        let data = "tree t1 = not_a_tree\n";
        std::fs::write(&path, data).unwrap();
        let trees = parse_nexus_trees(path.to_str().unwrap()).unwrap();
        assert_eq!(trees, vec!["not_a_tree".to_string()]);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn parse_revbayes_trees_invalid_row() {
        let dir = std::env::temp_dir();
        let path = dir.join("mcmc_invalid_row.trees");
        let data = "Iteration\tTop\n1\n";
        std::fs::write(&path, data).unwrap();
        let err = parse_revbayes_trees(path.to_str().unwrap()).unwrap_err();
        assert!(err.contains("Invalid trees file row"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn check_convergence_rejects_nan_log_values() {
        let dir = std::env::temp_dir();
        let log_path = dir.join("mcmc_nan.log");
        let data = "Iteration\talpha\n1\tNA\n2\t1.0\n";
        std::fs::write(&log_path, data).unwrap();
        let control = Control {
            emit_logs: false,
            ..Control::default()
        };
        let err = check_convergence(
            &[log_path.to_string_lossy().to_string()],
            "revbayes",
            &control,
        )
        .unwrap_err();
        assert!(err.contains("non-numeric or missing values"));
        let _ = std::fs::remove_file(log_path);
    }
}
