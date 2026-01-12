#![allow(nonstandard_style)]

pub mod core;

use extendr_api::prelude::*;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

/// 计算 MCMC 链的 ESS（Geyer IPS）
#[extendr]
fn ess_tracer(x: Vec<f64>) -> f64 {
    let samples = x.len();
    if samples < 2 {
        return samples as f64;
    }

    let max_lag_limit = 2000usize;
    let max_lag = std::cmp::min(samples - 1, max_lag_limit);

    let mean = x.iter().sum::<f64>() / samples as f64;

    let mut gamma_stat = vec![0.0; max_lag];
    let mut var_stat = 0.0;

    for lag in 0..max_lag {
        let mut acc = 0.0;
        for j in 0..(samples - lag) {
            let del1 = x[j] - mean;
            let del2 = x[j + lag] - mean;
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

#[extendr]
fn effective_size(x: Vec<f64>) -> f64 {
    ess_tracer(x)
}

#[derive(Default)]
struct Node {
    children: Vec<usize>,
    label: Option<String>,
}

fn is_delim(b: u8) -> bool {
    matches!(b, b'(' | b')' | b',' | b':' | b';' | b'[' | b']')
}

fn skip_ws(bytes: &[u8], idx: &mut usize) {
    while *idx < bytes.len() && bytes[*idx].is_ascii_whitespace() {
        *idx += 1;
    }
}

fn skip_annotation(bytes: &[u8], idx: &mut usize) {
    if *idx < bytes.len() && bytes[*idx] == b'[' {
        *idx += 1;
        while *idx < bytes.len() && bytes[*idx] != b']' {
            *idx += 1;
        }
        if *idx < bytes.len() {
            *idx += 1;
        }
    }
}

fn parse_label(bytes: &[u8], idx: &mut usize) -> Option<String> {
    skip_ws(bytes, idx);
    if *idx >= bytes.len() || is_delim(bytes[*idx]) {
        return None;
    }
    let start = *idx;
    while *idx < bytes.len() && !is_delim(bytes[*idx]) {
        *idx += 1;
    }
    if *idx > start {
        Some(String::from_utf8_lossy(&bytes[start..*idx]).trim().to_string())
    } else {
        None
    }
}

fn skip_branch_length(bytes: &[u8], idx: &mut usize) {
    skip_ws(bytes, idx);
    if *idx < bytes.len() && bytes[*idx] == b':' {
        *idx += 1;
        while *idx < bytes.len() && !is_delim(bytes[*idx]) {
            *idx += 1;
        }
    }
}

fn parse_subtree(bytes: &[u8], idx: &mut usize, nodes: &mut Vec<Node>) -> Result<usize> {
    skip_ws(bytes, idx);
    if *idx >= bytes.len() {
        return Err(Error::Other("Unexpected end of Newick".into()));
    }

    if bytes[*idx] == b'(' {
        *idx += 1;
        let mut children = Vec::new();
        loop {
            let child = parse_subtree(bytes, idx, nodes)?;
            children.push(child);
            skip_ws(bytes, idx);
            if *idx >= bytes.len() {
                return Err(Error::Other("Unterminated Newick group".into()));
            }
            if bytes[*idx] == b',' {
                *idx += 1;
                continue;
            }
            if bytes[*idx] == b')' {
                *idx += 1;
                break;
            }
            return Err(Error::Other("Invalid Newick group separator".into()));
        }

        let label = parse_label(bytes, idx);
        skip_ws(bytes, idx);
        while *idx < bytes.len() && bytes[*idx] == b'[' {
            skip_annotation(bytes, idx);
            skip_ws(bytes, idx);
        }
        skip_branch_length(bytes, idx);
        while *idx < bytes.len() && bytes[*idx] == b'[' {
            skip_annotation(bytes, idx);
            skip_ws(bytes, idx);
        }

        let node = Node { children, label };
        nodes.push(node);
        Ok(nodes.len() - 1)
    } else {
        let label = parse_label(bytes, idx).ok_or_else(|| {
            Error::Other("Expected leaf label in Newick".into())
        })?;
        skip_ws(bytes, idx);
        while *idx < bytes.len() && bytes[*idx] == b'[' {
            skip_annotation(bytes, idx);
            skip_ws(bytes, idx);
        }
        skip_branch_length(bytes, idx);
        while *idx < bytes.len() && bytes[*idx] == b'[' {
            skip_annotation(bytes, idx);
            skip_ws(bytes, idx);
        }

        let node = Node {
            children: Vec::new(),
            label: Some(label),
        };
        nodes.push(node);
        Ok(nodes.len() - 1)
    }
}

pub(crate) fn parse_newick(tree: &str) -> Result<(usize, Vec<Node>)> {
    let bytes = tree.as_bytes();
    let mut idx = 0usize;
    let mut nodes = Vec::new();
    let root = parse_subtree(bytes, &mut idx, &mut nodes)?;
    skip_ws(bytes, &mut idx);
    if idx < bytes.len() && bytes[idx] == b';' {
        // End-of-tree delimiter; no further parsing needed here.
    }
    Ok((root, nodes))
}

pub(crate) fn collect_clades(root: usize, nodes: &[Node]) -> Result<Vec<String>> {
    let mut sets: Vec<Vec<String>> = vec![Vec::new(); nodes.len()];
    fn fill(
        idx: usize,
        nodes: &[Node],
        sets: &mut [Vec<String>],
    ) -> Vec<String> {
        if nodes[idx].children.is_empty() {
            let label = nodes[idx].label.clone().unwrap_or_default();
            sets[idx] = vec![label.clone()];
            return vec![label];
        }
        let mut all = Vec::new();
        for child in &nodes[idx].children {
            let child_set = fill(*child, nodes, sets);
            all.extend(child_set);
        }
        all.sort();
        all.dedup();
        sets[idx] = all.clone();
        all
    }

    let root_set = fill(root, nodes, &mut sets);
    let total = root_set.len();
    if total == 0 {
        return Ok(Vec::new());
    }
    let mut clades = Vec::new();
    fn walk(
        idx: usize,
        nodes: &[Node],
        sets: &[Vec<String>],
        total: usize,
        clades: &mut Vec<String>,
    ) {
        if nodes[idx].children.is_empty() {
            return;
        }
        let size = sets[idx].len();
        if size > 1 && size <= total {
            clades.push(sets[idx].join(" "));
        }
        for child in &nodes[idx].children {
            walk(*child, nodes, sets, total, clades);
        }
    }
    walk(root, nodes, &sets, total, &mut clades);
    Ok(clades)
}

#[extendr]
fn tree_clades(tree: String) -> Result<Vec<String>> {
    let (root, nodes) = parse_newick(&tree)?;
    collect_clades(root, &nodes)
}

#[extendr]
fn tree_tips(tree: String) -> Result<Vec<String>> {
    let (_root, nodes) = parse_newick(&tree)?;
    let mut tips = Vec::new();
    for node in nodes.iter() {
        if node.children.is_empty() {
            if let Some(label) = &node.label {
                tips.push(label.clone());
            }
        }
    }
    tips.sort();
    tips.dedup();
    Ok(tips)
}

#[extendr]
fn clade_counts(trees: Vec<String>) -> Result<List> {
    let stats = core::clade_stats_ids(&trees);
    let mut names = Vec::with_capacity(stats.order.len());
    let mut vals = Vec::with_capacity(stats.order.len());
    for &id in stats.order.iter() {
        names.push(stats.names[id].clone());
        vals.push(stats.counts[id] as f64);
    }

    Ok(list!(cladenames = names, counts = vals))
}

#[extendr]
fn clade_sets(trees: Vec<String>) -> Result<List> {
    let mut out: Vec<Vec<String>> = Vec::with_capacity(trees.len());
    for tree in trees.iter() {
        let (root, nodes) = parse_newick(tree)?;
        let clades = collect_clades(root, &nodes)?;
        out.push(clades);
    }
    let list_vals: Vec<Robj> = out
        .into_iter()
        .map(|v| r!(v))
        .collect();
    Ok(List::from_values(list_vals))
}

#[extendr]
fn clade_sets_and_counts(trees: Vec<String>) -> Result<List> {
    let mut order = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut counts: HashMap<String, usize> = HashMap::new();
    let mut sets: Vec<Vec<String>> = Vec::with_capacity(trees.len());

    for tree in trees.iter() {
        let (root, nodes) = parse_newick(tree)?;
        let clades = collect_clades(root, &nodes)?;
        let unique: HashSet<String> = clades.iter().cloned().collect();
        for clade in clades.iter() {
            if seen.insert(clade.clone()) {
                order.push(clade.clone());
            }
        }
        for clade in unique {
            *counts.entry(clade).or_insert(0) += 1;
        }
        sets.push(clades);
    }

    let vals: Vec<f64> = order
        .iter()
        .map(|name| counts.get(name).copied().unwrap_or(0) as f64)
        .collect();

    let list_vals: Vec<Robj> = sets.into_iter().map(|v| r!(v)).collect();
    let set_list = List::from_values(list_vals);

    Ok(list!(cladenames = order, counts = vals, sets = set_list))
}

fn table_from_df(obj: &Robj) -> Result<core::Table> {
    if obj.is_null() {
        return Ok(core::Table {
            headers: Vec::new(),
            columns: Vec::new(),
        });
    }
    let list = obj
        .as_list()
        .ok_or_else(|| Error::Other("Expected data.frame or list".into()))?;
    if list.len() == 0 {
        return Ok(core::Table {
            headers: Vec::new(),
            columns: Vec::new(),
        });
    }
    let names: Vec<String> = list
        .names()
        .map(|iter| iter.map(|s| s.to_string()).collect())
        .unwrap_or_else(|| vec![String::new(); list.len()]);
    let mut headers = Vec::with_capacity(list.len());
    let mut columns = Vec::with_capacity(list.len());
    for (idx, (_name, col)) in list.iter().enumerate() {
        let col_vals = if let Some(vals) = col.as_real_vector() {
            vals
        } else if let Some(vals) = col.as_integer_vector() {
            vals.into_iter().map(|v| v as f64).collect()
        } else {
            Vec::new()
        };
        headers.push(names.get(idx).cloned().unwrap_or_default());
        columns.push(col_vals);
    }
    Ok(core::Table { headers, columns })
}

fn df_from_table(table: &core::Table) -> Result<Robj> {
    let mut cols = Vec::with_capacity(table.headers.len());
    for col in &table.columns {
        cols.push(r!(col.clone()));
    }
    let list = if table.headers.is_empty() {
        List::from_values(cols)
    } else {
        List::from_names_and_values(&table.headers, &cols)
            .map_err(|e| Error::Other(e.into()))?
    };
    let mut df = Robj::from(list);
    let nrows = table.nrows();
    if nrows == 0 {
        let _ = df.set_attrib("row.names", r!(Vec::<i32>::new()));
    } else {
        let row_names: Vec<i32> = (1..=nrows as i32).collect();
        let _ = df.set_attrib("row.names", r!(row_names));
    }
    let _ = df.set_attrib("class", r!(vec!["data.frame"]));
    Ok(df)
}

fn thin_table(table: &core::Table, step: usize) -> core::Table {
    if step <= 1 {
        return table.clone();
    }
    let mut columns = Vec::with_capacity(table.columns.len());
    for col in &table.columns {
        let mut out = Vec::new();
        for (idx, val) in col.iter().enumerate() {
            if idx % step == 0 {
                out.push(*val);
            }
        }
        columns.push(out);
    }
    core::Table {
        headers: table.headers.clone(),
        columns,
    }
}

fn split_table_windows(table: &core::Table) -> (core::Table, core::Table) {
    let len = table.nrows();
    let (first_end, start2) = core::window_bounds(len);
    if len == 0 {
        return (
            core::Table {
                headers: table.headers.clone(),
                columns: Vec::new(),
            },
            core::Table {
                headers: table.headers.clone(),
                columns: Vec::new(),
            },
        );
    }
    (
        table.slice_rows(0, first_end),
        table.slice_rows(start2, len),
    )
}

fn run_from_list(list: &List) -> Result<(core::Run, Option<f64>)> {
    let mut trees = Vec::new();
    let mut ptable = core::Table {
        headers: Vec::new(),
        columns: Vec::new(),
    };
    let mut gens = None;
    for (name, val) in list.iter() {
        match name {
            "trees" => {
                if let Some(vals) = val.as_str_vector() {
                    trees = vals.iter().map(|s| s.to_string()).collect();
                } else if let Some(inner) = val.as_list() {
                    if inner.len() == 1 {
                        if let Ok(item) = inner.elt(0) {
                            if let Some(vals) = item.as_str_vector() {
                                trees = vals.iter().map(|s| s.to_string()).collect();
                            }
                        }
                    }
                }
            }
            "ptable" => {
                ptable = table_from_df(&val)?;
            }
            "gens.per.tree" => {
                gens = val.as_real().or_else(|| val.as_integer().map(|v| v as f64));
            }
            _ => {}
        }
    }
    Ok((core::Run { trees, ptable }, gens))
}

fn runs_from_list(list: &List) -> Result<(Vec<core::Run>, Vec<Option<f64>>, Vec<String>)> {
    let mut runs = Vec::new();
    let mut gens = Vec::new();
    for (_name, val) in list.iter() {
        let run_list = val
            .as_list()
            .ok_or_else(|| Error::Other("Expected list of runs".into()))?;
        let (run, gen) = run_from_list(&run_list)?;
        runs.push(run);
        gens.push(gen);
    }
    let names = list
        .names()
        .map(|iter| iter.map(|s| s.to_string()).collect())
        .unwrap_or_default();
    Ok((runs, gens, names))
}

fn run_to_list(run: &core::Run, gens: Option<f64>) -> Result<Robj> {
    let trees = if run.trees.is_empty() {
        r!(NULL)
    } else {
        r!(run.trees.clone())
    };
    let ptable = df_from_table(&run.ptable)?;
    let gens_val = gens.map_or(r!(NULL), |v| r!(v));
    let list = List::from_names_and_values(
        ["trees", "ptable", "gens.per.tree"],
        [trees, ptable, gens_val],
    )
    .map_err(|e| Error::Other(e.into()))?;
    let mut obj = Robj::from(list);
    let _ = obj.set_attrib("class", r!(vec!["rwty.chain"]));
    Ok(obj)
}

fn runs_to_list(
    runs: &[core::Run],
    gens: &[Option<f64>],
    names: &[String],
) -> Result<Robj> {
    let mut vals = Vec::with_capacity(runs.len());
    for (idx, run) in runs.iter().enumerate() {
        vals.push(run_to_list(run, gens.get(idx).copied().unwrap_or(None))?);
    }
    let mut list = List::from_values(vals);
    if !names.is_empty() && names.len() == runs.len() {
        let _ = list.set_names(names);
    }
    Ok(Robj::from(list))
}

fn list_get(list: &List, key: &str) -> Option<Robj> {
    for (name, val) in list.iter() {
        if name == key {
            return Some(val);
        }
    }
    None
}

fn as_string_vec(val: &Robj) -> Vec<String> {
    val.as_str_vector()
        .map(|vals| vals.iter().map(|s| s.to_string()).collect())
        .unwrap_or_default()
}

fn build_df_rows(
    rows: &[Vec<f64>],
    row_names: &[String],
    col_names: &[String],
) -> Result<Robj> {
    if rows.is_empty() || col_names.is_empty() {
        return Ok(r!(NULL));
    }
    let nrow = rows.len();
    let ncol = col_names.len();
    let mut cols = Vec::with_capacity(ncol);
    for c in 0..ncol {
        let mut col = Vec::with_capacity(nrow);
        for r in 0..nrow {
            col.push(*rows[r].get(c).unwrap_or(&f64::NAN));
        }
        cols.push(r!(col));
    }
    let list = List::from_names_and_values(col_names, &cols)
        .map_err(|e| Error::Other(e.into()))?;
    let mut df = Robj::from(list);
    let _ = df.set_attrib("row.names", r!(row_names.to_vec()));
    let _ = df.set_attrib("class", r!(vec!["data.frame"]));
    Ok(df)
}

fn named_numeric(vec: &Vec<(String, f64)>) -> Robj {
    let names: Vec<String> = vec.iter().map(|v| v.0.clone()).collect();
    let values: Vec<f64> = vec.iter().map(|v| v.1).collect();
    let mut robj = r!(values);
    robj.set_names(names).ok();
    robj
}

fn list_of_named(vecs: &Vec<Vec<(String, f64)>>, names: Option<&Vec<String>>) -> Result<List> {
    let vals: Vec<Robj> = vecs.iter().map(named_numeric).collect();
    if let Some(list_names) = names {
        List::from_names_and_values(list_names, &vals).map_err(|e| Error::Other(e.into()))
    } else {
        Ok(List::from_values(vals))
    }
}

fn ordered_union_names(vecs: &Vec<Vec<(String, f64)>>) -> Vec<String> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for vec in vecs {
        for (name, _) in vec {
            if seen.insert(name.clone()) {
                out.push(name.clone());
            }
        }
    }
    out
}

fn build_df_from_vecs(vecs: &Vec<Vec<(String, f64)>>, col_names: &[String]) -> Result<Robj> {
    if vecs.is_empty() || col_names.is_empty() {
        return Ok(r!(NULL));
    }
    let row_names = ordered_union_names(vecs);
    if row_names.is_empty() {
        return Ok(r!(NULL));
    }
    let mut cols = Vec::new();
    for vec in vecs {
        let map: HashMap<String, f64> = vec.iter().cloned().collect();
        let mut values = Vec::with_capacity(row_names.len());
        for name in &row_names {
            values.push(*map.get(name).unwrap_or(&f64::NAN));
        }
        cols.push(r!(values));
    }
    let list = List::from_names_and_values(col_names, &cols)
        .map_err(|e| Error::Other(e.into()))?;
    let mut df = Robj::from(list);
    let _ = df.set_attrib("row.names", r!(row_names));
    let _ = df.set_attrib("class", r!(vec!["data.frame"]));
    Ok(df)
}

fn build_matrix_from_vecs(vecs: &Vec<Vec<(String, f64)>>, col_names: &[String]) -> Result<Robj> {
    if vecs.is_empty() || col_names.is_empty() {
        return Ok(r!(NULL));
    }
    let row_names = ordered_union_names(vecs);
    if row_names.is_empty() {
        return Ok(r!(NULL));
    }
    let nrow = row_names.len();
    let ncol = col_names.len();
    let mut data = Vec::with_capacity(nrow * ncol);
    for vec in vecs {
        let map: HashMap<String, f64> = vec.iter().cloned().collect();
        for name in &row_names {
            data.push(*map.get(name).unwrap_or(&f64::NAN));
        }
    }
    let mut mat = r!(data);
    let _ = mat.set_attrib("dim", r!(vec![nrow as i32, ncol as i32]));
    let dimnames = List::from_values(vec![r!(row_names), r!(col_names.to_vec())]);
    let _ = mat.set_attrib("dimnames", dimnames);
    Ok(mat)
}

#[extendr]
fn align_named_vectors(vec_list: List) -> Result<Robj> {
    let mut all_names: Vec<String> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut rows: Vec<(Vec<String>, Vec<f64>)> = Vec::new();
    for (_name, val) in vec_list.iter() {
        let names: Vec<String> = val
            .names()
            .map(|iter| iter.map(|s| s.to_string()).collect())
            .unwrap_or_default();
        for name in &names {
            if seen.insert(name.clone()) {
                all_names.push(name.clone());
            }
        }
        let values = if let Some(vals) = val.as_real_vector() {
            vals
        } else if let Some(vals) = val.as_integer_vector() {
            vals.into_iter().map(|v| v as f64).collect()
        } else {
            Vec::new()
        };
        rows.push((names, values));
    }
    let nrow = rows.len();
    let ncol = all_names.len();
    let mut data = vec![f64::NAN; nrow * ncol];
    let mut name_idx = HashMap::new();
    for (idx, name) in all_names.iter().enumerate() {
        name_idx.insert(name.clone(), idx);
    }
    for (row_idx, (names, values)) in rows.iter().enumerate() {
        for (i, name) in names.iter().enumerate() {
            if let Some(&col_idx) = name_idx.get(name) {
                let val = values.get(i).copied().unwrap_or(f64::NAN);
                data[col_idx * nrow + row_idx] = val;
            }
        }
    }
    let mut mat = r!(data);
    let _ = mat.set_attrib("dim", r!(vec![nrow as i32, ncol as i32]));
    let dimnames = List::from_values(vec![r!(NULL), r!(all_names)]);
    let _ = mat.set_attrib("dimnames", dimnames);
    Ok(mat)
}

#[extendr]
fn read_trace(
    paths: Vec<String>,
    format: String,
    delim: String,
    burnin: f64,
    _check_names: bool,
    skip: i32,
) -> Result<Robj> {
    if paths.is_empty() {
        return Err(Error::Other("paths must be character strings".into()));
    }
    for p in &paths {
        if !Path::new(p).exists() {
            return Err(Error::Other(format!("Some files do not exist:\n{}", p)));
        }
    }
    if format != "simple" {
        return Err(Error::Other(
            "Complex trace type currently not supported".into(),
        ));
    }
    let delim_char = delim
        .as_bytes()
        .get(0)
        .ok_or_else(|| Error::Other("delim must be a single character string".into()))?;
    if burnin < 0.0 {
        return Err(Error::Other(
            "burnin must be a positive numeric value".into(),
        ));
    }
    let skip = if skip < 0 { 0 } else { skip as usize };

    let mut tables = Vec::new();
    let mut headers: Option<Vec<String>> = None;
    for path in &paths {
        let table = core::parse_table(path, skip, *delim_char)?;
        if let Some(h) = &headers {
            let set1: HashSet<&str> = h.iter().map(|s| s.as_str()).collect();
            let set2: HashSet<&str> = table.headers.iter().map(|s| s.as_str()).collect();
            if set1 != set2 {
                return Err(Error::Other("Not all headers of trace files match".into()));
            }
        } else {
            headers = Some(table.headers.clone());
        }
        tables.push(table);
    }
    let mut table = tables
        .into_iter()
        .next()
        .unwrap_or(core::Table {
            headers: Vec::new(),
            columns: Vec::new(),
        });
    let nrows = table.nrows();
    if burnin >= nrows as f64 && nrows > 0 {
        return Err(Error::Other(
            "Burnin larger than iterations in file".into(),
        ));
    }
    let discard = if burnin >= 1.0 {
        ((burnin / 100.0) * nrows as f64).ceil() as usize
    } else if burnin > 0.0 {
        (burnin * nrows as f64).ceil() as usize
    } else {
        0
    };
    if nrows > 0 && discard > 0 {
        table = table.slice_rows(discard.min(nrows), nrows);
    }
    let ptable = df_from_table(&table)?;
    let list = List::from_names_and_values(
        ["trees", "ptable", "gens.per.tree"],
        [r!(NULL), ptable, r!(NULL)],
    )
    .map_err(|e| Error::Other(e.into()))?;
    let mut obj = Robj::from(list);
    let _ = obj.set_attrib("class", r!(vec!["rwty.chain"]));
    Ok(obj)
}

#[extendr]
fn load_trees(
    file: String,
    format: String,
    gens_per_tree: f64,
    trim: i32,
    logfile: Nullable<String>,
    skip: i32,
) -> Result<Robj> {
    if format.to_lowercase() != "revbayes" {
        return Err(Error::Other(
            "Only revbayes tree parsing is supported without external dependencies."
                .into(),
        ));
    }
    let trim = if trim <= 1 { 1 } else { trim as usize };
    let skip = if skip < 0 { 0 } else { skip as usize };
    core::print_log(true, "Reading trees...");
    let (mut trees, mut rb_ptable) = core::parse_revbayes_trees(&file)?;
    if trim > 1 {
        trees = trees.into_iter().step_by(trim).collect();
        rb_ptable = thin_table(&rb_ptable, trim);
    }
    let gens = if gens_per_tree.is_nan() {
        if let Some(idx) = rb_ptable.column_index("Iteration") {
            if rb_ptable.columns[idx].len() >= 2 {
                (rb_ptable.columns[idx][1] - rb_ptable.columns[idx][0]).round()
            } else {
                1.0
            }
        } else {
            1.0
        }
    } else {
        gens_per_tree
    };
    core::print_log(true, &format!("{} generations per tree...", gens as i64));
    core::print_log(true, "Unrooting, this may take a while...");

    let mut log_path = logfile.into_option();
    if log_path.is_none() {
        if file.ends_with(".trees") {
            let candidate = file.trim_end_matches(".trees").to_string() + ".log";
            if Path::new(&candidate).exists() {
                log_path = Some(candidate);
            }
        }
    }
    let mut ptable = rb_ptable;
    if let Some(logfile) = log_path {
        if !Path::new(&logfile).exists() {
            return Err(Error::Other(format!("Logfile not found at {}", logfile)));
        }
        core::print_log(
            true,
            &format!(
                "Reading parameter values from {}",
                Path::new(&logfile).file_name().unwrap().to_string_lossy()
            ),
        );
        let mut log_table = core::parse_table(&logfile, skip, b'\t')?;
        if trim > 1 {
            log_table = thin_table(&log_table, trim);
        }
        ptable = core::merge_tables(&ptable, &log_table);
    }
    core::print_log(true, "rerooting trees...");
    if let Some(first) = trees.get(0) {
        let tips = tree_tips(first.clone())?;
        if let Some(outgroup) = tips.get(0) {
            core::print_log(true, &format!("Outgroup {}", outgroup));
        }
    }
    let list = run_to_list(&core::Run { trees, ptable }, Some(gens))?;
    Ok(list)
}

#[extendr]
fn load_multi(
    path: Nullable<String>,
    tree_files: Nullable<Vec<String>>,
    log_files: Nullable<Vec<String>>,
    format: String,
    labels: Nullable<Vec<String>>,
    skip: i32,
) -> Result<Robj> {
    let format = format.to_lowercase();
    let (tree_ext, log_ext) = match format.as_str() {
        "revbayes" => (".trees", ".log"),
        "mb" | "mrbayes" => (".t", ".p"),
        "beast" => (".trees", ".log"),
        "*beast" => (".trees", ".log"),
        "phylobayes" => (".treelist", ".trace"),
        "pyrate" => (".trees", ".log"),
        _ => return Err(Error::Other("Provide format!".into())),
    };
    let mut tfiles = Vec::new();
    let mut pfiles = Vec::new();
    if let Some(dir) = path.into_option() {
        for entry in fs::read_dir(&dir).map_err(|e| Error::Other(e.to_string()))? {
            let entry = entry.map_err(|e| Error::Other(e.to_string()))?;
            let path = entry.path();
            if let Some(name) = path.file_name().and_then(|s| s.to_str()) {
                if name.ends_with(tree_ext) {
                    tfiles.push(path.to_string_lossy().to_string());
                }
            }
        }
        for t in &tfiles {
            let base = t.trim_end_matches(tree_ext);
            pfiles.push(format!("{}{}", base, log_ext));
        }
    } else if let Some(files) = tree_files.into_option() {
        tfiles = files;
        if let Some(logs) = log_files.into_option() {
            pfiles = logs;
        } else {
            for t in &tfiles {
                let base = t.trim_end_matches(tree_ext);
                pfiles.push(format!("{}{}", base, log_ext));
            }
        }
    }
    if tfiles.is_empty() {
        return Err(Error::Other("Couldn't find any tree files".into()));
    }
    let mut runs = Vec::new();
    let mut gens = Vec::new();
    for (idx, tree) in tfiles.iter().enumerate() {
        core::print_log(
            true,
            &Path::new(tree).file_name().unwrap().to_string_lossy(),
        );
        let log_path = pfiles.get(idx).filter(|p| Path::new(p).exists());
        let list = load_trees(
            tree.clone(),
            format.clone(),
            f64::NAN,
            1,
            Nullable::from(log_path.cloned()),
            skip,
        )?;
        let run_list = list
            .as_list()
            .ok_or_else(|| Error::Other("Expected rwty.chain".into()))?;
        let (run, gen) = run_from_list(&run_list)?;
        runs.push(run);
        gens.push(gen);
    }
    let names = if let Some(vals) = labels.into_option() {
        vals
    } else {
        tfiles
            .iter()
            .map(|t| Path::new(t).file_name().unwrap().to_string_lossy().to_string())
            .collect()
    };
    runs_to_list(&runs, &gens, &names)
}

#[extendr]
fn load_files(
    path: Nullable<String>,
    list_files: Nullable<Vec<String>>,
    format: String,
    _tree_name: String,
) -> Result<Robj> {
    let format = format.to_lowercase();
    let (log_ext, tree_ext, delim, skip) = match format.as_str() {
        "revbayes" => (".log", ".trees", b'\t', 0usize),
        "mb" | "mrbayes" => (".p", ".t", b'\t', 1usize),
        "beast" => (".log", ".trees", b'\t', 2usize),
        "*beast" => (".log", ".trees", b',', 0usize),
        "phylobayes" => (".trace", ".treelist", b'\t', 0usize),
        "pyrate" => (".log", ".trees", b'\t', 0usize),
        _ => return Err(Error::Other("Provide format!".into())),
    };
    let files: Vec<String> = if let Some(dir) = path.into_option() {
        fs::read_dir(&dir)
            .map_err(|e| Error::Other(e.to_string()))?
            .filter_map(|entry| entry.ok())
            .map(|entry| entry.path().to_string_lossy().to_string())
            .collect()
    } else if let Some(list) = list_files.into_option() {
        for f in &list {
            if !Path::new(f).exists() {
                return Err(Error::Other(format!(
                    "Some files do not exist:\n{}",
                    f
                )));
            }
        }
        list
    } else {
        return Err(Error::Other("Some paths are not character strings".into()));
    };
    let mut log_files = Vec::new();
    let mut tree_files = Vec::new();
    for f in &files {
        if f.ends_with(log_ext) {
            log_files.push(f.clone());
        } else if f.ends_with(tree_ext) {
            tree_files.push(f.clone());
        }
    }
    if log_files.is_empty() && tree_files.is_empty() {
        return Err(Error::Other("No files to read".into()));
    }
    if !tree_files.is_empty() {
        let mut list = load_multi(
            Nullable::from(None::<String>),
            Nullable::from(Some(tree_files.clone())),
            if log_files.is_empty() {
                Nullable::from(None::<Vec<String>>)
            } else {
                Nullable::from(Some(log_files.clone()))
            },
            format.clone(),
            Nullable::from(None::<Vec<String>>),
            skip as i32,
        )?;
        if log_files.is_empty() {
            if let Some(runs) = list.as_list() {
                let (mut parsed, gens, names) = runs_from_list(&runs)?;
                for run in parsed.iter_mut() {
                    run.ptable = core::Table {
                        headers: Vec::new(),
                        columns: Vec::new(),
                    };
                }
                list = runs_to_list(&parsed, &gens, &names)?;
            }
        }
        return Ok(list);
    }
    let mut runs = Vec::new();
    let mut gens = Vec::new();
    for log in &log_files {
        let table = core::parse_table(log, skip, delim)?;
        runs.push(core::Run {
            trees: Vec::new(),
            ptable: table,
        });
        gens.push(None);
    }
    runs_to_list(&runs, &gens, &Vec::new())
}

#[extendr]
fn get_info(
    all_runs: List,
    run: i32,
    names_to_exclude: String,
    trees: bool,
    split_windows: bool,
) -> Result<Robj> {
    let (runs, _gens, _names) = runs_from_list(&all_runs)?;
    let run_idx = (run - 1).max(0) as usize;
    let run = runs
        .get(run_idx)
        .ok_or_else(|| Error::Other("Invalid run index".into()))?;
    if trees {
        if split_windows {
            let (w1, w2) = core::split_windows(&run.trees);
            return Ok(Robj::from(List::from_values(vec![r!(w1), r!(w2)])));
        }
        return Ok(r!(run.trees.clone()));
    }
    let names_re = Regex::new(&names_to_exclude).map_err(|e| Error::Other(e.to_string()))?;
    let filtered = core::filter_table(&run.ptable, &names_re);
    if split_windows {
        if filtered.headers.is_empty() {
            let empty = df_from_table(&core::Table {
                headers: Vec::new(),
                columns: Vec::new(),
            })?;
            return Ok(Robj::from(List::from_values(vec![empty.clone(), empty])));
        }
        let (w1, w2) = split_table_windows(&filtered);
        return Ok(Robj::from(List::from_values(vec![
            df_from_table(&w1)?,
            df_from_table(&w2)?,
        ])));
    }
    df_from_table(&filtered)
}

#[extendr]
fn remove_burnin(output: List, burnin: f64) -> Result<Robj> {
    let (mut runs, gens, names) = runs_from_list(&output)?;
    core::remove_burnin(&mut runs, burnin)?;
    runs_to_list(&runs, &gens, &names)
}

#[extendr]
fn clade_freq_named(x: Robj, start: i32, end: i32) -> Result<Robj> {
    let trees = if let Some(list) = x.as_list() {
        if let Some(t) = list_get(&list, "trees") {
            as_string_vec(&t)
        } else {
            as_string_vec(&x)
        }
    } else {
        as_string_vec(&x)
    };
    let start = (start - 1).max(0) as usize;
    let end = end.max(0) as usize;
    let slice_end = end.min(trees.len());
    let slice = if start < slice_end {
        &trees[start..slice_end]
    } else {
        &[][..]
    };
    let stats = core::clade_stats_ids(slice);
    let mut names = Vec::with_capacity(stats.order.len());
    let mut freqs = Vec::new();
    for &id in stats.order.iter() {
        names.push(stats.names[id].clone());
        let count = stats.counts[id] as f64;
        freqs.push(if slice.is_empty() {
            f64::NAN
        } else {
            count / slice.len() as f64
        });
    }
    let df = List::from_names_and_values(
        ["cladenames_post", "cladefreqs_post"],
        [r!(names), r!(freqs)],
    )
    .map_err(|e| Error::Other(e.into()))?;
    let mut obj = Robj::from(df);
    let _ = obj.set_attrib("class", r!(vec!["data.frame"]));
    Ok(obj)
}

#[extendr]
fn clade_freq_tree(x: String) -> Result<Robj> {
    let clades = tree_clades(x)?;
    let freqs = vec![1.0; clades.len()];
    let df = List::from_names_and_values(
        ["cladenames", "cladefreqs"],
        [r!(clades), r!(freqs)],
    )
    .map_err(|e| Error::Other(e.into()))?;
    let mut obj = Robj::from(df);
    let _ = obj.set_attrib("class", r!(vec!["data.frame"]));
    Ok(obj)
}

#[extendr]
fn clade_freq_trees(x: Robj, start: i32, end: i32) -> Result<Robj> {
    let trees = if let Some(list) = x.as_list() {
        if let Some(t) = list_get(&list, "trees") {
            as_string_vec(&t)
        } else {
            as_string_vec(&x)
        }
    } else {
        as_string_vec(&x)
    };
    let start = (start - 1).max(0) as usize;
    let end = end.max(0) as usize;
    let slice_end = end.min(trees.len());
    let slice = if start < slice_end {
        &trees[start..slice_end]
    } else {
        &[][..]
    };
    let stats = core::clade_stats_ids(slice);
    let mut kept_names = Vec::new();
    let mut kept_counts = Vec::new();
    for &id in stats.order.iter() {
        let name = stats.names[id].clone();
        let count = stats.counts[id] as f64;
        let freq = if slice.is_empty() {
            0.0
        } else {
            count / slice.len() as f64
        };
        if freq <= 0.975 && freq >= 0.025 {
            kept_names.push(name);
            kept_counts.push(count);
        }
    }
    let df = List::from_names_and_values(
        ["cladenames_post", "cladefreqs_post"],
        [r!(kept_names), r!(kept_counts)],
    )
    .map_err(|e| Error::Other(e.into()))?;
    let mut obj = Robj::from(df);
    let _ = obj.set_attrib("class", r!(vec!["data.frame"]));
    Ok(obj)
}

#[extendr]
fn check_clades_freq(runs: List, freq: f64) -> Result<List> {
    let (parsed, _gens, _names) = runs_from_list(&runs)?;
    let mut out = Vec::new();
    for run in parsed {
        let stats = core::clade_stats_ids(&run.trees);
        let mut names = Vec::with_capacity(stats.order.len());
        let mut freqs = Vec::with_capacity(stats.order.len());
        for &id in stats.order.iter() {
            names.push(stats.names[id].clone());
            let count = stats.counts[id] as f64;
            freqs.push(if run.trees.is_empty() {
                f64::NAN
            } else {
                count / run.trees.len() as f64
            });
        }
        if let Some(first) = run.trees.get(0) {
            let tips = tree_tips(first.clone())?;
            if !tips.is_empty() {
                names.push(tips.join(" "));
                freqs.push(1.0);
            }
        }
        let mut selected = Vec::new();
        for (idx, name) in names.iter().enumerate() {
            if freq > 0.5 {
                if freqs.get(idx).copied().unwrap_or(0.0) > freq {
                    selected.push(name.clone());
                }
            } else if freqs.get(idx).copied().unwrap_or(0.0) < freq {
                selected.push(name.clone());
            }
        }
        out.push(r!(selected));
    }
    Ok(List::from_values(out))
}

#[extendr]
fn ess_cont_param(
    runs: List,
    windows: bool,
    names_to_exclude: String,
    tracer: bool,
) -> Result<Robj> {
    let (parsed, _gens, _names) = runs_from_list(&runs)?;
    let names_re = Regex::new(&names_to_exclude).map_err(|e| Error::Other(e.to_string()))?;
    if !windows {
        let mut per_run: Vec<Vec<(String, f64)>> = Vec::new();
        for run in &parsed {
            let filtered = core::filter_table(&run.ptable, &names_re);
            let mut vals = Vec::new();
            for (idx, name) in filtered.headers.iter().enumerate() {
                let col = &filtered.columns[idx];
                let ess = if tracer {
                    core::ess_tracer(col)
                } else {
                    core::ess_tracer(col)
                };
                vals.push((name.clone(), ess));
            }
            per_run.push(vals);
        }
        let col_names: Vec<String> = (1..=per_run.len())
            .map(|i| format!("ESS_run_{}", i))
            .collect();
        return build_df_from_vecs(&per_run, &col_names);
    }
    let mut rows = Vec::new();
    let mut row_names = Vec::new();
    let mut col_names = Vec::new();
    for (idx, run) in parsed.iter().enumerate() {
        let filtered = core::filter_table(&run.ptable, &names_re);
        let (w1, w2) = split_table_windows(&filtered);
        if col_names.is_empty() {
            col_names = filtered.headers.clone();
        }
        let mut row1 = Vec::new();
        let mut row2 = Vec::new();
        for col in &w1.columns {
            row1.push(core::ess_tracer(col));
        }
        for col in &w2.columns {
            row2.push(core::ess_tracer(col));
        }
        rows.push(row1);
        rows.push(row2);
        row_names.push(format!("Run_{}_window_1", idx + 1));
        row_names.push(format!("Run_{}window_2", idx + 1));
    }
    build_df_rows(&rows, &row_names, &col_names)
}

#[extendr]
fn ess_split_freq(runs: List, tracer: bool) -> Result<List> {
    let (parsed, _gens, _names) = runs_from_list(&runs)?;
    let mut out = Vec::new();
    let mut out_names = Vec::new();
    for (idx, run) in parsed.iter().enumerate() {
        let stats = core::clade_stats_ids(&run.trees);
        let mut names = Vec::with_capacity(stats.order.len());
        let mut freqs = Vec::with_capacity(stats.order.len());
        for &id in stats.order.iter() {
            names.push(stats.names[id].clone());
            let count = stats.counts[id] as f64;
            freqs.push(if run.trees.is_empty() {
                f64::NAN
            } else {
                count / run.trees.len() as f64
            });
        }
        let mut vec_names = Vec::new();
        for (i, name) in names.iter().enumerate() {
            let f = freqs.get(i).copied().unwrap_or(0.0);
            if f <= 0.975 && f >= 0.025 {
                vec_names.push(name.clone());
            }
        }
        let mut ess_vals = Vec::new();
        if !vec_names.is_empty() {
            for clade in &vec_names {
                let Some(&id) = stats.index.get(clade) else { continue };
                let mut is_split = Vec::with_capacity(stats.sets.len());
                for set in &stats.sets {
                    is_split.push(if set.contains(&id) { 1.0 } else { 0.0 });
                }
                let ess = if tracer {
                    core::ess_tracer(&is_split)
                } else {
                    core::ess_tracer(&is_split)
                };
                ess_vals.push((clade.clone(), ess));
            }
        }
        out.push(named_numeric(&ess_vals));
        out_names.push(format!("Run_{}", idx + 1));
    }
    List::from_names_and_values(&out_names, &out).map_err(|e| Error::Other(e.into()))
}

#[extendr]
fn split_freq(runs: List, windows: bool) -> Result<Robj> {
    let (parsed, _gens, _names) = runs_from_list(&runs)?;
    if !windows {
        let mut all_df = Vec::new();
        for run in &parsed {
            let stats = core::clade_stats_ids(&run.trees);
            let mut names = Vec::with_capacity(stats.order.len());
            let mut freqs_by_id = vec![f64::NAN; stats.names.len()];
            for id in 0..stats.names.len() {
                freqs_by_id[id] = if run.trees.is_empty() {
                    f64::NAN
                } else {
                    stats.counts[id] as f64 / run.trees.len() as f64
                };
            }
            for &id in stats.order.iter() {
                names.push(stats.names[id].clone());
            }
            all_df.push((stats, names, freqs_by_id));
        }
        let mut check_null = true;
        for (_stats, names, _freqs) in &all_df {
            if names.is_empty() {
                check_null = false;
            }
        }
        let mut list_diff: Vec<Robj> = Vec::new();
        let mut list_freqs: Vec<Robj> = Vec::new();
        if check_null {
            if all_df.len() == 1 {
                let stats = &all_df[0].0;
                let names = &all_df[0].1;
                let freqs_by_id = &all_df[0].2;
                let mut freqs = Vec::with_capacity(names.len());
                for name in names {
                    if let Some(&id) = stats.index.get(name) {
                        freqs.push(freqs_by_id[id]);
                    }
                }
                let mut data = Vec::with_capacity(names.len() * 2);
                for (name, freq) in names.iter().zip(freqs.iter()) {
                    data.push(name.clone());
                    data.push(freq.to_string());
                }
                let mut mat = r!(data);
                let _ = mat.set_attrib("dim", r!(vec![2, names.len() as i32]));
                let dimnames = List::from_values(vec![
                    r!(vec!["listDiffSplits", "listFrequencies"]),
                    r!(NULL),
                ]);
                let _ = mat.set_attrib("dimnames", dimnames);
                return Ok(mat);
            }
            for r1 in 0..(all_df.len() - 1) {
                for r2 in (r1 + 1)..all_df.len() {
                    let (stats1, _names1, freqs1) = &all_df[r1];
                    let (stats2, _names2, freqs2) = &all_df[r2];
                    let mut vec_splits = Vec::new();
                    let mut vec_names = Vec::new();
                    let mut vec_freqs = Vec::new();
                    for &id1 in stats1.order.iter() {
                        let name = &stats1.names[id1];
                        let Some(&id2) = stats2.index.get(name) else { continue };
                        let diff = (freqs1[id1] - freqs2[id2]).abs();
                        vec_splits.push(diff);
                        vec_names.push(name.clone());
                        vec_freqs.push((freqs1[id1] + freqs2[id2]) / 2.0);
                    }
                    let mut splits = r!(vec_splits);
                    let _ = splits.set_attrib("names", r!(vec_names.clone()));
                    let mut freqs = r!(vec_freqs);
                    let _ = freqs.set_attrib("names", r!(vec_names));
                    list_diff.push(splits);
                    list_freqs.push(freqs);
                }
            }
        }
        let mut interleaved = Vec::with_capacity(list_diff.len() * 2);
        for i in 0..list_diff.len() {
            interleaved.push(list_diff[i].clone());
            interleaved.push(list_freqs[i].clone());
        }
        let mut list = List::from_values(interleaved);
        let _ = list.set_attrib("dim", r!(vec![2, list_diff.len() as i32]));
        let dimnames = List::from_values(vec![
            r!(vec!["listDiffSplits", "listFrequencies"]),
            r!(NULL),
        ]);
        let _ = list.set_attrib("dimnames", dimnames);
        return Ok(Robj::from(list));
    }
    let mut list_diff: Vec<Robj> = Vec::new();
    let mut list_freqs: Vec<Robj> = Vec::new();
    for run in &parsed {
        let (w1, w2) = core::split_windows(&run.trees);
        let stats1 = core::clade_stats_ids(&w1);
        let stats2 = core::clade_stats_ids(&w2);
        let mut freqs1 = vec![0.0; stats1.names.len()];
        let mut freqs2 = vec![0.0; stats2.names.len()];
        for id in 0..stats1.names.len() {
            freqs1[id] = stats1.counts[id] as f64 / w1.len().max(1) as f64;
        }
        for id in 0..stats2.names.len() {
            freqs2[id] = stats2.counts[id] as f64 / w2.len().max(1) as f64;
        }
        let mut vec_splits = Vec::new();
        let mut vec_freqs = Vec::new();
        let mut vec_names = Vec::new();
        for &id1 in stats1.order.iter() {
            let name = &stats1.names[id1];
            let Some(&id2) = stats2.index.get(name) else { continue };
            let f1 = freqs1[id1];
            let f2 = freqs2[id2];
            vec_splits.push((f1 - f2).abs());
            vec_freqs.push((f1 + f2) / 2.0);
            vec_names.push(name.clone());
        }
        let mut splits = r!(vec_splits);
        let _ = splits.set_attrib("names", r!(vec_names.clone()));
        let mut freqs = r!(vec_freqs);
        let _ = freqs.set_attrib("names", r!(vec_names));
        list_diff.push(splits);
        list_freqs.push(freqs);
    }
    let mut interleaved = Vec::with_capacity(list_diff.len() * 2);
    for i in 0..list_diff.len() {
        interleaved.push(list_diff[i].clone());
        interleaved.push(list_freqs[i].clone());
    }
    let mut list = List::from_values(interleaved);
    let _ = list.set_attrib("dim", r!(vec![2, list_diff.len() as i32]));
    let dimnames = List::from_values(vec![
        r!(vec!["listDiffSplits", "listFrequencies"]),
        r!(NULL),
    ]);
    let _ = list.set_attrib("dimnames", dimnames);
    Ok(Robj::from(list))
}

#[extendr]
fn mean_cont_param(runs: List, names_to_exclude: String) -> Result<Robj> {
    let (parsed, _gens, _names) = runs_from_list(&runs)?;
    let names_re = Regex::new(&names_to_exclude).map_err(|e| Error::Other(e.to_string()))?;
    let mut tables = Vec::new();
    for run in &parsed {
        tables.push(core::filter_table(&run.ptable, &names_re));
    }
    if tables.is_empty() {
        return Ok(r!(NULL));
    }
    let names_parameters = tables[0].headers.clone();
    let mut rows = Vec::new();
    let mut count = 0;
    if tables[0].headers.is_empty() {
        return Ok(r!(NULL));
    }
    for df1 in 0..(tables.len().saturating_sub(1)) {
        for _df2 in (df1 + 1)..tables.len() {
            let mut vec_means = Vec::new();
            for col in &tables[df1].columns {
                let mean = if col.is_empty() {
                    f64::NAN
                } else {
                    col.iter().sum::<f64>() / col.len() as f64
                };
                vec_means.push(mean);
            }
            rows.push(vec_means);
            count += 1;
        }
    }
    if count == 0 {
        return Ok(r!(NULL));
    }
    build_df_rows(&rows, &(1..=count).map(|i| i.to_string()).collect::<Vec<_>>(), &names_parameters)
}

#[extendr]
fn ks_test(runs: List, windows: bool, names_to_exclude: String) -> Result<Robj> {
    let (parsed, _gens, _names) = runs_from_list(&runs)?;
    let names_re = Regex::new(&names_to_exclude).map_err(|e| Error::Other(e.to_string()))?;
    if !windows {
        let mut tables = Vec::new();
        for run in &parsed {
            tables.push(core::filter_table(&run.ptable, &names_re));
        }
        if tables.is_empty() {
            return Ok(r!(NULL));
        }
        let names_parameters = tables[0].headers.clone();
        let mut rows = Vec::new();
        let mut count = 0;
        for df1 in 0..(tables.len().saturating_sub(1)) {
            for df2 in (df1 + 1)..tables.len() {
                let mut vec_ks = Vec::new();
                for i in 0..tables[df1].columns.len() {
                    let d = core::ks_statistic(&tables[df1].columns[i], &tables[df2].columns[i]);
                    vec_ks.push(d);
                }
                rows.push(vec_ks);
                count += 1;
            }
        }
        if count == 0 {
            return Ok(r!(NULL));
        }
        return build_df_rows(
            &rows,
            &(1..=count).map(|i| i.to_string()).collect::<Vec<_>>(),
            &names_parameters,
        );
    }
    let mut rows = Vec::new();
    let mut row_names = Vec::new();
    let mut col_names = Vec::new();
    for (idx, run) in parsed.iter().enumerate() {
        let filtered = core::filter_table(&run.ptable, &names_re);
        let (w1, w2) = split_table_windows(&filtered);
        if col_names.is_empty() {
            col_names = filtered.headers.clone();
        }
        let mut vec_ks = Vec::new();
        for i in 0..w1.columns.len() {
            let d = core::ks_statistic(&w1.columns[i], &w2.columns[i]);
            vec_ks.push(d);
        }
        rows.push(vec_ks);
        row_names.push(format!("Run_{}", idx + 1));
    }
    build_df_rows(&rows, &row_names, &col_names)
}

#[extendr]
fn ks_threshold(alpha: f64, ess: f64) -> f64 {
    core::ks_threshold(alpha, ess)
}

#[extendr]
fn min_ess(per: f64) -> f64 {
    core::min_ess(per)
}

#[extendr]
fn expected_diff_splits(ess: f64) -> Result<Robj> {
    let ess = ess.round().max(1.0) as usize;
    let (probs, thresh) = core::expected_diff_splits_cached(ess);
    let mut data = Vec::with_capacity(probs.len() * 2);
    for (p, t) in probs.iter().zip(thresh.iter()) {
        data.push(*p);
        data.push(*t);
    }
    let mut mat = r!(data);
    let _ = mat.set_attrib("dim", r!(vec![2, probs.len() as i32]));
    let dimnames = List::from_values(vec![r!(vec!["prob", "threshold"]), r!(NULL)]);
    let _ = mat.set_attrib("dimnames", dimnames);
    Ok(mat)
}

#[extendr]
fn se(x: Vec<f64>) -> f64 {
    let n = x.len();
    if n < 2 {
        return f64::NAN;
    }
    let mean = x.iter().sum::<f64>() / n as f64;
    let mut var = 0.0;
    for v in &x {
        var += (v - mean).powi(2);
    }
    let sd = (var / (n as f64 - 1.0)).sqrt();
    let ess = ess_tracer(x);
    if ess <= 0.0 {
        f64::NAN
    } else {
        sd / ess.sqrt()
    }
}

#[extendr]
fn quants(x: Vec<f64>) -> Robj {
    if x.is_empty() {
        return r!(vec![f64::NAN, f64::NAN]);
    }
    let mut sorted = x;
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len() as f64;
    let probs = [0.025, 0.975];
    let mut out = Vec::new();
    for p in probs {
        let h = (n - 1.0) * p + 1.0;
        let j = h.floor();
        let g = h - j;
        let j_idx = (j as usize).saturating_sub(1);
        let j1_idx = (j_idx + 1).min(sorted.len() - 1);
        let val = if j1_idx == j_idx {
            sorted[j_idx]
        } else {
            sorted[j_idx] + g * (sorted[j1_idx] - sorted[j_idx])
        };
        out.push(val);
    }
    let mut robj = r!(out);
    let _ = robj.set_attrib("names", r!(vec!["2.5%", "97.5%"]));
    robj
}

#[extendr]
fn check_convergence_r(
    list_files: Vec<String>,
    format: String,
    tracer: Option<bool>,
    burnin: Option<f64>,
    precision: Option<f64>,
    names_to_exclude: Option<String>,
    emit_logs: Option<bool>,
) -> List {
    let mut control = core::Control::default();
    if let Some(val) = tracer {
        control.tracer = val;
    }
    if let Some(val) = burnin {
        control.burnin = val;
    }
    if let Some(val) = precision {
        control.precision = val;
    }
    if let Some(val) = names_to_exclude {
        control.names_to_exclude = val;
    }
    if let Some(val) = emit_logs {
        control.emit_logs = val;
    }

    let result = match core::check_convergence(&list_files, &format, &control) {
        Ok(r) => r,
        Err(err) => throw_r_error(err),
    };

    let minimum_ess = (1.0 / (control.precision * 4.0)).powi(2);

    let run_count = std::cmp::max(result.tree_ess.len(), result.cont_ess.len());
    let run_names: Vec<String> = (1..=run_count).map(|i| format!("Run_{}", i)).collect();
    let mut compar_names = Vec::new();
    if run_count > 1 {
        for r1 in 1..run_count {
            for r2 in (r1 + 1)..=run_count {
                compar_names.push(format!("Run_{}_Run_{}", r1, r2));
            }
        }
    }

    let mut tree_fail_names = Vec::new();
    for (idx, ess_vec) in result.tree_ess.iter().enumerate() {
        let fails: Vec<(String, f64)> = ess_vec
            .iter()
            .filter(|(_, val)| *val < minimum_ess)
            .map(|(name, val)| (name.clone(), *val))
            .collect();
        if !fails.is_empty() {
            tree_fail_names.push((format!("ESS_of_Run_{}", idx + 1), named_numeric(&fails)));
        }
    }
    for (idx, compare_vec) in result.tree_compare.iter().enumerate() {
        let fails: Vec<String> = compare_vec
            .iter()
            .filter(|(_, val)| *val <= 0.0)
            .map(|(name, _)| name.clone())
            .collect();
        if !fails.is_empty() {
            let name = format!(
                "Between_{}",
                compar_names.get(idx).cloned().unwrap_or_default()
            );
            tree_fail_names.push((name, r!(fails)));
        }
    }

    let mut cont_fail_names = Vec::new();
    for (idx, ess_vec) in result.cont_ess.iter().enumerate() {
        let fails: Vec<(String, f64)> = ess_vec
            .iter()
            .filter(|(_, val)| *val < minimum_ess)
            .map(|(name, val)| (name.clone(), *val))
            .collect();
        if !fails.is_empty() {
            cont_fail_names.push((format!("ESS_run_{}", idx + 1), named_numeric(&fails)));
        }
    }
    for (idx, compare_vec) in result.cont_compare.iter().enumerate() {
        let fails: Vec<String> = compare_vec
            .iter()
            .filter(|(_, val)| *val < 0.0)
            .map(|(name, _)| name.clone())
            .collect();
        if !fails.is_empty() {
            let name = compar_names.get(idx).cloned().unwrap_or_default();
            cont_fail_names.push((name, r!(fails)));
        }
    }

    let mut failed_names_keys = Vec::new();
    let mut failed_names_vals = Vec::new();
    if !tree_fail_names.is_empty() {
        let names: Vec<String> = tree_fail_names.iter().map(|(n, _)| n.clone()).collect();
        let vals: Vec<Robj> = tree_fail_names.iter().map(|(_, v)| v.clone()).collect();
        let tree_list = List::from_names_and_values(&names, &vals)
            .map_err(|e| Error::Other(e.into()))
            .unwrap_or_else(|e| throw_r_error(e.to_string()));
        failed_names_keys.push("tree_parameters".to_string());
        failed_names_vals.push(Robj::from(tree_list));
    }
    if !cont_fail_names.is_empty() {
        let names: Vec<String> = cont_fail_names.iter().map(|(n, _)| n.clone()).collect();
        let vals: Vec<Robj> = cont_fail_names.iter().map(|(_, v)| v.clone()).collect();
        let cont_list = List::from_names_and_values(&names, &vals)
            .map_err(|e| Error::Other(e.into()))
            .unwrap_or_else(|e| throw_r_error(e.to_string()));
        failed_names_keys.push("continuous_parameters".to_string());
        failed_names_vals.push(Robj::from(cont_list));
    }
    let failed_names_list = if failed_names_keys.is_empty() {
        List::new(0)
    } else {
        List::from_names_and_values(&failed_names_keys, &failed_names_vals)
            .map_err(|e| Error::Other(e.into()))
            .unwrap_or_else(|e| throw_r_error(e.to_string()))
    };

    let tree_exclude_high: Vec<Robj> = result
        .tree_exclude_high
        .iter()
        .map(|v| r!(v.clone()))
        .collect();
    let tree_exclude_low: Vec<Robj> = result
        .tree_exclude_low
        .iter()
        .map(|v| r!(v.clone()))
        .collect();
    let cont_exclude: Vec<Robj> = result
        .cont_exclude
        .iter()
        .map(|v| r!(v.clone()))
        .collect();

    let tree_ess = list_of_named(&result.tree_ess, Some(&run_names))
        .unwrap_or_else(|e| throw_r_error(e.to_string()));
    let tree_freqs = list_of_named(&result.tree_freqs, Some(&run_names))
        .unwrap_or_else(|e| throw_r_error(e.to_string()));
    let tree_compare_freqs = if result.tree_compare_freqs.is_empty() {
        list_of_named(&result.tree_compare_freqs, None)
            .unwrap_or_else(|e| throw_r_error(e.to_string()))
    } else {
        list_of_named(&result.tree_compare_freqs, Some(&compar_names))
            .unwrap_or_else(|e| throw_r_error(e.to_string()))
    };
    let tree_compare = if result.tree_compare.is_empty() {
        list_of_named(&result.tree_compare, None)
            .unwrap_or_else(|e| throw_r_error(e.to_string()))
    } else {
        list_of_named(&result.tree_compare, Some(&compar_names))
            .unwrap_or_else(|e| throw_r_error(e.to_string()))
    };
    let cont_ess = list_of_named(&result.cont_ess, Some(&run_names))
        .unwrap_or_else(|e| throw_r_error(e.to_string()));
    let cont_compare = if result.cont_compare.is_empty() {
        list_of_named(&result.cont_compare, None)
            .unwrap_or_else(|e| throw_r_error(e.to_string()))
    } else {
        list_of_named(&result.cont_compare, Some(&compar_names))
            .unwrap_or_else(|e| throw_r_error(e.to_string()))
    };
    let tree_ess_df = build_df_from_vecs(&result.tree_ess, &run_names)
        .unwrap_or_else(|e| throw_r_error(e.to_string()));
    let cont_ess_df = build_df_from_vecs(&result.cont_ess, &run_names)
        .unwrap_or_else(|e| throw_r_error(e.to_string()));
    let cont_compare_df = if compar_names.is_empty() {
        Ok(r!(NULL))
    } else {
        build_df_from_vecs(&result.cont_compare, &compar_names)
    }
    .unwrap_or_else(|e| throw_r_error(e.to_string()));
    let freq_per_run = build_matrix_from_vecs(&result.tree_freqs, &run_names)
        .unwrap_or_else(|e| throw_r_error(e.to_string()));

    list!(
        converged = result.converged,
        burnin = result.burnin,
        message = result.message,
        message_complete = result.message_complete,
        failed = result.fail_msgs,
        failed_names = failed_names_list,
        compar_names = compar_names,
        tree_exclude_high = List::from_values(tree_exclude_high),
        tree_exclude_low = List::from_values(tree_exclude_low),
        tree_ess = tree_ess,
        tree_freqs = tree_freqs,
        tree_compare_freqs = tree_compare_freqs,
        tree_compare = tree_compare,
        tree_ess_df = tree_ess_df,
        cont_ess_df = cont_ess_df,
        cont_compare_df = cont_compare_df,
        freq_per_run = freq_per_run,
        cont_means = named_numeric(&result.cont_means),
        cont_ess = cont_ess,
        cont_compare = cont_compare,
        cont_exclude = List::from_values(cont_exclude)
    )
}

// Module registration for extendr.
extendr_module! {
    mod mcmcCheckConvergence;
    fn ess_tracer;
    fn effective_size;
    fn tree_clades;
    fn tree_tips;
    fn clade_counts;
    fn clade_sets;
    fn clade_sets_and_counts;
    fn align_named_vectors;
    fn read_trace;
    fn load_trees;
    fn load_multi;
    fn load_files;
    fn get_info;
    fn remove_burnin;
    fn clade_freq_named;
    fn clade_freq_tree;
    fn clade_freq_trees;
    fn check_clades_freq;
    fn ess_cont_param;
    fn ess_split_freq;
    fn split_freq;
    fn mean_cont_param;
    fn ks_test;
    fn ks_threshold;
    fn min_ess;
    fn expected_diff_splits;
    fn se;
    fn quants;
    fn check_convergence_r;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tree_tips_sorted() {
        let tips = tree_tips("(B,A);".to_string()).unwrap();
        assert_eq!(tips, vec!["A".to_string(), "B".to_string()]);
    }

    #[test]
    fn tree_clades_basic() {
        let clades = tree_clades("((A,B),C);".to_string()).unwrap();
        assert!(clades.iter().any(|c| c == "A B"));
    }
}
