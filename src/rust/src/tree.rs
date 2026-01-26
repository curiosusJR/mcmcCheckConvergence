#[derive(Default)]
pub(crate) struct Node {
    pub(crate) children: Vec<usize>,
    pub(crate) label: Option<String>,
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
        Some(
            String::from_utf8_lossy(&bytes[start..*idx])
                .trim()
                .to_string(),
        )
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

fn parse_subtree(bytes: &[u8], idx: &mut usize, nodes: &mut Vec<Node>) -> Result<usize, String> {
    skip_ws(bytes, idx);
    if *idx >= bytes.len() {
        return Err("Unexpected end of Newick".to_string());
    }

    if bytes[*idx] == b'(' {
        *idx += 1;
        let mut children = Vec::new();
        loop {
            let child = parse_subtree(bytes, idx, nodes)?;
            children.push(child);
            skip_ws(bytes, idx);
            if *idx >= bytes.len() {
                return Err("Unterminated Newick group".to_string());
            }
            if bytes[*idx] == b',' {
                *idx += 1;
                continue;
            }
            if bytes[*idx] == b')' {
                *idx += 1;
                break;
            }
            return Err("Invalid Newick group separator".to_string());
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
        let label = parse_label(bytes, idx)
            .ok_or_else(|| "Expected leaf label in Newick".to_string())?;
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

pub(crate) fn parse_newick(tree: &str) -> Result<(usize, Vec<Node>), String> {
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

pub(crate) fn collect_clades(root: usize, nodes: &[Node]) -> Result<Vec<String>, String> {
    let mut sets: Vec<Vec<String>> = vec![Vec::new(); nodes.len()];
    fn fill(idx: usize, nodes: &[Node], sets: &mut [Vec<String>]) -> Vec<String> {
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

pub(crate) fn tree_tips(tree: &str) -> Result<Vec<String>, String> {
    let (_root, nodes) = parse_newick(tree)?;
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
