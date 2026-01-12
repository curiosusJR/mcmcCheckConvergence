#!/usr/bin/env bash
set -euo pipefail

data_dir="${1:-output}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
bin_path="${repo_root}/src/rust/target/release/convergence_cli"

required_files=(
  "${data_dir}/posterior_run_1.log"
  "${data_dir}/posterior_run_2.log"
  "${data_dir}/posterior_run_1.trees"
  "${data_dir}/posterior_run_2.trees"
)

missing=0
for f in "${required_files[@]}"; do
  if [[ ! -f "${f}" ]]; then
    missing=1
  fi
done

if [[ "${missing}" -ne 0 ]]; then
  echo "Missing required files in ${data_dir}" >&2
  exit 1
fi

cargo build --manifest-path "${repo_root}/src/rust/Cargo.toml" --bin convergence_cli --release >/dev/null 2>&1

"${bin_path}" --path "${data_dir}" --format revbayes --message-only
