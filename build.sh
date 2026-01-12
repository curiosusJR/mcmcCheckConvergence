#!/bin/sh

set -eu

usage() {
  cat <<'EOF'
Usage:
  ./build.sh [--cli-only] [--r-only] [--check] [--help]

Options:
  --cli-only   Build only the Rust CLI (release).
  --r-only     Build only the R package tarball.
  --check      Run R CMD check --no-manual on the tarball.
  --help       Show this help.

Defaults:
  Builds both the Rust CLI (release) and the R package tarball.
EOF
}

mode="all"
run_check="false"

while [ $# -gt 0 ]; do
  case "$1" in
    --cli-only)
      mode="cli"
      ;;
    --r-only)
      mode="r"
      ;;
    --check)
      run_check="true"
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
  shift
done

build_cli() {
  cargo build --manifest-path src/rust/Cargo.toml --bin convergence_cli --release
}

build_r() {
  R CMD build .
  tarball="$(ls -t mcmcCheckConvergence_*.tar.gz | head -n 1)"
  if [ -z "${tarball}" ]; then
    echo "No tarball produced." >&2
    exit 1
  fi
  if [ "${run_check}" = "true" ]; then
    R CMD check --no-manual "${tarball}"
  fi
}

case "${mode}" in
  cli)
    build_cli
    ;;
  r)
    build_r
    ;;
  all)
    build_cli
    build_r
    ;;
esac
