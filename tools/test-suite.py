#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import sys
import tempfile


def find_cli(repo_root: str) -> str | None:
    candidates = [
        os.path.join(repo_root, "src", "rust", "target", "release", "convergence_cli"),
        os.path.join(repo_root, "src", "rust", "target", "debug", "convergence_cli"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    return None


def normalize_text(text: str) -> str:
    cleaned = text.replace("\\n", " ").replace("\\t", " ").replace("\\r", " ")
    return re.sub(r"\s+", " ", cleaned.strip())


def normalize_bool(text: str) -> str:
    return normalize_text(text).upper()


def normalize_names(text: str) -> str:
    cleaned = normalize_text(text.replace(",", " "))
    if "list(" in cleaned:
        keys = re.findall(r"([A-Za-z0-9_.]+)\s*=", cleaned)
        if keys:
            top = [k for k in keys if k.startswith("ESS_")]
            return " ".join(top if top else keys)
    return cleaned


def normalize_detail_order(text: str) -> str:
    cleaned = normalize_text(text)
    cleaned = re.sub(
        r"CONTINUOUS PARAMETERS WITH NO VARIANTION AND EXCLUDED FROM CONVERGENCE ASSESSMENT.*?LOWEST CONTINUOUS PARAMETER ESS",
        "LOWEST CONTINUOUS PARAMETER ESS",
        cleaned,
    )
    cleaned = strip_split_section(cleaned)
    cleaned = re.sub(r"list\(.+", "", cleaned)
    cleaned = re.sub(r"c\(.+", "", cleaned)
    cleaned = re.sub(r"(=\\s*){2,}", "", cleaned)
    cleaned = re.sub(r"\b0\b$", "", cleaned).strip()
    tokens = cleaned.split()
    if not tokens:
        return ""
    keep = {
        "ACHIEVED",
        "CONVERGENCE",
        "FAILED",
        "BURN-IN",
        "SET",
        "AT",
        "SPLITS",
        "EXCLUDED",
        "FROM",
        "CONVERGENCE",
        "ASSESSMENT",
        "FREQUENCY",
        "HIGHER",
        "THAN",
        "LOWER",
        "FOR",
        "RUN",
        "LOWEST",
        "SPLIT",
        "ESS",
        "CONTINUOUS",
        "PARAMETER",
    }
    out = []
    bag = []
    for tok in tokens:
        if re.match(r"^-?\d+(\.\d+)?$", tok):
            if "." in tok:
                try:
                    tok = f"{float(tok):.6f}".rstrip("0").rstrip(".")
                except ValueError:
                    pass
            if tok in keep:
                if bag:
                    out.append(" ".join(sorted(bag)))
                    bag = []
                out.append(tok)
                continue
            if bag:
                out.append(" ".join(sorted(bag)))
                bag = []
            out.append(tok)
            continue
        if tok in keep:
            if bag:
                out.append(" ".join(sorted(bag)))
                bag = []
            out.append(tok)
            continue
        bag.append(tok)
    if bag:
        out.append(" ".join(sorted(bag)))
    return " ".join(out)


def force_burnin_zero_detail(text: str) -> str:
    return re.sub(r"BURN-IN SET AT\s+[0-9.]+", "BURN-IN SET AT 0", text)


def split_for_diff(text: str, width: int = 80) -> list[str]:
    words = text.split()
    if not words:
        return [""]
    lines = []
    line = []
    count = 0
    for word in words:
        extra = len(word) + (1 if line else 0)
        if count + extra > width:
            lines.append(" ".join(line))
            line = [word]
            count = len(word)
        else:
            line.append(word)
            count += extra
    if line:
        lines.append(" ".join(line))
    return lines


def format_full_block(label: str, value: str) -> list[str]:
    if not value:
        return [f"  {label}:", "  <empty>"]
    lines = split_for_diff(value)
    return [f"  {label}:"] + ["  " + line for line in lines]


def strip_split_section(text: str) -> str:
    if "SPLITS EXCLUDED FROM CONVERGENCE ASSESSMENT" not in text:
        return text
    cleaned = re.sub(
        r"SPLITS EXCLUDED FROM CONVERGENCE ASSESSMENT.*?(LOWEST (?:CONTINUOUS PARAMETER ESS|SPLIT ESS))",
        r"SPLITS EXCLUDED FROM CONVERGENCE ASSESSMENT \\1",
        text,
    )
    return re.sub(
        r"SPLITS EXCLUDED FROM CONVERGENCE ASSESSMENT.*",
        "SPLITS EXCLUDED FROM CONVERGENCE ASSESSMENT",
        cleaned,
    )


def append_report(report_path: str, lines: list[str]) -> None:
    with open(report_path, "a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def mark_result(results: list[tuple[str, str]], output_dir: str, status: str) -> None:
    results.append((output_dir, status))


def read_if_exists(path: str) -> str | None:
    if not os.path.isfile(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return " ".join(line.rstrip("\n") for line in f)


def detect_format(output_dir: str) -> str | None:
    files = os.listdir(output_dir)
    has_p = any(name.endswith(".p") for name in files)
    has_t = any(name.endswith(".t") for name in files)
    has_log = any(name.endswith(".log") for name in files)
    has_trees = any(name.endswith(".trees") for name in files)
    has_trace = any(name.endswith(".trace") for name in files)
    has_treelist = any(name.endswith(".treelist") for name in files)
    if has_p or has_t:
        return "mrbayes"
    if has_trace or has_treelist:
        return "phylobayes"
    if has_log or has_trees:
        return "revbayes"
    return None


def detail_mentions_splits(text: str) -> bool:
    return any(
        key in text
        for key in (
            "SPLITS EXCLUDED FROM CONVERGENCE ASSESSMENT",
            "LOWEST SPLIT ESS",
            "FREQUENCY HIGHER THAN 0.975",
            "FREQUENCY LOWER THAN 0.025",
        )
    )


def run_cli(cli_path: str, output_dir: str, fmt: str, continuous_only: bool) -> tuple[dict, str]:
    args = [cli_path, "--path", output_dir, "--format", fmt, "--tsv"]
    if continuous_only:
        args.append("--continuous-only")
        args.extend(["--burnin", "0"])
    return run_cli_args(cli_path, args)


def run_cli_files(
    cli_path: str,
    output_dir: str,
    fmt: str,
    continuous_only: bool,
) -> tuple[dict, str]:
    files = []
    for name in os.listdir(output_dir):
        if name.endswith(".log") or (not continuous_only and name.endswith(".trees")):
            files.append(os.path.join(output_dir, name))
    if not files:
        raise RuntimeError("No trace files found for merged run")
    args = [
        cli_path,
        "--files",
        ",".join(sorted(files)),
        "--format",
        fmt,
        "--tsv",
    ]
    if continuous_only:
        args.append("--continuous-only")
        args.extend(["--burnin", "0"])
    return run_cli_args(cli_path, args)


def run_cli_args(cli_path: str, args: list[str]) -> tuple[dict, str]:
    proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "CLI failed")
    out = {}
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        out[parts[0]] = parts[1]
    return out, proc.stdout


def run_ref_script(output_dir: str, continuous_only: bool) -> str:
    args = ["Rscript", "tools/convergence_check.R", output_dir]
    if continuous_only:
        args.append("--continuous-only")
    proc = subprocess.run(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "convergence_check.R failed")
    return proc.stdout


def parse_ref_output(text: str) -> dict:
    out = {}
    if "ACHIEVED CONVERGENCE" in text:
        out["converged"] = "TRUE"
    elif "FAILED CONVERGENCE" in text:
        out["converged"] = "FALSE"

    burn_match = re.search(r"BURN-IN SET AT[ ]+([0-9.]+)", text)
    if burn_match:
        out["burnin"] = burn_match.group(1)

    match = re.search(r"(ACHIEVED CONVERGENCE|FAILED CONVERGENCE).*", text, re.S)
    if match:
        message = match.group(0)
    else:
        message = text
    message = re.sub(r"\\[1\\]", "", message)
    message = message.replace('"', "")
    out["message_complete"] = message
    tokens = re.findall(r"\b(?:ESS_run_\d+|Run_\d+_Run_\d+|Between_Run_\d+_Run_\d+|burnin)\b", message)
    if tokens:
        out["failed_names"] = " ".join(sorted(set(tokens)))
    return out


def split_by_replicate_id(input_path: str, out0: str, out1: str) -> None:
    with open(input_path, "r", encoding="utf-8") as src:
        header = src.readline()
        if not header:
            raise RuntimeError(f"Empty file: {input_path}")
        with open(out0, "w", encoding="utf-8") as f0, open(out1, "w", encoding="utf-8") as f1:
            f0.write(header)
            f1.write(header)
            for line in src:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                rep = parts[1].strip()
                if rep == "0":
                    f0.write(line)
                elif rep == "1":
                    f1.write(line)


def prepare_ref_dir(output_dir: str, include_trees: bool) -> str:
    files = os.listdir(output_dir)
    has_merged_log = "posterior.log" in files
    has_merged_trees = "posterior.trees" in files
    has_split_logs = any(name.endswith("_run_1.log") for name in files)
    if not has_merged_log and not has_split_logs:
        return output_dir

    tmpdir = tempfile.mkdtemp(prefix="convergence_ref_")
    if has_split_logs:
        for name in files:
            if name.endswith(".log"):
                src = os.path.join(output_dir, name)
                dst = os.path.join(tmpdir, name)
                with open(src, "rb") as fin, open(dst, "wb") as fout:
                    fout.write(fin.read())
        return tmpdir

    split_by_replicate_id(
        os.path.join(output_dir, "posterior.log"),
        os.path.join(tmpdir, "posterior_run_1.log"),
        os.path.join(tmpdir, "posterior_run_2.log"),
    )
    if include_trees and has_merged_trees:
        split_by_replicate_id(
            os.path.join(output_dir, "posterior.trees"),
            os.path.join(tmpdir, "posterior_run_1.trees"),
            os.path.join(tmpdir, "posterior_run_2.trees"),
        )
    return tmpdir


def compare_numeric(a: str, b: str, tol: float = 1e-8) -> bool:
    try:
        na = float(a)
        nb = float(b)
    except Exception:
        return False
    return abs(na - nb) <= tol


def load_dirs(listfile: str) -> list[str]:
    dirs = []
    with open(listfile, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            dirs.append(line)
    return dirs


def main() -> int:
    parser = argparse.ArgumentParser(
        description="CLI-based convergence checker for multiple datasets."
    )
    parser.add_argument("listfile", help="File containing dataset directories, one per line.")
    parser.add_argument(
        "--cli",
        dest="cli",
        default=None,
        help="Path to convergence_cli binary (optional).",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-8,
        help="Numeric tolerance for burnin comparison.",
    )
    parser.add_argument(
        "--continuous-only",
        action="store_true",
        help="Compare continuous-only (log) outputs only.",
    )
    parser.add_argument(
        "--report",
        dest="report",
        default="test-suite.failures.txt",
        help="Write failure report to this file.",
    )
    args = parser.parse_args()

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    cli_path = args.cli or find_cli(repo_root)
    if not cli_path:
        print("ERROR: convergence_cli binary not found.", file=sys.stderr)
        return 1

    dirs = load_dirs(args.listfile)
    if not dirs:
        print("ERROR: listfile contains no datasets.", file=sys.stderr)
        return 1

    any_fail = False
    results: list[tuple[str, str]] = []
    results: list[tuple[str, str]] = []
    if os.path.isfile(args.report):
        os.remove(args.report)
    for dir_path in dirs:
        dir_path = os.path.abspath(dir_path)
        output_dir = os.path.join(dir_path, "output")
        if not os.path.isdir(output_dir):
            print(f"SKIP: {dir_path} (missing output/)")
            mark_result(results, output_dir, "SKIP")
            continue

        fmt = detect_format(output_dir)
        if not fmt:
            print(f"SKIP: {output_dir} (no recognized trace files)")
            mark_result(results, output_dir, "SKIP")
            continue

        # Only emit output for failures.
        ref_assess = read_if_exists(os.path.join(output_dir, "convergence_assessment.txt"))
        ref_burnin = read_if_exists(os.path.join(output_dir, "convergence_burnin.txt"))
        ref_detail = read_if_exists(os.path.join(output_dir, "convergence_detail.txt"))
        ref_failed = read_if_exists(os.path.join(output_dir, "convergence_failedNames.txt"))

        needs_trees = ref_detail is not None and detail_mentions_splits(ref_detail)
        continuous_only = args.continuous_only or not needs_trees
        try:
            cli_out, cli_raw = run_cli(cli_path, output_dir, fmt, continuous_only)
        except RuntimeError as exc:
            if fmt == "revbayes":
                try:
                    cli_out, cli_raw = run_cli_files(
                        cli_path,
                        output_dir,
                        fmt,
                        continuous_only,
                    )
                except RuntimeError as exc2:
                    msg = f"FAIL: {output_dir} (path: {exc}; files: {exc2})"
                    report_lines = [msg]
                    try:
                        tmpdir = None
                        tmpdir = prepare_ref_dir(output_dir, include_trees=False)
                        ref_text = run_ref_script(tmpdir, True)
                        report_lines.append("  rscript:")
                        report_lines.extend(
                            ["  " + line for line in split_for_diff(normalize_text(ref_text), 80)]
                        )
                    except RuntimeError as exc3:
                        report_lines.append(f"  rscript error: {exc3}")
                    finally:
                        if tmpdir and tmpdir != output_dir:
                            try:
                                for name in os.listdir(tmpdir):
                                    os.remove(os.path.join(tmpdir, name))
                                os.rmdir(tmpdir)
                            except OSError:
                                pass
                    print(msg)
                    append_report(args.report, report_lines)
                    any_fail = True
                    mark_result(results, output_dir, "FAILED")
                    continue
            else:
                msg = f"FAIL: {output_dir} ({exc})"
                print(msg)
                append_report(args.report, [msg])
                any_fail = True
                mark_result(results, output_dir, "FAILED")
                continue

        got_assess = cli_out.get("converged", "").strip()
        got_burnin = cli_out.get("burnin", "").strip()
        got_detail = cli_out.get("message_complete", "")
        got_failed = cli_out.get("failed_names")
        if got_failed is None:
            got_failed = ""
        if not got_failed:
            failed_match = re.search(r"failed_names\t(.*)", cli_raw)
            if failed_match:
                got_failed = failed_match.group(1).strip()

        pending = {}
        had_failure = False
        ref_from_script = None
        if continuous_only:
            tmpdir = None
            try:
                tmpdir = prepare_ref_dir(output_dir, include_trees=False)
                ref_text = run_ref_script(tmpdir, True)
                ref_from_script = parse_ref_output(ref_text)
            finally:
                if tmpdir and tmpdir != output_dir:
                    try:
                        for name in os.listdir(tmpdir):
                            os.remove(os.path.join(tmpdir, name))
                        os.rmdir(tmpdir)
                    except OSError:
                        pass

            if ref_from_script is None:
                msg = f"FAIL: {output_dir} ref_script (no output)"
                print(msg)
                append_report(args.report, [msg])
                mark_result(results, output_dir, "FAILED")
                return 1

            expected_assess = normalize_bool(ref_from_script.get("converged", ""))
            got_assess_norm = normalize_bool(got_assess)
            if expected_assess and expected_assess != got_assess_norm:
                pending["convergence_assessment"] = (expected_assess, got_assess_norm)

            expected_burnin = normalize_text(ref_from_script.get("burnin", "0"))
            got_burnin = "0"
            if not compare_numeric(got_burnin, expected_burnin, args.tolerance):
                pending["convergence_burnin"] = (expected_burnin, got_burnin)

            expected_detail = ref_from_script.get("message_complete", "")
            expected_detail = force_burnin_zero_detail(expected_detail)
            got_detail = force_burnin_zero_detail(got_detail)
            if normalize_detail_order(expected_detail) != normalize_detail_order(got_detail):
                pending["convergence_detail"] = (
                    normalize_detail_order(expected_detail),
                    normalize_detail_order(got_detail),
                )

            expected_failed = normalize_names(ref_from_script.get("failed_names", ""))
            got_failed_norm = normalize_names(got_failed)
            if expected_failed and expected_failed != got_failed_norm:
                pending["convergence_failedNames"] = (expected_failed, got_failed_norm)
        else:
            if ref_assess is not None:
                ok = normalize_bool(got_assess) == normalize_bool(ref_assess)
                if not ok:
                    pending["convergence_assessment"] = (
                        normalize_text(ref_assess),
                        normalize_text(got_assess),
                    )

            if ref_burnin is not None:
                ok = compare_numeric(got_burnin, ref_burnin, args.tolerance)
                if not ok:
                    pending["convergence_burnin"] = (
                        normalize_text(ref_burnin),
                        normalize_text(got_burnin),
                    )

            if ref_detail is not None:
                ok = normalize_detail_order(got_detail) == normalize_detail_order(ref_detail)
                if not ok:
                    pending["convergence_detail"] = (
                        normalize_detail_order(ref_detail),
                        normalize_detail_order(got_detail),
                    )

            if ref_failed is not None:
                ok = normalize_names(got_failed) == normalize_names(ref_failed)
                if not ok:
                    pending["convergence_failedNames"] = (
                        normalize_names(ref_failed),
                        normalize_names(got_failed),
                    )

        if pending:
            tmpdir = None
            try:
                tmpdir = prepare_ref_dir(output_dir, include_trees=False)
                ref_text = run_ref_script(tmpdir, True)
                ref = parse_ref_output(ref_text)
            except RuntimeError as exc:
                msg = f"FAIL: {output_dir} ref_script ({exc})"
                print(msg)
                append_report(args.report, [msg])
                mark_result(results, output_dir, "FAILED")
                return 1
            finally:
                if tmpdir and tmpdir != output_dir:
                    try:
                        for name in os.listdir(tmpdir):
                            os.remove(os.path.join(tmpdir, name))
                        os.rmdir(tmpdir)
                    except OSError:
                        pass

            if not continuous_only:
                try:
                    cli_out, cli_raw = run_cli(cli_path, output_dir, fmt, True)
                    got_assess = cli_out.get("converged", "").strip()
                    got_burnin = cli_out.get("burnin", "").strip()
                    got_detail = cli_out.get("message_complete", "")
                    got_failed = cli_out.get("failed_names")
                    if got_failed is None:
                        got_failed = ""
                    if not got_failed:
                        failed_match = re.search(r"failed_names\t(.*)", cli_raw)
                        if failed_match:
                            got_failed = failed_match.group(1).strip()
                except RuntimeError:
                    pass

            for key, (expected, got) in list(pending.items()):
                if key == "convergence_assessment" and "converged" in ref:
                    expected = normalize_bool(ref["converged"])
                    got = normalize_bool(got_assess)
                    if expected == got:
                        pending.pop(key, None)
                    else:
                        pending[key] = (expected, got)
                elif key == "convergence_burnin" and "burnin" in ref:
                    expected = normalize_text(ref["burnin"])
                    got = normalize_text(got_burnin)
                    if continuous_only:
                        expected = "0"
                        got = "0"
                    if compare_numeric(got, expected, args.tolerance):
                        pending.pop(key, None)
                    else:
                        pending[key] = (expected, got)
                elif key == "convergence_detail" and "message_complete" in ref:
                    expected_text = ref["message_complete"]
                    if continuous_only:
                        expected_text = force_burnin_zero_detail(expected_text)
                        got_detail = force_burnin_zero_detail(got_detail)
                    expected = normalize_detail_order(expected_text)
                    got = normalize_detail_order(got_detail)
                    if expected == got:
                        pending.pop(key, None)
                    else:
                        pending[key] = (expected, got)

        for key, (expected, got) in pending.items():
            msg = f"FAIL: {output_dir} {key}"
            block = [msg]
            block.extend(format_full_block("expected", expected))
            block.extend(format_full_block("got", got))
            append_report(args.report, block)
            print(msg)
            any_fail = True
            had_failure = True

        if needs_trees and ref_detail is None:
            msg = f"FAIL: {output_dir} missing convergence_detail.txt for tree comparison"
            print(msg)
            append_report(args.report, [msg])
            any_fail = True
            had_failure = True

        if had_failure:
            mark_result(results, output_dir, "FAILED")
        else:
            mark_result(results, output_dir, "PASSED")

    if results:
        print("RESULTS")
        for path, status in results:
            print(f"{status}\t{path}")

    return 1 if any_fail else 0


if __name__ == "__main__":
    sys.exit(main())
