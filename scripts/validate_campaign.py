#!/usr/bin/env python3
"""
validate_campaign.py DEST [--cross-tol FRAC] [--min-root-mb MB]

Post-run check of a finished run_process.py campaign. Prints one line per check, prefixed
"OK", "WARN" or "ERROR", and a final verdict. Exit status 1 if any ERROR.

Checks
  1. Every resbos log (resbos/resbos_output_*.log, latest job per piece) has a finite
     TOTAL CROSS SECTION, reached "finished with exit code 0", and contains no NaN.
  2. Each resbos piece produced a .root in DEST with a plausible size (stubs are ~60 KB,
     real outputs are ~1 GB; below --min-root-mb is an ERROR).
  3. The merged grids (*.out in legacy/, w_pert/, w_asym/, get_yk_new/, resbos/Resbos_grids/)
     contain no NaN/Inf tokens.
  4. If resNLO, nloasy and nlodsi are all present for a process: resNLO should agree with
     nloasy + nlodsi (WARN if the relative difference exceeds --cross-tol).
"""
import argparse
import glob
import math
import os
import re
import sys

RESBOS_LOG = re.compile(r"^resbos_output_(.+)_(resNLO|nloasy|nlodsi)_(\d+)\.log$")
XSEC = re.compile(r"TOTAL CROSS SECTION\s*=\s*(\S+)\s+PB\s*\+/-\s*(\S+)")
NAN_TOKEN = re.compile(r"(?<![A-Za-z0-9_])(NaN|Inf|Infinity)(?![A-Za-z0-9_])", re.IGNORECASE)

errors, warns = [], []


def ok(msg):
    print(f"OK    {msg}")


def warn(msg):
    warns.append(msg)
    print(f"WARN  {msg}")


def error(msg):
    errors.append(msg)
    print(f"ERROR {msg}")


def latest_logs(dest):
    """{(base, piece): (logpath, jobid)} keeping only the highest job id per (base, piece)."""
    found = {}
    for path in glob.glob(os.path.join(dest, "resbos", "resbos_output_*.log")):
        m = RESBOS_LOG.match(os.path.basename(path))
        if not m:
            continue
        key = (m.group(1), m.group(2))
        jid = int(m.group(3))
        if key not in found or jid > found[key][1]:
            found[key] = (path, jid)
    return found


def parse_log(path):
    with open(path, errors="replace") as f:
        text = f.read()
    m = None
    for m in XSEC.finditer(text):
        pass  # keep the last TOTAL CROSS SECTION line
    xs = (m.group(1), m.group(2)) if m else None
    nan_count = len(re.findall(r"\bNaN\b", text))
    finished = "finished with exit code 0" in text
    has_events = "EVENTS PASSING CUTS" in text
    return xs, nan_count, finished, has_events


def check_resbos(dest, logs):
    results = {}
    for (base, piece), (path, jid) in sorted(logs.items()):
        name = f"{base} {piece} (job {jid})"
        xs, nan_count, finished, has_events = parse_log(path)
        if xs is None:
            error(f"{name}: no 'TOTAL CROSS SECTION' line -- run did not finish")
            continue
        val, err = xs
        try:
            v, e = float(val), float(err)
        except ValueError:
            v, e = float("nan"), float("nan")
        if not (math.isfinite(v) and math.isfinite(e)):
            error(f"{name}: cross section is not finite ({val} +/- {err})")
            continue
        if nan_count:
            error(f"{name}: {nan_count} NaN tokens in the resbos log, although the final value is finite")
        if not finished:
            error(f"{name}: missing 'finished with exit code 0'")
        if not has_events:
            error(f"{name}: no 'EVENTS PASSING CUTS' line")
        if finished and has_events and not nan_count:
            ok(f"{name}: sigma = {v:.6g} +/- {e:.4g} pb, finite, finished")
        results[(base, piece)] = (v, e)
    return results


def check_roots(dest, logs, min_mb):
    pieces = sorted({piece for (_, piece) in logs})
    for piece in pieces:
        matches = sorted(glob.glob(os.path.join(dest, f"*_{piece}.root")))
        if not matches:
            error(f"no *_{piece}.root in {dest}")
            continue
        for path in matches:
            mb = os.path.getsize(path) / 1e6
            if mb < min_mb:
                error(f"{os.path.basename(path)}: {mb:.2f} MB -- stub (real outputs are ~1 GB)")
            else:
                ok(f"{os.path.basename(path)}: {mb:.0f} MB")


def scan_grid_nan(dest):
    patterns = ["legacy/*.out", "w_pert/*.out", "w_asym/*.out", "get_yk_new/*.out",
                "resbos/Resbos_grids/*.out"]
    files = []
    for pat in patterns:
        files += [p for p in glob.glob(os.path.join(dest, pat)) if "_shards" not in p]
    for path in sorted(files):
        bad = 0
        with open(path, errors="replace") as f:
            for line in f:
                if NAN_TOKEN.search(line):
                    bad += 1
                    if bad >= 1:
                        break
        rel = os.path.relpath(path, dest)
        if bad:
            error(f"{rel}: contains NaN/Inf values -- the grid is corrupt")
        else:
            ok(f"{rel}: no NaN/Inf")


def check_consistency(dest, results, tol):
    bases = sorted({b for (b, _) in results})
    for base in bases:
        need = ("resNLO", "nloasy", "nlodsi")
        if not all((base, p) in results for p in need):
            continue
        (rn, en) = results[(base, "resNLO")]
        (a, ea) = results[(base, "nloasy")]
        (d, ed) = results[(base, "nlodsi")]
        total = a + d
        rel = (rn - total) / total
        sigma = math.hypot(en, math.hypot(ea, ed))
        msg = (f"{base}: resNLO = {rn:.6g} +/- {en:.3g}, nloasy + nlodsi = {total:.6g} "
               f"+/- {sigma:.3g} pb, relative difference {100*rel:+.3f} %")
        if abs(rel) > tol:
            warn(msg + f"  (above tolerance {100*tol:.1f} %)")
        else:
            ok(msg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dest")
    ap.add_argument("--cross-tol", type=float, default=0.03,
                    help="relative tolerance for resNLO vs nloasy+nlodsi (default 3 percent)")
    ap.add_argument("--min-root-mb", type=float, default=100.0,
                    help="ROOT files smaller than this are reported as stubs (default 100 MB)")
    args = ap.parse_args()

    print(f"validate_campaign: {args.dest}")
    logs = latest_logs(args.dest)
    if not logs:
        error("no resbos logs found under resbos/ -- nothing to validate")
    else:
        results = check_resbos(args.dest, logs)
        check_roots(args.dest, logs, args.min_root_mb)
        check_consistency(args.dest, results, args.cross_tol)
    scan_grid_nan(args.dest)

    verdict = "PASSED" if not errors else f"FAILED ({len(errors)} error(s))"
    extra = f", {len(warns)} warning(s)" if warns else ""
    print(f"VALIDATION {verdict}{extra}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
