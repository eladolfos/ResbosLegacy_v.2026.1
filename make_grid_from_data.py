#!/usr/bin/env python3
"""
make_grid_from_data.py -- build q_grid.inp/y_grid.inp from an experimental data table,
so a run_process.py campaign's rectangular (Q,qT,y) grid actually lands on the (Q,y)
points a fixed-target experiment measured, instead of a generic dense grid.

This is the "convert the experimental table to inp/ format" step: it does NOT loop one
job per point (like Workspace/*/e605/run.sh does) -- it takes the UNIQUE Q values and the
UNIQUE y values in the table and writes them as two grid files, so run_process.py's normal
single rectangular-grid campaign (Q list x qT list x y list) covers every real point plus,
usually, only a handful of unmeasured (Q,y) combinations (e.g. for E605: 119 real points,
18 unique Q x 7 unique y = 126 grid points -- only 7 extra). Check "extra" below to see how
many that is for your table; if it's a large fraction of the grid, this approach wastes a
lot of compute on unmeasured points and a true per-point campaign (like run.sh) fits better.

The qT grid is untouched -- point this script's output at a template whose qt_grid.inp is
already the one you want (e.g. the real per-experiment qT binning), and only q_grid.inp/
y_grid.inp get replaced.

Usage:
    python3 make_grid_from_data.py <data_table> <out_q_grid.inp> <out_y_grid.inp> \\
        [--y-col N] [--q-col N]

    data_table    whitespace-separated table, one row per measured point (blank lines and
                  lines starting with '#' are skipped)
    out_*.inp     one value per line, sorted ascending -- the same format already used by
                  templates/*/legacy/inp/q_grid.inp etc.
    --y-col N     1-based column holding y (default 1, matching the e605/e866/e906 tables
                  under New_kFactorCT25/FixedTarget_pp830a016_yao_09182026/)
    --q-col N     1-based column holding Q (default 2)

Example (E605):
    python3 make_grid_from_data.py \\
        ".../Workspace/FixedTargetKFactor/e605/e605" \\
        templates/E605_full/legacy/inp/q_grid.inp \\
        templates/E605_full/legacy/inp/y_grid.inp
    # then copy the same two files into w_pert/inp/ and w_asym/inp/ so all three match by
    # md5 (run_process.py's grid-lockstep check requires that)
"""
import argparse
import sys


def unique_sorted(values):
    """De-dup by float value, keep first-seen string form, sort ascending by value."""
    seen = {}
    for v in values:
        key = float(v)
        if key not in seen:
            seen[key] = v
    return [seen[k] for k in sorted(seen)]


def read_column(path, col):
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < col:
                sys.exit(f"{path}: line has fewer than {col} columns: {line!r}")
            out.append(parts[col - 1])
    if not out:
        sys.exit(f"{path}: no data rows found")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("data_table")
    ap.add_argument("out_q_grid")
    ap.add_argument("out_y_grid")
    ap.add_argument("--y-col", type=int, default=1)
    ap.add_argument("--q-col", type=int, default=2)
    args = ap.parse_args()

    y_vals = unique_sorted(read_column(args.data_table, args.y_col))
    q_vals = unique_sorted(read_column(args.data_table, args.q_col))
    n_rows = sum(1 for line in open(args.data_table)
                 if line.strip() and not line.strip().startswith("#"))
    grid_pts = len(q_vals) * len(y_vals)

    for path, vals in ((args.out_q_grid, q_vals), (args.out_y_grid, y_vals)):
        with open(path, "w") as f:
            f.write("\n".join(vals) + "\n")

    print(f"{args.data_table}: {n_rows} measured points -> "
          f"{len(q_vals)} unique Q x {len(y_vals)} unique y = {grid_pts} grid points "
          f"({grid_pts - n_rows} not actually measured)")
    print(f"wrote {args.out_q_grid} ({len(q_vals)} values)")
    print(f"wrote {args.out_y_grid} ({len(y_vals)} values)")


if __name__ == "__main__":
    main()
