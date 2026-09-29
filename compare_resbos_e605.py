#!/usr/bin/env python3
"""
compare_resbos_e605.py -- compare a resbos .root result against the real E605 measurement
(read_hepdata_e605.py) in d^2(sigma)/dQ/dy [nb/GeV, per nucleon].

STATUS: verified against a real resbos .root (E201_e605_fine campaign, hpcc_result/
E201_e605n_A0_nocuts.root, 25M events). The tree is "h10" as guessed, but resbos writes Q/qT/y
as READY-MADE branches (M_B, pT_B, y_B -- boson mass/pT/rapidity), not a 4-momentum to
reconstruct from -- simpler than the original guess. Sanity check that passed: M_B in
[4,20] and pT_B up to 9.06 match the campaign's own grid range exactly, and sum(WT00)=101.5
(pb, at lumi=1) is the right order of magnitude for the HEPData cross sections integrated
over the same (Q,y) region.

Physics recap (see read_hepdata_e605.py's docstring for the HEPData side):
  - Branches used: M_B (=Q), pT_B (=qT), y_B (=y), WT00 (event weight, can be negative --
    normal for a matched NLO+resummed calculation with subtraction terms).
  - The measured observable is d^2(sigma)/dQ/dy, integrated over qT and over the muon decay
    angles -- so events are binned by (Q,y) only; qT is summed over, not cut on.
  - resbos's per-event weight becomes a cross section via the campaign's luminosity ([resbos]
    lumi in the .ini, pb^-1; E201_e605_fine.ini didn't override it, so it's the resbos
    template's own default, 1.0 pb^-1 -- check LUMI_PB_DEFAULT below matches what you actually
    ran with): dsigma/bin = sum(weight in bin) / LUMI_PB / (bin's DeltaQ * Deltay).
  - Bin edges are NOT in the HEPData table (only bin centers) -- estimated here as the
    midpoints between neighboring Q values within each y-slice (extended symmetrically at the
    edges) and +/-0.05 around each y (the y grid is uniform, step 0.1). This is an
    approximation Yao's own analysis likely refines somehow; revisit if the comparison looks
    systematically off in a way finer binning would explain.

Usage:
    python3 compare_resbos_e605.py <resbos_output.root> [--hepdata-dir DIR] [--lumi PB]
                                    [--plot OUTPUT.png]
"""
import argparse
import sys

# ------------------------------------------------------------------ CONFIG
TREE_NAME = "h10"                                    # resbos_root.f/froot.c: InitROOTNT -> TTree("h10","h10")
BRANCHES = {                                          # resbos_root.f ~5465-5480: ready-made boson kinematics
    "q": "M_B", "qt": "pT_B", "y": "y_B", "weight": "WT00",
}
LUMI_PB_DEFAULT = 1.0                                 # pb^-1; must match the campaign's [resbos] lumi


def load_events(root_path, tree_name=TREE_NAME, branches=BRANCHES):
    import uproot
    with uproot.open(root_path) as f:
        if tree_name not in [k.split(";")[0] for k in f.keys()]:
            sys.exit(f"{root_path}: no tree '{tree_name}' -- found {f.keys()}. "
                     f"Fix TREE_NAME in this script's CONFIG.")
        tree = f[tree_name]
        avail = set(tree.keys())
        missing = [b for b in branches.values() if b not in avail]
        if missing:
            sys.exit(f"{root_path}: tree '{tree_name}' is missing branch(es) {missing} -- "
                     f"available: {sorted(avail)}. Fix BRANCHES in this script's CONFIG.")
        arrays = tree.arrays(list(branches.values()), library="np")
    return arrays


def extract_qty(arrays, branches=BRANCHES):
    """Returns per-event (Q, qT, y, weight) numpy arrays -- resbos writes them ready-made,
    no reconstruction from a 4-momentum needed."""
    return (arrays[branches["q"]], arrays[branches["qt"]],
            arrays[branches["y"]], arrays[branches["weight"]])


def bin_edges(centers):
    """Midpoints between neighboring sorted values, extended symmetrically at both ends --
    same convention for the Q-slice (non-uniform) and the y grid (uniform, step 0.1)."""
    c = sorted(centers)
    if len(c) == 1:
        return [c[0] - 0.5, c[0] + 0.5]
    edges = [c[0] - (c[1] - c[0]) / 2]
    edges += [(c[i] + c[i + 1]) / 2 for i in range(len(c) - 1)]
    edges.append(c[-1] + (c[-1] - c[-2]) / 2)
    return edges


def theory_d2sig_dQdy(Q, qT, y, w, hep_rows, lumi_pb):
    """For every (y, Q) point in hep_rows, sum resbos event weights landing in that point's
    (Q,y) bin (qT integrated over, i.e. summed regardless of value) and normalize to
    d^2(sigma)/dQ/dy [nb/GeV/nucleon] -- the same quantity read_hepdata_e605.py reports."""
    import numpy as np
    by_y = {}
    for r in hep_rows:
        by_y.setdefault(r["y"], []).append(r["Q"])

    y_centers = sorted(by_y)
    y_edges = bin_edges(y_centers)
    q_edges_by_y = {yc: bin_edges(by_y[yc]) for yc in y_centers}

    out = []
    for r in hep_rows:
        yc, qc = r["y"], r["Q"]
        iy = y_centers.index(yc)
        ylo, yhi = y_edges[iy], y_edges[iy + 1]
        qs = sorted(by_y[yc])
        iq = qs.index(qc)
        qedges = q_edges_by_y[yc]
        qlo, qhi = qedges[iq], qedges[iq + 1]

        sel = (y >= ylo) & (y < yhi) & (Q >= qlo) & (Q < qhi)
        n = int(sel.sum())
        wsum = float(w[sel].sum()) if n else 0.0
        dQ, dy = qhi - qlo, yhi - ylo
        theory = wsum / lumi_pb / dQ / dy / 1000.0     # pb/GeV -> nb/GeV (1 nb = 1000 pb)
        out.append({**r, "n_events": n, "theory_d2sig_dQdy": theory})
    return out


def make_plot(rows, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    y_vals = sorted(set(r["y"] for r in rows))
    ncols = 4
    nrows = -(-len(y_vals) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.2 * nrows), squeeze=False)

    for ax, yc in zip(axes.flat, y_vals):
        rs = sorted((r for r in rows if r["y"] == yc), key=lambda r: r["Q"])
        Qv = [r["Q"] for r in rs]
        data = [r["d2sig_dQdy"] for r in rs]
        data_err = [(r["stat_err"] or 0) / S_TO_D2SIGDQDY for r in rs]   # same unit conversion as the value
        theory = [r["theory_d2sig_dQdy"] for r in rs]

        ax.errorbar(Qv, data, yerr=data_err, fmt="ko", ms=4, label="E605 data", capsize=2)
        ax.plot(Qv, theory, "r-o", ms=3, lw=1, label="resbos (resNLO)")
        ax.set_yscale("log")
        ax.set_title(f"y = {yc:.1f}")
        ax.set_xlabel("Q [GeV]")
        ax.set_ylabel(r"$d^2\sigma/dQ/dy$ [nb/GeV]")

    for ax in axes.flat[len(y_vals):]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")


S_TO_D2SIGDQDY = None   # set in main() once we know ECM (avoids importing read_hepdata twice)


def main():
    global S_TO_D2SIGDQDY
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("root_file")
    ap.add_argument("--hepdata-dir", default="templates/FixTargetTables/HEPData-e605")
    ap.add_argument("--lumi", type=float, default=LUMI_PB_DEFAULT, help="pb^-1")
    ap.add_argument("--plot", metavar="OUTPUT.png")
    args = ap.parse_args()

    import read_hepdata_e605 as hep
    S_TO_D2SIGDQDY = hep.S ** 1.5
    hep_rows = hep.read_table(args.hepdata_dir)

    arrays = load_events(args.root_file)
    Q, qT, y, w = extract_qty(arrays)
    print(f"{args.root_file}: {len(Q)} events, Q range [{Q.min():.2f}, {Q.max():.2f}], "
          f"y range [{y.min():.2f}, {y.max():.2f}], sum(weight)={w.sum():.4g} pb (at lumi=1)")

    rows = theory_d2sig_dQdy(Q, qT, y, w, hep_rows, args.lumi)

    print(f"\n{'y':>6} {'Q':>8} {'data':>12} {'theory':>12} {'theory/data':>12} {'n_ev':>8}")
    ratios = []
    for r in rows:
        data = r["d2sig_dQdy"]
        th = r["theory_d2sig_dQdy"]
        ratio = th / data if data else float("nan")
        if data:
            ratios.append(ratio)
        print(f"{r['y']:>6.1f} {r['Q']:>8.3f} {data:>12.6g} {th:>12.6g} {ratio:>12.4f} "
              f"{r['n_events']:>8d}")

    if ratios:
        ratios.sort()
        n = len(ratios)
        print(f"\ntheory/data: min={ratios[0]:.3f} median={ratios[n//2]:.3f} "
              f"max={ratios[-1]:.3f} mean={sum(ratios)/n:.3f}  ({n} points)")

    if args.plot:
        make_plot(rows, args.plot)


if __name__ == "__main__":
    main()
