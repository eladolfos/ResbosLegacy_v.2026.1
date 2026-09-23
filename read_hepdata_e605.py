#!/usr/bin/env python3
"""
read_hepdata_e605.py -- parse the real E605 measurement (HEPData-e605/Table1..7.yaml,
https://doi.org/10.17182/hepdata.22831.v1) into a clean (y, Q, cross_section, errors) table,
instead of guessing at the meaning of the extra, unexplained columns in Yao's own working
"e605" data file (New_kFactorCT25/.../Workspace/FixedTargetKFactor/e605/e605) -- checked: its
column 1/2 (y, Q) match these tables' 119 points exactly (to Q rounding), but its columns
4+ do NOT match the real cross section value here by any simple rescaling; don't use them.

Tables 1-7 hold the observable S*D2(SIG)/D(SQRT(TAU))/DYRAP, i.e. S * d^2(sigma)/d(sqrt(tau))/dy,
in nb*GeV^2/nucleon, at fixed y (Table i -> y = -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4 for i=1..7),
as a function of sqrt(tau) = Q/sqrt(s) (s = ECM^2 = 38.8^2 GeV^2). "/NUCLEON" means it's already
averaged over the Cu target's proton/neutron mix -- the same physical quantity FRACT_N=0.54 in
our Legacy .in files is meant to reproduce.

Unit conversion to the more familiar d^2(sigma)/dQ/dy (nb/GeV, per nucleon):
    d(sqrt(tau))/dQ = 1/sqrt(s)   (since sqrt(tau) = Q/sqrt(s))
    => d^2(sigma)/d(sqrt(tau))/dy = sqrt(s) * d^2(sigma)/dQ/dy
    => S * d^2(sigma)/d(sqrt(tau))/dy = s^(3/2) * d^2(sigma)/dQ/dy
    => d^2(sigma)/dQ/dy = (S*D2SIG/DSQRTTAU/DYRAP value) / s**1.5

Usage:
    python3 read_hepdata_e605.py [hepdata_dir]     # default: templates/FixTargetTables/HEPData-e605
    # prints the parsed table; import read_table() to use it from another script.
"""
import glob
import os
import sys

try:
    import yaml
except ImportError:
    sys.exit("needs PyYAML: pip install pyyaml")

ECM = 38.8
S = ECM ** 2


def read_table(hepdata_dir, table_range=range(1, 8)):
    """Returns a list of dicts: y, Q, sqrt_tau, s_d2sig_dsqrttau_dy (the raw HEPData value,
    nb*GeV^2/nucleon), d2sig_dQdy (converted, nb/GeV/nucleon), stat_err, sys_norm_pct,
    sys_pp_pct -- for every point that has a real measurement (HEPData's '-' placeholders,
    i.e. (Q,y) combinations in the table's own Q list with no actual measurement, are skipped).
    """
    rows = []
    for i in table_range:
        path = os.path.join(hepdata_dir, f"Table{i}.yaml")
        with open(path) as f:
            doc = yaml.safe_load(f)
        dv = doc["dependent_variables"][0]
        y = float(next(q["value"] for q in dv["qualifiers"] if q["name"] == "YRAP"))
        sqrt_taus = [v["value"] for v in doc["independent_variables"][0]["values"]]
        for st, entry in zip(sqrt_taus, dv["values"]):
            if entry.get("value") == "-":
                continue
            val = float(entry["value"])
            errs = entry.get("errors", [])
            stat = errs[0]["symerror"] if errs and "symerror" in errs[0] else None
            sys_norm = sys_pp = None
            for e in errs[1:]:
                pct = float(e["symerror"].rstrip("%")) if isinstance(e["symerror"], str) else None
                if "Normalisation" in e.get("label", ""):
                    sys_norm = pct
                elif "Point-to-point" in e.get("label", ""):
                    sys_pp = pct
            rows.append({
                "y": y, "Q": st * ECM, "sqrt_tau": st,
                "s_d2sig_dsqrttau_dy": val,
                "d2sig_dQdy": val / S ** 1.5,
                "stat_err": stat, "sys_norm_pct": sys_norm, "sys_pp_pct": sys_pp,
            })
    return rows


def main():
    hepdata_dir = sys.argv[1] if len(sys.argv) > 1 else "templates/FixTargetTables/HEPData-e605"
    rows = read_table(hepdata_dir)
    print(f"{len(rows)} measured points from {hepdata_dir} (Tables 1-7, ECM={ECM} GeV)")
    print(f"{'y':>6} {'Q':>8} {'sqrt_tau':>10} {'S*d2sig/dsqrttau/dy':>22} "
          f"{'d2sig/dQdy':>14} {'stat':>8} {'sys_norm%':>10} {'sys_pp%':>8}")
    for r in rows:
        print(f"{r['y']:>6.1f} {r['Q']:>8.3f} {r['sqrt_tau']:>10.4f} "
              f"{r['s_d2sig_dsqrttau_dy']:>22.4g} {r['d2sig_dQdy']:>14.6g} "
              f"{r['stat_err'] if r['stat_err'] is not None else float('nan'):>8.3g} "
              f"{r['sys_norm_pct'] if r['sys_norm_pct'] is not None else float('nan'):>10.1f} "
              f"{r['sys_pp_pct'] if r['sys_pp_pct'] is not None else float('nan'):>8.1f}")


if __name__ == "__main__":
    main()
