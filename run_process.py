#!/usr/bin/env python3
"""
run_process.py -- copy a template ResBos-Legacy folder to a new place, put the
right energy / PDF / boson in every .in file, and chain the jobs with shards.

    python3 run_process.py 7TeV_WpWm_example.ini            # prepare <dest>, write <dest>/submit_all.sh
    python3 run_process.py 7TeV_WpWm_example.ini --submit   # ... and run it (on the HPCC login node)

The .ini names a SOURCE folder that holds get_yk_new/ legacy/ resbos/ w_asym/ w_pert/
(e.g. templates/7TeV_WpWm) and a DEST folder.  Only what a run needs is copied to
DEST: the executables (the fresh build in bin/ made by setup_resbos_legacy.sb, or the
ones named in [executables] in the .ini; never the old ones inside the source), the grids (./inp/*.inp) and make_dummy_rai.py.  The source
does NOT need the outputs of earlier steps: get_yk_new needs the w_pert / w_asym /
legacy outputs and resbos needs the legacy and Yk grids, so each of those jobs
brings its inputs next to itself (symlink or copy) when it starts, once the upstream
jobs have produced them.  Everything else in the source (outputs, logs, .git, old
.in files, old shard scripts) is left behind.

The .in files are the source's W+ files used as templates: the script rewrites ECM,
PDF, Type_V / JWTYPE and the active grid range; for resbos only what the .ini
overrides (cuts, seed, luminosity, ...) is changed.

Per boson:
    w_pert  ----------------\\
    w_asym  ------------------> get_yk_new --> Yk grid --\\
    legacy LTO=3 (Y piece) --/                            +--> resbos  (never sharded:
    legacy LTO=0 (main grid) -----------------------------/             it is a VEGAS MC)

w_pert, w_asym and legacy are split along Q with the scripts of shrds_scripts/
(make_shards.py, merge_shards.py, run_*_array.sb, merge_*_array.sb, submit_*_array.sh),
copied into each folder with the executable name and the header length of that code.
"""
import argparse
import configparser
import hashlib
import os
import re
import shutil
import subprocess
import sys

SELF = os.path.abspath(__file__)
HERE = os.path.dirname(SELF)
HDR_RE = re.compile(rb"\s*Q\s*,\s*qT\s*,\s*y")
PROCS = {                                            # Type_V -> (file tag, JWTYPE)
    "W+": ("Wp", 1),
    "W-": ("Wm", -1),
    "Z0": ("Z0", 2),                                 # w_pert/w_asym: JWTYPE=2 ("neutral current"); Legacy Type_V
    "A0": ("A0", 2),                                 # distinguishes Z0 (resonance) from A0 (pure photon) -- both
}                                                     # use the SAME JWTYPE=2 in w_pert/w_asym (confirmed against
                                                      # a real working E605 w_pert.in in New_kFactorCT25/
                                                      # FixedTarget_pp830a016_yao_09182026/, NOT the naive
                                                      # DATA jwm,jwp,jz,jph,JHB /-1,1,2,3,4/ jph=3 reading).
                                                      # NOTE: get_yk_new.f only knows ZU/ZD/W+/W-/A0 (get_yk_new.f:
                                                      # 181-200) -- Z0 still STOPs ('convert is not assigned to
                                                      # this Boson yet!'); A0 works at NLO (its convert factor is
                                                      # provably unused there) but STOPs at NNLO (no real
                                                      # conversion constant implemented). Z0/NNLO-A0: known, not
                                                      # handled here.
# stage: (folder, executable, label used in the shard scripts, header lines, lines/point)
# legacy_dsi (LTO=1, NLO_Sig/DeltaSigma) and legacy_asy (LTO=2, Asymptotic) are the two
# extra pieces [campaign] compute = NLO needs on top of legacy_Y (LTO=3, shared with
# resNLO) -- confirmed against legacy_final_vesion/main.for: the ~14-line preamble
# (ECM/IBEAM/PDF/masses/TYPE_V/g1-g3/... echo, written unconditionally before the LTO
# branch) plus one LTO-specific column-header line = 15 total header lines for every LTO
# value (matches the already-validated LTO=0/3 count), and both LTO=1 and LTO=2 write
# exactly one data line per (Q,qT,y) point (same shape as LTO=0/legacy_main, not LTO=3's
# 3-lines/point) -- WRITE(22,122)/WRITE(22,102) each called once per point, no inner loop
# writes multiple lines. Not yet run for real (no gfortran here); spot-check against a
# real LTO=1/LTO=2 output on the HPCC before trusting this at scale.
STAGES = {
    "w_pert":      ("w_pert", "w_pert", "w_pert", 3, 1),
    "w_asym":      ("w_asym", "w_asym", "w_asym", 3, 1),
    "legacy_Y":    ("legacy", "main", "legacy", 15, 3),
    "legacy_main": ("legacy", "main", "legacy", 15, 1),
    "legacy_dsi":  ("legacy", "main", "legacy", 15, 1),
    "legacy_asy":  ("legacy", "main", "legacy", 15, 1),
}
LTO_BY_STAGE = {"legacy_Y": "3", "legacy_main": "0", "legacy_dsi": "1", "legacy_asy": "2"}
LEGACY_SUFFIX = {"legacy_Y": "Y", "legacy_main": "main", "legacy_dsi": "dsi", "legacy_asy": "asy"}
RESNLO_STAGES = ("w_pert", "w_asym", "legacy_Y", "legacy_main")
NLO_STAGES = ("legacy_Y", "legacy_dsi", "legacy_asy")
CHUNK = 16 << 20


def die(msg):
    sys.exit("ERROR: " + msg)


# ------------------------------------------------------------------ .in editing
def find(lines, marker, what):
    """Index of the first line whose trailing '> comment' contains marker."""
    for i, l in enumerate(lines):
        if ">" in l and marker in l.split(">", 1)[1]:
            return i
    die(f"template for {what} has no line commented '> ...{marker}...'")


def put(lines, i, text):
    val, com = lines[i].rstrip("\n").split(">", 1)
    lines[i] = text.ljust(max(len(val), len(text) + 1)) + ">" + com + "\n"


def set_token(lines, i, k, new):
    """Replace the k-th comma-separated token of line i, keeping its spacing."""
    parts = lines[i].split(">", 1)[0].split(",")
    m = re.match(r"(\s*)(.*?)(\s*)$", parts[k])
    if m.group(2) == new:
        return
    parts[k] = m.group(1) + new + m.group(3)
    put(lines, i, ",".join(parts).rstrip())


def set_active(lines, act, full):
    for marker, want in (("Active setting", act), ("Full range", full)):
        i = find(lines, marker, marker)
        if [int(t) for t in lines[i].split(">", 1)[0].split()] != list(want):
            put(lines, i, "  ".join(str(a) for a in want))


def grid_paths(lines):
    return [lines[find(lines, m, m)].split(">", 1)[0].strip()
            for m in ("Q grid file name", "qT grid file name", "y grid file name")]


def patch_stage(cfg, tlines, stage, tag, jw, vtype, ecm):
    """Apply the ECM/PDF/boson/[legacy] overrides common to a normal grid campaign
    (build()) and a per-experimental-point campaign (build_points()); caller still
    has to set_active() (build()) or write single-point grid files (build_points())."""
    lines = list(tlines[stage])
    if stage in ("w_pert", "w_asym"):
        set_token(lines, find(lines, "ECM,iBeam", "ECM"), 0, ecm)
        jwtype_line = find(lines, "JWTYPE", "JWTYPE")
        set_token(lines, jwtype_line, 0, str(jw))
        if jw == 2:
            # JWTYPE=2 (jz): w_pert.f/w_asym.f pick BOSON from JZ_TYPE (1=ZU,-1=ZD,0=Z0/general
            # couplings), a DIFFERENT token the W+/W- template this was copied from leaves at 1
            # (ZU) -- confirmed against a real working w_pert.in (New_kFactorCT25/
            # FixedTarget_pp830a016_yao_09182026/Workspace/FixedTargetKFactor/e605/w_pert.in:
            # "2,0,0") that it must be reset to 0 for both our Z0 and A0 (Legacy's Type_V is what
            # actually distinguishes them there), or w_pert/w_asym silently compute wrong (ZU)
            # physics instead. w_pert: JWTYPE,IDO_CBAR,JZ_TYPE; w_asym: JWTYPE,JZ_TYPE.
            set_token(lines, jwtype_line, 2 if stage == "w_pert" else 1, "0")
        put(lines, find(lines, "PDF file name", "PDF"), "lha_" + cfg.pdf)
    else:
        l2 = find(lines, "ECM,LTO", "ECM,LTO")
        set_token(lines, l2, 0, ecm if "." in ecm else ecm + ".0")
        set_token(lines, l2, 1, LTO_BY_STAGE[stage])
        if stage == "legacy_dsi":
            # res.for's SetC1_4: LTO=1 (DeltaSigma) hard-STOPs in main.for ("THIS ONLY WORKS
            # FOR IFLAG_C3 = 1") unless the scale choice is the canonical one (IFLAG_C3=1).
            # IFLAG_C3=99 (custom C1/B0,C2,C3/B0) is only safe to auto-convert when its values
            # are numerically identical to canonical (C1/B0=C2=C3/B0=1 -> C1=B0,C2=1,C3=B0,C4=1
            # either way) -- otherwise forcing IFLAG_C3=1 would silently give legacy_dsi a
            # different renormalization/factorization scale than legacy_asy/legacy_Y use,
            # breaking the qT_Sep cancellation between them (see README.md's NLO section).
            l3 = find(lines, "IFLAG_C3", "IFLAG_C3")
            toks = [t.strip() for t in lines[l3].split(">", 1)[0].split(",")]
            iflag_c3, c_vals = toks[1], toks[2:5]
            if iflag_c3 != "1":
                same_as_canonical = iflag_c3 == "99" and all(
                    float(t.lower().replace("d", "e")) == 1.0 for t in c_vals)
                if same_as_canonical:
                    set_token(lines, l3, 1, "1")
                else:
                    die(f"legacy_dsi (LTO=1) needs IFLAG_C3=1 (Legacy's res.for STOPs otherwise: "
                        f"'THIS ONLY WORKS FOR IFLAG_C3 = 1'), but the legacy_dsi template has "
                        f"IFLAG_C3={iflag_c3} with C1/B0,C2,C3/B0={','.join(c_vals)}, not "
                        f"numerically identical to canonical (C1/B0=C2=C3/B0=1) -- auto-forcing it "
                        f"would silently run DeltaSigma at a different scale than legacy_asy/"
                        f"legacy_Y, breaking the qT_Sep cancellation between them. Point "
                        f"[templates] legacy_dsi at a template using IFLAG_C3=1 (or IFLAG_C3=99 "
                        f"with C1/B0=C2=C3/B0=1), or drop 'NLO' from [campaign] compute.")
        put(lines, find(lines, "Type_V", "Type_V"), vtype)
        put(lines, find(lines, "Evolved PDF file", "PDF"), "lha_" + cfg.pdf)
        if cfg.get("legacy", "bmax"):
            put(lines, find(lines, "bMax", "bMax"), cfg.get("legacy", "bmax"))
        if cfg.get("legacy", "nonpert"):
            put(lines, find(lines, "g1, g2, g3, Q0, nG", "nonpert"), cfg.get("legacy", "nonpert"))
        if cfg.get("legacy", "ibeam"):
            set_token(lines, find(lines, "FRACT_N", "ibeam/fract_n"), 0, cfg.get("legacy", "ibeam"))
        if cfg.get("legacy", "fract_n"):
            set_token(lines, find(lines, "FRACT_N", "ibeam/fract_n"), 1, cfg.get("legacy", "fract_n"))
    return lines


def read_points(path, y_col, q_col):
    """Parse an experimental data table (whitespace-separated, blank/'#' lines skipped)
    into a list of (y, Q) string pairs, 1-based y_col/q_col -- same convention as Yao's
    run.sh scripts (awk '{print $1}'=y, '{print $2}'=Q) in New_kFactorCT25/
    FixedTarget_pp830a016_yao_09182026/Workspace/*/e605(etc.)/run.sh."""
    if not os.path.isfile(path):
        die(f"[grids] experimental: {path} not found")
    pts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            need = max(y_col, q_col)
            if len(parts) < need:
                die(f"{path}: line has fewer than {need} columns: {line!r}")
            pts.append((parts[y_col - 1], parts[q_col - 1]))
    if not pts:
        die(f"{path}: no data rows found")
    return pts


def read_lines(path):
    if not os.path.isfile(path):
        die(f"missing template {path}")
    with open(path) as f:
        return f.readlines()


def write(path, text, exe=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)
    if exe:
        os.chmod(path, 0o755)


def copy(src, dst, exe=False):
    if not os.path.isfile(src):
        die(f"missing {src}")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)
    os.chmod(dst, 0o755 if exe else 0o644)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(CHUNK), b""):
            h.update(b)
    return h.hexdigest()


def check_jobname(job):
    name = f"{job}_shards/{job}_shard00"      # Fortran job names are character*100
    if len(name) > 95:
        die(f"'{name}' is {len(name)} chars; the codes truncate at 100. Shorten [campaign] name")


# ------------------------------------------------------------------ config
class Cfg:
    def __init__(self, path):
        cp = configparser.ConfigParser(inline_comment_prefixes=("#",), interpolation=None)
        cp.optionxform = str
        if not cp.read(path):
            die(f"cannot read {path}")
        self.cp, self.path = cp, os.path.abspath(path)
        base = os.path.dirname(self.path)
        self.rel = lambda p: os.path.abspath(os.path.join(base, os.path.expanduser(p)))
        self.src = self.rel(self.req("paths", "source"))
        self.dest = self.rel(self.req("paths", "dest"))
        self.shard_scripts = self.rel(self.get("paths", "shard_scripts", os.path.join(HERE, "shrds_scripts")))
        self.stage = self.get("paths", "stage_outputs", "link")      # how a job receives upstream outputs
        if self.stage not in ("link", "copy"):
            die("[paths] stage_outputs must be 'link' or 'copy'")
        self.name = self.req("campaign", "name")
        self.ecm = self.req("campaign", "ecm")
        float(self.ecm)
        self.pdf = self.req("campaign", "pdf")
        self.order = self.req("campaign", "order").upper()
        if self.order not in ("NLO", "NNLO"):
            die("[campaign] order must be NLO or NNLO")
        self.procs = [p.strip() for p in self.req("campaign", "processes").split(",") if p.strip()]
        for p in self.procs:
            if p not in PROCS:
                die(f"process '{p}' not supported yet (supported: {', '.join(PROCS)})")
        # [campaign] compute: resNLO (default, unchanged behaviour), NLO (fixed order via
        # Legacy LTO=1+LTO=2+LTO=3 phase-space slicing + two resbos runs + hadd -- see
        # README.md "NLO (fixed order)"), or "NLO, resNLO" for both in one campaign (legacy_Y,
        # LTO=3, is computed once and shared between them).
        norm = {"nlo": "NLO", "resnlo": "resNLO"}
        self.compute = set()
        for c in self.get("campaign", "compute", "resNLO").split(","):
            c = c.strip()
            if not c:
                continue
            key = norm.get(c.lower())
            if not key:
                die(f"[campaign] compute: '{c}' not recognized (use NLO, resNLO, or 'NLO, resNLO')")
            self.compute.add(key)
        if not self.compute:
            die("[campaign] compute: at least one of NLO, resNLO required")
        # [grids] experimental: "ExpCustomGrid" mode -- one job per (y,Q) row of a real
        # experimental data table (like Yao's run.sh), instead of one shared rectangular
        # grid. Only worth it when the table's unique-Q x unique-y product is close to its
        # row count (e.g. E605: 119 points, 18x7=126); otherwise a normal template grid
        # (this key unset) or a grid built with make_grid_from_data.py fits better -- see
        # that script's docstring and E201_e605_full.ini vs the run.sh-loop experiments
        # (e866f/e906aF/e866ppxf) discussed when this feature was added.
        self.experimental = self.get("grids", "experimental")
        if self.experimental:
            self.experimental = self.rel(self.experimental)
        self.exp_y_col = int(self.get("grids", "exp_y_col", "1"))
        self.exp_q_col = int(self.get("grids", "exp_q_col", "2"))
        self.exp_throttle = self.get("grids", "throttle")   # e.g. "20" -> SLURM --array=1-N%20

        # [grids] generate: build a fine, dense Q x y rectangle (make_grid_from_data.py's
        # "generate" mode, called in-process) into the template's own legacy/w_pert/w_asym
        # inp/ before the normal build() runs -- so [templates] doesn't need a pre-baked
        # grid checked in. Mutually exclusive with [grids] experimental (that mode doesn't
        # use a shared rectangular grid at all).
        self.gen_grid = self.get("grids", "generate", "").strip().lower() in ("yes", "true", "1")
        if self.gen_grid:
            if self.experimental:
                die("[grids] generate and [grids] experimental are mutually exclusive")
            self.gen_n_q = int(self.req("grids", "n_q"))
            self.gen_n_y = int(self.req("grids", "n_y"))
            self.gen_q_min = self.get("grids", "q_min")
            self.gen_q_max = self.get("grids", "q_max")
            self.gen_q_from = self.get("grids", "q_from")
            if self.gen_q_from:
                self.gen_q_from = self.rel(self.gen_q_from)
            self.gen_q_col = int(self.get("grids", "q_col", "2"))
            self.gen_q_spacing = self.get("grids", "q_spacing", "linear")
            self.gen_y_min = self.get("grids", "y_min")
            self.gen_y_max = self.get("grids", "y_max")
            self.gen_y_from = self.get("grids", "y_from")
            if self.gen_y_from:
                self.gen_y_from = self.rel(self.gen_y_from)
            self.gen_y_col = int(self.get("grids", "y_col", "1"))
            self.gen_y_spacing = self.get("grids", "y_spacing", "linear")

    def get(self, sec, key, default=None):
        return self.cp.get(sec, key, fallback=default)

    def req(self, sec, key):
        v = self.get(sec, key)
        if v is None or not v.strip():
            die(f"missing [{sec}] {key}")
        return v.strip()

    def exe_src(self, folder, exe):
        """Executable to copy for <folder>: [executables] <folder> = file; else <dir>/<folder>/<exe> with
        [executables] dir (default: bin/ next to this script, where setup_resbos_legacy.sb puts the fresh build).
        The executables stored in the source template are never used unless dir points at it."""
        p = self.get("executables", folder, "").strip()
        if p:
            return self.rel(p)
        d = self.get("executables", "dir", "").strip()
        return os.path.join(self.rel(d) if d else os.path.join(HERE, "bin"), folder, exe)

    def shards(self, stage):
        return int(self.get("shards", stage.split("_")[0] if stage.startswith("legacy") else stage, "1"))


# ------------------------------------------------------------------ pieces
def instantiate_shard_scripts(cfg, folder, label, exe, header_lines, npts, clean_shards):
    """Copy shrds_scripts/ into dest/<folder>, renaming w_asym -> label (exe for the srun line)."""
    d = os.path.join(cfg.dest, folder)
    sd = cfg.shard_scripts

    def conv(text):
        text = text.replace("srun ./w_asym", "srun ./@EXE@").replace("w_asym", label)
        return text.replace("@EXE@", exe)

    def read(name):
        p = os.path.join(sd, name)
        if not os.path.isfile(p):
            die(f"{p} not found (set [paths] shard_scripts)")
        return open(p).read()

    write(os.path.join(d, "make_shards.py"), conv(read("make_shards.py")), exe=True)
    ms = conv(read("merge_shards.py"))
    ms, n = re.subn(r"(?m)^HEADER_LINES = \d+.*$", f"HEADER_LINES = {header_lines}", ms)
    if n != 1:
        die("could not find 'HEADER_LINES = N' in shrds_scripts/merge_shards.py")
    write(os.path.join(d, "merge_shards.py"), ms, exe=True)
    write(os.path.join(d, f"run_{label}_array.sb"), conv(read("run_w_asym_array.sb")), exe=True)
    # after merging, verify the row count (legacy Y-piece files have 3 lines per point, the rest 1)
    check = (f'\n# row-count check + timing log added by run_process.py\n'
             f'STAGE_T0="$2"   # submit_*_array.sh passes when the array job was submitted\n'
             f'case "$JOBNAME" in *_legacy_*_Y) K=3;; *) K=1;; esac\n'
             f'python3 "{SELF}" _check "${{JOBNAME}}.out" $K {npts}\n'
             f'CHECK_EXIT=$?\n'
             f'{timing_fn(cfg)}'
             f'log_elapsed "{label} merge ($JOBNAME, exit=$CHECK_EXIT)" "$STAGE_T0"\n'
             f'[ $CHECK_EXIT -eq 0 ] || exit 1\n')
    if clean_shards:
        # only reached if the check above passed (it exits 1 otherwise): the shards are kept for debugging on failure
        check += (f'echo "--clean-shards: removing ${{JOBNAME}}_shards/ (merged output verified above)"\n'
                  f'rm -rf "${{JOBNAME}}_shards"\n')
    write(os.path.join(d, f"merge_{label}_array.sb"), conv(read("merge_w_asym_array.sb")) + check, exe=True)
    write(os.path.join(d, f"submit_{label}_array.sh"), conv(read("submit_w_asym_array.sb")), exe=True)


def env_lines(cfg, root=False):
    g = lambda k, d: cfg.get("environment", k, d)
    out = ["module purge", "module load " + g("modules", "GCC/13.2.0 OpenMPI/4.1.6-GCC-13.2.0 powertools")]
    if root:
        out.append("source " + g("root_setup", "/mnt/home/lopezels/InstallSources/ROOT/bin/thisroot.sh"))
    for key, d in (("lhapdf", "/mnt/home/lopezels/InstallSources/LHAPDF"),
                   ("hoppet", "/mnt/home/lopezels/InstallSources/HOPPET1")):
        p = g(key, d)
        out += [f"export PATH={p}/bin:$PATH", f"export LD_LIBRARY_PATH={p}/lib:$LD_LIBRARY_PATH"]
    return "\n".join(out) + "\n"


def stage_fn(cfg):
    """Bash function `stage SRC DST`: bring the output of an earlier step next to the job that needs it."""
    if cfg.stage == "copy":
        op = 'cp -f "$1" "$2"'
    else:
        op = 'ln -sfn "$(realpath --relative-to="$(dirname "$2")" "$1")" "$2"'
    return ('stage() {   # stage SRC DST\n'
            '  [ -s "$1" ] || { echo "ERROR: missing input $1 (did the upstream job fail?)"; exit 1; }\n'
            f'  {op}\n'
            '}\n')


def timing_header_lines(cfg):
    """submit_all.sh: write dest/timing.log's first line (T0, the campaign's own start time,
    as an epoch timestamp) so every later stage's log_elapsed() can report elapsed-since-start."""
    log = os.path.join(cfg.dest, "timing.log")
    return [f'TIMING_LOG="{log}"',
            'echo "T0 $(date +%s)  ($(date \'+%F %T\'), campaign start)" > "$TIMING_LOG"', '']


def timing_fn(cfg):
    """Bash function `log_elapsed STAGE [STAGE_START_EPOCH]`: append a locked, timestamped
    line to dest/timing.log. Always reports elapsed-since-campaign-start (T0, written by
    submit_all.sh). Also reports that STAGE's own isolated duration, two ways:
    - STAGE_START_EPOCH given (array+merge stages: w_pert/w_asym/legacy_Y/legacy_main --
      the merge job doesn't know how long the preceding array took, so build()/build_points()
      pass down the epoch when the array was submitted -- the array itself has no
      dependency, so submission time is effectively its start time).
    - STAGE_START_EPOCH omitted (get_yk_new/resbos -- single, non-array jobs that DO have a
      --dependency on earlier stages, so bash's own $SECONDS, counting from when this
      script itself started running, i.e. only once SLURM actually dispatched it, already
      excludes time spent waiting on those dependencies -- no argument-passing needed)."""
    log = os.path.join(cfg.dest, "timing.log")
    return (f'TIMING_LOG="{log}"\n'
            'log_elapsed() {   # log_elapsed STAGE [STAGE_START_EPOCH]\n'
            '  local t0 now elapsed stage_t0 dur dur_str\n'
            '  t0=$(awk \'/^T0 /{print $2; exit}\' "$TIMING_LOG" 2>/dev/null)\n'
            '  now=$(date +%s)\n'
            '  elapsed=$([ -n "$t0" ] && echo $((now - t0)) || echo -1)\n'
            '  stage_t0="$2"\n'
            '  dur=$([ -n "$stage_t0" ] && echo $((now - stage_t0)) || echo $SECONDS)\n'
            '  dur_str=$(printf \'%ss (%dh%02dm%02ds)\' "$dur" $((dur/3600)) $(((dur%3600)/60)) $((dur%60)))\n'
            '  { flock -x 201\n'
            '    printf \'%s  %-40s  duration=%-18s elapsed_since_start=%ss (%dh%02dm%02ds)\\n\' \\\n'
            '      "$(date \'+%F %T\')" "$1" "$dur_str" "$elapsed" \\\n'
            '      $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60)) \\\n'
            '      >> "$TIMING_LOG"\n'
            '  } 201>>"$TIMING_LOG.lock"\n'
            '}\n')


def yk_script(cfg, tag, vtype, jobs, npts, rai_file, make_rai):
    job = f"{cfg.name}_Yk_{tag}_{cfg.order}"
    ins = [(f'../{folder}/{jobs[stage]}.out', f'{jobs[stage]}.out')
           for folder, stage in (("w_pert", "w_pert"), ("w_asym", "w_asym"), ("legacy", "legacy_Y"))]
    staging = "".join(f"stage {src} {dst}\n" for src, dst in ins)
    body = f"""#!/bin/bash --login
#SBATCH --job-name={job}
#SBATCH --output=slurm_yk_{tag}_%j.out
#SBATCH --error=slurm_yk_{tag}_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=general-long

cd ${{SLURM_SUBMIT_DIR}}
{env_lines(cfg)}
{stage_fn(cfg)}
# inputs = outputs of w_pert, w_asym and legacy (LTO=3)
{staging}{make_rai}srun ./get_yk_new {" ".join(dst for _, dst in ins)} {rai_file} {job}.out {vtype} {cfg.order} {cfg.pdf}
EXIT_CODE=$?
echo "get_yk_new finished with exit code $EXIT_CODE"
[ $EXIT_CODE -eq 0 ] && python3 "{SELF}" _check {job}.out 3 {npts} || EXIT_CODE=1
{timing_fn(cfg)}
log_elapsed "get_yk_new ({tag}, exit=$EXIT_CODE)"
exit $EXIT_CODE
"""
    return job, body


def resbos_script(cfg, rjob, main_src, y_src):
    """main_src/y_src are (folder, job) pairs naming the "Main data grid"/"Y piece grid"
    inputs -- for resNLO these are legacy(LTO=0)/get_yk_new's Yk grid; for NLO's two
    phase-space-slicing runs they're both legacy outputs (asy=LTO=2 main + Y=LTO=3 Y-piece
    directly, or dsi=LTO=1 main + no Y piece at all, y_src=None -> the .in gets "-")."""
    main_folder, main_job = main_src
    staging = f"stage ../{main_folder}/{main_job}.out Resbos_grids/{main_job}.out\n"
    if y_src:
        y_folder, y_job = y_src
        staging += f"stage ../{y_folder}/{y_job}.out Resbos_grids/{y_job}.out\n"
    return f"""#!/bin/bash --login
#SBATCH --job-name={rjob}
#SBATCH --output=slurm_{rjob}_%j.out
#SBATCH --error=slurm_{rjob}_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1      # resbos_root is serial
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=general-long

cd ${{SLURM_SUBMIT_DIR}}
{env_lines(cfg, root=True)}
{stage_fn(cfg)}
mkdir -p Resbos_grids
{staging}
LOG_FILE="resbos_output_{rjob}_${{SLURM_JOB_ID}}.log"
srun ./resbos_root {rjob} >> "$LOG_FILE" 2>&1
EXIT_CODE=$?
echo "resbos_root finished with exit code $EXIT_CODE" | tee -a "$LOG_FILE"
# resbos_root.f writes "<jobname>.root" in the CWD -- move the final result up to the
# campaign's top-level dest/ folder so it's easy to find/copy without digging into resbos/
if [ $EXIT_CODE -eq 0 ] && [ -f "{rjob}.root" ]; then
  mv -f "{rjob}.root" "{cfg.dest}/{rjob}.root"
  echo "moved {rjob}.root -> {cfg.dest}/{rjob}.root" | tee -a "$LOG_FILE"
fi
{timing_fn(cfg)}
log_elapsed "resbos ({rjob}, exit=$EXIT_CODE)"
exit $EXIT_CODE
"""


def hadd_script(cfg, hjob, asy_root, dsi_root, out_root):
    """NLO = phase-space-sliced fixed order: qT>qT_Sep (asy+Y piece resbos run) + qT<qT_Sep
    (dsi-alone resbos run), combined event-by-event with ROOT's hadd -- the qT_Sep
    dependence is designed to cancel between the two files (see README.md "NLO (fixed
    order)"). hadd needs ROOT on PATH, same as resbos_root (env_lines(cfg, root=True))."""
    return f"""#!/bin/bash --login
#SBATCH --job-name={hjob}
#SBATCH --output=slurm_{hjob}_%j.out
#SBATCH --error=slurm_{hjob}_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=general-long

cd ${{SLURM_SUBMIT_DIR}}
{env_lines(cfg, root=True)}
hadd -f "{cfg.dest}/{out_root}" "{cfg.dest}/{asy_root}" "{cfg.dest}/{dsi_root}"
EXIT_CODE=$?
echo "hadd finished with exit code $EXIT_CODE"
{timing_fn(cfg)}
log_elapsed "hadd NLO ({hjob}, exit=$EXIT_CODE)"
exit $EXIT_CODE
"""


def write_resbos_run(cfg, rdir, tag, run, main_src, y_src, resbos_lines, suffix=None):
    """Write resbos/<name>_<tag>_<run>[_<suffix>].in + its .sb from the resbos template +
    [cuts.<run>]; return the job name. main_src/y_src: see resbos_script(); y_src=None
    writes "-" (no Y piece) for the NLO dsi-alone run."""
    sec = "cuts." + run
    if run != "default" and not cfg.cp.has_section(sec):
        die(f"[resbos] runs lists '{run}' but there is no [{sec}] section")
    lines = list(resbos_lines)
    # ISEED and the VEGAS settings: only what the .ini overrides
    i = find(lines, "seed for random", "seed")
    toks = [t.strip() for t in lines[i].split(">", 1)[0].split(",")]
    if cfg.get("resbos", "vegas"):
        toks[:5] = [t.strip() for t in cfg.req("resbos", "vegas").split(",")]
    if cfg.get("resbos", "seed"):
        toks[5] = cfg.req("resbos", "seed")
    put(lines, i, ",".join(toks))
    put(lines, find(lines, "Main data grid", "main grid"), f"./Resbos_grids/{main_src[1]}.out")
    put(lines, find(lines, "Y piece grid", "Y grid"),
        f"./Resbos_grids/{y_src[1]}.out" if y_src else "-")
    if y_src and y_src[0] == "legacy":
        # NLO's asy+Y run feeds Legacy's raw LTO=3 output directly (10 columns, 15 header
        # lines: 'Q,qT,y, Singular (L0,A3), Pert. (L0,A3,A1,A2,A4)', main.for ~1872) --
        # resbos_root.f's iYGrid selects the Y-grid *format* it expects (read from resbos.in,
        # not auto-detected from the file): iYGrid=1 is exactly this raw 15-header/10-column
        # layout; iYGrid=2 (what every resbos.in template ships with, tuned for resNLO) is
        # get_yk_new's newer 17-header/11-column format with the extra R_Ai columns appended
        # -- feeding it a raw Legacy file makes resbos_root.f print "Must have 17 comment
        # lines in the Y-Grid file" and Call Exit (silently, exit code 0, no events). Force
        # iYGrid=1 only for this run; resNLO's own run (y_src=("get_yk_new", ...)) keeps
        # whatever the template already has.
        set_token(lines, find(lines, "iYG", "iYGrid"), 3, "1")
    for key, marker in (("lepton", "Cuts(1)"), ("mass_qt_y", "Cuts(2)"), ("mt_met", "Cuts(3)")):
        v = cfg.get(sec, key)
        if v:
            if "{qwindow}" in v:
                v = v.replace("{qwindow}", cfg.req("resbos", "qwindow"))
            put(lines, find(lines, marker, key), v)
    if cfg.get("resbos", "lumi"):
        set_token(lines, find(lines, "Luminosity", "luminosity"), 0, cfg.req("resbos", "lumi"))
    if cfg.get("resbos", "output"):
        put(lines, find(lines, "Output fromat", "output"), cfg.req("resbos", "output"))
    rjob = f"{cfg.name}_{tag}_{run}" + (f"_{suffix}" if suffix else "")
    write(os.path.join(rdir, rjob + ".in"), "".join(lines))
    write(os.path.join(rdir, f"run_{rjob}.sb"), resbos_script(cfg, rjob, main_src, y_src), exe=True)
    return rjob


def resbos_add(cfg, args):
    """--resbos-only: dest already has the w_pert/w_asym/legacy/get_yk_new outputs of an earlier run;
    only (re)generate and optionally submit the resbos job(s) for [resbos] runs, reusing them directly."""
    if "resNLO" not in cfg.compute:
        die("--resbos-only re-runs the resNLO resbos step (Legacy main grid + get_yk_new's Yk grid); "
            "add 'resNLO' to [campaign] compute")
    if not os.path.isdir(cfg.dest):
        die(f"{cfg.dest} does not exist: run the full campaign first (no --resbos-only)")
    resbos_lines = read_lines(os.path.join(cfg.src, cfg.req("templates", "resbos")))
    rdir = os.path.join(cfg.dest, "resbos")
    os.makedirs(os.path.join(rdir, "Resbos_grids"), exist_ok=True)
    exe_src = cfg.exe_src("resbos", "resbos_root")
    if not os.path.isfile(exe_src):
        die(f"executable not found: {exe_src}\n       build it first with:  sbatch setup_resbos_legacy.sb")
    copy(exe_src, os.path.join(rdir, "resbos_root"), exe=True)
    runs = [r.strip() for r in cfg.get("resbos", "runs", "default").split(",") if r.strip()]
    to_submit = []
    for vtype in cfg.procs:
        tag, _ = PROCS[vtype]
        main_job = f"{cfg.name}_legacy_{tag}_main"
        yk_job = f"{cfg.name}_Yk_{tag}_{cfg.order}"
        main_out = os.path.join(cfg.dest, "legacy", main_job + ".out")
        yk_out = os.path.join(cfg.dest, "get_yk_new", yk_job + ".out")
        missing = [p for p in (main_out, yk_out) if not os.path.isfile(p)]
        if missing:
            die("missing output(s) of the main campaign, run (and wait for) it first:\n       "
                + "\n       ".join(missing))
        print(f"\n== {vtype}")
        for run in runs:
            rjob = write_resbos_run(cfg, rdir, tag, run, ("legacy", main_job), ("get_yk_new", yk_job), resbos_lines)
            print(f"  resbos       {rjob}.in")
            to_submit.append((rdir, f"run_{rjob}.sb", rjob))
    if args.submit:
        if not shutil.which("sbatch"):
            die("sbatch not found: use --submit on the HPCC login node")
        for d, script, rjob in to_submit:
            out = subprocess.check_output(["sbatch", "--parsable", script], cwd=d, text=True).strip()
            print(f"submitted {rjob}: job {out}")
    else:
        print("\nnext:")
        for d, script, _ in to_submit:
            print(f'  (cd "{d}" && sbatch {script})')


# ------------------------------------------------------------------ ExpCustomGrid (ye. one job per point)
def points_array_script(cfg, exe, prefix, width, npts):
    """SLURM array job: task i cds into points/pt<i>/ (already holding <prefix>_pt<i>.in
    and its own single-point inp/) and runs the executable there, Yao-run.sh-style."""
    job = f"{prefix}_points"
    arr = f"1-{npts}" + (f"%{cfg.exp_throttle}" if cfg.exp_throttle else "")
    return job, f"""#!/bin/bash --login
#SBATCH --job-name={job}
#SBATCH --output=slurm_{job}_%A_%a.out
#SBATCH --error=slurm_{job}_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=general-long
#SBATCH --array={arr}

cd ${{SLURM_SUBMIT_DIR}}
{env_lines(cfg)}
PT=$(printf "%0{width}d" $SLURM_ARRAY_TASK_ID)
JOBNAME="{prefix}_pt${{PT}}"
cd "points/pt${{PT}}"
srun ../../{exe} "$JOBNAME"
EXIT_CODE=$?
echo "$JOBNAME finished with exit code $EXIT_CODE"
exit $EXIT_CODE
"""


def points_merge_script(cfg, prefix, width, npts, header_lines, per_point, nqt):
    """Strip each point's header (merge_points.py CLI below) and concatenate, then verify
    the combined row count -- same idea as merge_shards.py/the shard-array's row-count check."""
    job = f"{prefix}_merge"
    return job, f"""#!/bin/bash --login
#SBATCH --job-name={job}
#SBATCH --output=slurm_{job}_%j.out
#SBATCH --error=slurm_{job}_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=general-long

cd ${{SLURM_SUBMIT_DIR}}
STAGE_T0="$1"   # when the points array job was submitted, passed by build_points()
python3 "{SELF}" _merge_points {prefix} {width} {npts} {header_lines} {prefix}_combined.out
EXIT_CODE=$?
[ $EXIT_CODE -eq 0 ] && python3 "{SELF}" _check {prefix}_combined.out {per_point} {npts * nqt} || EXIT_CODE=1
{timing_fn(cfg)}
log_elapsed "{prefix} merge (exit=$EXIT_CODE)" "$STAGE_T0"
exit $EXIT_CODE
"""


def unique_sorted_vals(values):
    """De-dup a list of number strings by float value (keep first-seen string form),
    sorted ascending -- same rule as make_grid_from_data.py's unique_sorted()."""
    seen = {}
    for v in values:
        k = float(v)
        if k not in seen:
            seen[k] = v
    return [seen[k] for k in sorted(seen)]


def write_rect_rai(path, q_vals, y_vals, qt_vals, ecm, vtype, pdf):
    """Write a dummy (all R_Ai=1) R_Ai file as a genuine rectangular grid (Q outer, y
    middle, qT inner -- CheckSum's required nesting), covering every unique Q/y pair from
    an ExpCustomGrid campaign's points instead of just its (irregular) real point set.
    Exact for NLO (get_yk_new.f only applies R_Ai when iorder=="NNLO"): the combined
    w_pert/w_asym/legacy_Y files still only have the real points -- this file only has to
    satisfy get_yk_new.f's CheckSum (Sample()'s block-counting requires NQ*Ny*NqT ==
    total rows, which an irregular point set can't provide even with all-1 values)."""
    lines = ["ECM, TYPE_V, PDF\n", f"{ecm} {vtype} {pdf}\n", "  Q,qT,y R_A0 R_A1 R_A2 R_A3\n"]
    for q in q_vals:
        for y in y_vals:
            for qt in qt_vals:
                lines.append(f"{q} {qt} {y} 1 1 1 1\n")
    write(path, "".join(lines))


def build_points(cfg, args):
    """[grids] experimental set: one job per (y,Q) row of the data table (Yao's run.sh
    pattern) instead of one shared rectangular grid. w_pert/w_asym/legacy_Y/legacy_main
    run per point and are merged into <prefix>_combined.out each; get_yk_new and resbos
    then run once per boson against those combined files, exactly like build()'s normal
    campaign does against its single shared grid's outputs -- the merged files keep the
    same (Q,qT,y) row order across w_pert/w_asym/legacy_Y (same points, same qT grid,
    written identically for all four stages), which is all get_yk_new's row-by-row match
    check (eps=1e-8) needs. Only resNLO is supported here (NLO fixed order needs a genuine
    dense Q x y rectangle for its LTO=1/2/3 pieces, same reason resbos itself can't run on
    an irregular ExpCustomGrid point set -- see the WARNING below)."""
    if cfg.compute != {"resNLO"}:
        die("[grids] experimental only supports [campaign] compute = resNLO "
            "(NLO fixed order needs a dense Q x y rectangle -- use [grids] generate instead)")
    for d in ("get_yk_new", "legacy", "resbos", "w_asym", "w_pert"):
        if not os.path.isdir(os.path.join(cfg.src, d)):
            die(f"source folder {cfg.src} has no '{d}/'")
    if os.path.exists(cfg.dest) and not args.reuse:
        die(f"{cfg.dest} already exists (use --reuse to refresh its .in files and scripts)")
    if os.path.realpath(cfg.dest) == os.path.realpath(cfg.src):
        die("dest must differ from source")
    if args.clean_shards:
        print("note: --clean-shards has no effect in [grids] experimental mode (no shard dirs are made)")

    tpl = {k: cfg.req("templates", k) for k in ("legacy_Y", "legacy_main", "w_pert", "w_asym", "resbos")}
    tlines = {k: read_lines(os.path.join(cfg.src, v)) for k, v in tpl.items()}
    points = read_points(cfg.experimental, cfg.exp_y_col, cfg.exp_q_col)
    npts = len(points)
    width = len(str(npts))

    # ---- qT grid: every code must read the identical one. Q/y become single-value files
    # the driver writes itself (identically for all four stages), so no lockstep check is
    # needed for those the way build() needs one for its one shared Q/qT/y grid files.
    ref_qt = None
    for stage in ("w_pert", "w_asym", "legacy_Y", "legacy_main"):
        folder = STAGES[stage][0]
        qt_rel = grid_paths(tlines[stage])[1]
        p = os.path.join(cfg.src, folder, qt_rel)
        if not os.path.isfile(p):
            die(f"grid {p} not found")
        if ref_qt is None:
            ref_qt = (p, md5(p))
        elif md5(p) != ref_qt[1]:
            die(f"qT grid differs between {ref_qt[0]} and {p}: all codes must use the same qT grid")
    nqt = sum(1 for l in open(ref_qt[0]) if l.strip())
    total_pts = npts * nqt
    print(f"{cfg.experimental}: {npts} experimental points x {nqt} qT points each "
          f"({total_pts} rows/stage/boson)")

    exes = [(STAGES[k][0], STAGES[k][1]) for k in ("w_pert", "w_asym", "legacy_Y")]
    exes += [("get_yk_new", "get_yk_new"), ("resbos", "resbos_root")]
    missing = [cfg.exe_src(f, e) for f, e in exes if not os.path.isfile(cfg.exe_src(f, e))]
    if missing:
        die("executable(s) not found:\n       " + "\n       ".join(missing) + "\n"
            "       build them first with:  sbatch setup_resbos_legacy.sb   (writes bin/<code>/<exe>)\n"
            "       or point to them in the .ini:  [executables] dir = <folder>  (or one key per code)")

    print(f"copying to {cfg.dest}")
    for folder, exe in exes:
        src = cfg.exe_src(folder, exe)
        print(f"  {folder}/{exe} <- {src}")
        copy(src, os.path.join(cfg.dest, folder, exe), exe=True)
    # (make_dummy_rai.py isn't used here: it reads a real w_pert.out's own point set, which
    # in this mode is the irregular one that fails get_yk_new.f's CheckSum -- see
    # write_rect_rai() below, which builds the full rectangular Q x y closure instead)
    q_vals = unique_sorted_vals([q for _, q in points])
    y_vals = unique_sorted_vals([y for y, _ in points])
    qt_vals = [l.strip() for l in open(ref_qt[0]) if l.strip()]

    if npts != len(q_vals) * len(y_vals):
        print(f"WARNING: these {npts} points are not a perfect Q x y rectangle "
              f"({len(q_vals)} unique Q x {len(y_vals)} unique y = "
              f"{len(q_vals) * len(y_vals)}, {len(q_vals) * len(y_vals) - npts} missing). "
              f"get_yk_new's R_Ai dummy is padded to the full rectangle to pass its own\n"
              f"         CheckSum (harmless -- unused values), but resbos_root.f's OWN "
              f"CheckSum on the 'Main data grid' (Legacy's REAL, non-dummy output) has the\n"
              f"         same NQ*Ny*NqT==rows requirement and CANNOT be padded the same way "
              f"without faking real physics -- get_yk_new will likely still produce a Yk\n"
              f"         grid, but resbos will very likely fail with 'File checksum BAD. The "
              f"grid file is corrupt!' on legacy_..._main_combined.out. For a real resbos/\n"
              f"         .root result, use a rectangular grid instead (make_grid_from_data.py "
              f"+ a normal [templates] campaign, e.g. E201_e605_full.ini for E605).")

    sub = ["#!/bin/bash", "set -e", *timing_header_lines(cfg)]
    ecm = cfg.ecm
    for vtype in cfg.procs:
        tag, jw = PROCS[vtype]
        print(f"\n== {vtype}")
        combined, merge_ids = {}, {}
        for stage in ("w_pert", "w_asym", "legacy_Y", "legacy_main"):
            folder, exe, _, hdr, k = STAGES[stage]
            prefix = (f"{cfg.name}_{stage}_{tag}" if stage in ("w_pert", "w_asym") else
                      f"{cfg.name}_legacy_{tag}_{'Y' if stage == 'legacy_Y' else 'main'}")
            for i, (y_val, q_val) in enumerate(points, start=1):
                pt = str(i).zfill(width)
                job = f"{prefix}_pt{pt}"
                if len(job) > 95:
                    die(f"'{job}' is {len(job)} chars; the codes truncate at 100. "
                        f"Shorten [campaign] name")
                lines = patch_stage(cfg, tlines, stage, tag, jw, vtype, ecm)
                q_rel, qt_rel, y_rel = grid_paths(lines)
                set_active(lines, (1, nqt, 1, 1, 1, 1, 1, 1, 1), (1, nqt, 1, 1, 1, 1, 1, 1, 1))
                pdir = os.path.join(cfg.dest, folder, "points", f"pt{pt}")
                write(os.path.join(pdir, q_rel), q_val + "\n")
                write(os.path.join(pdir, y_rel), y_val + "\n")
                copy(os.path.join(cfg.src, folder, qt_rel), os.path.join(pdir, qt_rel))
                write(os.path.join(pdir, job + ".in"), "".join(lines))
            fdir = os.path.join(cfg.dest, folder)
            aj, abody = points_array_script(cfg, exe, prefix, width, npts)
            write(os.path.join(fdir, f"run_{aj}.sb"), abody, exe=True)
            mj, mbody = points_merge_script(cfg, prefix, width, npts, hdr, k, nqt)
            write(os.path.join(fdir, f"run_{mj}.sb"), mbody, exe=True)
            var = f"M_{stage.upper()}_{tag}"
            merge_ids[stage] = var
            combined[stage] = f"{prefix}_combined"       # basename (no .out) of the merged file
            sub.append('STAGE_T0=$(date +%s)')
            sub.append(f'ARRAY_ID=$(cd "{fdir}" && sbatch --parsable run_{aj}.sb)')
            sub.append(f'{var}=$(cd "{fdir}" && sbatch --parsable --dependency=afterok:$ARRAY_ID '
                       f'run_{mj}.sb "$STAGE_T0")')
            print(f"  {stage:<12} {npts} points -> {prefix}_combined.out")

        # ---- get_yk_new: same as build(), but fed the merged _combined.out files
        yk_dir = os.path.join(cfg.dest, "get_yk_new")
        rai = cfg.get("get_yk_new", f"r_ai_{tag}")
        if rai:
            rai_file = os.path.basename(rai)
            copy(cfg.rel(rai), os.path.join(yk_dir, rai_file))
            make_rai = ""
        else:
            if cfg.order == "NNLO":
                die(f"NNLO needs a real R_Ai table: set [get_yk_new] r_ai_{tag}")
            rai_file = f"{cfg.name}_dummy_R_Ai_{tag}.txt"           # all ones; exact for NLO
            write_rect_rai(os.path.join(yk_dir, rai_file), q_vals, y_vals, qt_vals,
                            int(float(ecm)), vtype, cfg.pdf)
            make_rai = ""
        yk_job, body = yk_script(cfg, tag, vtype, combined, total_pts, rai_file, make_rai)
        write(os.path.join(yk_dir, f"run_get_yk_new_{tag}.sb"), body, exe=True)
        sub.append(f'YK_{tag}=$(cd "{yk_dir}" && sbatch --parsable --dependency=afterok:'
                   f'${merge_ids["w_pert"]}:${merge_ids["w_asym"]}:${merge_ids["legacy_Y"]} '
                   f'--kill-on-invalid-dep=yes run_get_yk_new_{tag}.sb)')
        print(f"  get_yk_new   {yk_job}.out  ({cfg.order})")

        # ---- resbos: one job per [resbos] run, fed the merged legacy_main + the Yk grid above
        rdir = os.path.join(cfg.dest, "resbos")
        os.makedirs(os.path.join(rdir, "Resbos_grids"), exist_ok=True)
        runs = [r.strip() for r in cfg.get("resbos", "runs", "default").split(",") if r.strip()]
        for run in runs:
            rjob = write_resbos_run(cfg, rdir, tag, run, ("legacy", combined["legacy_main"]),
                                     ("get_yk_new", yk_job), tlines["resbos"])
            sub.append(f'(cd "{rdir}" && sbatch --parsable --dependency=afterok:$YK_{tag}:'
                       f'${merge_ids["legacy_main"]} --kill-on-invalid-dep=yes run_{rjob}.sb)')
            print(f"  resbos       {rjob}.in")
        sub.append("")

    sub.append('echo "all jobs submitted; watch with: squeue -u $USER"')
    write(os.path.join(cfg.dest, "submit_all.sh"), "\n".join(sub) + "\n", exe=True)
    shutil.copyfile(cfg.path, os.path.join(cfg.dest, "campaign.ini"))
    print(f"\nprepared {cfg.dest}\nnext: bash {os.path.join(cfg.dest, 'submit_all.sh')}")
    if args.submit:
        if not shutil.which("sbatch"):
            die("sbatch not found: use --submit on the HPCC login node")
        subprocess.check_call(["bash", os.path.join(cfg.dest, "submit_all.sh")])


# ------------------------------------------------------------------ main flow
def generate_grid_files(cfg, tlines):
    """[grids] generate: build a fine Q x y rectangle in-process (same math as
    make_grid_from_data.py's "generate" mode) and write it into the template's own
    legacy/w_pert/w_asym inp/, at whatever relative paths their .in files already name --
    so [templates] doesn't need a pre-baked grid checked in, just the .in structure.
    w_pert/w_asym are skipped when compute = NLO only (build() doesn't load their
    templates then -- see tlines' conditional keys in build()); legacy_dsi/legacy_asy
    don't need their own entry since their grid paths come from the same file as
    legacy_Y by default (patch_stage() only rewrites the LTO token, not the grid paths)."""
    import make_grid_from_data as mgfd
    from types import SimpleNamespace

    q_ns = SimpleNamespace(q_min=float(cfg.gen_q_min) if cfg.gen_q_min else None,
                            q_max=float(cfg.gen_q_max) if cfg.gen_q_max else None,
                            q_from=cfg.gen_q_from, q_col=cfg.gen_q_col)
    y_ns = SimpleNamespace(y_min=float(cfg.gen_y_min) if cfg.gen_y_min else None,
                            y_max=float(cfg.gen_y_max) if cfg.gen_y_max else None,
                            y_from=cfg.gen_y_from, y_col=cfg.gen_y_col)
    q_lo, q_hi = mgfd.resolve_range(q_ns, "q")
    y_lo, y_hi = mgfd.resolve_range(y_ns, "y")
    q_vals = [mgfd.fmt(v) for v in mgfd.spaced(q_lo, q_hi, cfg.gen_n_q, cfg.gen_q_spacing)]
    y_vals = [mgfd.fmt(v) for v in mgfd.spaced(y_lo, y_hi, cfg.gen_n_y, cfg.gen_y_spacing)]

    targets = [("legacy_Y", "legacy")]
    targets += [(k, k) for k in ("w_pert", "w_asym") if k in tlines]
    for key, folder in targets:
        q_rel, qt_rel, y_rel = grid_paths(tlines[key])
        write(os.path.join(cfg.src, folder, q_rel), "\n".join(q_vals) + "\n")
        write(os.path.join(cfg.src, folder, y_rel), "\n".join(y_vals) + "\n")

    print(f"[grids] generate: wrote {cfg.gen_n_q} Q ({cfg.gen_q_spacing}, [{mgfd.fmt(q_lo)}, "
          f"{mgfd.fmt(q_hi)}]) x {cfg.gen_n_y} y ({cfg.gen_y_spacing}, [{mgfd.fmt(y_lo)}, "
          f"{mgfd.fmt(y_hi)}]) into {cfg.src}'s {'/'.join(f for _, f in targets)} inp/")


def build(cfg, args):
    for d in ("get_yk_new", "legacy", "resbos", "w_asym", "w_pert"):
        if not os.path.isdir(os.path.join(cfg.src, d)):
            die(f"source folder {cfg.src} has no '{d}/'")
    if os.path.exists(cfg.dest) and not args.reuse:
        die(f"{cfg.dest} already exists (use --reuse to refresh its .in files and scripts)")
    if os.path.realpath(cfg.dest) == os.path.realpath(cfg.src):
        die("dest must differ from source")

    do_resnlo = "resNLO" in cfg.compute
    do_nlo = "NLO" in cfg.compute
    needed = set()
    if do_resnlo:
        needed |= set(RESNLO_STAGES)
    if do_nlo:
        needed |= set(NLO_STAGES)
    active_stages = [s for s in STAGES if s in needed]

    # ---- templates: legacy_dsi/legacy_asy (LTO=1/2, NLO only) default to the legacy_Y
    # template -- patch_stage() overwrites the LTO field regardless of what the template
    # started with, so no new template file is required unless [templates] overrides one.
    tpl = {"legacy_Y": cfg.req("templates", "legacy_Y"), "resbos": cfg.req("templates", "resbos")}
    if do_resnlo:
        tpl["legacy_main"] = cfg.req("templates", "legacy_main")
        tpl["w_pert"] = cfg.req("templates", "w_pert")
        tpl["w_asym"] = cfg.req("templates", "w_asym")
    if do_nlo:
        tpl["legacy_dsi"] = cfg.get("templates", "legacy_dsi", tpl["legacy_Y"])
        tpl["legacy_asy"] = cfg.get("templates", "legacy_asy", tpl["legacy_Y"])
    tlines = {k: read_lines(os.path.join(cfg.src, v)) for k, v in tpl.items()}

    if cfg.gen_grid:
        generate_grid_files(cfg, tlines)

    # ---- grids: every code that runs must read identical ones. With resNLO, w_pert/w_asym
    # are cross-checked against legacy_Y's grid; with NLO alone only Legacy itself runs, so
    # its own grid (shared automatically -- legacy_Y/legacy_main/legacy_dsi/legacy_asy all
    # read the same inp/ files) needs no cross-code check.
    ref = {}
    grid_stages = [("w_pert", "w_pert"), ("w_asym", "w_asym"), ("legacy", "legacy_Y")] if do_resnlo \
        else [("legacy", "legacy_Y")]
    for folder, key in grid_stages:
        for role, rel in zip(("Q", "qT", "y"), grid_paths(tlines[key])):
            p = os.path.join(cfg.src, folder, rel)
            if not os.path.isfile(p):
                die(f"grid {p} not found")
            ref.setdefault(role, (p, md5(p)))
            if md5(p) != ref[role][1]:
                die(f"{role} grid differs between {ref[role][0]} and {p}: all codes must use the same grid")
    counts = {r: sum(1 for l in open(p) if l.strip()) for r, (p, _) in ref.items()}
    full = (1, counts["qT"], 1, 1, counts["y"], 1, 1, counts["Q"], 1)
    act = tuple(int(t) for t in cfg.get("grids", "active", " ".join(map(str, full))).split())
    if len(act) != 9:
        die("[grids] active needs 9 integers")
    n_of = lambda a, b, s: len(range(a, b + 1, s))
    nqt, ny = n_of(*act[0:3]), n_of(*act[3:6])
    npts = nqt * ny * n_of(*act[6:9])
    print(f"compute: {', '.join(sorted(cfg.compute))}; grid: {counts['Q']} Q x {counts['qT']} qT x "
          f"{counts['y']} y; active {npts} points")

    # ---- copy only what a run needs
    exes = [("legacy", "main"), ("resbos", "resbos_root")]
    if do_resnlo:
        exes += [("w_pert", "w_pert"), ("w_asym", "w_asym"), ("get_yk_new", "get_yk_new")]
    missing = [cfg.exe_src(f, e) for f, e in exes if not os.path.isfile(cfg.exe_src(f, e))]
    if missing:                     # check all before copying anything: no half-made dest
        die("executable(s) not found:\n       " + "\n       ".join(missing) + "\n"
            "       build them first with:  sbatch setup_resbos_legacy.sb   (writes bin/<code>/<exe>)\n"
            "       or point to them in the .ini:  [executables] dir = <folder>  (or one key per code)")

    print(f"copying to {cfg.dest}")

    def copy_exe(folder, exe):
        src = cfg.exe_src(folder, exe)
        print(f"  {folder}/{exe} <- {src}")
        copy(src, os.path.join(cfg.dest, folder, exe), exe=True)

    for folder, exe in exes:
        copy_exe(folder, exe)
    copy(os.path.join(cfg.src, "legacy", grid_paths(tlines["legacy_Y"])[0]),
         os.path.join(cfg.dest, "legacy", grid_paths(tlines["legacy_Y"])[0]))
    copy(os.path.join(cfg.src, "legacy", grid_paths(tlines["legacy_Y"])[1]),
         os.path.join(cfg.dest, "legacy", grid_paths(tlines["legacy_Y"])[1]))
    copy(os.path.join(cfg.src, "legacy", grid_paths(tlines["legacy_Y"])[2]),
         os.path.join(cfg.dest, "legacy", grid_paths(tlines["legacy_Y"])[2]))
    if do_resnlo:
        for folder, key in (("w_pert", "w_pert"), ("w_asym", "w_asym")):
            for rel in grid_paths(tlines[key]):
                copy(os.path.join(cfg.src, folder, rel), os.path.join(cfg.dest, folder, rel))
        if cfg.order == "NLO":
            copy(os.path.join(cfg.src, "get_yk_new", "make_dummy_rai.py"),
                 os.path.join(cfg.dest, "get_yk_new", "make_dummy_rai.py"))
    for folder, exe, label, hdr in {(STAGES[s][0], STAGES[s][1], STAGES[s][2], STAGES[s][3]) for s in active_stages}:
        instantiate_shard_scripts(cfg, folder, label, exe, hdr, npts, args.clean_shards)

    # ---- per boson: .in files + submit lines. Stages whose merged output already exists in
    # dest and passes its own row-count check are reused (not resubmitted); only stages
    # actually (re)submitted this run feed into downstream --dependency= clauses, so a
    # dest folder that already has resNLO's outputs and now also asks for NLO only computes
    # the new legacy_dsi/legacy_asy/resbos/hadd pieces (legacy_Y is shared either way).
    sub = ["#!/bin/bash", "set -e", *timing_header_lines(cfg),
           "run_shards() {   # dir script job nshards -> prints the merge job id",
           '  local out; out=$(cd "$1" && bash "$2" "$3" "0-$(($4 - 1))") || return 1',
           '  echo "$out" >&2',
           "  echo \"$out\" | sed -n 's/^Submitted merge job \\([0-9]*\\).*/\\1/p'",
           "}", ""]
    ecm = cfg.ecm
    for vtype in cfg.procs:
        tag, jw = PROCS[vtype]
        print(f"\n== {vtype}")
        jobs, merge_ids, submitted = {}, {}, set()
        for stage in active_stages:
            folder, exe, label, hdr, k = STAGES[stage]
            lines = patch_stage(cfg, tlines, stage, tag, jw, vtype, ecm)
            job = (f"{cfg.name}_{stage}_{tag}" if stage in ("w_pert", "w_asym") else
                   f"{cfg.name}_legacy_{tag}_{LEGACY_SUFFIX[stage]}")
            set_active(lines, act, full)
            check_jobname(job)
            inpath = os.path.join(cfg.dest, folder, job + ".in")
            if os.path.isdir(os.path.join(cfg.dest, folder, job + "_shards")) and \
                    (not os.path.isfile(inpath) or open(inpath).read() != "".join(lines)):
                # submit_*_array.sh reuses existing shards, which would keep the OLD values
                die(f"{folder}/{job}_shards/ exists from an earlier input and {job}.in changed; "
                    f"remove that folder (or use a new dest) so the shards are regenerated")
            write(inpath, "".join(lines))
            jobs[stage] = job
            out_path = os.path.join(cfg.dest, folder, job + ".out")
            if valid_output(out_path, k, npts):
                print(f"  {stage:<12} {job}.out already exists ({k * npts} lines) -- reusing, not resubmitted")
                continue
            submitted.add(stage)
            n = cfg.shards(stage)
            var = f"M_{stage.upper()}_{tag}"
            merge_ids[stage] = var
            sub.append(f'{var}=$(run_shards "{os.path.join(cfg.dest, folder)}" submit_{label}_array.sh {job} {n})')
            print(f"  {stage:<12} {job}.in   {n} shard(s)")

        def dep_clause(dep_vars):
            """--dependency=... clause built only from stages actually (re)submitted this
            run (an unset bash var in --dependency=afterok:$VAR would break sbatch). Entries
            are stage names (resolved through merge_ids) or bare bash variable names
            (e.g. "YK_Wp") used as-is."""
            names = [merge_ids[s] if s in merge_ids else s for s in dep_vars]
            return (f'--dependency=afterok:{":".join("${" + n + "}" for n in names)} '
                    f'--kill-on-invalid-dep=yes ') if names else ""

        # ---- get_yk_new (resNLO only)
        yk_job = None
        if do_resnlo:
            yk_dir = os.path.join(cfg.dest, "get_yk_new")
            yk_job = f"{cfg.name}_Yk_{tag}_{cfg.order}"
            yk_out = os.path.join(yk_dir, yk_job + ".out")
            yk_upstream = [s for s in ("w_pert", "w_asym", "legacy_Y") if s in submitted]
            if not yk_upstream and valid_output(yk_out, 3, npts):
                print(f"  get_yk_new   {yk_job}.out already exists -- reusing, not resubmitted")
            else:
                rai = cfg.get("get_yk_new", f"r_ai_{tag}")
                if rai:
                    rai_file = os.path.basename(rai)
                    copy(cfg.rel(rai), os.path.join(yk_dir, rai_file))
                    make_rai = ""
                else:
                    if cfg.order == "NNLO":
                        die(f"NNLO needs a real R_Ai table: set [get_yk_new] r_ai_{tag}")
                    rai_file = f"{cfg.name}_dummy_R_Ai_{tag}.txt"           # all ones; exact for NLO
                    make_rai = (f'[ -f {rai_file} ] || python3 make_dummy_rai.py ../w_pert/{jobs["w_pert"]}.out '
                                f'{rai_file} {int(float(ecm))} {vtype} {cfg.pdf} || exit 1\n')
                _, body = yk_script(cfg, tag, vtype, jobs, npts, rai_file, make_rai)
                write(os.path.join(yk_dir, f"run_get_yk_new_{tag}.sb"), body, exe=True)
                sub.append(f'YK_{tag}=$(cd "{yk_dir}" && sbatch --parsable {dep_clause(yk_upstream)}'
                           f'run_get_yk_new_{tag}.sb)')
                print(f"  get_yk_new   {yk_job}.out  ({cfg.order})")
                submitted.add("get_yk_new")

        # ---- resbos (resNLO): one job per run, never sharded; its grids are staged by the job itself
        if do_resnlo:
            rdir = os.path.join(cfg.dest, "resbos")
            os.makedirs(os.path.join(rdir, "Resbos_grids"), exist_ok=True)
            runs = [r.strip() for r in cfg.get("resbos", "runs", "default").split(",") if r.strip()]
            for run in runs:
                rjob = write_resbos_run(cfg, rdir, tag, run, ("legacy", jobs["legacy_main"]),
                                         ("get_yk_new", yk_job), tlines["resbos"])
                root_out = os.path.join(cfg.dest, rjob + ".root")
                rdeps = [s for s in ("get_yk_new", "legacy_main") if s in submitted]
                if not rdeps and os.path.isfile(root_out):
                    print(f"  resbos       {rjob}.root already exists -- reusing, not resubmitted")
                    continue
                dc = dep_clause(["YK_" + tag if s == "get_yk_new" else s for s in rdeps])
                sub.append(f'(cd "{rdir}" && sbatch --parsable {dc}run_{rjob}.sb)')
                print(f"  resbos       {rjob}.in")
            sub.append("")

        # ---- NLO fixed order: Legacy LTO=2 (asy) + LTO=3 (Y, shared with resNLO) -> resbos,
        # and Legacy LTO=1 (dsi) alone -> resbos (Y piece = "-"), then hadd the two .root
        # files -- see README.md "NLO (fixed order)" for why this equals the real-emission
        # fixed-order cross section.
        if do_nlo:
            rdir = os.path.join(cfg.dest, "resbos")
            os.makedirs(os.path.join(rdir, "Resbos_grids"), exist_ok=True)
            runs = [r.strip() for r in cfg.get("resbos", "runs", "default").split(",") if r.strip()]
            for run in runs:
                asy_job = write_resbos_run(cfg, rdir, tag, run, ("legacy", jobs["legacy_asy"]),
                                            ("legacy", jobs["legacy_Y"]), tlines["resbos"], suffix="nloasy")
                dsi_job = write_resbos_run(cfg, rdir, tag, run, ("legacy", jobs["legacy_dsi"]),
                                            None, tlines["resbos"], suffix="nlodsi")
                asy_root = os.path.join(cfg.dest, asy_job + ".root")
                dsi_root = os.path.join(cfg.dest, dsi_job + ".root")
                asy_deps = [s for s in ("legacy_asy", "legacy_Y") if s in submitted]
                dsi_deps = [s for s in ("legacy_dsi",) if s in submitted]
                asy_reused = not asy_deps and os.path.isfile(asy_root)
                dsi_reused = not dsi_deps and os.path.isfile(dsi_root)
                if asy_reused:
                    print(f"  resbos       {asy_job}.root already exists -- reusing, not resubmitted")
                else:
                    sub.append(f'ASY_{tag}_{run}=$(cd "{rdir}" && sbatch --parsable {dep_clause(asy_deps)}'
                               f'run_{asy_job}.sb)')
                    print(f"  resbos       {asy_job}.in")
                if dsi_reused:
                    print(f"  resbos       {dsi_job}.root already exists -- reusing, not resubmitted")
                else:
                    sub.append(f'DSI_{tag}_{run}=$(cd "{rdir}" && sbatch --parsable {dep_clause(dsi_deps)}'
                               f'run_{dsi_job}.sb)')
                    print(f"  resbos       {dsi_job}.in")

                hjob = f"{cfg.name}_{tag}_{run}_NLO"
                hadd_out = os.path.join(cfg.dest, hjob + ".root")
                if asy_reused and dsi_reused and os.path.isfile(hadd_out):
                    print(f"  hadd         {hjob}.root already exists -- reusing, not resubmitted")
                else:
                    write(os.path.join(rdir, f"run_{hjob}.sb"),
                          hadd_script(cfg, hjob, f"{asy_job}.root", f"{dsi_job}.root", f"{hjob}.root"), exe=True)
                    hdeps = ([] if asy_reused else [f"${{ASY_{tag}_{run}}}"]) + \
                            ([] if dsi_reused else [f"${{DSI_{tag}_{run}}}"])
                    dc = f'--dependency=afterok:{":".join(hdeps)} --kill-on-invalid-dep=yes ' if hdeps else ""
                    sub.append(f'(cd "{rdir}" && sbatch --parsable {dc}run_{hjob}.sb)')
                    print(f"  hadd         {hjob}.root  ({asy_job}.root + {dsi_job}.root)")
            sub.append("")

    sub.append('echo "all jobs submitted; watch with: squeue -u $USER"')
    write(os.path.join(cfg.dest, "submit_all.sh"), "\n".join(sub) + "\n", exe=True)
    shutil.copyfile(cfg.path, os.path.join(cfg.dest, "campaign.ini"))
    print(f"\nprepared {cfg.dest}\nnext: bash {os.path.join(cfg.dest, 'submit_all.sh')}")
    if args.submit:
        if not shutil.which("sbatch"):
            die("sbatch not found: use --submit on the HPCC login node")
        subprocess.check_call(["bash", os.path.join(cfg.dest, "submit_all.sh")])


# ------------------------------------------------------------------ helper used inside jobs
def count_data_rows(path):
    """Rows after the first 'Q,qT,y' header line, or None if the file/header is missing --
    shared by cmd_check (job scripts' own row-count check) and valid_output (the driver's
    own reuse/skip check, [campaign] "don't repeat what's already there")."""
    if not os.path.isfile(path):
        return None
    with open(path, "rb") as f:
        while True:
            line = f.readline()
            if not line:
                return None
            if HDR_RE.match(line):
                break
        rows, last = 0, b"\n"
        while True:
            buf = f.read(CHUNK)
            if not buf:
                break
            rows += buf.count(b"\n")
            last = buf[-1:]
    return rows + (0 if last == b"\n" else 1)


def valid_output(path, per_point, npts):
    """True iff path exists, has a 'Q,qT,y' header, and its row count matches
    per_point*npts exactly -- used to decide whether a stage's output can be reused
    instead of resubmitted."""
    rows = count_data_rows(path)
    return rows is not None and rows == per_point * npts


def cmd_check(path, per_point, npts):
    """Exit 0 iff the data lines after the 'Q,qT,y' header line number per_point*npts."""
    rows = count_data_rows(path)
    if rows is None:
        print(f"{path}: file not found" if not os.path.isfile(path) else f"{path}: no 'Q,qT,y' header line")
        sys.exit(1)
    want = int(per_point) * int(npts)
    print(f"{path}: {rows} data lines, expected {want}")
    sys.exit(0 if rows == want else 1)


def cmd_merge_points(prefix, width, npts, header_lines, out_file):
    """Concatenate points/pt<i>/<prefix>_pt<i>.out (i=1..npts), stripping each one's own
    header_lines-line header and keeping only the first point's header -- combine.sh's job,
    run from run_process.py itself instead of a separate script."""
    width, npts, header_lines = int(width), int(npts), int(header_lines)
    header, body_chunks, total = None, [], 0
    for i in range(1, npts + 1):
        pt = str(i).zfill(width)
        path = os.path.join("points", f"pt{pt}", f"{prefix}_pt{pt}.out")
        if not os.path.isfile(path):
            print(f"{path}: file not found -- that point hasn't finished or failed")
            sys.exit(1)
        with open(path) as f:
            lines = f.readlines()
        if len(lines) <= header_lines:
            print(f"{path}: only {len(lines)} lines, expected a {header_lines}-line header plus data")
            sys.exit(1)
        if header is None:
            header = lines[:header_lines]
        body_chunks.append(lines[header_lines:])
        total += len(lines) - header_lines
    with open(out_file, "w") as f:
        f.writelines(header)
        for chunk in body_chunks:
            f.writelines(chunk)
    print(f"{out_file}: {header_lines} header lines + {total} data rows from {npts} points")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_check":
        return cmd_check(*sys.argv[2:5])
    if len(sys.argv) > 1 and sys.argv[1] == "_merge_points":
        return cmd_merge_points(*sys.argv[2:7])
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("campaign", help="campaign .ini (see 7TeV_WpWm_example.ini)")
    ap.add_argument("--submit", action="store_true", help="run submit_all.sh after preparing dest")
    ap.add_argument("--reuse", action="store_true", help="dest already exists: refresh .in files and scripts")
    ap.add_argument("--resbos-only", action="store_true",
                    help="dest already has finished w_pert/w_asym/legacy/get_yk_new outputs: only (re)generate "
                         "and, with --submit, run the resbos job(s) of [resbos] runs, reusing those outputs directly")
    ap.add_argument("--clean-shards", action="store_true",
                    help="w_pert/w_asym/legacy: delete each stage's <job>_shards/ folder from inside its merge job, "
                         "right after the row-count check passes (saves disk; kept on failure for debugging)")
    args = ap.parse_args()
    cfg = Cfg(args.campaign)
    if args.resbos_only:
        if cfg.experimental:
            die("--resbos-only doesn't apply in [grids] experimental mode "
                "(get_yk_new/resbos aren't run there at all)")
        return resbos_add(cfg, args)
    (build_points if cfg.experimental else build)(cfg, args)


if __name__ == "__main__":
    main()
