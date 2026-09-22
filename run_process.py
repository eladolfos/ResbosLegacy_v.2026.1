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
PROCS = {"W+": ("Wp", 1), "W-": ("Wm", -1)}          # Type_V -> (file tag, JWTYPE)
# stage: (folder, executable, label used in the shard scripts, header lines, lines/point)
STAGES = {
    "w_pert":      ("w_pert", "w_pert", "w_pert", 3, 1),
    "w_asym":      ("w_asym", "w_asym", "w_asym", 3, 1),
    "legacy_Y":    ("legacy", "main", "legacy", 15, 3),
    "legacy_main": ("legacy", "main", "legacy", 15, 1),
}
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
    check = (f'\n# row-count check added by run_process.py\n'
             f'case "$JOBNAME" in *_legacy_*_Y) K=3;; *) K=1;; esac\n'
             f'python3 "{SELF}" _check "${{JOBNAME}}.out" $K {npts} || exit 1\n')
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
exit $EXIT_CODE
"""
    return job, body


def resbos_script(cfg, rjob, main_job, yk_job):
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
# the two grids resbos reads: Legacy LTO=0 ("W321 main grid") and the Yk grid from get_yk_new
mkdir -p Resbos_grids
stage ../legacy/{main_job}.out Resbos_grids/{main_job}.out
stage ../get_yk_new/{yk_job}.out Resbos_grids/{yk_job}.out

LOG_FILE="resbos_output_{rjob}_${{SLURM_JOB_ID}}.log"
srun ./resbos_root {rjob} >> "$LOG_FILE" 2>&1
EXIT_CODE=$?
echo "resbos_root finished with exit code $EXIT_CODE" | tee -a "$LOG_FILE"
exit $EXIT_CODE
"""


def write_resbos_run(cfg, rdir, tag, run, main_job, yk_job, resbos_lines):
    """Write resbos/<name>_<tag>_<run>.in + its .sb from the resbos template + [cuts.<run>]; return the job name."""
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
    put(lines, find(lines, "Main data grid", "main grid"), f'./Resbos_grids/{main_job}.out')
    put(lines, find(lines, "Y piece grid", "Y grid"), f"./Resbos_grids/{yk_job}.out")
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
    rjob = f"{cfg.name}_{tag}_{run}"
    write(os.path.join(rdir, rjob + ".in"), "".join(lines))
    write(os.path.join(rdir, f"run_{rjob}.sb"), resbos_script(cfg, rjob, main_job, yk_job), exe=True)
    return rjob


def resbos_add(cfg, args):
    """--resbos-only: dest already has the w_pert/w_asym/legacy/get_yk_new outputs of an earlier run;
    only (re)generate and optionally submit the resbos job(s) for [resbos] runs, reusing them directly."""
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
            rjob = write_resbos_run(cfg, rdir, tag, run, main_job, yk_job, resbos_lines)
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


# ------------------------------------------------------------------ main flow
def build(cfg, args):
    for d in ("get_yk_new", "legacy", "resbos", "w_asym", "w_pert"):
        if not os.path.isdir(os.path.join(cfg.src, d)):
            die(f"source folder {cfg.src} has no '{d}/'")
    if os.path.exists(cfg.dest) and not args.reuse:
        die(f"{cfg.dest} already exists (use --reuse to refresh its .in files and scripts)")
    if os.path.realpath(cfg.dest) == os.path.realpath(cfg.src):
        die("dest must differ from source")

    tpl = {k: cfg.req("templates", k) for k in ("legacy_Y", "legacy_main", "w_pert", "w_asym", "resbos")}
    tlines = {k: read_lines(os.path.join(cfg.src, v)) for k, v in tpl.items()}

    # ---- grids: every code must read identical ones
    ref = {}
    for folder, key in (("w_pert", "w_pert"), ("w_asym", "w_asym"), ("legacy", "legacy_Y")):
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
    print(f"grid: {counts['Q']} Q x {counts['qT']} qT x {counts['y']} y; active {npts} points")

    # ---- copy only what a run needs
    exes = [(f, STAGES[k][1]) for f, k in (("w_pert", "w_pert"), ("w_asym", "w_asym"), ("legacy", "legacy_Y"))]
    exes += [("get_yk_new", "get_yk_new"), ("resbos", "resbos_root")]
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

    for folder, key in (("w_pert", "w_pert"), ("w_asym", "w_asym"), ("legacy", "legacy_Y")):
        copy_exe(folder, STAGES[key][1])
        for rel in grid_paths(tlines[key]):
            copy(os.path.join(cfg.src, folder, rel), os.path.join(cfg.dest, folder, rel))
    copy_exe("get_yk_new", "get_yk_new")
    copy_exe("resbos", "resbos_root")
    if cfg.order == "NLO":
        copy(os.path.join(cfg.src, "get_yk_new", "make_dummy_rai.py"),
             os.path.join(cfg.dest, "get_yk_new", "make_dummy_rai.py"))
    for folder, exe, label, hdr in {(f, e, l, h) for f, e, l, h, _ in STAGES.values()}:
        instantiate_shard_scripts(cfg, folder, label, exe, hdr, npts, args.clean_shards)

    # ---- per boson: .in files + submit lines
    sub = ["#!/bin/bash", "set -e",
           "run_shards() {   # dir script job nshards -> prints the merge job id",
           '  local out; out=$(cd "$1" && bash "$2" "$3" "0-$(($4 - 1))") || return 1',
           '  echo "$out" >&2',
           "  echo \"$out\" | sed -n 's/^Submitted merge job \\([0-9]*\\).*/\\1/p'",
           "}", ""]
    ecm = cfg.ecm
    for vtype in cfg.procs:
        tag, jw = PROCS[vtype]
        print(f"\n== {vtype}")
        jobs, merge_ids = {}, {}
        for stage, (folder, exe, label, hdr, k) in STAGES.items():
            lines = list(tlines[stage])
            if stage in ("w_pert", "w_asym"):
                set_token(lines, find(lines, "ECM,iBeam", "ECM"), 0, ecm)
                set_token(lines, find(lines, "JWTYPE", "JWTYPE"), 0, str(jw))
                put(lines, find(lines, "PDF file name", "PDF"), "lha_" + cfg.pdf)
                job = f"{cfg.name}_{stage}_{tag}"
            else:
                l2 = find(lines, "ECM,LTO", "ECM,LTO")
                set_token(lines, l2, 0, ecm if "." in ecm else ecm + ".0")
                set_token(lines, l2, 1, "3" if stage == "legacy_Y" else "0")
                put(lines, find(lines, "Type_V", "Type_V"), vtype)
                put(lines, find(lines, "Evolved PDF file", "PDF"), "lha_" + cfg.pdf)
                if cfg.get("legacy", "bmax"):
                    put(lines, find(lines, "bMax", "bMax"), cfg.get("legacy", "bmax"))
                if cfg.get("legacy", "nonpert"):
                    put(lines, find(lines, "g1, g2, g3, Q0, nG", "nonpert"), cfg.get("legacy", "nonpert"))
                job = f"{cfg.name}_legacy_{tag}_{'Y' if stage == 'legacy_Y' else 'main'}"
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
            n = cfg.shards(stage)
            var = f"M_{stage.upper()}_{tag}"
            merge_ids[stage] = var
            sub.append(f'{var}=$(run_shards "{os.path.join(cfg.dest, folder)}" submit_{label}_array.sh {job} {n})')
            print(f"  {stage:<12} {job}.in   {n} shard(s)")

        # ---- get_yk_new
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
            make_rai = (f'[ -f {rai_file} ] || python3 make_dummy_rai.py ../w_pert/{jobs["w_pert"]}.out '
                        f'{rai_file} {int(float(ecm))} {vtype} {cfg.pdf} || exit 1\n')
        yk_job, body = yk_script(cfg, tag, vtype, jobs, npts, rai_file, make_rai)
        write(os.path.join(yk_dir, f"run_get_yk_new_{tag}.sb"), body, exe=True)
        sub.append(f'YK_{tag}=$(cd "{yk_dir}" && sbatch --parsable --dependency=afterok:'
                   f'${merge_ids["w_pert"]}:${merge_ids["w_asym"]}:${merge_ids["legacy_Y"]} '
                   f'--kill-on-invalid-dep=yes run_get_yk_new_{tag}.sb)')
        print(f"  get_yk_new   {yk_job}.out  ({cfg.order})")

        # ---- resbos: one job per run, never sharded; its grids are staged by the job itself
        rdir = os.path.join(cfg.dest, "resbos")
        os.makedirs(os.path.join(rdir, "Resbos_grids"), exist_ok=True)
        runs = [r.strip() for r in cfg.get("resbos", "runs", "default").split(",") if r.strip()]
        for run in runs:
            rjob = write_resbos_run(cfg, rdir, tag, run, jobs["legacy_main"], yk_job, tlines["resbos"])
            sub.append(f'(cd "{rdir}" && sbatch --parsable --dependency=afterok:$YK_{tag}:${merge_ids["legacy_main"]} '
                       f'--kill-on-invalid-dep=yes run_{rjob}.sb)')
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


# ------------------------------------------------------------------ helper used inside jobs
def cmd_check(path, per_point, npts):
    """Exit 0 iff the data lines after the 'Q,qT,y' header line number per_point*npts."""
    if not os.path.isfile(path):
        print(f"{path}: file not found")
        sys.exit(1)
    with open(path, "rb") as f:
        while True:
            line = f.readline()
            if not line:
                print(f"{path}: no 'Q,qT,y' header line")
                sys.exit(1)
            if HDR_RE.match(line):
                break
        rows, last = 0, b"\n"
        while True:
            buf = f.read(CHUNK)
            if not buf:
                break
            rows += buf.count(b"\n")
            last = buf[-1:]
    rows += 0 if last == b"\n" else 1
    want = int(per_point) * int(npts)
    print(f"{path}: {rows} data lines, expected {want}")
    sys.exit(0 if rows == want else 1)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_check":
        return cmd_check(*sys.argv[2:5])
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
    (resbos_add if args.resbos_only else build)(cfg, args)


if __name__ == "__main__":
    main()
