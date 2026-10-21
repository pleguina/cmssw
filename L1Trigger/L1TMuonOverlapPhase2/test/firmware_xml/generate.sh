#!/bin/bash
# Generates the OMTF firmware golden XML set from the GEN-SIM-DIGI-L1 samples and writes
# provenance.json next to it: CMSSW release and commit, configuration, inputs and outputs, each
# with its sha256, and per sample the cmsRun arguments, the input path and every
# simOmtfPhase2Digis parameter. <tag>/config_dump.py is the full configuration as run.
#
#   generate.sh <input-dir> <output-dir> [first-debug-event last-debug-event]
#
# <input-dir> holds SingleMuFlatPt2To100_Run4D110_GENSIMDIGIL1_<tag>.root for the tags below.
# The debug events (default 1 100) are dumped as DetailedSingleEv XML for the pos and neg samples.
# Run from a CMSSW area set up with cmsenv; refuses to run if tracked files under src/ differ
# from the commit, so that the recorded commit is what produced the files.
set -euo pipefail

TAGS="pos neg extra1 extra2 extra3 pos_norpc"
DEBUG_TAGS="pos neg"
IN="$(cd "${1:?input dir}" && pwd)"
OUT="$(mkdir -p "${2:?output dir}" && cd "$2" && pwd)"
FIRST="${3:-1}"
LAST="${4:-100}"
HERE="$(cd "$(dirname "$0")" && pwd)"
CFG="$HERE/firmware_xml_cfg.py"

: "${CMSSW_BASE:?run cmsenv first}"
SRC="$CMSSW_BASE/src"
if [[ -n "$(git -C "$SRC" status --porcelain --untracked-files=no)" ]]; then
    echo "ERROR: $SRC has uncommitted changes:" >&2
    git -C "$SRC" status --short --untracked-files=no >&2
    exit 1
fi
# Conditions need a site configuration; on a machine without one (SITECONF/local empty) use CERN's.
if [[ ! -f "${SITECONFIG_PATH:-}/JobConfig/site-local-config.xml" ]]; then
    export SITECONFIG_PATH=/cvmfs/cms.cern.ch/SITECONF/T2_CH_CERN
fi

debug_events="$(seq -s, "$FIRST" "$LAST")"
pids=()
for tag in $TAGS; do
    input="$IN/SingleMuFlatPt2To100_Run4D110_GENSIMDIGIL1_$tag.root"
    [[ -f "$input" ]] || { echo "ERROR: missing $input" >&2; exit 1; }
    args=("inputFiles=file:$input" "outputName=TestEvents_swnew_p1scale_ntest9_real_$tag")
    [[ "$tag" == *norpc ]] && args+=("noRPC=1")
    [[ " $DEBUG_TAGS " == *" $tag "* ]] && args+=("debugEvents=$debug_events")
    mkdir -p "$OUT/$tag"
    printf '%s\n' "${args[@]}" > "$OUT/$tag/cmsRun.args"
    edmConfigDump "$CFG" "${args[@]}" > "$OUT/$tag/config_dump.py"
    (cd "$OUT/$tag" && cmsRun "$CFG" "${args[@]}" > cmsRun.log 2>&1) &
    pids+=($!)
done
rc=0
for pid in "${pids[@]}"; do wait "$pid" || rc=1; done
if [[ $rc -ne 0 ]]; then
    echo "ERROR: cmsRun failed, see $OUT/<tag>/cmsRun.log" >&2
    exit 1
fi

python3 - "$OUT" "$SRC" "$CFG" "$IN" "$TAGS" "$FIRST" "$LAST" <<'PY'
import hashlib, json, os, platform, re, subprocess, sys, time
from pathlib import Path

out, src, cfg, inp, tags, first, last = sys.argv[1:]
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
def git(*a):
    return subprocess.run(["git", "-C", src, *a], capture_output=True, text=True, check=True).stdout.strip()

files = []
for p in sorted(Path(out).glob("*/*.xml")):
    text = p.read_text(errors="replace")
    files.append({
        "path": str(p.relative_to(out)),
        "sha256": sha(p),
        "events": len(re.findall(r"<Event ", text)),
        "has_processor_data": "<Processor" in text,
    })
def emulator_parameters(dump):
    import FWCore.ParameterSet.Config as cms
    scope = {}
    exec(compile(Path(dump).read_text(), str(dump), "exec"), scope)
    module = scope["process"].simOmtfPhase2Digis
    def plain(v):
        v = v.value() if hasattr(v, "value") else v
        if isinstance(v, (list, tuple)):
            return [plain(x) for x in v]
        return v if isinstance(v, (bool, int, float, str)) or v is None else str(v)
    return {name: plain(getattr(module, name)) for name in sorted(module.parameterNames_())}

samples = {}
for t in tags.split():
    d = Path(out, t)
    samples[t] = {
        "input": str(Path(inp) / f"SingleMuFlatPt2To100_Run4D110_GENSIMDIGIL1_{t}.root"),
        "cmsrun_args": d.joinpath("cmsRun.args").read_text().split("\n")[:-1],
        "config_dump": {"file": f"{t}/config_dump.py", "sha256": sha(d / "config_dump.py")},
        "simOmtfPhase2Digis": emulator_parameters(d / "config_dump.py"),
    }
omtf_config = Path(src) / "L1Trigger/L1TMuon/data/omtf_config"
prov = {
    "schema_version": 1,
    "cmssw_release": os.environ.get("CMSSW_VERSION"),
    "scram_arch": os.environ.get("SCRAM_ARCH"),
    "cmssw_commit": git("rev-parse", "HEAD"),
    "cmssw_branch": git("rev-parse", "--abbrev-ref", "HEAD"),
    "cmssw_remote": git("remote", "get-url", git("remote").split()[0]) if git("remote") else "",
    "config": {"file": os.path.relpath(cfg, src), "sha256": sha(cfg)},
    "omtf_config": {p.name: sha(p) for p in sorted(omtf_config.glob("*")) if p.is_file()},
    "inputs": {f"SingleMuFlatPt2To100_Run4D110_GENSIMDIGIL1_{t}.root":
               sha(Path(inp) / f"SingleMuFlatPt2To100_Run4D110_GENSIMDIGIL1_{t}.root") for t in tags.split()},
    "debug_events": [int(first), int(last)],
    "samples": samples,
    "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "host": platform.node(),
    "files": files,
}
Path(out, "provenance.json").write_text(json.dumps(prov, indent=1) + "\n")
print(f"provenance.json: {len(files)} files, commit {prov['cmssw_commit'][:12]}")
PY
