#!/usr/bin/env python3
import os
import json
import csv
import argparse
from pathlib import Path
from typing import Dict, List, Tuple

# ============================================================
# Helpers
'''
ls toy*/*sh | xargs -I {} -P 2 bash -c "{}"


sample_id	hap1	hap2	bam
toy1	/CDATA/pangenome/flagger_test/WDL/toy/KPPD019.hap1.fa	/CDATA/pangenome/flagger_test/WDL/toy/KPPD019.hap2.fa	/CDATA/pangenome/flagger_test/toy/test/tets.sorted.bam
toy2	/CDATA/pangenome/flagger_test/WDL/toy/KPPD019.hap1.fa	/CDATA/pangenome/flagger_test/WDL/toy/KPPD019.hap2.fa	/CDATA/pangenome/flagger_test/toy/test/tets.sorted.bam
'''
'''
  python3 mk.sh.py \
  --sample-sheet sample.tsv \
  --cromwell-jar /CDATA/pangenome/flagger_test/cromwell-85.jar \
  --wdl /CDATA/pangenome/flagger_test/flagger-1.1.0/wdls/workflows/hmm_flagger_end_to_end.wdl \
  --out-root /CDATA/pangenome/flagger_test/toy/wdl_test \
  --max-threads 128 \
  --max-mem-gb 256 \
  --ref-dir /CDATA/pangenome/flagger_test/test/ref \
  --projection-ref-fasta /CDATA/pangenome/flagger_test/test/ref/chm13v2.0_maskedY_rCRS.fa \
  --alpha-tsv /CDATA/pangenome/flagger_test/test/ref/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv \
  --enable-call-cache
'''
# ============================================================

def die(msg: str, code: int = 2) -> None:
    raise SystemExit(f"[ERROR] {msg}")

def read_table(path: str) -> List[Dict[str, str]]:
    """
    Read TSV/CSV into list of dict rows.
    Auto-detect delimiter by file extension (.tsv/.csv) but also
    accepts tabs in csv if user provides.
    Required headers:
      sample_id, hap1, hap2, bam
    """
    p = Path(path)
    if not p.exists():
        die(f"Sample sheet not found: {path}")

    # Detect delimiter
    suffix = p.suffix.lower()
    if suffix == ".tsv":
        delim = "\t"
    elif suffix == ".csv":
        delim = ","
    else:
        # fallback: try tab first, then comma
        delim = "\t"

    rows: List[Dict[str, str]] = []
    with p.open(newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if reader.fieldnames is None:
            die(f"Could not read header from sample sheet: {path}")

        # If delimiter guess was wrong, try comma fallback (common)
        if set(reader.fieldnames) == {p.name} or len(reader.fieldnames) == 1:
            f.seek(0)
            reader = csv.DictReader(f, delimiter=",")
            if reader.fieldnames is None:
                die(f"Could not read header from sample sheet: {path}")

        for r in reader:
            # Normalize keys: strip spaces
            rr = {k.strip(): (v.strip() if isinstance(v, str) else v) for k, v in r.items()}
            # skip empty lines
            if not any(rr.values()):
                continue
            rows.append(rr)

    return rows

def validate_rows(rows: List[Dict[str, str]]) -> None:
    required = ["sample_id", "hap1", "hap2", "bam"]
    if not rows:
        die("Sample sheet is empty (no rows).")

    missing_cols = [c for c in required if c not in rows[0]]
    if missing_cols:
        die(f"Sample sheet missing required columns: {missing_cols}\n"
            f"Required: {required}\n"
            f"Found: {list(rows[0].keys())}")

    # Validate each row
    for i, r in enumerate(rows, start=2):  # header is line 1
        for c in required:
            if not r.get(c):
                die(f"Missing value in row {i} for column '{c}'")

def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

# ============================================================
# Build JSON exactly like input.winnowmap.ori.json (all fields)
# ============================================================

def build_input_json_full(
    sample_id: str,
    hap1_fa: str,
    hap2_fa: str,
    bam: str,
    ref_dir: str,
    projection_ref_fasta: str,
    alpha_tsv: str,
    max_threads: int,
    max_mem_gb: int,
    enable_decompose_cntr: bool = False,
    enable_bigwig: bool = False,
    enable_secphase: bool = True,
    secphase_docker: str = "mobinasri/secphase:v0.4.4",
    secphase_options: str = "--hifi",
    secphase_version: str = "v0.4.4",
) -> Dict:

    bai = bam + ".bai"

    # Paths derived from ref_dir (same as your example json)
    sex_bed = f"{ref_dir}/stratifications/sex/chm13v2.0_sex.bed"
    sd_bed = f"{ref_dir}/stratifications/sd/chm13v2.0_SD.all.bed"
    cntr_bed = f"{ref_dir}/stratifications/censat/chm13v2.0_no_ct.bed"
    cntr_ct_bed = f"{ref_dir}/stratifications/censat/chm13v2.0_only_ct.bed"

    bias_beds = [
        f"{ref_dir}/potential_biases/chm13v2.0_bsat.bed",
        f"{ref_dir}/potential_biases/chm13v2.0_hor.bed",
        f"{ref_dir}/potential_biases/chm13v2.0_hsat1A.bed",
        f"{ref_dir}/potential_biases/chm13v2.0_hsat1B.bed",
        f"{ref_dir}/potential_biases/chm13v2.0_hsat2.bed",
        f"{ref_dir}/potential_biases/chm13v2.0_hsat3.bed",
    ]

    # suffix/trackName with sampleID prefix
    suffix_base = "hmm_flagger_v1.1.0"
    track_base = "hmm_flagger_v1.1.0"
    suffix = f"{sample_id}.{suffix_base}"
    track_name = f"{sample_id}.{track_base}"

    # IMPORTANT: keep types consistent with WDL expectations.
    # Your example used strings for many numeric fields; Cromwell accepts numbers too.
    # We'll use numbers for numeric fields (more standard), but if you want exact
    # string-typed JSON like your example, set stringify=True (not requested).
    return {
        # --- sample naming ---
        "HMMFlaggerEndToEnd.sampleName": sample_id,
        "HMMFlaggerEndToEnd.suffix": suffix,
        "HMMFlaggerEndToEnd.trackName": track_name,

        # --- assemblies ---
        "HMMFlaggerEndToEnd.hap1AssemblyFasta": hap1_fa,
        "HMMFlaggerEndToEnd.hap2AssemblyFasta": hap2_fa,

        # --- input alignments ---
        "HMMFlaggerEndToEnd.readAlignmentBam": bam,
        "HMMFlaggerEndToEnd.readAlignmentBai": bai,

        # --- references ---
        "HMMFlaggerEndToEnd.alphaTsv": alpha_tsv,
        "HMMFlaggerEndToEnd.projectionReferenceFasta": projection_ref_fasta,

        # --- projection beds ---
        "HMMFlaggerEndToEnd.sexBedToBeProjected": sex_bed,
        "HMMFlaggerEndToEnd.SDBedToBeProjected": sd_bed,
        "HMMFlaggerEndToEnd.cntrBedToBeProjected": cntr_bed,
        "HMMFlaggerEndToEnd.cntrCtBedToBeProjected": cntr_ct_bed,

        # --- bias beds ---
        "HMMFlaggerEndToEnd.biasAnnotationsBedArrayToBeProjected": bias_beds,

        # --- core model params (from your example) ---
        "HMMFlaggerEndToEnd.windowLen": 4000,
        "HMMFlaggerEndToEnd.chunkLen": 20000000,
        "HMMFlaggerEndToEnd.modelType": "trunc_exp_gaussian",
        "HMMFlaggerEndToEnd.numberOfIterations": 100,
        "HMMFlaggerEndToEnd.convergenceTolerance": 0.001,
        "HMMFlaggerEndToEnd.downSamplingRate": 1.0,

        # --- filters (from your example) ---
        "HMMFlaggerEndToEnd.minReadLength": 5000,
        "HMMFlaggerEndToEnd.minAlignmentLength": 5000,
        "HMMFlaggerEndToEnd.bam2cov.minAlignmentLength": 5000,
        "HMMFlaggerEndToEnd.maxReadDivergence": 0.1,
        "HMMFlaggerEndToEnd.minHighMapqRatio": 0.75,
        "HMMFlaggerEndToEnd.maxHighMapqRatio": 0.25,

        # --- labels ---
        "HMMFlaggerEndToEnd.labelNames": "Err,Dup,Hap,Col",

        # --- output toggles ---
        "HMMFlaggerEndToEnd.enableDecomposingCntrBed": bool(enable_decompose_cntr),
        "HMMFlaggerEndToEnd.enableOutputtingBigWig": bool(enable_bigwig),

        # --- task resources / docker / threads (from your example) ---
        "HMMFlaggerEndToEnd.produceFai.dockerImage": "mobinasri/bio_base:v0.4.0",
        "HMMFlaggerEndToEnd.produceFai.memSize": max_mem_gb,
        "HMMFlaggerEndToEnd.produceFai.threadCount": max_threads,

        "HMMFlaggerEndToEnd.createDipAsm.dockerImage": "mobinasri/bio_base:v0.4.0",
        "HMMFlaggerEndToEnd.createDipAsm.memSize": max_mem_gb,
        "HMMFlaggerEndToEnd.createDipAsm.threadCount": max_threads,

        "HMMFlaggerEndToEnd.collectAnnotations.memSize": max_mem_gb,
        "HMMFlaggerEndToEnd.collectAnnotations.threadCount": max_threads,

        "HMMFlaggerEndToEnd.augmentCoverageByLabels.threadCount": 8,

        "HMMFlaggerEndToEnd.cov2bigwig.memSize": 32,
        "HMMFlaggerEndToEnd.cov2bigwig.threadCount": 8,

        # Flagger core resources (you explicitly want CLI control)
        "HMMFlaggerEndToEnd.flaggerMemSize": max_mem_gb,
        "HMMFlaggerEndToEnd.flaggerThreadCount": max_threads,

        "HMMFlaggerEndToEnd.labelPrediction.dockerImage": "mobinasri/flagger:v1.1.0",
        "HMMFlaggerEndToEnd.labelPrediction.memSize": max_mem_gb,
        "HMMFlaggerEndToEnd.labelPrediction.threadCount": max_threads,

        "HMMFlaggerEndToEnd.labelTruth.dockerImage": "mobinasri/flagger:v1.1.0",
        "HMMFlaggerEndToEnd.labelTruth.memSize": max_mem_gb,
        "HMMFlaggerEndToEnd.labelTruth.threadCount": max_threads,

        "HMMFlaggerEndToEnd.getFinalBed.memSize": max_mem_gb,
        "HMMFlaggerEndToEnd.getFinalBed.threadCount": max_threads,

        "HMMFlaggerEndToEnd.makeSummaryTable.memSize": 32,
        "HMMFlaggerEndToEnd.makeSummaryTable.threadCount": 8,

        # hap-to-ref split flags (as in your example)
        "HMMFlaggerEndToEnd.hap1ToRef.splitAssembly": False,
        "HMMFlaggerEndToEnd.hap2ToRef.splitAssembly": False,

        # Secphase block (as in your example)
        "HMMFlaggerEndToEnd.enableRunningSecphase": bool(enable_secphase),
        "HMMFlaggerEndToEnd.secphaseDockerImage": secphase_docker,
        "HMMFlaggerEndToEnd.secphaseOptions": secphase_options,
        "HMMFlaggerEndToEnd.secphaseVersion": secphase_version,
    }

# ============================================================
# Cromwell config & run script
# ============================================================

def write_conf(conf_path: Path, workflow_root: Path, enable_call_cache: bool = True) -> None:
    conf = f"""
backend {{
  default = "Local"
  providers {{
    Local {{
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {{
        root = "{workflow_root}"
      }}
    }}
  }}
}}

call-caching {{
  enabled = {"true" if enable_call_cache else "false"}
}}
"""
    ensure_parent(conf_path)
    conf_path.write_text(conf.strip() + "\n")

def write_run_sh(
    sh_path: Path,
    cromwell_jar: str,
    wdl_path: str,
    conf_path: Path,
    workflow_root: Path,
    input_json_path: Path,
    metadata_json_path: Path,
) -> None:
    # Cromwell v85: config file is passed via -Dconfig.file or --config only works for some older patterns.
    script = f"""#!/usr/bin/env bash
set -euo pipefail

java -Dconfig.file="{conf_path}" \\
  -jar "{cromwell_jar}" run \\
  --workflow-root "{workflow_root}" \\
  "{wdl_path}" \\
  -i "{input_json_path}" \\
  -m "{metadata_json_path}"
"""
    #> run.log 2>&1 &
    ensure_parent(sh_path)
    sh_path.write_text(script)
    os.chmod(sh_path, 0o755)

# ============================================================
# Main
# ============================================================

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Prepare per-sample Cromwell config + full HMMFlaggerEndToEnd input JSON from TSV/CSV."
    )

    ap.add_argument("--sample-sheet", required=True, help="TSV/CSV with columns: sample_id, hap1, hap2, bam")
    ap.add_argument("--out-root", required=True, help="Output root directory (will create per-sample dirs)")
    ap.add_argument("--cromwell-jar", required=True, help="Path to cromwell-85.jar")
    ap.add_argument("--wdl", required=True, help="Path to hmm_flagger_end_to_end.wdl")

    # You requested these become arguments
    ap.add_argument("--ref-dir", required=True, help="Reference directory containing stratifications/ and potential_biases/")
    ap.add_argument("--projection-ref-fasta", required=True, help="Projection reference fasta (e.g., chm13v2.0_maskedY_rCRS.fa)")
    ap.add_argument("--alpha-tsv", required=True, help="Alpha TSV path")

    # You explicitly requested these be added
    ap.add_argument("--max-threads", type=int, required=True, help="Max threads to use for Flagger tasks")
    ap.add_argument("--max-mem-gb", type=int, required=True, help="Max memory (GB) to use for Flagger tasks")

    # Optional toggles (default match your example json)
    ap.add_argument("--enable-call-cache", action="store_true", help="Enable call-caching in conf (default: off unless set)")
    ap.add_argument("--disable-secphase", action="store_true", help="If set, write enableRunningSecphase=false")
    ap.add_argument("--enable-bigwig", action="store_true", help="If set, enableOutputtingBigWig=true")
    ap.add_argument("--enable-decompose-cntr", action="store_true", help="If set, enableDecomposingCntrBed=true")

    args = ap.parse_args()

    rows = read_table(args.sample_sheet)
    validate_rows(rows)

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    for r in rows:
        sample_id = r["sample_id"]
        hap1 = r["hap1"]
        hap2 = r["hap2"]
        bam = r["bam"]

        sample_dir = out_root / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        workflow_root = sample_dir / "cromwell_root"
        workflow_root.mkdir(exist_ok=True)

        conf_path = sample_dir / f"{sample_id}.conf"
        input_json_path = sample_dir / f"{sample_id}.input.winnowmap.hmm_flagger_end_to_end_WDL.json"
        meta_json_path = sample_dir / f"{sample_id}.output.winnowmap.hmm_flagger_end_to_end_WDL.json"
        run_sh_path = sample_dir / "run.wdl.sh"

        # conf (per-sample root to avoid collisions)
        write_conf(conf_path, workflow_root, enable_call_cache=bool(args.enable_call_cache))

        # input json (FULL keys from your ori json)
        inp = build_input_json_full(
            sample_id=sample_id,
            hap1_fa=hap1,
            hap2_fa=hap2,
            bam=bam,
            ref_dir=args.ref_dir,
            projection_ref_fasta=args.projection_ref_fasta,
            alpha_tsv=args.alpha_tsv,
            max_threads=args.max_threads,
            max_mem_gb=args.max_mem_gb,
            enable_decompose_cntr=bool(args.enable_decompose_cntr),
            enable_bigwig=bool(args.enable_bigwig),
            enable_secphase=(not args.disable_secphase),
        )
        input_json_path.write_text(json.dumps(inp, indent=2) + "\n")

        # run script
        write_run_sh(
            sh_path=run_sh_path,
            cromwell_jar=args.cromwell_jar,
            wdl_path=args.wdl,
            conf_path=conf_path,
            workflow_root=workflow_root,
            input_json_path=input_json_path,
            metadata_json_path=meta_json_path,
        )

        print(f"[OK] {sample_id}")
        print(f"  {conf_path}")
        print(f"  {input_json_path}")
        print(f"  {run_sh_path}")
        print(f"  Run: cd {sample_dir} && ./run.wdl.sh")
        print()

if __name__ == "__main__":
    main()
