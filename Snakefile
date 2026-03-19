"""
HOTAIR nanopore m6A analysis pipeline

Workflow:
  1. WarpDemuX demultiplexing (pod5 → per-sample pod5)
  2. Dorado basecalling with modification calling (all mods)
  3. minimap2 alignment to GENCODE transcriptome
  4. modkit modification pileup and per-read extraction

Run from aa-tRNA-seq-pipeline directory:
  pixi run snakemake -s ../Snakefile --configfile=../config/config-hotair.yml -n
  pixi run snakemake -s ../Snakefile --configfile=../config/config-hotair.yml --cores 8
"""

import os
import sys
import glob
import yaml

# Project root (where this Snakefile lives)
PROJECT_DIR = os.path.dirname(workflow.snakefile)

# Pipeline directory (submodule)
PIPELINE_DIR = os.path.join(PROJECT_DIR, "aa-tRNA-seq-pipeline")
SCRIPT_DIR = os.path.join(PIPELINE_DIR, "workflow", "scripts")

# Include base config defaults
configfile: os.path.join(PIPELINE_DIR, "config", "config-base.yml")

outdir = config["output_directory"]
if not os.path.isabs(outdir):
    outdir = os.path.join(PROJECT_DIR, outdir)
fasta = config["fasta"]
if not os.path.isabs(fasta):
    fasta = os.path.join(PROJECT_DIR, fasta)
minimap2_opts = config.get("opts", {}).get("minimap2", " -ax map-ont ")

# Resolve base_calling_model path
_bcm = config.get("base_calling_model", "")
if _bcm and not os.path.isabs(_bcm):
    config["base_calling_model"] = os.path.join(PIPELINE_DIR, _bcm)
dorado_opts = config.get("opts", {}).get("dorado", " --modified-bases m5C_2OmeC inosine_m6A_2OmeA pseU_2OmeU 2OmeG ")


# --- Sample parsing ---

def parse_samples_yaml(fl):
    """Parse YAML-format sample file with barcode assignments."""
    samples = {}
    with open(fl) as f:
        data = yaml.safe_load(f)
    for run_idx, run in enumerate(data["runs"]):
        run_path = run["path"]
        run_id = os.path.basename(run_path.rstrip("/"))
        barcode_kit = run.get("barcode_kit", config.get("warpdemux", {}).get("barcode_kit"))
        for sample_name, sample_val in run["samples"].items():
            if isinstance(sample_val, str) or sample_val is None:
                barcode = sample_val
            elif isinstance(sample_val, dict):
                barcode = sample_val.get("wdx")
            else:
                sys.exit(f"Invalid sample value for '{sample_name}': {sample_val}")
            samples[sample_name] = {
                "path": {run_path},
                "barcode": barcode,
                "run_id": run_id,
                "barcode_kit": barcode_kit,
            }
    return samples


def find_raw_inputs(sample_dict):
    """Find pod5 files for each sample."""
    POD5_DIRS = ["pod5_pass", "pod5_fail", "pod5"]
    for sample, info in sample_dict.items():
        raw_fls = []
        for path in info["path"]:
            for subdir in POD5_DIRS:
                data_path = os.path.join(path, subdir, "*.pod5")
                raw_fls += glob.glob(data_path)
        if not raw_fls:
            sys.exit(f"No pod5 files found for sample: {sample}")
        sample_dict[sample]["raw_files"] = raw_fls
    return sample_dict


# Resolve samples path relative to project root if relative
_samples_path = config["samples"]
if not os.path.isabs(_samples_path):
    _samples_path = os.path.join(PROJECT_DIR, _samples_path)
samples = parse_samples_yaml(_samples_path)
samples = find_raw_inputs(samples)

wildcard_constraints:
    sample="|".join(samples.keys()),


# --- Helper functions ---

def get_run_ids():
    run_ids = set()
    for sample, info in samples.items():
        if info.get("barcode") and info.get("run_id"):
            run_ids.add(info["run_id"])
    return list(run_ids)


def get_run_path(run_id):
    for sample, info in samples.items():
        if info.get("run_id") == run_id:
            return list(info["path"])[0]
    return None


def get_run_pod5_dirs(run_id):
    POD5_DIRS = ["pod5_pass", "pod5_fail", "pod5"]
    run_path = get_run_path(run_id)
    return [
        os.path.join(run_path, d)
        for d in POD5_DIRS
        if os.path.isdir(os.path.join(run_path, d))
    ]


def get_run_raw_inputs(wildcards):
    POD5_DIRS = ["pod5_pass", "pod5_fail", "pod5"]
    run_path = get_run_path(wildcards.run_id)
    raw_fls = []
    for subdir in POD5_DIRS:
        raw_fls += glob.glob(os.path.join(run_path, subdir, "*.pod5"))
    if not raw_fls:
        sys.exit(f"No pod5 files found for run: {wildcards.run_id}")
    return raw_fls


def get_sample_barcode_mapping(wildcards):
    run_id = samples[wildcards.sample]["run_id"]
    return os.path.join(outdir, "demux", "read_ids", run_id, "barcode_mapping.tsv.gz")


def get_sample_run_raw_inputs(wildcards):
    run_id = samples[wildcards.sample]["run_id"]
    POD5_DIRS = ["pod5_pass", "pod5_fail", "pod5"]
    run_path = get_run_path(run_id)
    raw_fls = []
    for subdir in POD5_DIRS:
        raw_fls += glob.glob(os.path.join(run_path, subdir, "*.pod5"))
    if not raw_fls:
        sys.exit(f"No pod5 files found for run: {run_id}")
    return raw_fls


def get_barcode_kit_for_run(run_id):
    for sample, info in samples.items():
        if info.get("run_id") == run_id:
            return info.get("barcode_kit")
    return config.get("warpdemux", {}).get("barcode_kit")


def get_modkit_threshold_opts():
    """Build modkit threshold CLI options from config."""
    modkit_cfg = config.get("modkit", {})
    opts = []
    ft = modkit_cfg.get("filter_threshold")
    if ft is not None:
        opts.append(f"--filter-threshold {ft}")
    for mod_code, thresh in modkit_cfg.get("mod_thresholds", {}).items():
        opts.append(f"--mod-thresholds {mod_code}:{thresh}")
    return " ".join(opts)


# --- Target rule ---

rule all:
    input:
        # Final aligned BAMs with modification tags
        expand(
            os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam"),
            sample=samples.keys(),
        ),
        # Modkit pileup
        expand(
            os.path.join(outdir, "summary", "modkit", "{sample}", "{sample}.pileup.bed.gz"),
            sample=samples.keys(),
        ),
        # Modkit per-read calls
        expand(
            os.path.join(outdir, "summary", "modkit", "{sample}", "{sample}.mod_calls.tsv.gz"),
            sample=samples.keys(),
        ),


# ===== WarpDemuX Demultiplexing =====

rule warpdemux:
    """Run WarpDemuX barcode demultiplexing on raw POD5 files."""
    input:
        get_run_raw_inputs,
    output:
        outdir=directory(os.path.join(outdir, "demux", "warpdemux_output", "{run_id}")),
        done=os.path.join(outdir, "demux", "warpdemux_output", "{run_id}", ".done"),
    log:
        os.path.join(outdir, "logs", "warpdemux", "{run_id}"),
    params:
        model=lambda wildcards: get_barcode_kit_for_run(wildcards.run_id),
        save_boundaries=lambda wildcards: (
            "true"
            if config.get("warpdemux", {}).get("save_boundaries", True)
            else "false"
        ),
        pod5_dirs=lambda wildcards: " ".join(get_run_pod5_dirs(wildcards.run_id)),
    threads: config.get("warpdemux", {}).get("threads", 8)
    shell:
        """
        warpdemux demux \
            -i {params.pod5_dirs} \
            -o {output.outdir} \
            -m {params.model} \
            -j {threads} \
            --save_boundaries {params.save_boundaries} 2>&1 | tee {log}
        touch {output.done}
        """


rule parse_warpdemux:
    """Parse WarpDemuX predictions and create barcode mapping per run."""
    input:
        demux_done=os.path.join(outdir, "demux", "warpdemux_output", "{run_id}", ".done"),
        demux_dir=os.path.join(outdir, "demux", "warpdemux_output", "{run_id}"),
    output:
        mapping=os.path.join(outdir, "demux", "read_ids", "{run_id}", "barcode_mapping.tsv.gz"),
        summary=os.path.join(outdir, "demux", "read_ids", "{run_id}", "demux_summary.tsv.gz"),
    log:
        os.path.join(outdir, "logs", "parse_warpdemux", "{run_id}"),
    run:
        import gzip
        import pandas as pd
        from pathlib import Path

        demux_base = Path(input.demux_dir)
        subdirs = [d for d in demux_base.iterdir() if d.is_dir() and d.name.startswith("warpdemux_")]
        pred_dir = subdirs[0] / "predictions" if subdirs else demux_base / "predictions"

        all_predictions = []
        for pred_file in pred_dir.glob("*.csv.gz"):
            with gzip.open(pred_file, "rt") as f:
                df = pd.read_csv(f, comment=None)
                df.columns = [c.lstrip("#") for c in df.columns]
                all_predictions.append(df)

        predictions = (
            pd.concat(all_predictions, ignore_index=True)
            if all_predictions
            else pd.DataFrame(columns=["read_id", "predicted_barcode"])
        )
        predictions["predicted_barcode"] = predictions["predicted_barcode"].apply(
            lambda x: f"barcode{int(x):02d}" if x != -1 else "unclassified"
        )
        predictions[["read_id", "predicted_barcode"]].to_csv(
            output.mapping, sep="\t", index=False, compression="gzip"
        )
        summary_data = predictions.groupby("predicted_barcode").size().reset_index(name="n_reads")
        summary_data.to_csv(output.summary, sep="\t", index=False, compression="gzip")


rule extract_sample_reads:
    """Extract read IDs for a specific sample based on barcode assignment."""
    input:
        mapping=get_sample_barcode_mapping,
    output:
        read_ids=os.path.join(outdir, "demux", "read_ids", "{sample}", "{sample}.txt"),
    params:
        barcode=lambda wildcards: samples[wildcards.sample]["barcode"],
    run:
        import pandas as pd
        mapping = pd.read_csv(input.mapping, sep="\t", compression="gzip")
        sample_reads = mapping[mapping["predicted_barcode"] == params.barcode]["read_id"]
        with open(output.read_ids, "w") as f:
            for read_id in sample_reads:
                f.write(f"{read_id}\n")


rule split_pod5:
    """Filter raw POD5 files by sample using read IDs from demultiplexing."""
    input:
        pod5=get_sample_run_raw_inputs,
        read_ids=rules.extract_sample_reads.output.read_ids,
    output:
        os.path.join(outdir, "demux", "pod5", "{sample}", "{sample}.pod5"),
    log:
        os.path.join(outdir, "logs", "split_pod5", "{sample}"),
    params:
        pod5_dirs=lambda wildcards: " ".join(
            get_run_pod5_dirs(samples[wildcards.sample]["run_id"])
        ),
    shell:
        """
        pod5 filter {params.pod5_dirs} --ids {input.read_ids} --missing-ok --output {output} 2>&1 | tee {log}
        """


# ===== Dorado Basecalling =====

rule download_mod_models:
    """Download dorado modification models if not already present."""
    output:
        sentinel=os.path.join(PIPELINE_DIR, "resources", "models", ".mod_models_ready"),
    params:
        models_dir=os.path.join(PIPELINE_DIR, "resources", "models"),
        base_model=config["dorado_model"],
        mod_bases=" ".join(dorado_opts.replace("--modified-bases", "").replace("--emit-moves", "").split()),
    shell:
        """
        for mod in {params.mod_bases}; do
            model="{params.base_model}_${{mod}}@v1"
            model_path="{params.models_dir}/$model"
            if [ -f "$model_path/.downloaded" ]; then
                echo "Model $model already downloaded"
            else
                echo "Downloading $model..."
                dorado download --model "$model" --models-directory {params.models_dir}
                touch "$model_path/.downloaded"
                echo "Downloaded $model"
            fi
        done
        touch {output.sentinel}
        """


rule rebasecall:
    """Basecall with dorado including modification detection."""
    input:
        pod5=os.path.join(outdir, "demux", "pod5", "{sample}", "{sample}.pod5"),
        mod_models=rules.download_mod_models.output.sentinel,
    output:
        os.path.join(outdir, "bam", "rebasecall", "{sample}", "{sample}.rbc.bam"),
    log:
        os.path.join(outdir, "logs", "rebasecall", "{sample}"),
    params:
        model=config["base_calling_model"],
        dorado_opts=dorado_opts,
        models_dir=os.path.join(PIPELINE_DIR, "resources", "models"),
    shell:
        """
        if [[ "${{CUDA_VISIBLE_DEVICES:-}}" ]]; then
            echo "CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES"
            export CUDA_VISIBLE_DEVICES
        fi

        dorado basecaller \
            --models-directory {params.models_dir} \
            {params.dorado_opts} \
            {params.model} {input.pod5} > {output} 2> {log}
        """


# ===== minimap2 Alignment to GENCODE =====

rule ubam_to_fastq:
    """Extract reads from uBAM to FASTQ for alignment."""
    input:
        rules.rebasecall.output,
    output:
        os.path.join(outdir, "fq", "{sample}", "{sample}.fq.gz"),
    shell:
        """
        samtools fastq {input} | gzip > {output}
        """


rule minimap2_align:
    """Align reads to GENCODE transcriptome with minimap2."""
    input:
        reads=rules.ubam_to_fastq.output,
    output:
        bam=os.path.join(outdir, "bam", "aln", "{sample}", "{sample}.aln.bam"),
        bai=os.path.join(outdir, "bam", "aln", "{sample}", "{sample}.aln.bam.bai"),
    params:
        ref=fasta,
        mm2_opts=minimap2_opts,
    log:
        os.path.join(outdir, "logs", "minimap2_align", "{sample}"),
    threads: 16
    shell:
        """
        minimap2 {params.mm2_opts} -t {threads} {params.ref} {input.reads} \
            | samtools sort -@ 4 -m 2G -o {output.bam} 2> {log}

        samtools index {output.bam}
        """


rule inject_ubam_tags:
    """Transfer modification tags (MM/ML) from dorado uBAM to aligned BAM."""
    input:
        source_bam=rules.rebasecall.output,
        target_bam=rules.minimap2_align.output.bam,
        target_bai=rules.minimap2_align.output.bai,
    output:
        bam=os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam"),
        bai=os.path.join(outdir, "bam", "final", "{sample}", "{sample}.bam.bai"),
    params:
        src=SCRIPT_DIR,
    log:
        os.path.join(outdir, "logs", "inject_ubam_tags", "{sample}"),
    shell:
        """
        python {params.src}/transfer_tags.py \
            --all-tags \
            --source {input.source_bam} \
            --target {input.target_bam} \
            --output {output.bam} 2> {log}

        samtools index {output.bam}
        """


# ===== Modification Calling =====

rule modkit_pileup:
    """Run modkit pileup for modification frequency at each position."""
    input:
        bam=rules.inject_ubam_tags.output.bam,
        bai=rules.inject_ubam_tags.output.bai,
    output:
        bed=os.path.join(outdir, "summary", "modkit", "{sample}", "{sample}.pileup.bed.gz"),
    log:
        os.path.join(outdir, "logs", "modkit", "pileup", "{sample}"),
    params:
        ref=fasta,
        threshold_opts=get_modkit_threshold_opts(),
    shell:
        """
        modkit pileup \
            --log-filepath {log} \
            --ref {params.ref} \
            {params.threshold_opts} \
            {input.bam} - \
            | gzip > {output.bed}
        """


rule modkit_extract_calls:
    """Extract per-read modification calls."""
    input:
        bam=rules.inject_ubam_tags.output.bam,
        bai=rules.inject_ubam_tags.output.bai,
    output:
        tsv=os.path.join(outdir, "summary", "modkit", "{sample}", "{sample}.mod_calls.tsv.gz"),
    log:
        os.path.join(outdir, "logs", "modkit", "extract_calls", "{sample}"),
    params:
        ref=fasta,
        threshold_opts=get_modkit_threshold_opts(),
    shell:
        """
        modkit extract calls \
            --reference {params.ref} \
            --log-filepath {log} \
            --mapped --pass \
            {params.threshold_opts} \
            {input.bam} - \
            | gzip > {output.tsv}
        """
