from dataclasses import dataclass
from pathlib import Path

FASTQ_SUFFIXES = ["_downsampled", "_filtered", ".porechop_abi", ".trimmed.trim", ""]


@dataclass
class TranaOutputs:
    rel_abundance: Path | None = None
    read_assignment: Path | None = None
    alignment: Path | None = None
    emu_log: Path | None = None
    multiqc_data: Path | None = None
    software_versions: Path | None = None


def existing(path):
    return path if path.is_file() else None


def find_emu_file(results_dir, sample_name, suffix):
    for fastq_suffix in FASTQ_SUFFIXES:
        path = Path(results_dir) / f"{sample_name}{fastq_suffix}.fastq{suffix}"
        if path.is_file():
            return path
    return None


def locate_outputs(outdir, sample_name):
    outdir = Path(outdir)
    results = outdir / "results"
    return TranaOutputs(
        rel_abundance=find_emu_file(results, sample_name, "_rel-abundance.tsv"),
        read_assignment=find_emu_file(results, sample_name, "_read-assignment-distributions.tsv"),
        alignment=find_emu_file(results, sample_name, "_emu_alignments.sam"),
        emu_log=existing(results / "emu_logs" / f"{sample_name}_emu_log.log"),
        multiqc_data=existing(outdir / "multiqc" / "multiqc_data" / "multiqc_data.json"),
        software_versions=existing(outdir / "pipeline_info" / "software_versions.yml"),
    )
