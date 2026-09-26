from dataclasses import asdict, dataclass
import json
import re

import yaml

MULTIQC_KEYS = {
    "number_of_reads": "Number of reads_fastq",
    "mean_read_length": "Mean read length_fastq",
    "mean_read_quality": "Mean read quality_fastq",
    "median_read_length": "Median read length_fastq",
    "median_read_quality": "Median read quality_fastq",
    "read_length_n50": "Read length N50_fastq",
    "stdev_read_length": "STDEV read length_fastq",
    "total_bases": "Total bases_fastq",
}


@dataclass
class SummaryStats:
    number_of_reads: float | None = None
    mean_read_length: float | None = None
    mean_read_quality: float | None = None
    median_read_length: float | None = None
    median_read_quality: float | None = None
    read_length_n50: float | None = None
    stdev_read_length: float | None = None
    total_bases: float | None = None
    mapped_reads: int | None = None

    def to_dict(self):
        return asdict(self)


def nanoplot_stats(multiqc_data, sample_name):
    key = f"{sample_name}_nanoplot_processed"
    for section in multiqc_data.get("report_general_stats_data", []):
        if key in section:
            return section[key]
    return {}


def parse_mapped_reads(emu_log_text):
    match = re.search(r"mapped (\d+) sequences", emu_log_text)
    return int(match.group(1)) if match else None


def summary_stats(multiqc_json_path, sample_name, emu_log_path=None):
    with open(multiqc_json_path) as f:
        stats = nanoplot_stats(json.load(f), sample_name)

    mapped_reads = None
    if emu_log_path is not None:
        with open(emu_log_path) as f:
            mapped_reads = parse_mapped_reads(f.read())

    values = {field: stats.get(key) for field, key in MULTIQC_KEYS.items()}
    return SummaryStats(**values, mapped_reads=mapped_reads)


def quality_class(quality):
    if quality is None:
        return None
    if quality < 15:
        return "not_acceptable"
    if quality < 20:
        return "good"
    return "excellent"


def trana_version(software_versions_path):
    with open(software_versions_path) as f:
        versions = yaml.safe_load(f) or {}
    return versions.get("Workflow", {}).get("genomic-medicine-sweden/TRANA")
