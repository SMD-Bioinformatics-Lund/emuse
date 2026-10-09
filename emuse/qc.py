from dataclasses import asdict, dataclass
import json

import yaml

PROCESSED_KEYS = {
    "reads_after_filtering": "Number of reads_fastq",
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
    reads_before_filtering: float | None = None
    reads_after_filtering: float | None = None
    mean_read_length: float | None = None
    mean_read_quality: float | None = None
    median_read_length: float | None = None
    median_read_quality: float | None = None
    read_length_n50: float | None = None
    stdev_read_length: float | None = None
    total_bases: float | None = None

    def to_dict(self):
        return asdict(self)


def nanoplot_stats(multiqc_data, sample_name, stage):
    key = f"{sample_name}_nanoplot_{stage}"
    for section in multiqc_data.get("report_general_stats_data", []):
        if key in section:
            return section[key]
    return {}


def summary_stats(multiqc_json_path, sample_name):
    with open(multiqc_json_path) as f:
        multiqc_data = json.load(f)
    unprocessed = nanoplot_stats(multiqc_data, sample_name, "unprocessed")
    processed = nanoplot_stats(multiqc_data, sample_name, "processed")
    values = {field: processed.get(key) for field, key in PROCESSED_KEYS.items()}
    return SummaryStats(reads_before_filtering=unprocessed.get("Number of reads_fastq"), **values)


def quality_class(quality):
    if quality is None:
        return None
    if quality < 15:
        return "not_acceptable"
    if quality < 20:
        return "good"
    return "excellent"


def trana_version(software_versions_path):
    with open(software_versions_path) as v:
        software_versions = yaml.safe_load(v)
    return software_versions["Workflow"]["genomic-medicine-sweden/TRANA"]
