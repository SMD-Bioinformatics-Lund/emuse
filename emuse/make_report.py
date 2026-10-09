import argparse
from importlib import resources
from pathlib import Path
import tomllib

from emuse.report import render_report

# Bundled package data (templates, CSS, taxonomy mapping, default config)
DATA_DIR = resources.files("emuse") / "data"
DEFAULT_CONFIG = DATA_DIR / "configs" / "config.toml"


def main():
    argp = argparse.ArgumentParser()
    argp.add_argument("-i", "--input-dir", type=str, required=True, help="Path to the input directory containing results")
    argp.add_argument("-o", "--output-file", type=str, required=True, help="Path to the output report file")
    argp.add_argument("-s", "--sample-name", type=str, required=True, help="Name of the sample")
    argp.add_argument("-n", "--neg-control", type=str, required=True, help="Name of the negative control")
    argp.add_argument("-c", "--config", type=Path, default=DEFAULT_CONFIG, help="Path to config file. Default is the bundled config.toml")
    argp.add_argument("-p", "--prob-score", action="store_true", help="Include probability score in the report")
    argp.add_argument("-m", "--alignment-metrics", action="store_true", help="Include metrics based on the raw alignment of reads to the database (percent identity and percent coverage)")

    args = argp.parse_args()

    # Spike species
    with open(args.config, "rb") as f:
        config = tomllib.load(f)

    html = render_report(
        args.input_dir,
        args.sample_name,
        args.neg_control,
        spike_species=config.get("spike_species", []),
        normalising_spike_species=config.get("normalising_spike_species"),
        prob_score=args.prob_score,
        alignment_metrics=args.alignment_metrics,
    )

    with open(args.output_file, "w") as f:
        f.write(html)


if __name__ == "__main__":
    main()
