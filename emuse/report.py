from datetime import date
from importlib import resources
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

from emuse.abundance import abundance_table, read_rel_abundance, to_records
from emuse.alignment import alignment_metrics as compute_alignment_metrics
from emuse.files import find_emu_file
from emuse.negative_control import compare_to_negative_controls
from emuse.qc import summary_stats, trana_version
from emuse.read_assignment import read_assignment_stats

DATA_DIR = resources.files("emuse") / "data"

DISPLAY_COLUMNS = {
    "tax_id": "tax id",
    "estimated_counts": "estimated read counts",
    "median_probability": "median probability*",
    "mean_probability": "mean probability*",
    "median_identity": "median aligned identity",
    "median_coverage": "median aligned coverage",
}


def row_style(flags, attribute, style):
    def apply(row):
        return [style if getattr(flags[row.name], attribute) else ""] * len(row)
    return apply


def render_report(input_dir, sample_name, neg_control, spike_species=(), normalising_spike_species=None,
                  prob_score=False, alignment_metrics=False):
    # Read CSS file content
    with open(DATA_DIR / "static" / "style.css", "r") as f:
        css_content = f.read()

    env = Environment(loader=FileSystemLoader(str(DATA_DIR / "templates")))
    template = env.get_template("report.html.j2")

    results_dir = f"{input_dir}/results"

    # Load sample read assignment table
    assignment_summary = read_assignment_stats(find_emu_file(results_dir, sample_name, "_read-assignment-distributions.tsv")).rename(columns=DISPLAY_COLUMNS)

    # Load neg control abundance table
    neg_control_ordered = abundance_table(read_rel_abundance(find_emu_file(results_dir, neg_control, "_rel-abundance.tsv"))).rename(columns=DISPLAY_COLUMNS)
    # Create fake index column for styling purposes (need it to start from 1 instead of 0)
    neg_control_ordered.insert(0, "row", range(1, len(neg_control_ordered) + 1))

    # Load sample abundance table
    abundance_ordered = abundance_table(read_rel_abundance(find_emu_file(results_dir, sample_name, "_rel-abundance.tsv"))).rename(columns=DISPLAY_COLUMNS)

    # Merge abundance and assignment if prob_score is given
    if prob_score:
        abundance_assignment = abundance_ordered.merge(assignment_summary, on="tax id", how="left")
    else:
        abundance_assignment = abundance_ordered

    # Create fake index column for styling purposes (need it to start from 1 instead of 0)
    abundance_assignment.insert(0, "row", range(1, len(abundance_ordered) + 1))

    # Merge abundance and alignment based metrics
    if alignment_metrics:
        sample_alignment_metrics = compute_alignment_metrics(find_emu_file(results_dir, sample_name, "_emu_alignments.sam"))
        abundance_assignment = abundance_assignment.merge(sample_alignment_metrics.rename(columns=DISPLAY_COLUMNS), on="tax id", how="left")
        neg_control_alignment_metrics = compute_alignment_metrics(find_emu_file(results_dir, neg_control, "_emu_alignments.sam"))
        neg_control_ordered = neg_control_ordered.merge(neg_control_alignment_metrics.rename(columns=DISPLAY_COLUMNS), on="tax id", how="left")

    highlight = set(spike_species)

    sample_flags = compare_to_negative_controls(
        to_records(abundance_ordered), {neg_control: to_records(neg_control_ordered)}, highlight, normalising_spike_species
    )
    neg_control_flags = compare_to_negative_controls(to_records(neg_control_ordered), {}, highlight)

    # Apply functions for spike species and unique species
    styled_abundance = (abundance_assignment.style
        .apply(row_style(sample_flags, "enriched", "background-color: #dcfce7"), axis=1)
        .apply(row_style(sample_flags, "absent_in_controls", "background-color: #dcfce7"), axis=1)
        .apply(row_style(sample_flags, "spike", "background-color: #ddd6fe"), axis=1)
        .apply(row_style(sample_flags, "low_abundance", "color: #9ca3af"), axis=1)
        .format({
            "estimated read counts": "{:.0f}",
            "abundance": "{:.2%}",
            "median aligned identity": "{:.2%}",
            "median aligned coverage": "{:.2%}",
            })
        .set_properties(subset=["row"], **{
            "background-color": "#f2f2f2",
            "font-weight": "bold"
        })
        .hide(axis="index")
    )

    # Convert to html table
    html_table = styled_abundance.to_html(index=False, border=0, escape=False)

    # Apply function for spike species
    styled_neg_control = (neg_control_ordered.style
        .apply(row_style(neg_control_flags, "spike", "background-color: #ddd6fe"), axis=1)
        .apply(row_style(neg_control_flags, "low_abundance", "color: #9ca3af"), axis=1)
        .format({
            "estimated read counts": "{:.0f}",
            "abundance": "{:.2%}",
            "median aligned identity": "{:.2%}",
            "median aligned coverage": "{:.2%}",
            })
        .set_properties(subset=["row"], **{
            "background-color": "#f2f2f2",
            "font-weight": "bold"
        })
        .hide(axis="index")            
    )
    # Convert to html table
    neg_control_html_table = styled_neg_control.to_html(index=False, border=0, escape=False)

    legend_lines = [
        '<span class="inline-block w-4 h-4 bg-purple-200 mr-2 border"></span> Purple rows indicate spike species<br>',
        '<span class="inline-block w-4 h-4 bg-green-100 mr-2 border"></span> Green rows indicate species not found in negative control or with a normalised abundance ratio &ge; 25.<br>'
    ]
    # Optional line for probability score
    if prob_score:
        legend_lines.append(
            '<span class="not-italic text-sm"><span class="font-bold text-black">median/mean probability*</span>: Median/Mean probability of the assigned taxon across reads</span>'
        )

    legend_html = '<p class="text-gray-700 italic mt-2">'+ "".join(legend_lines) + '</p>'

    legend_neg_html = """
    <p class="text-gray-700 italic mt-2">
        <span class="inline-block w-4 h-4 bg-purple-200 mr-2 border"></span>
        Purple rows indicate spike species
    </p>
    """

    # Save date
    today = date.today().strftime("%Y-%m-%d")

    stats = summary_stats(f"{input_dir}/multiqc/multiqc_data/multiqc_data.json", sample_name)

    pipeline_version = trana_version(f"{input_dir}/pipeline_info/software_versions.yml")

    return template.render(
        css = css_content,
        table = html_table,
        neg_control_table = neg_control_html_table,
        legend = legend_html,
        legend_neg = legend_neg_html,
        today = today,
        pipeline_version = pipeline_version,
        stats = stats,
        input_dir = Path(input_dir).name,
        sample_name = sample_name,
        neg_control = neg_control
    )
