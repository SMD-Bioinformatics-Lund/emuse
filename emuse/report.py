from datetime import date
from importlib import resources
from pathlib import Path

from jinja2 import Environment, PackageLoader

from emuse.abundance import abundance_table, read_rel_abundance, to_records
from emuse.alignment import alignment_metrics
from emuse.files import locate_outputs
from emuse.negative_control import compare_to_negative_controls
from emuse.qc import summary_stats, trana_version
from emuse.read_assignment import read_assignment_stats

SPIKE_STYLE = "background-color: #ddd6fe"
NOT_IN_CONTROL_STYLE = "background-color: #dcfce7"
LOW_ABUNDANCE_STYLE = "color: #9ca3af"

DISPLAY_COLUMNS = {
    "tax_id": "tax id",
    "estimated_counts": "estimated read counts",
    "median_probability": "median probability*",
    "mean_probability": "mean probability*",
    "median_identity": "median aligned identity",
    "median_coverage": "median aligned coverage",
}

TABLE_FORMAT = {
    "estimated read counts": "{:.0f}",
    "abundance": "{:.2%}",
    "median aligned identity": "{:.2%}",
    "median aligned coverage": "{:.2%}",
}

LEGEND_NEG_HTML = """
    <p class="text-gray-700 italic mt-2">
        <span class="inline-block w-4 h-4 bg-purple-200 mr-2 border"></span>
        Purple rows indicate spike species
    </p>
    """


def row_styles(flags, attribute, style):
    def apply(row):
        return [style if getattr(flags[row.name], attribute) else ""] * len(row)
    return apply


def style_table(df, flags, compare_to_control):
    styler = df.style
    if compare_to_control:
        styler = (styler
            .apply(row_styles(flags, "enriched", NOT_IN_CONTROL_STYLE), axis=1)
            .apply(row_styles(flags, "absent_in_controls", NOT_IN_CONTROL_STYLE), axis=1))
    return (styler
        .apply(row_styles(flags, "spike", SPIKE_STYLE), axis=1)
        .apply(row_styles(flags, "low_abundance", LOW_ABUNDANCE_STYLE), axis=1)
        .format(TABLE_FORMAT)
        .set_properties(subset=["row"], **{
            "background-color": "#f2f2f2",
            "font-weight": "bold"
        })
        .hide(axis="index")
        .to_html(index=False, border=0, escape=False))


def build_legend(prob_score):
    legend_lines = [
        '<span class="inline-block w-4 h-4 bg-purple-200 mr-2 border"></span> Purple rows indicate spike species<br>',
        '<span class="inline-block w-4 h-4 bg-green-100 mr-2 border"></span> Green rows indicate species not found in negative control or with a normalised abundance ratio &ge; 25.<br>'
    ]
    if prob_score:
        legend_lines.append(
            '<span class="not-italic text-sm"><span class="font-bold text-black">median/mean probability*</span>: Median/Mean probability of the assigned taxon across reads</span>'
        )
    return '<p class="text-gray-700 italic mt-2">' + "".join(legend_lines) + '</p>'


def render_report(input_dir, sample_name, neg_control, spike_species=(),
                  normalising_spike_species=None, prob_score=False, include_alignment_metrics=False):
    sample_files = locate_outputs(input_dir, sample_name)
    control_files = locate_outputs(input_dir, neg_control)

    sample = abundance_table(read_rel_abundance(sample_files.rel_abundance))
    control = abundance_table(read_rel_abundance(control_files.rel_abundance))

    if prob_score:
        sample = sample.merge(read_assignment_stats(sample_files.read_assignment), on="tax_id", how="left")

    if include_alignment_metrics:
        sample = sample.merge(alignment_metrics(sample_files.alignment), on="tax_id", how="left")
        control = control.merge(alignment_metrics(control_files.alignment), on="tax_id", how="left")

    sample_flags = compare_to_negative_controls(
        to_records(sample), {neg_control: to_records(control)}, spike_species, normalising_spike_species
    )
    control_flags = compare_to_negative_controls(to_records(control), {}, spike_species)

    for df in (sample, control):
        df.insert(0, "row", range(1, len(df) + 1))
        df.rename(columns=DISPLAY_COLUMNS, inplace=True)

    stats = summary_stats(sample_files.multiqc_data, sample_name, sample_files.emu_log)
    template = Environment(loader=PackageLoader("emuse", "data/templates")).get_template("report.html.j2")

    return template.render(
        css=(resources.files("emuse") / "data" / "static" / "style.css").read_text(),
        table=style_table(sample, sample_flags, compare_to_control=True),
        neg_control_table=style_table(control, control_flags, compare_to_control=False),
        legend=build_legend(prob_score),
        legend_neg=LEGEND_NEG_HTML,
        today=date.today().strftime("%Y-%m-%d"),
        pipeline_version=trana_version(sample_files.software_versions),
        stats=stats,
        mapped_value=stats.mapped_reads if stats.mapped_reads is not None else "N/A",
        input_dir=Path(input_dir).name,
        sample_name=sample_name,
        neg_control=neg_control,
    )
