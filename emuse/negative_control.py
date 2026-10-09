from dataclasses import asdict, dataclass, field

ENRICHMENT_FOLD = 25
LOW_ABUNDANCE_CUTOFF = 0.005


@dataclass
class ControlComparison:
    control_name: str
    present_in_control: bool
    sample_ratio: float | None = None
    control_ratio: float | None = None
    enriched: bool = False


@dataclass
class TaxonFlags:
    species: str | None
    spike: bool = False
    absent_in_controls: bool = False
    enriched: bool = False
    low_abundance: bool = False
    comparisons: list[ControlComparison] = field(default_factory=list)

    @property
    def highlight(self):
        if self.spike:
            return "spike"
        if self.absent_in_controls or self.enriched:
            return "not_in_control"
        return None

    def to_dict(self):
        return dict(asdict(self), highlight=self.highlight)


def first_abundance(rows, species):
    for row in rows:
        if row.get("species") == species:
            return row.get("abundance")
    return None


def compare_to_control(sample_rows, control_rows, control_name, normalising_species=None, fold=ENRICHMENT_FOLD):
    control_species = {row.get("species") for row in control_rows}
    sample_spike = first_abundance(sample_rows, normalising_species)
    control_spike = first_abundance(control_rows, normalising_species)

    comparisons = []
    for row in sample_rows:
        species = row.get("species")
        comparison = ControlComparison(control_name, species in control_species)
        control_abundance = first_abundance(control_rows, species)
        if normalising_species and sample_spike and control_spike and control_abundance is not None:
            comparison.sample_ratio = row["abundance"] / sample_spike
            comparison.control_ratio = control_abundance / control_spike
            comparison.enriched = 0 < comparison.control_ratio * fold < comparison.sample_ratio
        comparisons.append(comparison)
    return comparisons


def compare_to_negative_controls(sample_rows, controls, spike_species=(), normalising_species=None,
                                 fold=ENRICHMENT_FOLD, low_abundance_cutoff=LOW_ABUNDANCE_CUTOFF):
    per_control = [
        compare_to_control(sample_rows, rows, name, normalising_species, fold)
        for name, rows in controls.items()
    ]

    flags = []
    for i, row in enumerate(sample_rows):
        comparisons = [results[i] for results in per_control]
        present = [c for c in comparisons if c.present_in_control]
        flags.append(TaxonFlags(
            species=row.get("species"),
            spike=row.get("species") in spike_species,
            absent_in_controls=bool(comparisons) and not present,
            enriched=bool(present) and all(c.enriched for c in present),
            low_abundance=row["abundance"] < low_abundance_cutoff,
            comparisons=comparisons,
        ))
    return flags
