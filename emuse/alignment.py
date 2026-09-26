import statistics

import pandas as pd
import pysam

METRICS_COLUMNS = ["tax_id", "median_identity", "median_coverage"]


def alignment_metrics(alignment_path):
    identities = {}
    coverages = {}
    with pysam.AlignmentFile(str(alignment_path)) as align_file:
        for aln in align_file:
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            taxid = aln.reference_name.split(":")[0]
            identity, coverage = get_align_stats(aln, align_file)
            identities.setdefault(taxid, []).append(identity)
            coverages.setdefault(taxid, []).append(coverage)

    rows = [
        [taxid, statistics.median(identities[taxid]), statistics.median(coverages[taxid])]
        for taxid in identities
    ]
    return pd.DataFrame(rows, columns=METRICS_COLUMNS)


def write_alignment_metrics(metrics, path):
    metrics[METRICS_COLUMNS].to_csv(path, sep="\t", index=False)


def read_alignment_metrics(path):
    return pd.read_csv(path, sep="\t", dtype={"tax_id": str})[METRICS_COLUMNS]


def get_align_stats(alignment, alignment_file):
    CIGAROP_MATCH = 0
    CIGAROP_INS = 1
    CIGAROP_DEL = 2
    CIGAROP_REF_SKIP = 3
    CIGAROP_SOFTCLIP = 4
    CIGAROP_HARDCLIP = 5
    CIGAROP_PAD = 6
    CIGAROP_EQUAL = 7
    CIGAROP_DIFF = 8
    CIGAROP_BACK = 9

    cigar_stats = alignment.get_cigar_stats()[0]
    edit_distance = alignment.get_tag("NM") if alignment.has_tag("NM") else None

    cigar_stats = alignment.get_cigar_stats()[0]

    insertions = cigar_stats[CIGAROP_INS]
    deletions = cigar_stats[CIGAROP_DEL]
    soft_clips = cigar_stats[CIGAROP_SOFTCLIP]

    query_len = alignment.query_length
    query_alignment_len = alignment.query_alignment_length

    reference_len = alignment_file.get_reference_length(alignment.reference_name)
    reference_alignment_len = alignment.reference_length

    identity, ref_coverage = calculate_align_stats(query_len,
                                                   query_alignment_len,
                                                   reference_len,
                                                   reference_alignment_len,
                                                   insertions,
                                                   deletions,
                                                   edit_distance,
                                                   soft_clips)

    return identity, ref_coverage



def calculate_align_stats(query_len,
                          query_alignment_len,
                          reference_len,
                          reference_alignment_len,
                          insertions,
                          deletions,
                          edit_distance,
                          soft_clips):
    """
    Return identity and coverage against the reference sequence

    Note that the identity in this calculation differs from the one in BLAST.
    While BLAST's identity measure does not count indels and instead divides
    the number of identical bases over the number of aligned bases, the
    identity metric below includes each aligned column corresponding to a base
    pair position in the reference or read, and counts the number of matches
    (identical bases) across the total aligned length, which is the length
    where gaps in both the reference and query are included.

    Args:
        query_len: The length in base pairs of the query sequence.
        query_alignment_len: The length in base pairs of the part of the query
            sequence covered by the alignment.
        reference_len: The length in base pairs of the reference sequence.
        reference_alignment_len: The length in base pairs of the part of the
            reference sequence covered by the alignment.
        insertions: Number of insertions measured in base pairs.
        deletions: Number of deletions measured in base pairs.
        edit_distance: Edit distance corresponding to the NM tag in SAM files.
        soft_clips: Bases in the ends of the alignment which do not align.

    Returns:
        identity: Identity as matching bases across the alignment length which
            contains both insertions and deletions.
        coverage: Percent of bases of the reference length covered by the
            alignment.
    """

    # ------------------------------------------------
    # Calculate identity
    # ------------------------------------------------
    matches = 0
    mismatches = 0
    if edit_distance is not None:
        mismatches = edit_distance - insertions - deletions
        matches = query_alignment_len - insertions - mismatches

    # We can not normalize over query or reference length, as these
    # won't contain either insertions or deletions
    divisor = matches + mismatches + insertions + deletions

    identity = 0
    if divisor > 0:
        identity = matches / divisor

    # ------------------------------------------------
    # Calculate coverage
    # ------------------------------------------------
    coverage = 0
    if reference_len > 0:
        coverage = reference_alignment_len / reference_len

    return identity, coverage
