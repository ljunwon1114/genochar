from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict

from .utils import canonicalize_name, find_column, open_maybe_gzip, parse_numeric, sniff_delimiter


def parse_metadata_table(path: Path | str) -> Dict[str, dict]:
    delimiter = sniff_delimiter(path)
    with open_maybe_gzip(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            return {}

        name_col = find_column(reader.fieldnames, ["Strain", "Strain name", "Sample", "Genome", "Bin", "Name"])
        coverage_col = find_column(reader.fieldnames, ["Sequencing coverage (×)", "Coverage (×)", "Coverage"])
        platforms_col = find_column(reader.fieldnames, ["Sequencing platforms", "Sequencing platform", "Platforms", "Platform"])
        assembly_col = find_column(reader.fieldnames, ["Assembly method", "Assembler", "Assembly"])
        genus_col = find_column(reader.fieldnames, ["Genus", "Genus name"])
        species_col = find_column(reader.fieldnames, ["Species", "Species name", "Specific epithet"])
        accession_col = find_column(reader.fieldnames, ["Accession", "Assembly accession", "Genome accession", "NCBI accession", "RefSeq accession", "GenBank accession"])
        repeat_col = find_column(reader.fieldnames, ["Repeat regions", "Repeat region", "Repeat count"])

        if not name_col:
            raise ValueError(f"Could not identify the strain/name column in metadata file: {path}")

        out: Dict[str, dict] = {}
        for row in reader:
            name = canonicalize_name(str(row.get(name_col, "")).strip())
            if not name:
                continue

            coverage = parse_numeric(row.get(coverage_col)) if coverage_col else None
            repeat_regions = parse_numeric(row.get(repeat_col)) if repeat_col else None

            out[name] = {
                "Sequencing coverage (×)": round(coverage, 2) if coverage is not None else None,
                "Sequencing platforms": row.get(platforms_col) if platforms_col else None,
                "Assembly method": row.get(assembly_col) if assembly_col else None,
                "Genus": row.get(genus_col).strip() if genus_col and row.get(genus_col) else None,
                "Species": row.get(species_col).strip() if species_col and row.get(species_col) else None,
                "Accession": row.get(accession_col).strip() if accession_col and row.get(accession_col) else None,
                "Repeat regions": int(repeat_regions) if repeat_regions is not None else None,
            }
        return out
