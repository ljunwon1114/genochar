from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Sequence

import pandas as pd

from . import __version__
from .assembly_stats import assembly_stats_to_row, compute_assembly_stats
from .checkm2 import parse_checkm2_report
from .coverage import parse_coverage_table
from .gff_stats import gff_stats_to_row, parse_gff_stats
from .metadata import parse_metadata_table
from .pipeline import (
    PipelineError,
    find_existing_gffs,
    run_bakta,
    run_checkm2,
    run_prokka,
)
from .utils import (
    canonicalize_name,
    ensure_parent,
    looks_like_fasta,
    looks_like_gff,
    pretty_float,
    pretty_int,
    resolve_inputs,
)


WIDE_COLUMNS = [
    "Strain",
    "Strain name",
    "Genome size (bp)",
    "GC content (%)",
    "No. of contigs",
    "N50 (bp)",
    "N90 (bp)",
    "L50",
    "L90",
    "Longest contig (bp)",
    "Gaps (N per 100 kb)",
    "Sequencing coverage (×)",
    "Sequencing platforms",
    "Assembly method",
    "CDSs",
    "tRNAs",
    "rRNAs",
    "tmRNA",
    "Repeat regions",
    "16S rRNA count",
    "16S rRNA length (bp)",
    "16S rRNA contig",
    "16S rRNA sequence",
    "Completeness (%)",
    "Contamination (%)",
]

FEATURE_ORDER = [
    "Strain name",
    "Genome size (bp)",
    "GC content (%)",
    "No. of contigs",
    "N50 (bp)",
    "N90 (bp)",
    "L50",
    "L90",
    "Longest contig (bp)",
    "Gaps (N per 100 kb)",
    "Sequencing coverage (×)",
    "Sequencing platforms",
    "Assembly method",
    "CDSs",
    "tRNAs",
    "rRNAs",
    "tmRNA",
    "Repeat regions",
    "16S rRNA count",
    "16S rRNA length (bp)",
    "16S rRNA contig",
    "16S rRNA sequence",
    "Completeness (%)",
    "Contamination (%)",
]


LEGACY_SUBCOMMANDS = {"summarize", "pipeline"}


def resolve_assemblies(paths: Sequence[str]) -> List[Path]:
    found = [p for p in resolve_inputs(paths, kind="fasta") if looks_like_fasta(p)]
    if not found:
        raise SystemExit("No assembly FASTA files found.")
    return found


def resolve_gffs(paths: Sequence[str] | None) -> List[Path]:
    if not paths:
        return []
    return [p for p in resolve_inputs(paths, kind="gff") if looks_like_gff(p)]


def _apply_optional(df: pd.DataFrame, other: Dict[str, dict]) -> pd.DataFrame:
    if not other:
        return df
    lookup = {canonicalize_name(key): val for key, val in other.items()}
    for idx, row in df.iterrows():
        can = canonicalize_name(str(row["Strain"]))
        rec = lookup.get(can)
        if not rec:
            continue
        for col, value in rec.items():
            if col not in df.columns:
                df[col] = pd.NA
            if value is not None and not (isinstance(value, str) and value == ""):
                df.at[idx, col] = value
    return df


def build_wide_dataframe(
    assemblies: Sequence[Path],
    gffs: Sequence[Path] | None = None,
    check_report_path: Path | None = None,
    coverage_path: Path | None = None,
    metadata_path: Path | None = None,
    sequencing_platforms: str | None = None,
    assembly_method: str | None = None,
) -> pd.DataFrame:
    assembly_rows = [assembly_stats_to_row(compute_assembly_stats(p)) for p in assemblies]
    df = pd.DataFrame(assembly_rows)
    if df.empty:
        return pd.DataFrame(columns=WIDE_COLUMNS)

    df["Strain name"] = df["Strain"]

    assembly_map = {
        canonicalize_name(str(row["Strain"])): Path(str(row["_assembly_path"]))
        for _, row in df.iterrows()
    }
    assembly_size_map = {
        canonicalize_name(str(row["Strain"])): int(row["Genome size (bp)"])
        for _, row in df.iterrows()
    }

    if gffs:
        gff_rows = []
        for gff in gffs:
            strain_key = canonicalize_name(Path(gff).stem)
            asm = assembly_map.get(strain_key)
            stats = parse_gff_stats(gff, assembly_path=asm) if asm is not None else parse_gff_stats(gff)
            gff_rows.append(gff_stats_to_row(stats))
        gff_map = {
            canonicalize_name(r["Strain"]): {k: v for k, v in r.items() if k != "Strain"}
            for r in gff_rows
        }
        df = _apply_optional(df, gff_map)

    if check_report_path:
        df = _apply_optional(df, parse_checkm2_report(check_report_path))

    if coverage_path:
        df = _apply_optional(df, parse_coverage_table(coverage_path, assembly_size_map=assembly_size_map))

    if metadata_path:
        df = _apply_optional(df, parse_metadata_table(metadata_path))

    if sequencing_platforms:
        df["Sequencing platforms"] = sequencing_platforms
    if assembly_method:
        df["Assembly method"] = assembly_method

    for col in WIDE_COLUMNS:
        if col not in df.columns:
            df[col] = pd.NA

    df = df[WIDE_COLUMNS].sort_values("Strain").reset_index(drop=True)
    return df


def _format_feature_value(feature: str, value) -> str:
    if pd.isna(value):
        return ""

    if feature == "16S rRNA sequence":
        return str(value)

    if feature in {
        "Genome size (bp)", "No. of contigs", "N50 (bp)", "N90 (bp)", "L50", "L90",
        "Longest contig (bp)", "CDSs", "tRNAs", "rRNAs", "tmRNA",
        "Repeat regions", "16S rRNA count", "16S rRNA length (bp)",
    }:
        try:
            return pretty_int(float(value)) or ""
        except Exception:
            return str(value)

    if feature in {"GC content (%)", "Gaps (N per 100 kb)", "Sequencing coverage (×)", "Completeness (%)", "Contamination (%)"}:
        try:
            return pretty_float(float(value), digits=2) or ""
        except Exception:
            return str(value)

    return str(value)


def build_feature_dataframe(wide_df: pd.DataFrame) -> pd.DataFrame:
    rows: List[dict] = []
    if wide_df.empty:
        return pd.DataFrame(columns=["Strain", "Feature", "Description"])

    for _, rec in wide_df.iterrows():
        strain = rec["Strain"]
        for feature in FEATURE_ORDER:
            rows.append(
                {
                    "Strain": strain,
                    "Feature": feature,
                    "Description": _format_feature_value(feature, rec.get(feature)),
                }
            )
    return pd.DataFrame(rows, columns=["Strain", "Feature", "Description"])


def write_table(df: pd.DataFrame, path: Path) -> None:
    path = ensure_parent(path)
    if path.suffix.lower() == ".csv":
        df.to_csv(path, index=False)
    else:
        df.to_csv(path, sep="	", index=False)


def write_outputs(
    wide_df: pd.DataFrame,
    output: Path,
    feature_df: pd.DataFrame | None = None,
    feature_output: Path | None = None,
    xlsx_path: Path | None = None,
) -> None:
    write_table(wide_df, output)

    if feature_output and feature_df is not None:
        write_table(feature_df, feature_output)

    if xlsx_path:
        xlsx_path = ensure_parent(xlsx_path)
        with pd.ExcelWriter(xlsx_path) as writer:
            wide_df.to_excel(writer, index=False, sheet_name="wide_table")
            if feature_df is not None:
                feature_df.to_excel(writer, index=False, sheet_name="feature_table")


def _run_summary(
    assemblies: Sequence[Path],
    gffs: Sequence[Path] | None,
    check_report_path: Path | None,
    coverage_path: Path | None,
    metadata_path: Path | None,
    sequencing_platforms: str | None,
    assembly_method: str | None,
    output: Path,
    feature_output: Path | None,
    xlsx_path: Path | None,
) -> int:
    wide_df = build_wide_dataframe(
        assemblies=assemblies,
        gffs=gffs,
        check_report_path=check_report_path,
        coverage_path=coverage_path,
        metadata_path=metadata_path,
        sequencing_platforms=sequencing_platforms,
        assembly_method=assembly_method,
    )
    feature_df = None
    if feature_output or xlsx_path:
        feature_df = build_feature_dataframe(wide_df)

    write_outputs(
        wide_df=wide_df,
        output=output,
        feature_df=feature_df,
        feature_output=feature_output,
        xlsx_path=xlsx_path,
    )
    print(f"Wrote wide table: {output}")
    if feature_output:
        print(f"Wrote feature table: {feature_output}")
    if xlsx_path:
        print(f"Wrote workbook: {xlsx_path}")
    return 0


def validate_args(args: argparse.Namespace) -> None:
    if args.check and args.check_report:
        raise SystemExit("Use either --check or --check-report, not both.")

    if args.annotate in {"prokka", "bakta"} and args.gff:
        raise SystemExit(
            "--gff is only used for existing annotation files. Remove --gff or use --annotate existing/none."
        )

    if args.annotate not in {"prokka", "bakta"} and args.annotate_args:
        raise SystemExit("--annotate-args can only be used with --annotate prokka or --annotate bakta.")

    if args.annotate != "prokka" and args.kingdom != "auto":
        raise SystemExit("--kingdom is only used with --annotate prokka.")


def run_workflow(args: argparse.Namespace) -> int:
    validate_args(args)

    assemblies = resolve_assemblies(args.input)
    explicit_gffs = resolve_gffs(args.gff)

    if args.annotate == "prokka":
        gffs = run_prokka(
            assemblies=assemblies,
            outdir=Path(args.workdir).resolve() / "prokka",
            threads=args.threads,
            kingdom=args.kingdom,
            extra_args=args.annotate_args,
            verbose=not args.quiet,
        )
    elif args.annotate == "bakta":
        gffs = run_bakta(
            assemblies=assemblies,
            outdir=Path(args.workdir).resolve() / "bakta",
            threads=args.threads,
            extra_args=args.annotate_args,
            verbose=not args.quiet,
        )
    elif args.annotate == "existing":
        gffs = find_existing_gffs(assemblies=assemblies, gff_inputs=explicit_gffs)
    else:
        gffs = explicit_gffs

    if args.check:
        check_report_path = run_checkm2(
            assemblies=assemblies,
            outdir=Path(args.workdir).resolve() / "checkm2",
            threads=args.threads,
            verbose=not args.quiet,
        )
    else:
        check_report_path = Path(args.check_report).resolve() if args.check_report else None

    coverage_path = Path(args.coverage).resolve() if args.coverage else None
    metadata_path = Path(args.metadata).resolve() if args.metadata else None

    return _run_summary(
        assemblies=assemblies,
        gffs=gffs or None,
        check_report_path=check_report_path,
        coverage_path=coverage_path,
        metadata_path=metadata_path,
        sequencing_platforms=args.sequencing_platforms,
        assembly_method=args.assembly_method,
        output=Path(args.output),
        feature_output=Path(args.feature_output) if args.feature_output else None,
        xlsx_path=Path(args.xlsx) if args.xlsx else None,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genochar",
        description=(
            "Generate publication-ready genome characterization tables from assembly FASTA files, "
            "with optional annotation and CheckM2 execution."
        ),
    )
    parser.add_argument("-V", "--version", action="version", version=f"GenoChar {__version__}")

    primary = parser.add_argument_group("Primary inputs")
    primary.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="Assembly FASTA files, directories, or glob patterns.",
    )
    primary.add_argument(
        "--gff",
        nargs="+",
        help="Existing GFF/GFF3 files, directories, or glob patterns.",
    )
    primary.add_argument(
        "--check-report",
        help="Existing CheckM2 quality_report.tsv file.",
    )
    primary.add_argument("--coverage", help="Coverage TSV/CSV path.")
    primary.add_argument("--metadata", help="Metadata TSV/CSV path.")

    workflow = parser.add_argument_group("Optional workflow")
    workflow.add_argument(
        "--annotate",
        choices=["none", "prokka", "bakta", "existing"],
        default="none",
        help=(
            "How to obtain GFF files: none (default), existing, prokka, or bakta. "
            "Use 'existing' to reuse nearby GFFs or explicitly supplied --gff files."
        ),
    )
    workflow.add_argument(
        "--annotate-args",
        help="Extra arguments passed to Prokka or Bakta as a single quoted string.",
    )
    workflow.add_argument(
        "--check",
        action="store_true",
        help="Run CheckM2 on the input assemblies and use the resulting quality_report.tsv.",
    )

    execution = parser.add_argument_group("Execution controls")
    execution.add_argument(
        "-k",
        "--kingdom",
        choices=["auto", "Bacteria", "Archaea"],
        default="auto",
        help="Passed to Prokka when --annotate prokka is used.",
    )
    execution.add_argument(
        "-t",
        "--threads",
        type=int,
        default=8,
        help="Threads used for Prokka, Bakta, and CheckM2.",
    )
    execution.add_argument(
        "-w",
        "--workdir",
        default="genochar_work",
        help="Working directory for annotation and CheckM2 outputs.",
    )
    execution.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Reduce logging from external tools.",
    )

    metadata_group = parser.add_argument_group("Optional fixed metadata")
    metadata_group.add_argument(
        "--sequencing-platforms",
        help="Set one sequencing platform string for all strains.",
    )
    metadata_group.add_argument(
        "--assembly-method",
        help="Set one assembly method string for all strains.",
    )

    outputs = parser.add_argument_group("Outputs")
    outputs.add_argument(
        "-o",
        "--output",
        default="genome_characterization.tsv",
        help="Main wide-table output (.tsv or .csv).",
    )
    outputs.add_argument(
        "-f",
        "--feature-output",
        help="Optional feature-table output (.tsv or .csv).",
    )
    outputs.add_argument(
        "-x",
        "--xlsx",
        help="Optional XLSX workbook with wide_table and feature_table sheets.",
    )

    return parser


def _normalize_legacy_argv(argv: Sequence[str]) -> List[str]:
    items = list(argv)
    if items and items[0] in LEGACY_SUBCOMMANDS:
        print(
            "Note: legacy subcommands are deprecated in v0.6.1. Use 'genochar' directly.",
            file=sys.stderr,
        )
        return items[1:]
    return items


def main(argv: Sequence[str] | None = None) -> int:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    normalized_argv = _normalize_legacy_argv(raw_argv)
    parser = build_parser()
    args = parser.parse_args(normalized_argv)
    try:
        return run_workflow(args)
    except PipelineError as exc:
        raise SystemExit(str(exc))


if __name__ == "__main__":
    raise SystemExit(main())
