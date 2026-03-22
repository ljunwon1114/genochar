from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

from .utils import infer_strain_name, open_maybe_gzip, read_fasta_dict, reverse_complement


@dataclass
class GFFStats:
    Strain: str
    CDSs: int
    tRNAs: int
    rRNAs: int
    tmRNA: int
    misc_rna: int
    repeat_regions: int
    rrna16s_count: int
    rrna16s_length_bp: int | None
    rrna16s_sequence: str | None
    rrna16s_contig: str | None


_RRNA16S_PATTERNS = (
    "16s",
    "16s rrna",
    "16s ribosomal rna",
    "small subunit ribosomal rna",
    "small subunit rrna",
    "ssu ribosomal rna",
    "ssu rrna",
)


def parse_gff_attributes(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    if not text or text == ".":
        return attrs
    for field in text.strip().split(";"):
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
        else:
            key, value = field, ""
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def _blob(feature_type: str, attrs: Dict[str, str]) -> str:
    return " ".join([feature_type, *attrs.keys(), *attrs.values()]).lower()


def _looks_like_16s(feature_type: str, attrs: Dict[str, str]) -> bool:
    feature_type_l = feature_type.lower()
    if feature_type_l not in {"rrna", "ribosomal_rna", "gene", "ncrna", "transcript"}:
        return False
    blob = _blob(feature_type, attrs)
    return any(pat in blob for pat in _RRNA16S_PATTERNS)


def _looks_like_tmrna(feature_type: str, attrs: Dict[str, str]) -> bool:
    feature_type_l = feature_type.lower()
    blob = _blob(feature_type, attrs)
    return (
        feature_type_l in {"tmrna", "tm_rna", "transfer_messenger_rna"}
        or "tmrna" in blob
        or "transfer-messenger rna" in blob
        or "ssra" in blob
    )


def _looks_like_misc_rna(feature_type: str, attrs: Dict[str, str]) -> bool:
    feature_type_l = feature_type.lower().replace("-", "_")
    blob = _blob(feature_type, attrs)

    if _looks_like_16s(feature_type, attrs) or _looks_like_tmrna(feature_type, attrs):
        return False

    if feature_type_l in {"misc_rna", "miscrna"}:
        return True
    if "misc_rna" in blob or "misc rna" in blob:
        return True

    if feature_type_l in {"ncrna", "ncrna_gene", "non_coding_rna", "ncrna_region"}:
        if any(pat in blob for pat in [
            "ribosomal rna",
            " rrna",
            "rrna ",
            "transfer rna",
            "trna",
            "tmrna",
            "transfer-messenger rna",
            "ssra",
            "16s",
        ]):
            return False
        return True

    return False


def _looks_like_repeat(feature_type: str, attrs: Dict[str, str]) -> bool:
    feature_type_l = feature_type.lower()
    blob = _blob(feature_type, attrs)
    return "repeat" in feature_type_l or "repeat" in blob


def _extract_feature_sequence(
    fasta_by_id: Dict[str, str],
    seqid: str,
    start: int,
    end: int,
    strand: str,
) -> str | None:
    seq = fasta_by_id.get(seqid)
    if seq is None:
        return None
    if start < 1 or end > len(seq) or start > end:
        return None
    subseq = seq[start - 1:end]
    if strand == "-":
        subseq = reverse_complement(subseq)
    return subseq.upper()


def parse_gff_stats(path: Path | str, assembly_path: Path | str | None = None) -> GFFStats:
    strain = infer_strain_name(path)
    cds = 0
    trna = 0
    rrna = 0
    tmrna = 0
    misc_rna = 0
    repeat_regions = 0
    rrna16s_records: List[Tuple[int, str | None, str | None]] = []

    fasta_by_id: Dict[str, str] = {}
    if assembly_path:
        fasta_by_id = read_fasta_dict(assembly_path)

    with open_maybe_gzip(path, "rt") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            seqid, _, feature, start_s, end_s, _, strand, _, attrs_text = parts[:9]
            feature_l = feature.strip().lower()
            attrs = parse_gff_attributes(attrs_text)

            if feature_l == "cds":
                cds += 1
            elif feature_l in {"trna", "transfer_rna"}:
                trna += 1
            elif feature_l in {"rrna", "ribosomal_rna"}:
                rrna += 1

            if _looks_like_tmrna(feature_l, attrs):
                tmrna += 1

            if _looks_like_misc_rna(feature_l, attrs):
                misc_rna += 1

            if _looks_like_repeat(feature_l, attrs):
                repeat_regions += 1

            if _looks_like_16s(feature_l, attrs):
                try:
                    start = int(start_s)
                    end = int(end_s)
                except ValueError:
                    continue
                seq = _extract_feature_sequence(fasta_by_id, seqid, start, end, strand)
                length = len(seq) if seq else (end - start + 1)
                rrna16s_records.append((length, seqid, seq))

    longest_len: int | None = None
    longest_seq: str | None = None
    longest_contig: str | None = None
    if rrna16s_records:
        rrna16s_records.sort(key=lambda x: x[0], reverse=True)
        longest_len, longest_contig, longest_seq = rrna16s_records[0]

    return GFFStats(
        Strain=strain,
        CDSs=cds,
        tRNAs=trna,
        rRNAs=rrna,
        tmRNA=tmrna,
        misc_rna=misc_rna,
        repeat_regions=repeat_regions,
        rrna16s_count=len(rrna16s_records),
        rrna16s_length_bp=longest_len,
        rrna16s_sequence=longest_seq,
        rrna16s_contig=longest_contig,
    )


def gff_stats_to_row(stats: GFFStats) -> dict:
    return {
        "Strain": stats.Strain,
        "CDSs": stats.CDSs,
        "tRNAs": stats.tRNAs,
        "rRNAs": stats.rRNAs,
        "tmRNA": stats.tmRNA,
        "misc RNA": stats.misc_rna,
        "Repeat regions": stats.repeat_regions,
        "16S rRNA count": stats.rrna16s_count,
        "16S rRNA length (bp)": stats.rrna16s_length_bp,
        "16S rRNA contig": stats.rrna16s_contig,
        "16S rRNA sequence": stats.rrna16s_sequence,
    }
