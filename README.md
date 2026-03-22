# GenoChar

**GenoChar** generates publication-ready genome characterization tables for bacterial and archaeal genomes.

Version **0.6.3.2** keeps the **single-command, summarize-first interface**, but adds a practical solution for external-tool dependency conflicts:

- **`genochar`** builds the final genome characterization table
- **`genochar setup`** prepares **managed Prokka and CheckM2 conda environments** under `~/.genochar`
- later, a normal `genochar ... --check --annotate prokka` run automatically calls those managed environments with `conda run -p ...`

That means users do **not** need to manually solve a shared Prokka/CheckM2 environment.

The main output is a wide table with one row per strain. Optional outputs include a feature-style table and an Excel workbook.

**Naming note:** the project/display name is **GenoChar**, while the Python package name and command-line executable remain lowercase as `genochar`.

## What the tool computes

### FASTA-derived fields
These work even if you provide only FASTA files:

- `Strain`
- `Strain name`
- `Genus`
- `Species`
- `Accession`
- `Genome size (bp)`
- `GC content (%)`
- `No. of contigs`
- `N50 (bp)`
- `N90 (bp)`
- `L50`
- `L90`
- `Longest contig (bp)`
- `Gaps (N per 100 kb)`

### GFF-derived fields
These are added when you provide GFF files or ask GenoChar to create annotation files:

- `CDSs`
- `tRNAs`
- `rRNAs`
- `tmRNA`
- `misc RNA`
- `Repeat regions`
- `16S rRNA count`
- `16S rRNA length (bp)`
- `16S rRNA contig`
- `16S rRNA sequence`

### CheckM2-derived fields
These are added when you provide or generate a CheckM2 report:

- `Completeness (%)`
- `Contamination (%)`

### User-supplied metadata
Optional input tables can add:

- `Sequencing coverage (×)`
- `Sequencing platforms`
- `Assembly method`
- `Genus`
- `Species`
- `Accession`
- `Repeat regions`

## Default output columns

The main output table contains:

- `Strain`
- `Strain name`
- `Genus`
- `Species`
- `Accession`
- `Genome size (bp)`
- `GC content (%)`
- `No. of contigs`
- `N50 (bp)`
- `N90 (bp)`
- `L50`
- `L90`
- `Longest contig (bp)`
- `Gaps (N per 100 kb)`
- `Sequencing coverage (×)`
- `Sequencing platforms`
- `Assembly method`
- `CDSs`
- `tRNAs`
- `rRNAs`
- `tmRNA`
- `misc RNA`
- `Repeat regions`
- `16S rRNA count`
- `16S rRNA length (bp)`
- `16S rRNA contig`
- `16S rRNA sequence`
- `Completeness (%)`
- `Contamination (%)`

## Installation

Clone the repository and install GenoChar into the active Python environment:

```bash
git clone https://github.com/ljunwon1114/GenoChar.git
cd GenoChar
pip install -e .
```

If you just want the latest GitHub version without cloning for development, you can also use:

```bash
pip install git+https://github.com/ljunwon1114/GenoChar.git
```

### Core Python requirement

GenoChar itself is lightweight. The **core package** now declares:

- `requires-python >=3.10`

The external tools are the difficult part, so v0.6.3.2 no longer assumes that Prokka and CheckM2 should share the same environment.

## One-time managed setup (recommended)

Run this once:

```bash
genochar setup
```

This creates managed environments under `~/.genochar`, typically:

```text
~/.genochar/
  config.json
  envs/
    prokka/
    checkm2/
  databases/
    CheckM2_database/
```

After that, normal workflow commands automatically use those managed environments when `--annotate prokka` and/or `--check` are requested.

When `--check` is used, GenoChar now passes the resolved input FASTA files directly to `checkm2 predict --input ...`, matching the official CheckM2 interface that accepts either a folder of bins or a list of FASTA files.

### Reuse an existing CheckM2 database

If you already downloaded the CheckM2 database, point setup at it directly:

```bash
genochar setup --checkm2-db /home/jwlee/databases/CheckM2_database/uniref100.KO.1.dmnd
```

You can also pass a directory that contains the `.dmnd` file.

### Optional setup flags

```bash
genochar setup --skip-prokka
genochar setup --skip-checkm2
genochar setup --force
```

- `--skip-prokka`: only prepare CheckM2
- `--skip-checkm2`: only prepare Prokka
- `--force`: recreate managed environments even if they already exist

## Command overview

### A. FASTA only

```bash
genochar -i "assemblies/*.fasta" -o genome_characterization.tsv
```

### B. FASTA + existing GFF + existing CheckM2 report

```bash
genochar -i "assemblies/*.fasta" --gff "annotations/*.gff*" --check-report checkm2_out/quality_report.tsv -o genome_characterization.tsv
```

### C. FASTA + managed CheckM2 first + managed Prokka annotation

```bash
genochar -i "assemblies/*.fasta" --check --annotate prokka -k Archaea -t 8 -w genochar_work -o genome_characterization.tsv
```

### D. Reuse existing GFF files automatically

```bash
genochar -i "assemblies/*.fasta" --annotate existing --check-report checkm2_out/quality_report.tsv -o genome_characterization.tsv
```

### E. Reuse explicitly supplied GFF files in existing-annotation mode

```bash
genochar -i "assemblies/*.fasta" --annotate existing --gff "annotations/*.gff*" --check-report checkm2_out/quality_report.tsv -o genome_characterization.tsv
```

## Optional extra outputs

### Feature-style table

```bash
genochar -i "assemblies/*.fasta" -f genome_characterization_feature.tsv -o genome_characterization.tsv
```

### Excel workbook

```bash
genochar -i "assemblies/*.fasta" -x genome_characterization.xlsx -o genome_characterization.tsv
```

## Coverage input

Coverage cannot be derived from FASTA alone. If you want to fill `Sequencing coverage (×)`, provide a coverage table.

Example:

```text
Strain	Coverage
IOH03	55.7
IOH05	50.3
```

or

```text
Strain	Total bases
IOH03	110.8 Mbp
IOH05	107.6 Mbp
```

If `Total bases` is provided, GenoChar computes:

```text
Sequencing coverage (×) = Total bases / Genome size
```

## Metadata input

Optional metadata columns include:

- `Strain`
- `Sequencing platforms`
- `Assembly method`
- `Genus`
- `Species`
- `Accession`
- `Repeat regions`
- `Sequencing coverage (×)`

Example:

```text
Strain	Genus	Species	Accession	Sequencing platforms	Assembly method
IOH03	Thermococcus	waiotapuensis	GCF_032304395	Illumina iSeq 100	Unicycler (short-read assembly)
IOH05	Thermococcus	sp.	GCA_000000000	Illumina iSeq 100	Unicycler (short-read assembly)
```

## Notes

- **GenoChar** is summarize-first by default. If you only pass FASTA, GFF, CheckM2, coverage, and metadata inputs, it behaves like a direct summarization tool.
- **`genochar setup`** is the recommended way to prepare Prokka and CheckM2 without forcing them into one shared environment.
- `--annotate prokka` tells GenoChar to create annotation files before building the final table.
- `--annotate existing` tells GenoChar to reuse nearby GFF files or explicitly supplied `--gff` inputs.
- `--check` runs CheckM2 internally **before annotation** and automatically integrates the resulting `quality_report.tsv` into the final table.
- `--check-report` reuses an existing CheckM2 `quality_report.tsv` file.
- `--check` and `--check-report` are mutually exclusive.
- `--gff` is intended for existing annotation files and should not be combined with `--annotate prokka`.
- If more than one 16S rRNA feature is found, GenoChar stores the longest 16S sequence in the main table.
- Legacy `genochar summarize ...` and `genochar pipeline ...` calls are still accepted in v0.6.3, but the preferred interface is the single-command form shown above.

## License

MIT
