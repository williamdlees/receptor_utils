# Copilot Instructions for receptor_utils

## Project Overview

This is a Python package for working with immunoglobulin/T-cell receptor (IG/TR) sequences, focused on sequence alignment, novel allele naming, and creating custom IgBlast databases. The project follows IMGT alignment standards and integrates with external databases like OGRDB.

## Architecture & Core Components

### Core Modules (`src/receptor_utils/`)
- **`simple_bio_seq.py`** - Simple wrappers around BioPython for FASTA/CSV I/O, sequence manipulation
- **`number_v.py`** - IMGT-gapping algorithms for V-gene sequences, conserved residue checking
- **`novel_allele_name.py`** - VDJbase-compliant novel allele naming with SNP notation
- **`aux_formats.py`** - File format converters for IgBlast (.ndm, .aux) and MiXCR from AIRR-C JSON/FASTA

### Command-Line Scripts (`src/scripts/`)
All scripts are entry points defined in `setup.cfg`. Key patterns:
- Use `argparse` with `get_parser()` function for consistency
- Support both FASTA and AIRR-C JSON input formats where applicable
- Handle species-specific alignment quirks (e.g., macaque insertions)

## Critical Data Formats & Standards

### IMGT Alignment Convention
- **Gapped sequences**: Use `.` for gaps (IMGT standard)
- **V-gene coordinates**: Default CDR positions at 79,114,166,195,313 (1-based)
- **Conserved residues**: C at position 23, W around position 41, specific motifs for quality checks
- **Macaque exception**: Contains inserted codons that break standard tools - use `fix_macaque_gaps`

### File Format Support
- **AIRR-C JSON**: Preferred for IgBlast file generation (contains precise coordinates)
- **FASTA**: Requires IMGT-gapped V-sequences for proper CDR delineation
- **IgBlast formats**: `.ndm` (V-gene regions), `.aux` (J-gene annotations)

### Sequence Naming Conventions
- Follow VDJbase scheme: `GENE*ALLELE_snp1_snp2_insertion` format
- Handle truncated sequences with ambiguous reference sets
- Support both gapped and ungapped input for V-genes

## Development Patterns

### Error Handling & Validation
- Always check for conserved residues in V-sequences (`number_v.check_conserved_residues()`)
- Validate file formats early (`.json` vs `.fasta` detection)
- Handle missing alleles gracefully (fallback to `*01` variants)

### Testing Approach
- Test files in `tests/` include reference data from IMGT and OGRDB
- Comparison scripts: `check_aux.py`, `check_ndm.py`, `check_fastas.py` for validation
- Expected mismatches documented in `download_germline_set_test_summaries.txt`

### External Integration Points
- **OGRDB API**: `https://ogrdb.airr-community.org/api_v2` for germline downloads
- **IMGT downloads**: Fixed URL pattern for reference sequences
- **IgBlast compatibility**: Strict coordinate requirements for downstream analysis

## Common Workflows

### Creating Custom IgBlast Databases
1. Download with `download_germline_set` (AIRR-C JSON preferred)
2. Generate `.ndm` with `make_igblast_ndm`
3. Generate `.aux` with `annotate_j`
4. Validate with checker scripts

### Novel Allele Analysis
1. Use `name_allele` with appropriate reference set (single locus only)
2. For V-genes: handles both gapped/ungapped input, outputs IMGT-gapped result
3. Reports SNPs, insertions, and sequence notes

### Sequence Processing Pipeline
1. Extract references: `extract_refs` → species-specific FASTA files
2. Gap sequences: `gap_sequences` using closest reference templates
3. Quality control: Check conserved residues, handle truncations

## Key Implementation Details

- **Deprecation pattern**: `number_ighv` → `number_v` with warning system
- **File I/O consistency**: `simple_bio_seq` functions use `(filename, data)` parameter order
- **Species handling**: Three-letter locus prefixes (IGH, TRB, etc.) for validation
- **Memory efficiency**: Use iterators, dict-based sequence storage as strings
- **Coordinate systems**: Mix of 0-based (internal) and 1-based (IMGT standard) - be explicit

## Other points to note
- Bioinformatic functions in `simple_bio_seq` should be used by preference. The functions simple_bio_seq.read_csv and simple_bio_seq.write_csv should be used for CSV I/O.
- The scripts in src/scripts/ should be used as command line entry points only. Their functions should be imported from the relevant modules in src/receptor_utils/ for use in other code.
- The scripts in src/scripts/ should use argparse for command line argument parsing, and should define a get_parser() function that returns an ArgumentParser object. This allows for consistent command line interfaces across scripts.
- A suitable conda environment is provided with the name digby_genomics310. All testing should be done in this environment. DO not use pip to install packages into this environment or create other venvs or conda envs.