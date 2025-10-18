# Procurator Pipeline

A bioinformatics pipeline for protein analysis and annotation.

## Structure

```
procurator-pipeline/
├── .gitignore
├── README.md
├── environment.yml
├── run_procurator.py       # Main entry point
├── db/                     # Database files (Pfam, etc.)
│   └── Pfam-A.hmm
└── procurator/             # Main package
    ├── __init__.py
    ├── main.py             # Core pipeline logic
    ├── io.py               # Input/output utilities
    ├── analysis.py         # Analysis functions
    ├── orf_finding.py      # ORF detection
    ├── annotation.py       # Sequence annotation
    └── protein.py          # Protein class and utilities
```

## Installation

Create a conda environment using the provided `environment.yml`:

```bash
conda env create -f environment.yml
conda activate procurator
```

## Usage

Run the pipeline:

```bash
python run_procurator.py <input_fasta> [options]
```

## Requirements

- Python 3.8+
- HMMER for protein domain detection
- Dependencies listed in `environment.yml`

## Output

The pipeline generates:
- Annotated sequence files
- ORF predictions
- Domain annotations
