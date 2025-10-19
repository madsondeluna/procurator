# Test Directory

This folder contains test data, test scripts, and results from testing the Procurator pipeline.

## Contents

- `generate_test_genome.py` - Script to generate synthetic prokaryotic test genome
- `*.csv`, `*.gff`, `*.faa` - Generated output files from running tests

## Quick Start

```bash
# Generate test genome
python test/generate_test_genome.py

# Run stats on test data
python run_procurator.py stats -i data/test_genome.fasta -o test/stats.csv

# Run find_orfs
python run_procurator.py find_orfs \
    -i data/test_genome.fasta \
    -o_gff test/genes.gff \
    -o_faa test/proteins.faa
```

## Test Data

The test genome contains 3 synthetic contigs totaling 12,000 bp with embedded genes for realistic testing.

## Cleaning Up

To remove test files:

```bash
rm -f test/*.csv test/*.gff test/*.faa test/*.tsv
```

Note: `.gitkeep` and this README should remain in the directory.
