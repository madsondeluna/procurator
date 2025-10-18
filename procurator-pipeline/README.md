# Procurator: Prokaryotic Genomic Analysis and Annotation Pipeline

Procurator is a modular, command-line orchestrated (CLI) bioinformatics pipeline for comprehensive prokaryotic genome analysis. Starting from a FASTA file of a sequenced genome, Procurator automates gene prediction steps, functional protein annotation, and preparation of targets for structural modeling.

The pipeline is designed to be robust, replacing manual implementations with industry-standard tools such as Prodigal, HMMER, and ColabFold, all accessible through a unified interface.

## Features

- **Statistical Analysis**: Computes essential genomic statistics (GC content, N50, contig count, etc.).

- **Gene Prediction**: Uses Prodigal (via pyrodigal-finder) for fast and accurate gene (ORF) prediction in prokaryotic genomes.

- **Functional Annotation**: Annotates predicted proteins against domain databases (e.g., Pfam) using HMMER (via pyhmmer).

- **Target Selection**: Filters and selects proteins of interest based on annotation keywords, preparing them for downstream analyses.

- **Structural Modeling**: Integrates with ColabFold to perform de novo protein structure modeling directly from the command line.

- **Modularity**: Each step is an independent sub-command, allowing flexibility in execution and integration with workflow managers like Snakemake or Nextflow.

- **Dependency Management**: Handles interactive installation of critical dependencies like ColabFold.

## Installation and Setup

### Prerequisites

- **Conda/Mamba**: It is highly recommended to manage the environment with Conda or, preferably, Mamba for faster dependency resolution.
- **Git**: To clone the repository.
- **System Requirements**: At least 4GB RAM for basic analysis; 20GB+ free disk space recommended (especially for ColabFold database setup).

### 1. Clone the Repository

```bash
git clone https://github.com/madsondeluna/procurator.git
cd procurator/procurator-pipeline
```

### 2. Create Conda Environment

The `environment.yml` file contains all necessary Python dependencies. Create and activate the environment with the following command:

```bash
# Recommended: use mamba for speed
mamba env create -f environment.yml

# Alternative: use conda
# conda env create -f environment.yml

# Activate the environment
conda activate procurator
```

### 3. ColabFold Setup (Critical Manual Step)

Protein modeling requires extensive databases. ColabFold can use a remote server (slower, requires internet) or local databases (much faster, recommended).

This step should be executed manually ONCE. It will download hundreds of gigabytes.

```bash
# Run this command in the terminal with the 'procurator' environment activated
colabfold_setup databases

# Follow the instructions. Choose a directory with sufficient disk space.
```

## Typical Workflow

The pipeline is designed to be executed in a logical sequence. A complete workflow, from raw genome to protein models, would be:

1. **stats**: (Optional) Assess the quality and characteristics of the genome.
2. **find_orfs**: Predict all genes/proteins from the input genome.
3. **annotate**: Functionally annotate the predicted proteins to understand what they do.
4. **select_proteins**: Select a subset of proteins of interest (e.g., all kinases) based on annotations.
5. **run_modeling**: Generate 3D structural models for selected proteins.

## Command-Line Usage (CLI)

The pipeline is executed through the `run_procurator.py` script, followed by a specific sub-command for each task.

### stats

Computes basic statistics from a genomic FASTA file.

**Usage:**

```bash
python run_procurator.py stats -i <genome.fna> -o <statistics.csv>
```

**Example:**

```bash
python run_procurator.py stats -i data/bac-seqs.fasta -o genome_stats.csv
```

**Output columns:**

- Total_Sequencias: Number of sequences in FASTA file
- Total_Bases_bp: Total base pairs
- Comprimento_Medio_bp: Average sequence length
- Comprimento_Max_bp: Maximum sequence length
- Comprimento_Min_bp: Minimum sequence length
- N50_bp: N50 statistic (weighted sequence length)
- GC_Medio_Perc: Average GC percentage

### find_orfs

Finds genes/ORFs using Prodigal.

**Outputs:**

- A GFF3 file with gene coordinates.
- A FASTA file (.faa) with translated protein sequences.

**Usage:**

```bash
python run_procurator.py find_orfs \
    -i <genome.fna> \
    -o_gff <genes.gff> \
    -o_faa <proteins.faa> \
    -t 11
```

**Arguments:**

- `-i, --input`: Input genome FASTA file (required)
- `-o_gff, --output_gff`: Output GFF3 file with gene coordinates (required)
- `-o_faa, --output_protein`: Output FASTA file with protein sequences (required)
- `-t, --translation_table`: Genetic code table (default: 11 for prokaryotes)

**Example:**

```bash
python run_procurator.py find_orfs \
    -i data/bac-seqs.fasta \
    -o_gff predicted_genes.gff \
    -o_faa predicted_proteins.faa \
    -t 11
```

### annotate

Annotates proteins against an HMM profile database (e.g., Pfam).

**Usage:**

```bash
python run_procurator.py annotate \
    -i <proteins.faa> \
    -db <path/to/Pfam-A.hmm> \
    -o <annotations.tsv> \
    --evalue 1e-5
```

**Arguments:**

- `-i, --input_faa`: Input protein FASTA file (required)
- `-db, --hmm_database`: Path to HMM database file (e.g., Pfam-A.hmm) (required)
- `-o, --output_tsv`: Output TSV file with HMMER hits (required)
- `--evalue`: Maximum e-value threshold (default: 1e-5)

**Example:**

```bash
# Prerequisite: Download and prepare Pfam-A.hmm
# wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
# gunzip Pfam-A.hmm.gz
# hmmpress Pfam-A.hmm

python run_procurator.py annotate \
    -i predicted_proteins.faa \
    -db ./db/Pfam-A.hmm \
    -o protein_annotations.tsv \
    --evalue 1e-5
```

**Output columns (TSV):**

- Query_ID: Protein identifier
- Target_HMM: HMM/Domain name
- Accession: Domain accession number
- Description: Domain description
- E-value: Statistical significance
- Score: HMMER bit score

### select_proteins

Selects proteins from a FASTA file based on keywords found in the annotation file.

**Usage:**

```bash
python run_procurator.py select_proteins \
    -i_faa <proteins.faa> \
    -i_ann <annotations.tsv> \
    -k "keyword1" "keyword2" \
    -o <selected_proteins.faa>
```

**Arguments:**

- `-i_faa, --input_faa`: Input protein FASTA file (required)
- `-i_ann, --input_annotations`: Input annotations TSV file from `annotate` command (required)
- `-k, --keywords`: Keywords to search in annotations; use multiple keywords separated by spaces (required)
- `-o, --output_fasta`: Output FASTA file with selected proteins (required)

**Example (select all kinases and transferases):**

```bash
python run_procurator.py select_proteins \
    -i_faa predicted_proteins.faa \
    -i_ann protein_annotations.tsv \
    -k "kinase" "transferase" \
    -o selected_kinases_transferases.faa
```

### run_modeling

Executes structural modeling using ColabFold for a given set of proteins.

**Usage:**

```bash
python run_procurator.py run_modeling \
    -i <selected_proteins.faa> \
    -o <output_directory> \
    --model_type multimer \
    --num_models 5
```

**Arguments:**

- `-i, --input_fasta`: Input protein FASTA file (required)
- `-o, --output_dir`: Output directory for results (required)
- `--model_type`: ColabFold model type - "multimer" (default), "monomer", or "monomer_v2"
- `--num_models`: Number of models to predict (default: 5)

**Example:**

```bash
python run_procurator.py run_modeling \
    -i selected_kinases_transferases.faa \
    -o ./structural_models \
    --model_type multimer \
    --num_models 5
```

**Output files (per protein):**

- `{protein_id}_relaxed_rank_1_model_1.pdb`: Final 3D structure in PDB format
- `{protein_id}_pae.json`: Predicted aligned error (confidence metric)
- `{protein_id}_scores.json`: Raw prediction scores
- `{protein_id}.log`: ColabFold execution log

**Important Notes:**

- ⚠️ **Database Setup**: For best performance and offline capability, run `colabfold_setup databases` once. This downloads ~100GB and takes several hours.
- If local databases are not found, ColabFold will attempt to use remote MMseqs2 server (requires internet, slower).
- ColabFold will be installed automatically if missing when this command is first run.
- Modeling time depends on protein size and available GPU. Large proteins may take hours.

## Project Structure

```
procurator-pipeline/
├── db/                     # Directory for databases (e.g., Pfam HMMs)
├── environment.yml         # Conda environment file with dependencies
├── README.md               # This file
├── run_procurator.py       # Executable CLI entry point
└── procurator/             # Python package source code
    ├── __init__.py
    ├── main.py             # Main argparse logic and sub-commands
    ├── io.py               # Functions to read/write files
    ├── analysis.py         # Module for 'stats' sub-command
    ├── orf_finding.py      # Module for 'find_orfs' sub-command
    ├── annotation.py       # Module for 'annotate' sub-command
    └── protein.py          # Modules for 'select_proteins' and 'run_modeling'
```

---

## Troubleshooting

### Common Issues

#### `command not found: python run_procurator.py`

- Make sure you are in the `procurator-pipeline/` directory
- Check that the Conda environment is activated: `conda activate procurator`
- Try using `python3` instead of `python`

#### `ModuleNotFoundError: No module named 'procurator'`

- Ensure the environment is activated: `conda activate procurator`
- Verify you are running the script from the correct directory

#### HMMER/Pfam database not found

- Download Pfam database:

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
mv Pfam-A.hmm procurator-pipeline/db/
```

#### ColabFold runs very slowly or fails

- Run `colabfold_setup databases` to download local databases (100GB, ~2-4 hours)
- Ensure sufficient free disk space (20GB minimum for intermediate files)
- Check internet connection if using remote server fallback

#### Out of Memory (OOM) errors

- Process fewer proteins at once
- Use a machine with more RAM (16GB+ recommended for large proteins)
- Split input FASTA into smaller batches

## System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| Python | 3.9 | 3.10+ |
| RAM | 4GB | 16GB+ (32GB for modeling) |
| Disk | 20GB | 100GB+ (if using ColabFold databases) |
| GPU | Optional | NVIDIA GPU with CUDA (10-100x speedup) |

## Dependencies

All dependencies are managed by Conda via `environment.yml`:

- **Python 3.9+**: Core language
- **Biopython**: Sequence I/O and manipulation
- **Pandas**: Data frame processing for statistics and annotations
- **Pyrodigal-finder**: Python wrapper for Prodigal gene prediction
- **Pyhmmer**: Python wrapper for HMMER protein domain detection
- **BCBio-GFF**: GFF3 file reading/writing

Optional (installed on demand):

- **ColabFold**: Protein structure prediction (installed automatically for `run_modeling`)

## Workflow Example: Complete E. coli Analysis

```bash
# 1. Activate environment
conda activate procurator

# 2. Compute genome statistics
python run_procurator.py stats \
    -i data/bac-seqs.fasta \
    -o ecoli_stats.csv

# 3. Predict all genes/proteins
python run_procurator.py find_orfs \
    -i data/bac-seqs.fasta \
    -o_gff ecoli_genes.gff \
    -o_faa ecoli_proteome.faa

# 4. Annotate with Pfam domains (download Pfam first if needed)
python run_procurator.py annotate \
    -i ecoli_proteome.faa \
    -db db/Pfam-A.hmm \
    -o ecoli_pfam_hits.tsv

# 5. Select proteins of interest
python run_procurator.py select_proteins \
    -i_faa ecoli_proteome.faa \
    -i_ann ecoli_pfam_hits.tsv \
    -k "kinase" "peptidase" "transport" \
    -o ecoli_selected_proteins.faa

# 6. Model selected protein structures
python run_procurator.py run_modeling \
    -i ecoli_selected_proteins.faa \
    -o ./ecoli_structures
```

## License

This project is provided as-is for research and bioinformatics applications.

## Citation

If you use Procurator in your research, please cite the tools it depends on:

- **Prodigal**: Hyatt et al. (2010)
- **HMMER**: Eddy (2011)
- **ColabFold**: Mirdita et al. (2022)
- **Biopython**: Cock et al. (2009)

## Contributing

Contributions are welcome! Please submit issues and pull requests to improve the pipeline.

## Contact

For questions or support, please open an issue on the repository:
[github.com/madsondeluna/procurator](https://github.com/madsondeluna/procurator)
