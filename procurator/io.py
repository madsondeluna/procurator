# procurator/io.py

from Bio import SeqIO
from BCBio import GFF
import pandas as pd
import sys

def read_fasta_sequences(filepath, as_dict=False):
    """Lê um arquivo FASTA e retorna um iterador ou dicionário de SeqRecords."""
    try:
        if as_dict:
            return SeqIO.to_dict(SeqIO.parse(filepath, "fasta"))
        else:
            return list(SeqIO.parse(filepath, "fasta"))
    except FileNotFoundError:
        print(f"Erro: Arquivo FASTA de entrada '{filepath}' não encontrado.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Erro ao ler o arquivo FASTA '{filepath}': {e}", file=sys.stderr)
        sys.exit(1)

def write_fasta(records, filepath):
    """Escreve uma lista de SeqRecords em um arquivo FASTA."""
    try:
        SeqIO.write(records, filepath, "fasta")
    except Exception as e:
        print(f"Erro ao escrever o arquivo FASTA '{filepath}': {e}", file=sys.stderr)
        sys.exit(1)

def write_gff(records, filepath):
    """Escreve uma lista de SeqRecords (com features) em um arquivo GFF3."""
    try:
        with open(filepath, "w") as out_handle:
            GFF.write(records, out_handle)
    except Exception as e:
        print(f"Erro ao escrever o arquivo GFF '{filepath}': {e}", file=sys.stderr)
        sys.exit(1)

def write_dataframe_to_csv(df, filepath):
    """Escreve um DataFrame do pandas em um arquivo CSV."""
    try:
        df.to_csv(filepath, index=False)
    except Exception as e:
        print(f"Erro ao escrever o arquivo CSV '{filepath}': {e}", file=sys.stderr)
        sys.exit(1)