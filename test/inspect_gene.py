#!/usr/bin/env python3
import pyrodigal
from Bio import SeqIO

record = next(SeqIO.parse("data/test_genome.fasta", "fasta"))
finder = pyrodigal.GeneFinder(meta=True)
seq_bytes = bytes(str(record.seq), 'utf-8')
genes = finder.find_genes(seq_bytes)

# Ver atributos do primeiro gene
if len(genes) > 0:
    gene = genes[0]
    print(f"Type: {type(gene)}")
    print(f"Dir: {[x for x in dir(gene) if not x.startswith('_')]}")
    print(f"Repr: {repr(gene)}")
    print(f"Str: {str(gene)}")
