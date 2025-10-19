#!/usr/bin/env python3
import pyrodigal
from Bio import SeqIO

# Ler sequência
record = next(SeqIO.parse("data/test_genome.fasta", "fasta"))
print(f"Testando com: {record.id} ({len(record.seq)} bp)")

# Criar finder em modo metagenomic (pré-treinado)
finder = pyrodigal.GeneFinder(meta=True)

# Encontrar genes
seq_bytes = bytes(str(record.seq), 'utf-8')
genes = finder.find_genes(seq_bytes)

print(f"✓ Encontrados {len(genes)} genes")
for i in range(min(3, len(genes))):  # Mostrar primeiro 3
    gene = genes[i]
    print(f"  Gene {i+1}: {gene.start}-{gene.end} ({gene.strand})")
    print(f"    Partial: {gene.partial_begin}, {gene.partial_end}")
