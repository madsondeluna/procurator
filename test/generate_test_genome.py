#!/usr/bin/env python3
"""
Generate a synthetic prokaryotic genome for testing Procurator pipeline.
Creates a realistic sequence with embedded genes for testing.
"""

import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write

# Set seed for reproducibility
random.seed(42)

def generate_dna_sequence(length, gc_content=0.5):
    """Generate random DNA sequence with specified GC content."""
    gc_bases = ['G', 'C']
    at_bases = ['A', 'T']
    
    num_gc = int(length * gc_content)
    num_at = length - num_gc
    
    bases = gc_bases * (num_gc // 2) + at_bases * (num_at // 2)
    random.shuffle(bases)
    return ''.join(bases[:length])

def generate_codon_sequence(length_bp=300):
    """Generate a coding sequence (start codon + codons + stop codon)."""
    codons = {
        'Met': ['ATG'],  # Start
        'Lys': ['AAA', 'AAG'],
        'Phe': ['TTT', 'TTC'],
        'Leu': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'Ser': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'Asp': ['GAT', 'GAC'],
        'Asn': ['AAT', 'AAC'],
        'Thr': ['ACT', 'ACC', 'ACA', 'ACG'],
        'Trp': ['TGG'],
        'Cys': ['TGT', 'TGC'],
        'His': ['CAT', 'CAC'],
        'Gln': ['CAA', 'CAG'],
        'Arg': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'Gly': ['GGT', 'GGC', 'GGA', 'GGG'],
        'Glu': ['GAA', 'GAG'],
        'Ala': ['GCT', 'GCC', 'GCA', 'GCG'],
        'Val': ['GTT', 'GTC', 'GTA', 'GTG'],
        'Pro': ['CCT', 'CCC', 'CCA', 'CCG'],
        'Ile': ['ATT', 'ATC', 'ATA'],
        'Stop': ['TAA', 'TAG', 'TGA'],
    }
    
    # Build sequence: ATG (start) + coding region + stop codon
    sequence = 'ATG'
    while len(sequence) < length_bp - 3:
        # Pick random amino acid
        amino_acid = random.choice(list(codons.keys()))
        if amino_acid != 'Stop':
            sequence += random.choice(codons[amino_acid])
    
    # Add stop codon
    sequence += random.choice(codons['Stop'])
    return sequence[:length_bp]

def create_synthetic_genome():
    """Create a synthetic prokaryotic genome with genes."""
    
    records = []
    
    # Contig 1: ~5 kb with genes
    contig1_seq = generate_dna_sequence(5000, gc_content=0.55)
    # Insert genes at specific positions
    gene1 = generate_codon_sequence(900)  # 300 aa
    gene2 = generate_codon_sequence(1200)  # 400 aa
    gene3 = generate_codon_sequence(600)   # 200 aa
    
    contig1_seq = (contig1_seq[:500] + gene1 + 
                   contig1_seq[1400:2500] + gene2 + 
                   contig1_seq[3700:4400] + gene3 + 
                   contig1_seq[5000:])[:5000]
    
    records.append(SeqRecord(
        Seq(contig1_seq),
        id="NC_000001",
        description="Synthetic prokaryotic chromosome 1, ~5 kb [synthetic]"
    ))
    
    # Contig 2: ~3 kb with genes
    contig2_seq = generate_dna_sequence(3000, gc_content=0.52)
    gene4 = generate_codon_sequence(750)   # 250 aa
    gene5 = generate_codon_sequence(450)   # 150 aa
    
    contig2_seq = (contig2_seq[:600] + gene4 + 
                   contig2_seq[1350:2400] + gene5 + 
                   contig2_seq[2850:])[:3000]
    
    records.append(SeqRecord(
        Seq(contig2_seq),
        id="NC_000002",
        description="Synthetic prokaryotic chromosome 2, ~3 kb [synthetic]"
    ))
    
    # Contig 3: ~4 kb with genes (higher GC)
    contig3_seq = generate_dna_sequence(4000, gc_content=0.60)
    gene6 = generate_codon_sequence(1050)  # 350 aa
    gene7 = generate_codon_sequence(825)   # 275 aa
    
    contig3_seq = (contig3_seq[:800] + gene6 + 
                   contig3_seq[1850:2800] + gene7 + 
                   contig3_seq[3625:])[:4000]
    
    records.append(SeqRecord(
        Seq(contig3_seq),
        id="NC_000003",
        description="Synthetic prokaryotic chromosome 3, ~4 kb [synthetic, high GC]"
    ))
    
    return records

if __name__ == "__main__":
    print("Generating synthetic prokaryotic genome for testing...")
    
    records = create_synthetic_genome()
    
    output_file = "data/test_genome.fasta"
    write(records, output_file, "fasta")
    
    total_length = sum(len(rec.seq) for rec in records)
    print(f"[OK] Generated {len(records)} contigs, {total_length} bp total")
    print(f"[OK] Saved to: {output_file}")
    
    for rec in records:
        print(f"  - {rec.id}: {len(rec.seq)} bp - {rec.description}")
