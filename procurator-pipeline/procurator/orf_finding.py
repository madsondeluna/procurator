# procurator/orf_finding.py

import pyrodigal
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from . import io

def run_orf_finding(args):
    """
    Executa a predição de genes com pyrodigal (Prodigal).
    Salva arquivos GFF3 e FASTA de proteínas.
    """
    print(f"Iniciando predição de genes em: {args.input}")
    
    sequences = io.read_fasta_sequences(args.input)
    finder = pyrodigal.Pyrodigal(meta=False) # meta=False para genomas isolados

    protein_records = []
    gff_records = [] # SeqRecord com features

    for record in sequences:
        seq_str = str(record.seq)
        genes = finder.find_genes(seq_str.encode("utf-8"))

        gff_features = []
        
        for i, gene in enumerate(genes):
            # 1. Preparar o GFF
            start = gene.begin # Base 1
            end = gene.end
            strand = 1 if gene.strand == '+' else -1
            
            feature = SeqFeature(
                FeatureLocation(start - 1, end, strand=strand), # Biopython é base 0
                type="CDS",
                qualifiers={
                    "ID": f"cds-{record.id}-{i+1}",
                    "locus_tag": f"{record.id}_gene_{i+1}",
                    "translation": gene.translate(args.translation_table),
                }
            )
            gff_features.append(feature)

            # 2. Preparar o FASTA de proteínas
            protein_id = f"{record.id}_gene_{i+1}"
            protein_desc = f"partial={gene.partial_begin or gene.partial_end} [contig={record.id}] [pos={start}-{end}({gene.strand})]"
            
            protein_records.append(
                SeqRecord(
                    Seq(gene.translate(args.translation_table)),
                    id=protein_id,
                    description=protein_desc
                )
            )
        
        record.features = gff_features
        gff_records.append(record)

    # 3. Salvar os resultados
    io.write_gff(gff_records, args.output_gff)
    io.write_fasta(protein_records, args.output_protein)
    
    print(f"Predição concluída.")
    print(f"  - {len(protein_records)} proteínas salvas em: {args.output_protein}")
    print(f"  - {len(gff_records)} contigs com features salvos em: {args.output_gff}")