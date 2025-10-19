# procurator/orf_finding.py

import pyrodigal
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from . import io
from . import ui

def run_orf_finding(args):
    """
    Executa a predição de genes com pyrodigal (Prodigal).
    Salva arquivos GFF3 e FASTA de proteínas.
    """
    tracker = ui.StepTracker()
    tracker.start_step("Loading sequences")
    
    sequences = io.read_fasta_sequences(args.input)
    if not sequences:
        ui.print_error("No sequences found in the input file.")
        return
    
    tracker.end_step()
    tracker.start_step("Initializing gene finder (Prodigal)")
    
    try:
        import pyrodigal
        finder = pyrodigal.GeneFinder(meta=True)  # Use meta mode (pré-treinado)
    except Exception as e:
        ui.print_error(f"Failed to initialize Prodigal: {e}")
        return
    
    tracker.end_step()

    protein_records = []
    gff_records = []

    # Create progress bar
    progress = ui.ProgressBar(len(sequences), title="Predicting genes")
    
    for record in sequences:
        seq_str = str(record.seq)
        
        try:
            seq_bytes = bytes(seq_str, 'utf-8')
            genes = finder.find_genes(seq_bytes)
        except Exception as e:
            ui.print_error(f"Error processing {record.id}: {e}")
            genes = []

        gff_features = []
        
        for i in range(len(genes)):
            gene = genes[i]
            
            # 1. Preparar o GFF
            start = gene.begin  # Base 1
            end = gene.end
            strand = 1 if gene.strand == '+' else -1
            
            feature = SeqFeature(
                FeatureLocation(start - 1, end, strand=strand),  # Biopython usa base 0
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
            protein_desc = f"[contig={record.id}] [pos={start}-{end}({gene.strand})] [GC={gene.gc_cont:.1f}%]"
            
            protein_records.append(
                SeqRecord(
                    Seq(gene.translate(args.translation_table)),
                    id=protein_id,
                    description=protein_desc
                )
            )
        
        record.features = gff_features
        gff_records.append(record)
        progress.update(1)
    
    progress.finish()

    # 3. Salvar os resultados
    tracker.start_step("Writing GFF3 file")
    io.write_gff(gff_records, args.output_gff)
    tracker.end_step()
    
    tracker.start_step("Writing protein FASTA file")
    io.write_fasta(protein_records, args.output_protein)
    tracker.end_step()
    
    # Display results
    print(ui.DataFormatter.format_results_summary(
        "GENE PREDICTION COMPLETE",
        {
            "Sequences processed": len(sequences),
            "Total genes predicted": len(protein_records),
            "Proteins saved to": args.output_protein,
            "Annotations saved to": args.output_gff,
        }
    ))
    
    tracker.print_summary()