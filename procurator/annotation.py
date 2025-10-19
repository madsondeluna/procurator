# procurator/annotation.py

import pyhmmer
import pandas as pd
import os
from . import io
from . import ui

def run_annotation(args):
    """
    Anota sequências de proteína usando pyhmmer (HMMER)
    contra um banco de dados de perfis HMM.
    """
    tracker = ui.StepTracker()
    tracker.start_step("Loading protein sequences")
    
    # Carregar as sequências de proteína
    sequences = list(io.read_fasta_sequences(args.input_faa))
    if not sequences:
        ui.print_error("No sequences found in the FASTA input file.")
        return
    
    tracker.end_step()
    tracker.start_step("Validating HMM database")

    # Converter para o formato que pyhmmer espera
    protein_alphabet = pyhmmer.easel.Alphabet.amino()
    digital_sequences = [
        pyhmmer.easel.TextSequence(name=seq.id.encode(), sequence=str(seq.seq)).digitize(protein_alphabet)
        for seq in sequences
    ]
    
    results = []
    
    # Detectar número de CPUs
    cpus_to_use = os.cpu_count()
    if cpus_to_use is None:
        cpus_to_use = 1 # Padrão
    
    tracker.end_step()
    tracker.start_step(f"Running hmmscan ({cpus_to_use} CPUs)")

    # Abrir o banco de dados HMM e executar o hmmscan
    try:
        with pyhmmer.plan7.HMMFile(args.hmm_database) as hmm_file:
            all_hits = list(pyhmmer.hmmscan(digital_sequences, hmm_file, cpus=cpus_to_use))
    except FileNotFoundError:
        ui.print_error(f"HMM database '{args.hmm_database}' not found.")
        return
    except Exception as e:
        ui.print_error(f"Error running hmmscan: {e}")
        return

    tracker.end_step()
    tracker.start_step("Processing results")

    # Processar os resultados
    for i, hits in enumerate(all_hits):
        query_name = hits.query_name.decode()
        
        for hit in hits:
            if hit.evalue <= args.evalue:
                results.append({
                    "Query_ID": query_name,
                    "Target_HMM": hit.name.decode(),
                    "Accession": hit.accession.decode() if hit.accession else "N/A",
                    "Description": hit.description.decode() if hit.description else "N/A",
                    "E-value": f"{hit.evalue:e}",
                    "Score": hit.score,
                })

    if not results:
        ui.print_warning("No significant hits found.")
        # Salvar um arquivo vazio com cabeçalho
        df = pd.DataFrame(columns=["Query_ID", "Target_HMM", "Accession", "Description", "E-value", "Score"])
    else:
        df = pd.DataFrame(results)

    tracker.end_step()
    tracker.start_step("Saving results")
    
    # Salvar em um TSV
    df.to_csv(args.output_tsv, sep="\t", index=False)
    
    tracker.end_step()
    
    # Display results
    print(ui.DataFormatter.format_results_summary(
        "ANNOTATION COMPLETE",
        {
            "Sequences analyzed": len(sequences),
            "Significant hits": len(df),
            "E-value threshold": f"{args.evalue}",
        }
    ))
    print(ui.DataFormatter.format_files_list(
        "OUTPUT FILES",
        [args.output_tsv]
    ))
    
    tracker.print_summary()