# procurator/annotation.py

import pyhmmer
import pandas as pd
import os
from . import io

def run_annotation(args):
    """
    Anota sequências de proteína usando pyhmmer (HMMER)
    contra um banco de dados de perfis HMM.
    """
    print(f"Iniciando anotação HMM com o banco: {args.hmm_database}")
    
    # Carregar as sequências de proteína
    sequences = list(io.read_fasta_sequences(args.input_faa))
    if not sequences:
        print("Nenhuma sequência encontrada no arquivo FASTA de entrada.")
        return

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
    print(f"Executando hmmscan com {cpus_to_use} CPUs...")

    # Abrir o banco de dados HMM e executar o hmmscan
    try:
        with pyhmmer.plan7.HMMFile(args.hmm_database) as hmm_file:
            all_hits = list(pyhmmer.hmmscan(digital_sequences, hmm_file, cpus=cpus_to_use))
    except FileNotFoundError:
        print(f"Erro: Banco de dados HMM '{args.hmm_database}' não encontrado.")
        return
    except Exception as e:
        print(f"Erro ao executar hmmscan: {e}")
        return

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
        print("Nenhum hit significativo encontrado.")
        # Salvar um arquivo vazio com cabeçalho
        df = pd.DataFrame(columns=["Query_ID", "Target_HMM", "Accession", "Description", "E-value", "Score"])
    else:
        df = pd.DataFrame(results)

    # Salvar em um TSV
    df.to_csv(args.output_tsv, sep="\t", index=False)
    
    print(f"Anotação concluída. {len(df)} hits salvos em: {args.output_tsv}")