# procurator/analysis.py

import pandas as pd
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
import Bio.SeqUtils
from . import io

def calculate_gc(seq):
    """Calculate GC content as percentage."""
    seq_str = str(seq).upper()
    gc_count = seq_str.count('G') + seq_str.count('C')
    return (gc_count / len(seq_str) * 100) if seq_str else 0

def calculate_n50(lengths):
    """Calculate N50 statistic."""
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= total / 2:
            return length
    return 0

def run_stats(args):
    """
    Calcula estatísticas básicas do genoma (GC, N50, etc.).
    """
    print(f"Calculando estatísticas para: {args.input}")
    
    records = io.read_fasta_sequences(args.input)
    
    if not records:
        print("Nenhuma sequência encontrada no arquivo.")
        return

    lengths = [len(rec) for rec in records]
    total_length = sum(lengths)
    
    # Calcular GC geral (ponderado pelo comprimento)
    gc_content_total = sum(calculate_gc(rec.seq) * len(rec) for rec in records)
    avg_gc = (gc_content_total / total_length) if total_length > 0 else 0
    
    # Calcular N50
    n50_val = calculate_n50(lengths)
    
    stats_data = {
        "Total_Sequencias": [len(records)],
        "Total_Bases_bp": [total_length],
        "Comprimento_Medio_bp": [f"{total_length / len(records):.2f}" if records else 0],
        "Comprimento_Max_bp": [max(lengths) if lengths else 0],
        "Comprimento_Min_bp": [min(lengths) if lengths else 0],
        "N50_bp": [n50_val],
        "GC_Medio_Perc": [f"{avg_gc:.2f}"],
    }
    
    df = pd.DataFrame(stats_data)
    
    # Salvar em CSV
    io.write_dataframe_to_csv(df, args.output)
    print(f"Estatísticas salvas em: {args.output}")
    print("\n--- Sumário das Estatísticas ---")
    print(df.transpose().to_string(header=False))
    print("------------------------------\n")