# procurator/analysis.py
# (Este módulo não foi detalhado antes, mas aqui está uma implementação)

import pandas as pd
from Bio.SeqUtils import GC, N50
from . import io

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
    gc_content_total = sum(GC(rec.seq) * len(rec) for rec in records)
    avg_gc = (gc_content_total / total_length) if total_length > 0 else 0
    
    # Calcular N50
    n50_val = N50(lengths)
    
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