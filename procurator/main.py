# procurator/main.py
import argparse
import sys

# Importar as funções dos seus módulos
from . import analysis
from . import orf_finding
from . import annotation
from . import protein
from . import io

def setup_argparse():
    parser = argparse.ArgumentParser(
        description="Procurator: Pipeline de Análise e Anotação Genômica de Procariotos."
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Sub-comando a ser executado.")

    # --- Sub-comando 'stats' ---
    stats_parser = subparsers.add_parser(
        "stats", 
        help="Calcula estatísticas básicas do genoma (GC, N50, etc.)."
    )
    stats_parser.add_argument(
        "-i", "--input", required=True, 
        help="Arquivo FASTA de entrada."
    )
    stats_parser.add_argument(
        "-o", "--output", required=True, 
        help="Arquivo CSV de saída para as estatísticas."
    )
    # Define a função que este sub-comando deve chamar
    stats_parser.set_defaults(func=analysis.run_stats) 

    # --- Sub-comando 'find_orfs' ---
    orf_parser = subparsers.add_parser(
        "find_orfs", 
        help="Encontra ORFs/genes usando Prodigal."
    )
    orf_parser.add_argument(
        "-i", "--input", required=True, 
        help="Arquivo FASTA de entrada (genoma completo)."
    )
    orf_parser.add_argument(
        "-o_gff", "--output_gff", required=True, 
        help="Arquivo GFF3 de saída com as coordenadas dos genes."
    )
    orf_parser.add_argument(
        "-o_faa", "--output_protein", required=True, 
        help="Arquivo FASTA de saída com as sequências de proteínas."
    )
    orf_parser.add_argument(
        "-t", "--translation_table", type=int, default=11, 
        help="Tabela de tradução genética (Padrão: 11 para procariotos)."
    )
    orf_parser.set_defaults(func=orf_finding.run_orf_finding)

    # --- Sub-comando 'annotate' (HMMER) ---
    ann_parser = subparsers.add_parser(
        "annotate", 
        help="Anota proteínas usando HMMER contra um banco (ex: Pfam)."
    )
    ann_parser.add_argument(
        "-i", "--input_faa", required=True, 
        help="Arquivo FASTA de proteínas (saída do 'find_orfs')."
    )
    ann_parser.add_argument(
        "-db", "--hmm_database", required=True, 
        help="Caminho para o banco de dados HMM (ex: Pfam-A.hmm)."
    )
    ann_parser.add_argument(
        "-o", "--output_tsv", required=True, 
        help="Arquivo TSV de saída com os 'hits' da anotação."
    )
    ann_parser.add_argument(
        "--evalue", type=float, default=1e-5, 
        help="E-value máximo para os hits."
    )
    ann_parser.set_defaults(func=annotation.run_annotation) # Função alvo original

    # --- Sub-comando 'select_proteins' ---
    af_parser = subparsers.add_parser(
        "select_proteins", 
        help="Seleciona proteínas para modelagem (ex: AlphaFold 3)."
    )
    af_parser.add_argument(
        "-i_faa", "--input_faa", required=True, 
        help="Arquivo FASTA de proteínas."
    )
    af_parser.add_argument(
        "-i_ann", "--input_annotations", required=True, 
        help="Arquivo TSV de anotações (saída do 'annotate')."
    )
    af_parser.add_argument(
        "-k", "--keywords", nargs='+', required=True, 
        help="Palavras-chave para buscar nas anotações (ex: 'kinase', 'transferase')."
    )
    af_parser.add_argument(
        "-o", "--output_fasta", required=True, 
        help="Arquivo FASTA de saída com as sequências selecionadas."
    )
    af_parser.set_defaults(func=protein.select_proteins)

    return parser

def main():
    parser = setup_argparse()
    args = parser.parse_args()
    
    # Chama a função associada ao sub-comando
    if hasattr(args, 'func'):
        try:
            args.func(args)
        except Exception as e:
            print(f"Erro ao executar o comando '{args.command}': {e}", file=sys.stderr)
            return 1
    return 0

if __name__ == "__main__":
    main()