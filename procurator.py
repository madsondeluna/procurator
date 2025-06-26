# =============================================================
# SEÇÃO DE VERIFICAÇÃO E INSTALAÇÃO DE DEPENDÊNCIAS
# =============================================================
import sys
import subprocess
import importlib.metadata
import statistics

pacotes_necessarios = ['biopython'] 
for pacote in pacotes_necessarios:
    try:
        importlib.metadata.distribution(pacote)
    except importlib.metadata.PackageNotFoundError:
        prompt = f"A biblioteca '{pacote}' é necessária, mas não está instalada. Deseja instalá-la agora? (s/n): "
        resposta = input(prompt)
        if resposta.lower().strip() == 's':
            print(f"Instalando '{pacote}'...")
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install", pacote])
                print(f"'{pacote}' instalado com sucesso.")
            except subprocess.CalledProcessError as e:
                print(f"Falha ao instalar '{pacote}'. Erro: {e}")
                sys.exit(1)
        else:
            print("Instalação cancelada. O script não pode continuar.")
            sys.exit(1)

# =============================================================
# SEÇÃO DE IMPORTAÇÕES PRINCIPAIS
# =============================================================
import argparse
import re
import csv
import os
from Bio import SeqIO

# =============================================================
# SEÇÃO 0: CONSTANTES
# =============================================================
CODING_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S',
    'TCG':'S', 'TCT':'S', 'AGA':'R', 'AGG':'R', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGA':'*', 'CAC':'H', 'CAT':'H', 'CAA':'Q',
    'CAG':'Q', 'TGC':'C', 'TGT':'C', 'TGG':'W', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
}
AA_RESIDUE_WEIGHTS = {
    'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388, 'E': 129.1155,
    'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
    'M': 131.1926, 'F': 147.1766, 'P': 97.1167, 'S': 87.0782, 'T': 101.1051, 'W': 186.2132,
    'Y': 163.1760, 'V': 99.1326
}
HYDROPATHY_SCALE = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7,
    'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5, 'Q': -3.5, 'D': -3.5,
    'N': -3.5, 'K': -3.9, 'R': -4.5
}
CHOU_FASMAN_PARAMS = {
    'A': {'P(a)': 142, 'P(b)':  83}, 'L': {'P(a)': 121, 'P(b)': 130}, 'R': {'P(a)':  98, 'P(b)':  93},
    'K': {'P(a)': 114, 'P(b)':  74}, 'N': {'P(a)':  67, 'P(b)':  89}, 'M': {'P(a)': 145, 'P(b)': 105},
    'D': {'P(a)': 101, 'P(b)':  54}, 'F': {'P(a)': 113, 'P(b)': 138}, 'C': {'P(a)':  70, 'P(b)': 119},
    'P': {'P(a)':  57, 'P(b)':  55}, 'Q': {'P(a)': 111, 'P(b)': 110}, 'S': {'P(a)':  77, 'P(b)':  75},
    'E': {'P(a)': 151, 'P(b)':  37}, 'T': {'P(a)':  83, 'P(b)': 119}, 'G': {'P(a)':  57, 'P(b)':  75},
    'W': {'P(a)': 108, 'P(b)': 137}, 'H': {'P(a)': 100, 'P(b)':  87}, 'Y': {'P(a)':  69, 'P(b)': 147},
    'I': {'P(a)': 108, 'P(b)': 160}, 'V': {'P(a)': 106, 'P(b)': 170}
}
PKA_VALUES = {
    'C_TERMINUS': 3.65, 'N_TERMINUS': 8.6, 'D': 3.65, 'E': 4.25, 'C': 8.33, 'Y': 10.07,
    'H': 6.0, 'K': 10.53, 'R': 12.48
}

SITIOS_REGEX_PROCARIOTOS = {
    "Minus_35_Box": r'TTG[AC]CA',
    "Minus_10_Box": r'TATAAT',
    "RBS": r'A?GGAGG',
    "Promotor_Completo": r'(TTG[AC]CA)[ATGC]{16,19}(TATAAT)'
}

DESCRICOES_ELEMENTOS = {
    "Minus_35_Box": "Localizada a ~35 pares de base ANTES do início da transcrição. É o primeiro ponto de reconhecimento para o Fator Sigma da RNA Polimerase, ajudando a 'ancorar' a enzima no local correto.",
    "Minus_10_Box": "Também conhecida como 'Pribnow Box', localizada a ~10 pares de base do início da transcrição. Rica em Adenina e Timina, é a região onde a dupla hélice de DNA efetivamente começa a se ABRIR para permitir que a transcrição comece.",
    "RBS": "Sítio de Ligação do Ribossomo (Sequência de Shine-Dalgarno). Localizado no mRNA, um pouco ANTES do códon de início (ATG). É crucial para iniciar a TRADUÇÃO, pois recruta o ribossomo para a síntese da proteína.",
    "Promotor_Completo": "Representa a estrutura funcional do promotor principal. A presença combinada da -35 Box e -10 Box, com o espaçamento correto, é um forte indicativo de um ponto de início de TRANSCRIÇÃO de um gene."
}

# <--- ALTERADO: A função agora contém a sua arte final
def exibir_splash_screen():
    """Exibe o título estilizado no início do programa."""
    
    # CORREÇÃO APLICADA AQUI: Trocado r""" por r'''
    splash_art = r'''
···········································································
:                                                                         :
:                                                                         :
:   8""""8 8"""8  8"""88 8""""8 8   8 8"""8  8""""8 ""8"" 8"""88 8"""8    :
:   8    8 8   8  8    8 8    " 8   8 8   8  8    8   8   8    8 8   8    :
:   8eeee8 8eee8e 8    8 8e     8e  8 8eee8e 8eeee8   8e  8    8 8eee8e   :
:   88     88   8 8    8 88     88  8 88   8 88   8   88  8    8 88   8   :
:   88     88   8 8    8 88   e 88  8 88   8 88   8   88  8    8 88   8   :
:   88     88   8 8eeee8 88eee8 88ee8 88   8 88   8   88  8eeee8 88   8   :
:                                                                         :
:                                                                         :
···········································································
'''
    print(splash_art)
    # A linha abaixo foi removida pois o novo design já é completo
    print("UM ANALISADOR DE SEQUÊNCIAS PROCARIÓTICAS")

# =============================================================
# SEÇÃO 2: FUNÇÕES DE ANÁLISE
# =============================================================
def calculate_gc_manually(sequence):
    seq_str = str(sequence).upper()
    total_length = len(seq_str)
    if total_length == 0:
        return 0.0
    g_count = seq_str.count('G')
    c_count = seq_str.count('C')
    gc_percent = ((g_count + c_count) / total_length) * 100
    return gc_percent

def get_reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, 'N') for base in dna_sequence)[::-1]

def translate_cds(cds):
    if len(cds) % 3 != 0: return "Erro: Sequência de CDS inválida."
    protein_sequence = ""
    for i in range(0, len(cds), 3):
        codon = cds[i:i+3].upper()
        amino_acid = CODING_TABLE.get(codon, '?')
        if amino_acid == '*':
            protein_sequence += '*'
            break
        protein_sequence += amino_acid
    return protein_sequence

def find_orfs_with_regex(sequence, min_len_bp=50):
    orfs_found = []
    orf_pattern = re.compile(r'(?=(ATG(?:[ATGC]{3})*?(?:TAA|TAG|TGA)))')
    for strand, dna_seq in [('+', sequence), ('-', get_reverse_complement(sequence))]:
        for match in orf_pattern.finditer(str(dna_seq).upper()):
            orf_seq = match.group(1)
            if len(orf_seq) >= min_len_bp:
                orfs_found.append({
                    "start": match.start() + 1, "end": match.start() + len(orf_seq),
                    "strand": strand, "length_bp": len(orf_seq), "sequence": orf_seq
                })
    return orfs_found

def predict_secondary_structure(protein_sequence):
    protein = protein_sequence.replace('*', '')
    structure = ""
    for aa in protein:
        params = CHOU_FASMAN_PARAMS.get(aa.upper())
        if not params:
            structure += "?"
            continue
        if params['P(a)'] > params['P(b)']: structure += "H"
        elif params['P(b)'] > params['P(a)']: structure += "E"
        else: structure += "C"
    return structure

def _calculate_net_charge_at_ph(protein, ph):
    charge = 1 / (1 + 10**(ph - PKA_VALUES['N_TERMINUS']))
    charge += -1 / (1 + 10**(PKA_VALUES['C_TERMINUS'] - ph))
    for aa in ['D', 'E', 'C', 'Y']:
        charge += protein.count(aa) * (-1 / (1 + 10**(PKA_VALUES[aa] - ph)))
    for aa in ['H', 'K', 'R']:
        charge += protein.count(aa) * (1 / (1 + 10**(ph - PKA_VALUES[aa])))
    return charge

def calculate_pi(protein_sequence):
    protein = protein_sequence.replace('*', '')
    if not protein: return 0.0
    low_ph, high_ph, precision = 0.0, 14.0, 0.001
    while (high_ph - low_ph) > precision:
        mid_ph = (low_ph + high_ph) / 2
        charge = _calculate_net_charge_at_ph(protein, mid_ph)
        if charge > 0: low_ph = mid_ph
        else: high_ph = mid_ph
    return (low_ph + high_ph) / 2

def calculate_physchem_properties(protein_sequence):
    protein = protein_sequence.replace('*', '')
    if not protein: return {}
    mw = sum(AA_RESIDUE_WEIGHTS.get(aa, 0) for aa in protein) - (len(protein) - 1) * 18.01528
    try: gravy_score = sum(HYDROPATHY_SCALE.get(aa, 0) for aa in protein) / len(protein)
    except ZeroDivisionError: gravy_score = 0
    isoelectric_point = calculate_pi(protein)
    secondary_structure = predict_secondary_structure(protein)
    return {
        "length_aa": len(protein), "molecular_weight_da": f"{mw:.2f}",
        "isoelectric_point": f"{isoelectric_point:.2f}", "gravy_score": f"{gravy_score:.3f}",
        "secondary_structure": secondary_structure
    }

def buscar_sitios_regex(sequencia, padrao_regex):
    posicoes = []
    for match in re.finditer(padrao_regex, str(sequencia), re.IGNORECASE):
        posicoes.append((match.start() + 1, match.end(), '+'))
    rev_comp = get_reverse_complement(str(sequencia))
    for match in re.finditer(padrao_regex, rev_comp, re.IGNORECASE):
        start_rev = match.start()
        end_rev = match.end()
        start_fwd = len(sequencia) - end_rev + 1
        end_fwd = len(sequencia) - start_rev
        posicoes.append((start_fwd, end_fwd, '-'))
    return posicoes

def analisar_fasta(caminho_arquivo, buscar_elementos=False, min_orf_len=None):
    try:
        sequencias = list(SeqIO.parse(caminho_arquivo, "fasta"))
    except FileNotFoundError:
        print(f"Erro: O arquivo '{caminho_arquivo}' não foi encontrado.")
        return None
    except Exception as e:
        print(f"Ocorreu um erro ao ler o arquivo: {e}")
        return None
    if not sequencias:
        print("Nenhuma sequência encontrada no arquivo.")
        return None
    
    analise_individual = {}
    comprimentos = []
    gc_contents = []
    
    for seq in sequencias:
        length = len(seq.seq)
        comprimentos.append(length)
        gc_val = calculate_gc_manually(seq.seq)
        gc_contents.append(gc_val)
        analise_individual[seq.id] = {
            "comprimento": length,
            "gc_percent": gc_val,
            "n_count": seq.seq.upper().count('N')
        }

    total_n_bases = sum(stats['n_count'] for stats in analise_individual.values())
    seq_max_len = max(sequencias, key=len)
    seq_min_len = min(sequencias, key=len)
    
    resultados = {
        "analise_individual": analise_individual,
        "sumario_geral": {
            "total_sequencias": len(sequencias),
            "total_bases": sum(comprimentos),
            "total_n_bases": total_n_bases,
            "comprimento_maximo": len(seq_max_len.seq),
            "max_len_id": seq_max_len.id,
            "comprimento_minimo": len(seq_min_len.seq),
            "min_len_id": seq_min_len.id,
            "comprimento_medio": f"{statistics.mean(comprimentos):.2f}" if comprimentos else "N/A",
            "stdev_comprimento": f"{statistics.stdev(comprimentos):.2f}" if len(comprimentos) > 1 else "N/A",
            "avg_gc": f"{statistics.mean(gc_contents):.2f}" if gc_contents else "N/A",
            "min_gc": f"{min(gc_contents):.2f}" if gc_contents else "N/A",
            "max_gc": f"{max(gc_contents):.2f}" if gc_contents else "N/A",
        },
        "resultados_busca": None,
        "analise_orfs": None
    }
    
    if buscar_elementos:
        resultados_busca_geral = {}
        for seq_record in sequencias:
            encontrados_nesta_seq = {}
            for elemento, padrao_regex in SITIOS_REGEX_PROCARIOTOS.items():
                posicoes = buscar_sitios_regex(seq_record.seq, padrao_regex)
                if posicoes:
                    encontrados_nesta_seq[elemento] = posicoes
            if encontrados_nesta_seq:
                resultados_busca_geral[seq_record.id] = encontrados_nesta_seq
        resultados["resultados_busca"] = resultados_busca_geral
    
    if min_orf_len is not None:
        resultados_orfs = {}
        for seq_record in sequencias:
            orfs = find_orfs_with_regex(seq_record.seq, min_orf_len)
            if orfs:
                for orf_data in orfs:
                    protein_seq = translate_cds(orf_data["sequence"])
                    orf_data["protein_sequence"] = protein_seq
                    orf_data["physchem_properties"] = calculate_physchem_properties(protein_seq)
                resultados_orfs[seq_record.id] = orfs
        resultados["analise_orfs"] = resultados_orfs
        
    return resultados

# =============================================================
# SEÇÃO 3: CÓDIGO PRINCIPAL, IMPRESSÃO E EXPORTAÇÃO
# =============================================================
def imprimir_resultados(resultados):
    if resultados is None: return

    sumario = resultados["sumario_geral"]
    print("\n--- Relatório Geral do Arquivo FASTA (Análise Procariótica) ---")
    print("SUMÁRIO GERAL:")
    print(f"  - Total de Sequências: {sumario['total_sequencias']}")
    print(f"  - Total de Bases (bp): {sumario['total_bases']:,}".replace(",", "."))
    print(f"  - Bases Ambíguas ('N'): {sumario['total_n_bases']:,}".replace(",", "."))
    
    print("\nESTATÍSTICAS DE COMPRIMENTO:")
    print(f"  - Médio: {sumario['comprimento_medio']} bp")
    print(f"  - Desvio Padrão: {sumario['stdev_comprimento']} bp")
    print(f"  - Máximo: {sumario['comprimento_maximo']} bp (ID: {sumario['max_len_id']})")
    print(f"  - Mínimo: {sumario['comprimento_minimo']} bp (ID: {sumario['min_len_id']})")

    print("\nCONTEÚDO GC (%):")
    print(f"  - Médio: {sumario['avg_gc']}")
    print(f"  - Máximo: {sumario['max_gc']}")
    print(f"  - Mínimo: {sumario['min_gc']}")

    if resultados.get("analise_individual"):
        print("\n\n--- Análise Individual por Sequência ---")
        for seq_id, stats in resultados["analise_individual"].items():
            print(f"\n-> Sequência ID: {seq_id}")
            print(f"  - Comprimento: {stats['comprimento']} bp")
            print(f"  - Conteúdo GC: {stats['gc_percent']:.2f} %")
            print(f"  - Bases 'N': {stats['n_count']}")

    if resultados.get("resultados_busca") and resultados["resultados_busca"]:
        busca = resultados["resultados_busca"]
        print("\n\n--- Resultados da Busca por Elementos Regulatórios Procarióticos ---")
        
        print("\nLEGENDA DOS ELEMENTOS ENCONTRADOS:")
        for elemento, descricao in DESCRICOES_ELEMENTOS.items():
            print(f"  -> {elemento}: {descricao}")
            
        for seq_id, elementos_encontrados in busca.items():
            print(f"\n-> Sítios encontrados em '{seq_id}':")
            for elemento, posicoes in elementos_encontrados.items():
                pos_str = ", ".join([f"[{p[0]}-{p[1]} Fita:{p[2]}]" for p in posicoes])
                print(f"  - {elemento}: encontrado na(s) posição(ões) {pos_str}")

    elif resultados.get("resultados_busca") is not None:
        print("\n\n--- Resultados da Busca por Elementos Regulatórios Procarióticos ---")
        print("Nenhum dos elementos definidos foi encontrado nas sequências.")
                  
    if resultados.get("analise_orfs") and resultados["analise_orfs"]:
        analise_orfs = resultados["analise_orfs"]
        print("\n\n--- Resultados da Análise de ORFs ---")
        for seq_id, orfs in analise_orfs.items():
            print(f"\n-> Em '{seq_id}': Foram encontradas {len(orfs)} ORF(s).")
            for i, orf in enumerate(orfs, 1):
                print(f"  [ORF {i}]")
                print(f"    Posição..: {orf['start']}-{orf['end']} (Fita: {orf['strand']})")
                print(f"    Tamanho..: {orf['length_bp']} bp")
                print(f"    Proteína.: {orf['protein_sequence']}")
                if orf.get("physchem_properties"):
                    props = orf['physchem_properties']
                    print(f"    Propriedades da Proteína ({props['length_aa']} aa):")
                    print(f"      - Peso Molecular: {props['molecular_weight_da']} Da")
                    print(f"      - Ponto Isoelétrico (pI): {props['isoelectric_point']}")
                    print(f"      - GRAVY Score: {props['gravy_score']}")
                    print(f"      - Estrutura Secundária Predita: {props['secondary_structure']}")
    elif resultados.get("analise_orfs") is not None:
        print("\n\n--- Resultados da Análise de ORFs ---")
        print("Nenhuma ORF (com o comprimento mínimo especificado) foi encontrada.")
    
    print("\n---------------------------------------------\n")

def gerar_csvs(resultados, nome_base):
    if not resultados:
        return
    
    print(f"Gerando arquivos CSV com o nome base '{nome_base}'...")

    sumario = resultados.get("sumario_geral")
    if sumario:
        nome_arquivo = f"{nome_base}_sumario.csv"
        with open(nome_arquivo, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['Metrica', 'Valor'])
            for chave, valor in sumario.items():
                writer.writerow([chave, valor])
        print(f"  - Arquivo de sumário gerado: {nome_arquivo}")

    individual = resultados.get("analise_individual")
    if individual:
        nome_arquivo = f"{nome_base}_analise_individual.csv"
        with open(nome_arquivo, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['ID_Sequencia', 'Comprimento_bp', 'GC_Percent', 'N_Count'])
            for seq_id, stats in individual.items():
                writer.writerow([seq_id, stats['comprimento'], f"{stats['gc_percent']:.2f}", stats['n_count']])
        print(f"  - Arquivo de análise individual gerado: {nome_arquivo}")

    busca = resultados.get("resultados_busca")
    if busca:
        nome_arquivo = f"{nome_base}_elementos_regulatorios.csv"
        with open(nome_arquivo, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['ID_Sequencia', 'Elemento', 'Posicao_Inicio', 'Posicao_Fim', 'Fita'])
            for seq_id, elementos in busca.items():
                for elemento, posicoes in elementos.items():
                    for pos in posicoes:
                        writer.writerow([seq_id, elemento, pos[0], pos[1], pos[2]])
        print(f"  - Arquivo de elementos regulatórios gerado: {nome_arquivo}")

    analise_orfs = resultados.get("analise_orfs")
    if analise_orfs:
        nome_arquivo = f"{nome_base}_analise_orfs.csv"
        with open(nome_arquivo, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            header = [
                'ID_Sequencia_Origem', 'ORF_Num', 'Posicao_Inicio', 'Posicao_Fim', 'Fita', 
                'Tamanho_bp', 'Tamanho_aa', 'Peso_Molecular_Da', 'Ponto_Isoeletrico_pI', 
                'GRAVY_Score', 'Estrutura_Secundaria_Predita', 'Sequencia_Proteica'
            ]
            writer.writerow(header)
            for seq_id, orfs in analise_orfs.items():
                for i, orf in enumerate(orfs, 1):
                    props = orf.get("physchem_properties", {})
                    row = [
                        seq_id, i, orf['start'], orf['end'], orf['strand'], orf['length_bp'],
                        props.get('length_aa', 'N/A'),
                        props.get('molecular_weight_da', 'N/A'),
                        props.get('isoelectric_point', 'N/A'),
                        props.get('gravy_score', 'N/A'),
                        props.get('secondary_structure', 'N/A'),
                        orf['protein_sequence']
                    ]
                    writer.writerow(row)
        print(f"  - Arquivo de análise de ORFs gerado: {nome_arquivo}")
    
    print("Exportação para CSV concluída.")


def main():
    exibir_splash_screen()

    parser = argparse.ArgumentParser(
        description="Analisa um arquivo FASTA de procariotos, buscando por elementos regulatórios e/ou ORFs.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("arquivo", help="O caminho para o arquivo FASTA de entrada.")
    parser.add_argument(
        "-e", "--elementos", action="store_true",
        help="Ativa a busca por elementos regulatórios procarióticos (-10, -35, RBS, etc.)."
    )
    parser.add_argument(
        "--orfs", type=int, nargs='?', const=50, default=None, metavar='MIN_LEN',
        help="Ativa a busca por ORFs. Opcionalmente, especifique um comprimento mínimo em bp.\nPadrão: 50 bp."
    )
    args = parser.parse_args()
    
    resultados_analise = analisar_fasta(
        caminho_arquivo=args.arquivo,
        buscar_elementos=args.elementos,
        min_orf_len=args.orfs
    )
    
    imprimir_resultados(resultados_analise)

    if resultados_analise:
        resposta = input("Deseja salvar estes resultados em arquivos CSV? (s/n): ").lower().strip()
        if resposta == 's':
            base_nome = os.path.splitext(os.path.basename(args.arquivo))[0]
            nome_saida = input(f"Digite o nome base para os arquivos de saída [padrão: {base_nome}]: ").strip()
            
            if not nome_saida:
                nome_saida = base_nome
            
            gerar_csvs(resultados_analise, nome_saida)

if __name__ == "__main__":
    main()