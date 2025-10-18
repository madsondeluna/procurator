# procurator/protein.py

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import subprocess
import shutil
import sys
import os
import importlib.metadata # <--- Importar
from . import io

# ... (função 'calculate_physchem_properties' permanece a mesma) ...
# ... (função 'select_proteins_for_modeling' permanece a mesma) ...


# --- FUNÇÃO ATUALIZADA ---
def _check_colabfold_installed():
    """
    Verifica se 'colabfold_batch' está instalado.
    Se não estiver, pergunta ao usuário se deseja instalar via pip.
    """
    try:
        # Verifica se o pacote python 'colabfold' está disponível
        importlib.metadata.distribution("colabfold")
        
        # Verifica se o executável 'colabfold_batch' está no PATH
        if shutil.which("colabfold_batch") is None:
            print("Aviso: O pacote 'colabfold' está instalado, mas o executável 'colabfold_batch' não foi encontrado no PATH.", file=sys.stderr)
            print("Por favor, verifique a instalação do seu ambiente (pode ser necessário reativá-lo).", file=sys.stderr)
            return False
        
        return True # Tudo certo
        
    except importlib.metadata.PackageNotFoundError:
        # O pacote não está instalado
        print("Aviso: A dependência 'colabfold' é necessária para a modelagem, mas não está instalada.", file=sys.stderr)
        
        # Pergunta ao usuário se deseja instalar (como no seu script original)
        prompt = "Deseja instalá-la agora (via 'pip install colabfold')? (s/n): "
        resposta = input(prompt)
        
        if resposta.lower().strip() == 's':
            print("Instalando 'colabfold'...")
            try:
                # Usa sys.executable para garantir que está instalando no pip correto
                subprocess.check_call([sys.executable, "-m", "pip", "install", "colabfold"])
                print("'colabfold' instalado com sucesso.")
                
                # Verifica novamente se o executável está no PATH
                if shutil.which("colabfold_batch") is None:
                     print("Erro: 'colabfold' foi instalado, mas 'colabfold_batch' não está no PATH.", file=sys.stderr)
                     print("Por favor, reinicie o script ou ajuste seu PATH.", file=sys.stderr)
                     return False
                return True
                
            except subprocess.CalledProcessError as e:
                print(f"Falha ao instalar 'colabfold'. Erro: {e}", file=sys.stderr)
                return False
        else:
            print("Instalação cancelada. O comando 'run_modeling' não pode continuar.", file=sys.stderr)
            return False

# --- FUNÇÃO ATUALIZADA ---
def run_modeling(args):
    """
    Executa a modelagem de proteínas usando 'colabfold_batch' (subprocess).
    """
    print("Iniciando a modelagem de proteínas com ColabFold...")

    # 1. Verificar (e tentar instalar) o colabfold_batch
    if not _check_colabfold_installed():
        sys.exit(1)

    # 2. Validar o FASTA de entrada
    if not os.path.exists(args.input_fasta):
        print(f"Erro: Arquivo FASTA de entrada '{args.input_fasta}' não encontrado.", file=sys.stderr)
        sys.exit(1)
    
    # 3. Criar o diretório de saída
    os.makedirs(args.output_dir, exist_ok=True)

    # 4. Montar o comando
    cmd = [
        "colabfold_batch",
        args.input_fasta,
        args.output_dir,
        "--model-type", args.model_type,
        "--num-models", str(args.num_models),
        "--stop-at-score", "80.0"
    ]

    print(f"Executando comando: {' '.join(cmd)}")
    
    # --- AVISO SOBRE BANCO DE DADOS ---
    print("\n--- AVISO IMPORTANTE SOBRE BANCOS DE DADOS ---")
    print("O 'colabfold_batch' tentará usar bancos de dados locais.")
    print("Se os bancos de dados NÃO estiverem instalados, ele tentará usar o servidor MMseqs2 remoto (requer internet).")
    print("Para melhor performance e uso offline, execute o setup de bancos de dados UMA VEZ:")
    print("  colabfold_setup databases")
    print("(Este comando baixa muitos GB e pode levar horas. Execute-o manualmente.)")
    print("------------------------------------------------\n")
    print("Iniciando a execução... Isso pode levar muito tempo.\n")
    
    # 5. Executar o subprocesso
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            encoding='utf-8'
        )
        
        print("--- Log do ColabFold (stdout) ---")
        print(result.stdout)
        print("---------------------------------")
        if result.stderr:
            print("--- Log do ColabFold (stderr) ---")
            print(result.stderr)
            print("---------------------------------")
            
        print("\nModelagem concluída com sucesso.")
        print(f"Modelos PDB e arquivos relacionados salvos em: {args.output_dir}")

    except subprocess.CalledProcessError as e:
        print(f"Erro ao executar 'colabfold_batch'.", file=sys.stderr)
        print("\n--- Erro (stdout) ---", file=sys.stderr)
        print(e.stdout, file=sys.stderr)
        print("\n--- Erro (stderr) ---", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        
        # --- ORIENTAÇÃO DE ERRO COMUM ---
        if "No such file or directory" in e.stderr or "Cannot find database" in e.stderr:
             print("\n--- POSSÍVEL SOLUÇÃO ---", file=sys.stderr)
             print("O erro pode indicar que os bancos de dados locais não foram encontrados E", file=sys.stderr)
             print("a conexão com o servidor remoto falhou ou foi desativada.", file=sys.stderr)
             print("Por favor, execute 'colabfold_setup databases' manualmente em seu terminal.", file=sys.stderr)
        
        sys.exit(1)