#!/bin/bash
# Script para salvar screenshot do Procurator em execução

echo "Este script ajuda a capturar uma screenshot de execução do Procurator"
echo ""
echo "Passos:"
echo "1. O comando abaixo será executado"
echo "2. Quando terminar, tire uma screenshot do terminal"
echo "3. Salve como: img/procurator-run.png"
echo ""
echo "Pressione ENTER para continuar..."
read

# Ativar ambiente
conda activate procurator

# Executar comando que gera boa saída visual
echo "Executando find_orfs para gerar saída visual..."
python run_procurator.py find_orfs \
    -i data/test_genome.fasta \
    -o_gff test/genes.gff \
    -o_faa test/proteins.faa

echo ""
echo "Comando completo! Agora tire uma screenshot:"
echo ""
echo "macOS:  Cmd+Shift+5 ou use Screenshot app"
echo "Linux:  gnome-screenshot ou PrintScreen"
echo "Windows: Win+Shift+S"
echo ""
echo "Salve como: img/procurator-run.png"
