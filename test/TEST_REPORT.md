# Procurator Pipeline - Teste com Pseudo-Genoma Procarioto

## Resumo dos Testes

Teste realizado em **18 de outubro de 2025** com um genoma procarioto sintético.

### Genoma Sintético Gerado

- **Total de sequências**: 3 contigs
- **Total de bases**: 12.000 bp
- **Comprimento médio**: 4.000 bp
- **Comprimento máximo**: 5.000 bp
- **N50**: 4.000 bp
- **GC médio**: 52,69%

### Testes Executados

#### ✅ 1. Comando `stats` - PASSOU

**Comando:**
```bash
python run_procurator.py stats -i data/test_genome.fasta -o test_results/stats.csv
```

**Resultado:**
- Estatísticas genômicas calculadas corretamente
- Arquivo CSV gerado: `test_results/stats.csv`
- Output esperado: ✓ Confirmado

**Saída:**
```
Calculando estatísticas para: data/test_genome.fasta
Estatísticas salvas em: test_results/stats.csv

--- Sumário das Estatísticas ---
Total_Sequencias            3
Total_Bases_bp          12000
Comprimento_Medio_bp  4000.00
Comprimento_Max_bp       5000
Comprimento_Min_bp       3000
N50_bp                   4000
GC_Medio_Perc           52.69
```

#### ✅ 2. Comando `find_orfs` - PASSOU

**Comando:**
```bash
python run_procurator.py find_orfs \
    -i data/test_genome.fasta \
    -o_gff test_results/genes.gff \
    -o_faa test_results/proteins.faa
```

**Resultado:**
- **13 proteínas preditas** (produtos de genes ORF)
- **3 contigs processados** com features GFF3
- Arquivo GFF3: `test_results/genes.gff`
- Arquivo FASTA: `test_results/proteins.faa`

**Saída:**
```
Iniciando predição de genes em: data/test_genome.fasta
Predição concluída.
  - 13 proteínas salvas em: test_results/proteins.faa
  - 3 contigs com features salvos em: test_results/genes.gff
```

**Exemplo de proteína predita:**
```
>NC_000001_gene_1 [contig=NC_000001] [pos=501-1400(+1)] [GC=0.5%]
MSIMEKHCIQHDVGFTHQNPLQTKRRPMVACHLPFNLPFMHFKWIWRRHRLMIFVCFERM...
```

## Arquivo de Dados de Entrada

- **Path**: `data/test_genome.fasta`
- **Gerado por**: `generate_test_genome.py`
- **Função**: Script para gerar genomas procariotos sintéticos para testes

## Ambiente

- **Python**: 3.13
- **Conda/Mamba**: Ativado
- **Ambiente**: `procurator`
- **Dependências principais**:
  - Biopython
  - Pandas
  - Pyrodigal (GeneFinder com meta=True)
  - BCBio-GFF

## Conclusões

✅ **Pipeline operacional para etapas testadas:**
1. **Leitura de FASTA genômico** → Funcionando
2. **Cálculo de estatísticas** (`stats`) → Funcionando
3. **Predição de genes** (`find_orfs`) → Funcionando com Prodigal
4. **Saída em GFF3 e proteínas FASTA** → Funcionando

⚠️ **Funcionalidades ainda não testadas:**
- `annotate` (necessário Pfam-A.hmm)
- `select_proteins` (depende de anotações)
- `run_modeling` (ColabFold - não instalado)

## Próximas Etapas (Opcionais)

1. Testar com `annotate` (download Pfam-A.hmm para `db/`)
2. Testar `select_proteins` com keywords
3. Setup do ColabFold para `run_modeling`
4. Testes com genomas reais

---

**Status**: ✅ **PRONTO PARA PRODUÇÃO** (para etapas testadas)
