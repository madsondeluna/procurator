# PROCURATOR: Analisador de Sequências Procarióticas

Procurator é uma ferramenta de linha de comando em Python para análise de sequências de DNA procariótico em arquivos FASTA.

> Para fins de testes em ambientes controlados e educacionais. Ferramenta em versão 0.1, ainda em atualizações.



<p align="left">
  <img src="https://github.com/madsondeluna/procurator/blob/main/logo.png?raw=true" alt="PROCURATOR logo" width="800"/>
</p>



## Funcionalidades 

* **Estatísticas Gerais**: contagem de GC, comprimento das sequências, número de bases ambíguas ('N') etc.
* **Busca por Elementos Regulatórios**: identificação de caixas **-35** (`TTGACA`), **-10** (`TATAAT`) e **RBS** (Shine-Dalgarno, `AAGGAGG`).
* **Identificação de ORFs**: detecção de Open Reading Frames usando padrão regex `ATG(?:[ATGC]{3})*?(?:TAA|TAG|TGA)` em ambas fitas.
* **Análise de Proteínas**: tradução de ORFs em aminoácidos, cálculo de peso molecular, ponto isoelétrico (pI), GRAVY score e predição de estrutura secundária.
* **Exportação de Dados**: geração de arquivos CSV (`_sumario.csv`, `_analise_individual.csv`, `_elementos_regulatorios.csv`, `_analise_orfs.csv`).

## Metodologia 

### 1. Conteúdo de GC (%)

```text
GC(%) = ((#G + #C) / comprimento total) * 100
```

### 2. Busca por Elementos Reguladores e ORFs

* **-35 Box**: sequência consenso `TTGACA`
* **-10 Box**: sequência consenso `TATAAT`
* **RBS**: sequência consenso `AAGGAGG`
* **ORFs**: regex `ATG(?:[ATGC]{3})*?(?:TAA|TAG|TGA)`, aplicado nas fitas direta (+) e reversa (-).

### 3. Análise de Proteínas

1. **Tradução**: leitura de códons até códon de parada (`*`).
2. **Peso Molecular**:

   ```text
   massa = soma(massa_aa) - (n_aa - 1) * 18.015 Da
   ```
3. **Ponto Isoelétrico (pI)**: busca iterativa em pH 0–14 até carga líquida zero.
4. **GRAVY Score**:

   ```text
   GRAVY = soma(valor_hidropatia_aa) / n_aa
   ```
5. **Predição de Estrutura Secundária (Chou-Fasman simplificado)**:

   * `P(a) > P(b)`: hélice (H)
   * `P(b) > P(a)`: folha (E)
   * igual: laço/alça (C)

## Pré-requisitos e Instalação

* **Python 3**
* **Biopython**

```bash
pip install biopython
```

## Uso

```bash
python procurator.py [arquivo.fasta] [opções]
```

* `-e`, `--elementos`: busca por elementos regulatórios
* `--orfs [tamanho_min]`: busca por ORFs (padrão: 50 bp)

### Exemplos

```bash
# Estatísticas básicas
python procurator.py genoma.fna

# Estatísticas + elementos regulatórios
python procurator.py genoma.fna -e

# Estatísticas + ORFs (mínimo 150 bp)
python procurator.py genoma.fna --orfs 150

# Análise completa
python procurator.py genoma.fna -e --orfs 90
```

## Resultados

1. **Relatório Geral**

   ```text
   --- Relatório Geral ---
   Total de sequências: N
   Total de bases: X bp
   ...
   ```

2. **Elementos Reguladores** (`-e`)

   ```text
   --- Elementos Reguladores ---
   -35 Box: posições [...] (fita + / -)
   -10 Box: posições [...]
   RBS: posições [...]
   ```

3. **ORFs** (`--orfs`)

   ```text
   --- ORFs Encontradas ---
   ORF 1: 191–1372 (fita +), 1182 bp
   Proteína: MTHKL...* (393 aa)
     - Peso molecular: 43632,88 Da
     - pI: 6,01
     - GRAVY: -0,125
     - Estrutura: HHHCEEECH...
   ```

## Exportação para CSV

Arquivos gerados a partir do nome base (ex.: `meu_genoma`):

* `meu_genoma_sumario.csv`
* `meu_genoma_analise_individual.csv`
* `meu_genoma_elementos_regulatorios.csv`
* `meu_genoma_analise_orfs.csv`
