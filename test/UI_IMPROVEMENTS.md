# UI/UX Improvements - Procurator CLI

## Resumo

O Procurator agora tem uma interface de linha de comando (CLI) muito mais amigável com:

- ✨ **Barra de progresso em tempo real** com estimativa de tempo restante
- ⏱️ **Rastreamento de etapas** mostrando timing individual para cada operação
- 📊 **Tabelas formatadas** com melhor legibilidade
- 🎯 **Mensagens claras** com ícones (✅, ❌, ⚠️, ℹ️)
- 📈 **Sumário de execução** com gráfico de distribuição de tempo

## Exemplos de Output

### Comando: `stats`

```
======================================================================
  PROCURATOR - Prokaryotic Genome Analysis Pipeline
======================================================================

🔹 Loading sequences...
  ✓ ✓ completed (17ms)

🔹 Computing statistics...
  ✓ ✓ completed (3ms)

🔹 Saving results...
  ✓ ✓ completed (2ms)

📋 Statistics Summary:
────────────────────────────────────────────────────────────
  Total Sequences       :                    3
  Total Bases (Bp)      :               12,000
  Avg Length (Bp)       :                4,000
  Max Length (Bp)       :                5,000
  Min Length (Bp)       :                3,000
  N50 (Bp)              :                4,000
  Gc Content (%)        :                52.69
────────────────────────────────────────────────────────────

✨ RESULTS SAVED
────────────────────────────────────────────────────────────
  📄 Output file                    : test/stats_new.csv
────────────────────────────────────────────────────────────

============================================================
📊 EXECUTION SUMMARY
============================================================
1. Loading sequences              ✓ completed           17ms ( 77.5%)
   └─ ███████████████████████
2. Computing statistics           ✓ completed            3ms ( 13.2%)
   └─ ███
3. Saving results                 ✓ completed            2ms (  9.3%)
   └─ ██
============================================================
Total Time: 22ms
============================================================
```

### Comando: `find_orfs` (com barra de progresso)

```
======================================================================
  PROCURATOR - Prokaryotic Genome Analysis Pipeline
======================================================================

🔹 Loading sequences...
  ✓ ✓ completed (17ms)

🔹 Initializing gene finder (Prodigal)...
  ✓ ✓ completed (0ms)

Predicting genes: [████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░]   33.3% |     1/3     | ⏱ 9ms | ⏰ 18ms
Predicting genes: [██████████████████████████░░░░░░░░░░░░░░]   66.7% |     2/3     | ⏱ 11ms | ⏰ 6ms
Predicting genes: [████████████████████████████████████████]  100.0% |     3/3     | ⏱ 20ms | ⏰ 0ms
✓ Predicting genes completed in 20ms

🔹 Writing GFF3 file...
  ✓ ✓ completed (1ms)

🔹 Writing protein FASTA file...
  ✓ ✓ completed (0ms)

✨ GENE PREDICTION COMPLETE
────────────────────────────────────────────────────────────
  🔢 Sequences processed            : 3
  🔢 Total genes predicted          : 13
  🔢 Proteins saved to              : test/proteins_new.faa
  🔢 Annotations saved to           : test/genes_new.gff
────────────────────────────────────────────────────────────

============================================================
📊 EXECUTION SUMMARY
============================================================
1. Loading sequences              ✓ completed           17ms ( 96.0%)
   └─ ████████████████████████████
2. Initializing gene finder (Prodigal) ✓ completed            0ms (  0.1%)
   └─ 
3. Writing GFF3 file              ✓ completed            1ms (  3.1%)
   └─ 
4. Writing protein FASTA file     ✓ completed            0ms (  0.8%)
   └─ 
============================================================
Total Time: 18ms
============================================================
```

## Componentes da UI

### 1. **ProgressBar** (`ui.ProgressBar`)

Barra de progresso em tempo real com:
- Porcentagem completada
- Contador (atual/total)
- Tempo decorrido (⏱)
- Tempo estimado restante (⏰)

Exemplo:
```python
progress = ui.ProgressBar(100, title="Processing")
for i in range(100):
    # ... do work ...
    progress.update(1)
progress.finish()
```

Output:
```
Processing: [██████████░░░░░░░░░░░░░░░░░░░░░░] 40.0% | 40/100 | ⏱ 2m 15s | ⏰ 3m 22s
```

### 2. **StepTracker** (`ui.StepTracker`)

Rastreamento de etapas com timing individual:
- Marca início e fim de etapas
- Calcula duração de cada etapa
- Exibe resumo ao final

Exemplo:
```python
tracker = ui.StepTracker()
tracker.start_step("Loading data")
# ... do work ...
tracker.end_step()

tracker.start_step("Processing")
# ... do work ...
tracker.end_step()

tracker.print_summary()
```

### 3. **DataFormatter** (`ui.DataFormatter`)

Formata dados para exibição em tabelas:
- `format_stats_table()` - Tabela de estatísticas
- `format_results_summary()` - Resumo de resultados

### 4. **Funções Auxiliares**

- `print_header(title)` - Cabeçalho formatado
- `print_info(message)` - Mensagem informativa
- `print_success(message)` - Mensagem de sucesso (✅)
- `print_error(message)` - Mensagem de erro (❌)
- `print_warning(message)` - Mensagem de aviso (⚠️)

## Integração em Todos os Comandos

Todos os cinco comandos foram atualizados:

1. **`stats`** - Mostra estatísticas em tabela formatada
2. **`find_orfs`** - Barra de progresso + tabela de resultados
3. **`annotate`** - Tracking de etapas + tabela de anotações
4. **`select_proteins`** - Tracking de etapas + sumário de seleção
5. **`run_modeling`** - Tracking de etapas + mensagens informativas

## Formatação de Tempo

Tempos são automaticamente formatados:
- `< 1s`: Milissegundos (ex: `500ms`)
- `1-60s`: Segundos (ex: `45s`)
- `1-60m`: Minutos e segundos (ex: `2m 30s`)
- `> 60m`: Horas, minutos e segundos (ex: `1h 30m 45s`)

## Benefícios

✅ **Melhor feedback do usuário** - Sabe exatamente o que está acontecendo
✅ **Estimativa de tempo** - Quando vai terminar
✅ **Rastreamento de performance** - Quais etapas são lentas
✅ **Mensagens claras** - Fácil de entender
✅ **Profissional** - Aparência mais polida

---

**Criado em**: 18 de outubro de 2025
**Versão**: Beta
**Commit**: 4ab9460
