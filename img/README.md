# Images Directory

This folder contains images and screenshots for the Procurator documentation.

## Screenshots

### procurator-run.png

Screenshot showing Procurator in execution with:
- ASCII logo display
- Step-by-step progress tracking
- Statistics output
- Results summary with generated files
- Total execution time

**To capture this screenshot:**

```bash
# Run the find_orfs command
conda activate procurator
python run_procurator.py find_orfs \
    -i data/test_genome.fasta \
    -o_gff test/genes.gff \
    -o_faa test/proteins.faa

# Take a screenshot of the terminal output
# On macOS: Cmd+Shift+5 or use Screenshot app
# On Linux: gnome-screenshot or similar
# On Windows: Win+Shift+S

# Save as: procurator-run.png
```

## Guidelines

- Keep images at a reasonable resolution (1200x800 or similar)
- Ensure terminal background is dark for better contrast
- Include the full output from command start to finish
- Use clear, readable fonts in the terminal
