# procurator/protein.py

import sys
from . import ui
from . import io

def select_proteins(args):
    """
    Selects proteins from a FASTA file based on keywords in annotations.
    """
    tracker = ui.StepTracker()
    tracker.start_step("Loading sequences")
    
    sequences = io.read_fasta_sequences(args.input_faa)
    if not sequences:
        ui.print_error("No sequences found in the input file.")
        return
    
    tracker.end_step()
    tracker.start_step("Filtering by keywords")
    
    keywords = [kw.strip().lower() for kw in args.keywords.split(',')]
    
    selected = []
    for seq in sequences:
        desc = (seq.description or "").lower()
        if any(kw in desc for kw in keywords):
            selected.append(seq)
    
    tracker.end_step()
    tracker.start_step("Writing selected proteins")
    
    if selected:
        io.write_fasta(selected, args.output_fasta)
        ui.print_success(f"Selected {len(selected)} out of {len(sequences)} proteins")
    else:
        ui.print_warning(f"No proteins matched keywords: {', '.join(keywords)}")
    
    tracker.end_step()
    
    print(ui.DataFormatter.format_results_summary(
        "PROTEIN SELECTION COMPLETE",
        {
            "Input sequences": len(sequences),
            "Selected sequences": len(selected),
            "Keywords": ", ".join(keywords),
            "Output file": args.output_fasta,
        }
    ))
    
    tracker.print_summary()

def run_modeling(args):
    """
    Executes structural modeling using ColabFold.
    """
    tracker = ui.StepTracker()
    tracker.start_step("Validating input")
    
    sequences = io.read_fasta_sequences(args.input_fasta)
    if not sequences:
        ui.print_error("No sequences found in the input file.")
        return
    
    ui.print_warning("ColabFold structural modeling not yet implemented in this environment.")
    ui.print_info(f"To use: install ColabFold and configure the database path")
    ui.print_info(f"Input file: {args.input_fasta}")
    ui.print_info(f"Output directory: {args.output_dir}")
    ui.print_info(f"Number of sequences: {len(sequences)}")
    
    tracker.end_step("skipped (ColabFold not available)")
    tracker.print_summary()
