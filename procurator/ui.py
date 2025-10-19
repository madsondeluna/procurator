# procurator/ui.py

"""
User Interface module for minimal, clean output.
Simple progress tracking and data display without excessive visuals.
"""

import sys
import time
from datetime import timedelta
from typing import Optional
from . import logo


class StepTracker:
    """Track and display execution steps with timing."""
    
    BOX_WIDTH = 76
    
    def __init__(self):
        self.steps = []
        self.current_step = None
        self.current_start = None
        
    def start_step(self, name: str):
        """Start a new step."""
        if self.current_step:
            self.end_step()
        
        self.current_step = name
        self.current_start = time.time()
        print(f"\n  > {name}...")
    
    def end_step(self, status: str = "[OK]"):
        """End current step and record timing."""
        if self.current_step:
            elapsed = time.time() - self.current_start
            self.steps.append({
                'name': self.current_step,
                'duration': elapsed,
                'status': status
            })
            print(f"    {status} {_format_time(elapsed)}")
        
        self.current_step = None
        self.current_start = None
    
    def print_summary(self):
        """Print simple summary of all steps."""
        if not self.steps:
            return
        
        total_time = sum(step['duration'] for step in self.steps)
        
        print("\n" + "=" * self.BOX_WIDTH)
        total_line = f"  Total execution time: {_format_time(total_time)}"
        print(total_line)
        print("=" * self.BOX_WIDTH + "\n")


class DataFormatter:
    """Format data for clean terminal display."""
    
    BOX_WIDTH = 76
    
    @staticmethod
    def format_stats_table(stats_dict: dict) -> str:
        """Format statistics as simple key-value pairs."""
        lines = []
        width = DataFormatter.BOX_WIDTH
        
        lines.append("\n" + "=" * width)
        lines.append("  STATISTICS")
        lines.append("=" * width)
        
        for key, value in stats_dict.items():
            label = key.replace('_', ' ').title()
            line = f"  {label:<50} {value:>22}"
            lines.append(line)
        
        lines.append("=" * width + "\n")
        return "\n".join(lines)
    
    @staticmethod
    def format_results_summary(title: str, items: dict) -> str:
        """Format results as simple key-value pairs."""
        lines = []
        width = DataFormatter.BOX_WIDTH
        
        lines.append("\n" + "=" * width)
        lines.append(f"  {title.upper()}")
        lines.append("=" * width)
        
        for key, value in items.items():
            line = f"  {key:<50} {value:>22}"
            lines.append(line)
        
        lines.append("=" * width + "\n")
        return "\n".join(lines)
    
    @staticmethod
    def format_files_list(title: str, files: list) -> str:
        """Format a list of generated files."""
        lines = []
        width = DataFormatter.BOX_WIDTH
        
        lines.append("\n" + "=" * width)
        lines.append(f"  {title.upper()}")
        lines.append("=" * width)
        
        for file_path in files:
            line = f"  • {file_path}"
            lines.append(line)
        
        lines.append("=" * width + "\n")
        return "\n".join(lines)


class ProgressBar:
    """Simple progress tracking - removed bars for minimal output."""
    
    def __init__(self, total: int, title: str = "Processing", width: int = 40):
        self.total = total
        self.title = title
        self.width = width
        self.current = 0
        self.start_time = time.time()
        
    def update(self, increment: int = 1):
        """Update progress by increment."""
        self.current = min(self.current + increment, self.total)
        
    def set_current(self, value: int):
        """Set current progress to specific value."""
        self.current = min(value, self.total)
    
    def finish(self):
        """Mark progress as complete."""
        self.current = self.total
        elapsed = time.time() - self.start_time
        print(f"  {self.title:<30} [OK] {_format_time(elapsed)}")



def _format_time(seconds: float) -> str:
    """Format seconds into human-readable time string."""
    if seconds < 1:
        return f"{seconds*1000:.0f}ms"
    
    td = timedelta(seconds=int(seconds))
    hours, remainder = divmod(int(td.total_seconds()), 3600)
    minutes, secs = divmod(remainder, 60)
    
    if hours > 0:
        return f"{hours}h {minutes}m {secs}s"
    elif minutes > 0:
        return f"{minutes}m {secs}s"
    else:
        return f"{secs}s"


def print_header(title: str = ""):
    """Print header with logo."""
    logo.print_header_modern()


def print_info(message: str):
    """Print info message."""
    print(f"  INFO: {message}")


def print_success(message: str):
    """Print success message."""
    print(f"  SUCCESS: {message}")


def print_error(message: str):
    """Print error message."""
    print(f"  ERROR: {message}")


def print_warning(message: str):
    """Print warning message."""
    print(f"  WARNING: {message}")

