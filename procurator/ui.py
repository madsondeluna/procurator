# procurator/ui.py

"""
User Interface module for better visualization and progress tracking.
Provides formatted output, progress bars, and timing information.
"""

import sys
import time
from datetime import timedelta
from typing import Optional, Callable


class ProgressBar:
    """Simple progress bar with timing information."""
    
    def __init__(self, total: int, title: str = "Processing", width: int = 40):
        self.total = total
        self.title = title
        self.width = width
        self.current = 0
        self.start_time = time.time()
        self.last_update = self.start_time
        
    def update(self, increment: int = 1):
        """Update progress by increment."""
        self.current = min(self.current + increment, self.total)
        self._draw()
        
    def set_current(self, value: int):
        """Set current progress to specific value."""
        self.current = min(value, self.total)
        self._draw()
    
    def _draw(self):
        """Draw progress bar to stdout."""
        elapsed = time.time() - self.start_time
        
        # Calculate remaining time estimate
        if self.current > 0 and self.current < self.total:
            rate = elapsed / self.current
            remaining = rate * (self.total - self.current)
        else:
            remaining = 0
        
        # Format times
        elapsed_str = _format_time(elapsed)
        remaining_str = _format_time(remaining)
        
        # Calculate percentage
        percent = (self.current / self.total * 100) if self.total > 0 else 0
        filled = int(self.width * self.current / self.total) if self.total > 0 else 0
        bar = '█' * filled + '░' * (self.width - filled)
        
        # Print progress
        line = f"\r{self.title}: [{bar}] {percent:6.1f}% | {self.current:>5}/{self.total:<5} | ⏱ {elapsed_str} | ⏰ {remaining_str}"
        sys.stdout.write(line)
        sys.stdout.flush()
    
    def finish(self):
        """Mark progress as complete."""
        self.current = self.total
        self._draw()
        elapsed = time.time() - self.start_time
        sys.stdout.write(f"\n✓ {self.title} completed in {_format_time(elapsed)}\n")
        sys.stdout.flush()


class StepTracker:
    """Track and display execution steps with timing."""
    
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
        _print_step_start(name)
    
    def end_step(self, status: str = "✓ completed"):
        """End current step and record timing."""
        if self.current_step:
            elapsed = time.time() - self.current_start
            self.steps.append({
                'name': self.current_step,
                'duration': elapsed,
                'status': status
            })
            _print_step_end(self.current_step, elapsed, status)
        
        self.current_step = None
        self.current_start = None
    
    def print_summary(self):
        """Print summary of all steps."""
        if not self.steps:
            return
        
        total_time = sum(step['duration'] for step in self.steps)
        
        print("\n" + "=" * 60)
        print("📊 EXECUTION SUMMARY")
        print("=" * 60)
        
        for i, step in enumerate(self.steps, 1):
            percent = (step['duration'] / total_time * 100) if total_time > 0 else 0
            bar_width = int(30 * step['duration'] / max(total_time, 0.001))
            bar = '█' * bar_width
            
            print(f"{i}. {step['name']:<30} {step['status']:<15} {_format_time(step['duration']):>10} ({percent:>5.1f}%)")
            print(f"   └─ {bar}")
        
        print("=" * 60)
        print(f"Total Time: {_format_time(total_time)}")
        print("=" * 60 + "\n")


class DataFormatter:
    """Format data for nice terminal display."""
    
    @staticmethod
    def format_stats_table(stats_dict: dict) -> str:
        """Format statistics dictionary as aligned table."""
        lines = []
        max_key_len = max(len(k) for k in stats_dict.keys())
        
        lines.append("📋 Statistics Summary:")
        lines.append("─" * 60)
        
        for key, value in stats_dict.items():
            formatted_key = key.replace('_', ' ').title()
            line = f"  {formatted_key:<{max_key_len+5}} : {value:>20}"
            lines.append(line)
        
        lines.append("─" * 60)
        return "\n".join(lines)
    
    @staticmethod
    def format_results_summary(title: str, items: dict) -> str:
        """Format results summary."""
        lines = []
        lines.append(f"\n✨ {title}")
        lines.append("─" * 60)
        
        for key, value in items.items():
            icon = "📄" if "file" in key.lower() else "🔢"
            lines.append(f"  {icon} {key:<30} : {value}")
        
        lines.append("─" * 60)
        return "\n".join(lines)


def _format_time(seconds: float) -> str:
    """Format seconds into human-readable time string."""
    if seconds < 1:
        return f"{seconds*1000:.0f}ms"
    
    td = timedelta(seconds=int(seconds))
    hours, remainder = divmod(int(td.total_seconds()), 3600)
    minutes, seconds = divmod(remainder, 60)
    
    if hours > 0:
        return f"{hours}h {minutes}m {seconds}s"
    elif minutes > 0:
        return f"{minutes}m {seconds}s"
    else:
        return f"{seconds}s"


def _print_step_start(name: str):
    """Print step start message."""
    print(f"\n🔹 {name}...")


def _print_step_end(name: str, elapsed: float, status: str):
    """Print step end message."""
    print(f"  ✓ {status} ({_format_time(elapsed)})")


def print_header(title: str, width: int = 60):
    """Print a formatted header."""
    print("\n" + "=" * width)
    print(f"  {title}")
    print("=" * width + "\n")


def print_info(message: str, icon: str = "ℹ"):
    """Print info message."""
    print(f"{icon} {message}")


def print_success(message: str):
    """Print success message."""
    print(f"✅ {message}")


def print_error(message: str):
    """Print error message."""
    print(f"❌ {message}")


def print_warning(message: str):
    """Print warning message."""
    print(f"⚠️  {message}")
