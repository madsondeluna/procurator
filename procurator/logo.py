# procurator/logo.py

"""
Logo and branding for Procurator CLI
"""

# Compact ASCII logo (single line friendly)
LOGO_COMPACT = """
╔═══════════════════════════════════════════════════════════════════╗
║  PROCURATOR - Prokaryotic Genome Analysis Pipeline               ║
║  Version 1.0.0 (Beta)                                            ║
╚═══════════════════════════════════════════════════════════════════╝
"""

# Alternative ultra-compact logo
LOGO_INLINE = "[PROCURATOR] Prokaryotic Genome Analysis Pipeline v1.0.0"

# ASCII art banner (still compact)
LOGO_ART = r"""

 ______   ______     ______     ______     __  __   
/\  == \ /\  == \   /\  __ \   /\  ___\   /\ \/\ \  
\ \  _-/ \ \  __<   \ \ \/\ \  \ \ \____  \ \ \_\ \ 
 \ \_\    \ \_\ \_\  \ \_____\  \ \_____\  \ \_____\
  \/_/     \/_/ /_/   \/_____/   \/_____/   \/_____/   
 ______     ______     ______   ______     ______   
/\  == \   /\  __ \   /\__  _\ /\  __ \   /\  == \  
\ \  __<   \ \  __ \  \/_/\ \/ \ \ \/\ \  \ \  __<  
 \ \_\ \_\  \ \_\ \_\    \ \_\  \ \_____\  \ \_\ \_\
  \/_/ /_/   \/_/\/_/     \/_/   \/_____/   \/_/ /_/

"""

def print_header_modern():
    """Print modern header with logo."""
    print("\n" + LOGO_ART)
#    print("=" * 76)
    print("  PROCURATOR - Prokaryotic Genome Annotator Pipeline")
    print("  Version 1.0.0 (Beta)")
#    print("=" * 76 + "\n")


def print_header_compact():
    """Print compact inline header."""
    print("\n[PROCURATOR] Prokaryotic Genome Analysis Pipeline v1.0.0\n")


def get_section_header(title: str) -> str:
    """Get formatted section header."""
    return f"\n{'=' * 70}\n  {title}\n{'=' * 70}\n"


def get_subsection_header(title: str) -> str:
    """Get formatted subsection header."""
    return f"\n  {title}"


def get_box_line(text: str, width: int = 70) -> str:
    """Get a line within a box."""
    return f"  {text:<{width-4}}  "


def print_section_box(title: str, lines: list):
    """Print a formatted box with content."""
    print("\n" + "┌" + "─" * 68 + "┐")
    print(f"│  {title:<64}  │")
    print("├" + "─" * 68 + "┤")
    for line in lines:
        print(f"│  {line:<64}  │")
    print("└" + "─" * 68 + "┘\n")


def print_status_line(status: str, message: str, elapsed: str = ""):
    """Print a status line with formatting."""
    if elapsed:
        return f"  [{status}] {message:<45} ({elapsed})"
    else:
        return f"  [{status}] {message}"
