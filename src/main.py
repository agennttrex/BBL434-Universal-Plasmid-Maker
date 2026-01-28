"""Main entry point for Plasmid Maker."""

import sys
import os

# Add project root to path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

from src.plasmid_maker import PlasmidMaker, main

if __name__ == '__main__':
    # If run as script, use command-line interface
    if len(sys.argv) > 1:
        main()
    else:
        # Example usage
        print("Plasmid Maker - Universal Plasmid Construction Tool")
        print("\nUsage:")
        print("  python src/main.py <input.fa> <design.txt> <output.fa> [markers.tab]")
        print("\nExample:")
        print("  python src/main.py data/pUC19.fa data/Design_pUC19.txt Output.fa data/markers.tab")
