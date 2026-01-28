#!/bin/bash
# Quick test script for Plasmid Maker

echo "Running Plasmid Maker test with pUC19..."
python src/main.py data/pUC19.fa data/Design_pUC19.txt Output.fa data/markers.tab

echo ""
echo "Checking if EcoRI sites were deleted..."
python -c "
import sys
sys.path.insert(0, '.')
from src.fasta_parser import read_fasta
from src.markers_db import parse_markers_tab, get_restriction_site_sequence
from src.restriction_handler import find_restriction_sites

markers_db = parse_markers_tab('data/markers.tab')
header, seq = read_fasta('Output.fa')
ecoRI_sites = find_restriction_sites(seq, 'EcoRI', markers_db)

if len(ecoRI_sites) == 0:
    print('✓ SUCCESS: EcoRI sites successfully deleted from output')
else:
    print(f'✗ WARNING: Found {len(ecoRI_sites)} EcoRI site(s) in output')
    print(f'  Positions: {ecoRI_sites}')
"

echo ""
echo "Test complete. Check Output.fa for the generated plasmid."
