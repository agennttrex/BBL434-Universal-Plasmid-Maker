"""Handler for restriction enzyme sites and their deletion."""

from typing import List, Set, Dict
import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.markers_db import get_restriction_site_sequence


def find_restriction_sites(sequence: str, enzyme_name: str, markers_db: Dict) -> List[int]:
    """
    Find all occurrences of a restriction site in a sequence.
    
    Args:
        sequence: DNA sequence to search
        enzyme_name: Name of restriction enzyme
        markers_db: Markers database
        
    Returns:
        List of start positions where the site occurs
    """
    site_sequence = get_restriction_site_sequence(enzyme_name, markers_db)
    if not site_sequence:
        return []
    
    positions = []
    seq_upper = sequence.upper()
    site_upper = site_sequence.upper()
    
    start = 0
    while True:
        pos = seq_upper.find(site_upper, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    
    return positions


def delete_restriction_sites(sequence: str, enzymes_to_delete: List[str], markers_db: Dict) -> str:
    """
    Delete restriction sites from sequence by mutating them.
    
    Strategy: Mutate one base in the recognition sequence to prevent recognition
    while maintaining sequence integrity. We mutate the last base of the site.
    
    Args:
        sequence: DNA sequence
        enzymes_to_delete: List of enzyme names whose sites should be deleted
        markers_db: Markers database
        
    Returns:
        Modified sequence with restriction sites mutated
    """
    result = list(sequence)
    sites_to_mutate = set()
    
    for enzyme in enzymes_to_delete:
        positions = find_restriction_sites(sequence, enzyme, markers_db)
        site_seq = get_restriction_site_sequence(enzyme, markers_db)
        if not site_seq:
            continue
        
        for pos in positions:
            # Mark all positions in this site for mutation
            for i in range(len(site_seq)):
                sites_to_mutate.add(pos + i)
    
    # Mutate sites: change last base to a different nucleotide
    for pos in sorted(sites_to_mutate, reverse=True):
        if pos >= len(result):
            continue
        
        current_base = result[pos].upper()
        # Mutate to a different base
        if current_base == 'A':
            result[pos] = 'G'
        elif current_base == 'T':
            result[pos] = 'C'
        elif current_base == 'G':
            result[pos] = 'A'
        elif current_base == 'C':
            result[pos] = 'T'
    
    return ''.join(result)


def verify_site_deletion(sequence: str, enzyme_name: str, markers_db: Dict) -> bool:
    """
    Verify that a restriction site has been deleted from the sequence.
    
    Args:
        sequence: DNA sequence to check
        enzyme_name: Name of restriction enzyme
        markers_db: Markers database
        
    Returns:
        True if site is not found, False if still present
    """
    positions = find_restriction_sites(sequence, enzyme_name, markers_db)
    return len(positions) == 0
