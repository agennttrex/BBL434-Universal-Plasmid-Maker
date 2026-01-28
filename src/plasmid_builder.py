"""Plasmid DNA sequence builder."""

from typing import List, Tuple, Dict, Optional
import sys
import os

# Add src to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.markers_db import get_restriction_site_sequence, get_marker_sequence


# Known marker sequences (these would ideally come from a database)
MARKER_SEQUENCES = {
    'AmpR': 'ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAACTGTCAGACCAAGTTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAG',
    'lacZ_alpha': 'ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC',
    # Minimal ORI sequences for common plasmid types
    'ori_pMB1': 'TTATCCACA',  # Minimal pMB1 ORI (simplified)
    'ori_pUC': 'TTATCCACA',
}


# Default replication genes (simplified - in practice these would be full sequences)
DEFAULT_REPLICATION_GENES = {
    'repA': 'ATG',  # Placeholder - would be full repA gene sequence
    'repB': 'ATG',  # Placeholder
}


def build_mcs_sequence(mcs_sites: List[Tuple[str, str]], markers_db: Dict) -> str:
    """
    Build multiple cloning site (MCS) sequence from specified sites.
    
    Args:
        mcs_sites: List of (site_name, enzyme_name) tuples
        markers_db: Markers database
        
    Returns:
        MCS sequence with restriction sites concatenated
    """
    mcs_sequence = []
    
    for site_name, enzyme_name in mcs_sites:
        site_seq = get_restriction_site_sequence(enzyme_name, markers_db)
        if site_seq:
            mcs_sequence.append(site_seq)
        else:
            # If enzyme not found, skip with warning (handled gracefully)
            print(f"Warning: Restriction site for {enzyme_name} not found, skipping")
    
    # Join sites with a short spacer (typically 2-4bp)
    spacer = 'AA'  # 2bp spacer between sites
    return spacer.join(mcs_sequence)


def get_marker_sequence_safe(marker_name: str, markers_db: Dict) -> Optional[str]:
    """
    Get marker sequence, handling missing markers gracefully.
    
    Args:
        marker_name: Name of marker
        markers_db: Markers database
        
    Returns:
        Marker sequence if available, None otherwise
    """
    # Try database first
    seq = get_marker_sequence(marker_name, markers_db)
    if seq:
        return seq
    
    # Try known sequences
    marker_key = marker_name.upper().replace('_', '').replace('-', '')
    for key, value in MARKER_SEQUENCES.items():
        if key.upper() == marker_key or marker_key in key.upper():
            return value
    
    # Try partial matches
    if 'AMP' in marker_key or 'AMPR' in marker_key:
        return MARKER_SEQUENCES.get('AmpR')
    if 'LACZ' in marker_key:
        return MARKER_SEQUENCES.get('lacZ_alpha')
    if 'ORI' in marker_key:
        return MARKER_SEQUENCES.get('ori_pMB1')
    
    return None


def build_plasmid_sequence(
    ori_sequence: str,
    mcs_sites: List[Tuple[str, str]],
    markers: List[Tuple[str, str]],
    markers_db: Dict,
    include_default_genes: bool = True
) -> str:
    """
    Build complete plasmid sequence from components.
    
    Structure:
    1. ORI sequence
    2. Default replication genes (if requested)
    3. MCS
    4. Markers
    
    Args:
        ori_sequence: Origin of replication sequence
        mcs_sites: List of MCS sites
        markers: List of markers
        markers_db: Markers database
        include_default_genes: Whether to include default replication genes
        
    Returns:
        Complete plasmid sequence
    """
    components = []
    
    # 1. ORI
    components.append(ori_sequence)
    
    # 2. Default replication genes (minimal sequences for plasmid maintenance)
    if include_default_genes:
        # Add minimal replication initiator sequences
        # In practice, these would be full gene sequences
        default_rep = 'ATG'  # Placeholder - minimal initiator
        components.append(default_rep)
    
    # 3. MCS
    mcs_seq = build_mcs_sequence(mcs_sites, markers_db)
    if mcs_seq:
        components.append(mcs_seq)
    
    # 4. Markers (skip ORI markers since we already use the found ORI)
    for marker_name, marker_type in markers:
        # Skip ORI markers - we already have ORI from genomic DNA
        if 'ori' in marker_name.lower() or 'replication' in marker_type.lower():
            print(f"Info: Skipping {marker_name} - using ORI from genomic DNA")
            continue
        
        marker_seq = get_marker_sequence_safe(marker_name, markers_db)
        if marker_seq:
            components.append(marker_seq)
        else:
            print(f"Warning: Marker sequence for {marker_name} not found, skipping")
    
    # Join components with short spacers
    spacer = 'AA'  # 2bp spacer between components
    plasmid_seq = spacer.join(components)
    
    return plasmid_seq
