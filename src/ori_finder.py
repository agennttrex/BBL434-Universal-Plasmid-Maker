"""ORI (Origin of Replication) finder using bioinformatics methods."""

from typing import Optional, Tuple
import re


def find_dnaa_boxes(sequence: str, min_matches: int = 2) -> Optional[Tuple[int, int]]:
    """
    Find DnaA box patterns in the sequence.
    
    DnaA boxes are typically: TTATCCACA or similar variants (TTATNCACA)
    These are binding sites for DnaA protein that initiates replication.
    
    Args:
        sequence: DNA sequence to search
        min_matches: Minimum number of DnaA boxes to find in a region
        
    Returns:
        Tuple of (start, end) positions if found, None otherwise
    """
    # DnaA box consensus: TTATCCACA or TTATNCACA (N = any nucleotide)
    # Also look for common variants and related sequences
    dnaa_patterns = [
        r'TTATCCACA',      # Exact consensus
        r'TTAT[ATCG]CACA',  # N variant
        r'TTATCC[ATCG]CA',  # N variant
        r'T[ATCG]ATCCACA',  # N variant
        r'TTAT[ATCG]{2}CACA',  # More flexible
        r'TTATCC',         # Partial match (common in ORIs)
    ]
    
    matches = []
    for pattern in dnaa_patterns:
        for match in re.finditer(pattern, sequence, re.IGNORECASE):
            matches.append(match.start())
    
    # Remove duplicates and sort
    matches = sorted(set(matches))
    
    if len(matches) >= min_matches:
        # Find region with multiple DnaA boxes
        # Look for clusters of DnaA boxes within 500bp
        for i in range(len(matches) - min_matches + 1):
            cluster = matches[i:i+min_matches]
            if cluster[-1] - cluster[0] < 500:
                # Extend region to include surrounding area (typical ORI is ~200-300bp)
                start = max(0, cluster[0] - 100)
                end = min(len(sequence), cluster[-1] + 200)
                return (start, end)
    
    # If we have at least one match, use it with extended region
    if len(matches) >= 1:
        match_pos = matches[0]
        start = max(0, match_pos - 150)
        end = min(len(sequence), match_pos + 250)
        return (start, end)
    
    return None


def find_at_rich_region(sequence: str, window_size: int = 200, at_threshold: float = 0.65) -> Optional[Tuple[int, int]]:
    """
    Find AT-rich regions which are characteristic of replication origins.
    
    Args:
        sequence: DNA sequence to search
        window_size: Size of sliding window
        at_threshold: Minimum AT content (0.0-1.0)
        
    Returns:
        Tuple of (start, end) positions if found, None otherwise
    """
    best_region = None
    best_at_content = 0.0
    
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        at_count = window.count('A') + window.count('T')
        at_content = at_count / len(window)
        
        if at_content >= at_threshold and at_content > best_at_content:
            best_at_content = at_content
            best_region = (i, i + window_size)
    
    if best_region:
        return best_region
    return None


def find_ori(sequence: str) -> Tuple[int, int, str]:
    """
    Find origin of replication in genomic DNA sequence.
    
    Uses multiple methods:
    1. DnaA box patterns (primary method)
    2. AT-rich regions (fallback)
    3. Default region if nothing found
    
    Args:
        sequence: Genomic DNA sequence
        
    Returns:
        Tuple of (start, end, method_used)
    """
    # Method 1: Look for DnaA boxes
    dnaa_result = find_dnaa_boxes(sequence)
    if dnaa_result:
        return (*dnaa_result, "dnaa_boxes")
    
    # Method 2: Look for AT-rich regions
    at_rich_result = find_at_rich_region(sequence)
    if at_rich_result:
        return (*at_rich_result, "at_rich")
    
    # Method 3: Default - use first 300bp (conservative fallback)
    # In practice, this would require more sophisticated methods
    default_size = min(300, len(sequence))
    return (0, default_size, "default")


def extract_ori_sequence(sequence: str, start: int, end: int) -> str:
    """
    Extract ORI sequence from genomic DNA.
    
    Args:
        sequence: Full genomic sequence
        start: Start position
        end: End position
        
    Returns:
        ORI sequence
    """
    return sequence[start:end]
