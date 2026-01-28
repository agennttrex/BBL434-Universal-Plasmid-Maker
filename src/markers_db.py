"""Parser and database for markers from markers.tab file."""

from typing import Dict, Optional
import re


def parse_markers_tab(file_path: str) -> Dict[str, Dict[str, str]]:
    """
    Parse markers.tab file to extract marker information.
    
    Args:
        file_path: Path to markers.tab file
        
    Returns:
        Dictionary mapping marker short names to their metadata
    """
    markers = {}
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        # Skip header line (first line with |)
        for line in lines[1:]:
            line = line.strip()
            if not line or not line.startswith('|'):
                continue
            
            # Parse table row: | Category | Name (short) | Recognition / Role | Typical use |
            parts = [p.strip() for p in line.split('|')]
            if len(parts) < 4:
                continue
            
            category = parts[1].strip()
            name_short = parts[2].strip()
            recognition = parts[3].strip()
            
            if name_short:
                markers[name_short] = {
                    'category': category,
                    'name': name_short,
                    'recognition': recognition
                }
    
    return markers


def get_restriction_site_sequence(enzyme_name: str, markers_db: Dict[str, Dict[str, str]]) -> Optional[str]:
    """
    Extract restriction site sequence from enzyme name using markers database.
    
    Args:
        enzyme_name: Name of restriction enzyme (e.g., 'EcoRI', 'BamHI')
        markers_db: Parsed markers database
        
    Returns:
        Recognition sequence if found, None otherwise
    """
    # Try exact match first
    if enzyme_name in markers_db:
        recognition = markers_db[enzyme_name]['recognition']
        # Extract sequence from recognition string (e.g., "Recognizes GAATTC, 5′ overhangs")
        match = re.search(r'([ATCG]{4,})', recognition.upper())
        if match:
            return match.group(1)
    
    # Try case-insensitive match
    for key, value in markers_db.items():
        if key.upper() == enzyme_name.upper():
            recognition = value['recognition']
            match = re.search(r'([ATCG]{4,})', recognition.upper())
            if match:
                return match.group(1)
    
    # Fallback: known restriction sites
    known_sites = {
        'ECORI': 'GAATTC',
        'BAMHI': 'GGATCC',
        'HINDIII': 'AAGCTT',
        'PSTI': 'CTGCAG',
        'KPRI': 'GGTACC',
        'SACI': 'GAGCTC',
        'SALI': 'GTCGAC',
        'XBAI': 'TCTAGA',
        'NOTI': 'GCGGCCGC',
        'SMAI': 'CCCGGG',
        'SPHI': 'GCATGC'
    }
    
    return known_sites.get(enzyme_name.upper())


def get_marker_sequence(marker_name: str, markers_db: Dict[str, Dict[str, str]]) -> Optional[str]:
    """
    Get marker sequence if available in database.
    For now, returns None as sequences are not in the tab file.
    This can be extended with a sequence database.
    
    Args:
        marker_name: Name of marker (e.g., 'AmpR', 'lacZα')
        markers_db: Parsed markers database
        
    Returns:
        Marker sequence if available, None otherwise
    """
    # Marker sequences would need to be in a separate database
    # For now, we'll handle this in the plasmid builder with known sequences
    return None
