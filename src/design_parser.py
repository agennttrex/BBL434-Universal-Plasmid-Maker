"""Parser for design files specifying MCS sites and markers."""

from typing import List, Tuple, Dict


def parse_design_file(file_path: str) -> Dict[str, List[Tuple[str, str]]]:
    """
    Parse a design file to extract MCS sites and markers.
    
    Format:
    Multiple_Cloning_Site1, RestrictionEnzyme1
    ...
    Antibiotic_marker1, Antibiotic Name 1
    ...
    
    Args:
        file_path: Path to the design file
        
    Returns:
        Dictionary with keys:
        - 'mcs_sites': List of (site_name, enzyme_name) tuples
        - 'markers': List of (marker_name, marker_type) tuples
    """
    mcs_sites = []
    markers = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Split by comma
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 2:
                continue
            
            site_name = parts[0]
            enzyme_or_type = parts[1]
            
            # Determine if it's an MCS site or marker based on naming convention
            # MCS sites typically end with "_site" or contain enzyme names
            # Markers typically end with "_gene", "_alpha", or contain antibiotic names
            if '_site' in site_name.lower() or any(
                enzyme in site_name.lower() 
                for enzyme in ['bamhi', 'hindiii', 'ecori', 'psti', 'kpn', 'sac', 'sal', 'xba', 'sma', 'sph', 'not']
            ):
                mcs_sites.append((site_name, enzyme_or_type))
            else:
                markers.append((site_name, enzyme_or_type))
    
    return {
        'mcs_sites': mcs_sites,
        'markers': markers
    }
