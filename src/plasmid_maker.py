"""Main orchestrator for plasmid construction."""

import sys
import os
from typing import Dict, List, Tuple

# Add src to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.fasta_parser import read_fasta, write_fasta
from src.design_parser import parse_design_file
from src.markers_db import parse_markers_tab
from src.ori_finder import find_ori, extract_ori_sequence
from src.plasmid_builder import build_plasmid_sequence
from src.restriction_handler import delete_restriction_sites, verify_site_deletion


class PlasmidMaker:
    """Main class for constructing plasmids from genomic DNA and design files."""
    
    def __init__(self, markers_db_path: str):
        """
        Initialize PlasmidMaker with markers database.
        
        Args:
            markers_db_path: Path to markers.tab file
        """
        self.markers_db = parse_markers_tab(markers_db_path)
    
    def make_plasmid(
        self,
        input_fasta: str,
        design_file: str,
        output_fasta: str,
        delete_sites: bool = True
    ) -> str:
        """
        Generate a plasmid sequence from inputs.
        
        Args:
            input_fasta: Path to input genomic DNA FASTA file
            design_file: Path to design file
            output_fasta: Path to output plasmid FASTA file
            delete_sites: Whether to delete restriction sites specified in design
            
        Returns:
            Generated plasmid sequence
        """
        # 1. Read genomic DNA
        print(f"Reading genomic DNA from {input_fasta}...")
        header, genomic_seq = read_fasta(input_fasta)
        print(f"Genomic sequence length: {len(genomic_seq)} bp")
        
        # 2. Parse design file
        print(f"Parsing design file {design_file}...")
        design = parse_design_file(design_file)
        mcs_sites = design['mcs_sites']
        markers = design['markers']
        print(f"Found {len(mcs_sites)} MCS sites and {len(markers)} markers")
        
        # 3. Find ORI in genomic DNA
        print("Identifying origin of replication...")
        ori_start, ori_end, method = find_ori(genomic_seq)
        ori_sequence = extract_ori_sequence(genomic_seq, ori_start, ori_end)
        print(f"ORI found at positions {ori_start}-{ori_end} using method: {method}")
        print(f"ORI sequence length: {len(ori_sequence)} bp")
        
        # 4. Build plasmid sequence
        print("Building plasmid sequence...")
        plasmid_seq = build_plasmid_sequence(
            ori_sequence=ori_sequence,
            mcs_sites=mcs_sites,
            markers=markers,
            markers_db=self.markers_db,
            include_default_genes=True
        )
        
        # 5. Delete restriction sites if specified
        if delete_sites:
            print("Deleting restriction sites not in design...")
            # Get all known restriction enzymes from markers database
            all_enzymes = [name for name, info in self.markers_db.items() 
                          if 'Restriction enzyme' in info.get('category', '')]
            
            # Enzymes in the design file should be kept (they're in the MCS)
            enzymes_in_design = {enzyme.upper() for _, enzyme in mcs_sites}
            
            # Find enzymes to delete (those NOT in the design)
            enzymes_to_delete = []
            for enzyme in all_enzymes:
                if enzyme.upper() not in enzymes_in_design:
                    enzymes_to_delete.append(enzyme)
            
            if enzymes_to_delete:
                print(f"Deleting sites for enzymes not in design: {', '.join(enzymes_to_delete)}")
                plasmid_seq = delete_restriction_sites(plasmid_seq, enzymes_to_delete, self.markers_db)
                
                # Verify deletions for key enzymes (like EcoRI for pUC19 test)
                for enzyme in enzymes_to_delete:
                    if not verify_site_deletion(plasmid_seq, enzyme, self.markers_db):
                        print(f"Warning: Could not verify deletion of {enzyme} sites")
            else:
                print("No restriction sites to delete")
        
        print(f"Final plasmid sequence length: {len(plasmid_seq)} bp")
        
        # 6. Write output
        print(f"Writing plasmid to {output_fasta}...")
        output_header = f"Plasmid constructed from {header} using design {design_file}"
        write_fasta(output_fasta, output_header, plasmid_seq)
        
        return plasmid_seq


def main():
    """Command-line interface."""
    if len(sys.argv) < 4:
        print("Usage: python -m src.plasmid_maker <input.fa> <design.txt> <output.fa> [markers.tab]")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    design_file = sys.argv[2]
    output_fasta = sys.argv[3]
    markers_db = sys.argv[4] if len(sys.argv) > 4 else 'data/markers.tab'
    
    if not os.path.exists(input_fasta):
        print(f"Error: Input FASTA file not found: {input_fasta}")
        sys.exit(1)
    
    if not os.path.exists(design_file):
        print(f"Error: Design file not found: {design_file}")
        sys.exit(1)
    
    if not os.path.exists(markers_db):
        print(f"Error: Markers database not found: {markers_db}")
        sys.exit(1)
    
    maker = PlasmidMaker(markers_db)
    maker.make_plasmid(input_fasta, design_file, output_fasta)


if __name__ == '__main__':
    main()
