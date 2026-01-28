# Universal Plasmid Maker

A bioinformatics tool for constructing functional plasmid DNA sequences from genomic DNA and design specifications.

## Overview

This tool takes a genomic DNA sequence from an unknown organism and a design file, then generates a functional plasmid sequence that can replicate in that organism. The solution identifies the origin of replication (ORI) from the genomic sequence and assembles a complete plasmid with specified multiple cloning sites (MCS) and antibiotic resistance markers.

## Features

- **ORI Identification**: Automatically identifies origin of replication using multiple bioinformatics methods:
  - DnaA box pattern recognition
  - AT-rich region detection
  - Fallback default region selection

- **Modular Design**: Clean, modular codebase with separate components for:
  - FASTA file parsing
  - Design file parsing
  - Markers database management
  - ORI finding
  - Plasmid construction
  - Restriction site handling

- **Graceful Error Handling**: Handles missing markers and invalid inputs gracefully

- **Restriction Site Management**: Automatically deletes restriction sites not specified in the design file

- **Comprehensive Testing**: Includes unit tests and integration tests, including pUC19 test case

## Project Structure

```
CC1/
├── src/
│   ├── __init__.py
│   ├── main.py                 # Main entry point
│   ├── fasta_parser.py         # FASTA file I/O
│   ├── design_parser.py        # Design file parsing
│   ├── markers_db.py           # Markers database parser
│   ├── ori_finder.py           # ORI identification methods
│   ├── plasmid_builder.py      # Plasmid sequence assembly
│   ├── restriction_handler.py  # Restriction site management
│   └── plasmid_maker.py        # Main orchestrator
├── tests/
│   ├── __init__.py
│   └── test_plasmid_maker.py   # Comprehensive test suite
├── data/
│   ├── pUC19.fa                # Test case input
│   ├── Design_pUC19.txt        # Test case design file
│   └── markers.tab             # Markers database
├── FEMS_Final_2013.pdf         # Reference paper
└── README.md                   # This file
```

## Installation

No external dependencies required - uses only Python standard library.

```bash
# Clone or navigate to the project directory
cd CC1
```

## Usage

### Command Line Interface

```bash
python src/main.py <input.fa> <design.txt> <output.fa> [markers.tab]
```

**Arguments:**
- `input.fa`: Path to input genomic DNA FASTA file
- `design.txt`: Path to design file specifying MCS sites and markers
- `output.fa`: Path to output plasmid FASTA file
- `markers.tab`: (Optional) Path to markers database (default: `data/markers.tab`)

### Example

```bash
python src/main.py data/pUC19.fa data/Design_pUC19.txt Output.fa data/markers.tab
```

### Programmatic Usage

```python
from src.plasmid_maker import PlasmidMaker

maker = PlasmidMaker('data/markers.tab')
plasmid_seq = maker.make_plasmid(
    input_fasta='data/pUC19.fa',
    design_file='data/Design_pUC19.txt',
    output_fasta='Output.fa',
    delete_sites=True
)
```

## Design File Format

The design file specifies multiple cloning sites and markers:

```
Multiple_Cloning_Site1, RestrictionEnzyme1
Multiple_Cloning_Site2, RestrictionEnzyme2
...
Antibiotic_marker1, Antibiotic Name 1
Antibiotic_marker2, Antibiotic Name 2
...
```

**Example (`Design_pUC19.txt`):**
```
BamHI_site, BamHI
HindIII_site, HindIII
PstI_site, PstI
SphI_site, SphI
SalI_site, SalI
XbaI_site, XbaI
KpnI_site, KpnI
SacI_site, SacI
SmaI_site, SmaI
AmpR_gene, Ampicillin
lacZ_alpha, Blue_White_Selection
ori_pMB1, High_Copy_Replication
```

## How It Works

1. **Read Inputs**: Parses genomic DNA FASTA file and design file
2. **Identify ORI**: Uses bioinformatics methods to find origin of replication:
   - Searches for DnaA box patterns (TTATCCACA variants)
   - Identifies AT-rich regions (>65% AT content)
   - Falls back to default region if no patterns found
3. **Build Plasmid**: Assembles components in order:
   - ORI sequence (from genomic DNA)
   - Default replication genes (minimal sequences for maintenance)
   - Multiple cloning site (from design file)
   - Antibiotic resistance markers (from design file)
4. **Delete Sites**: Removes restriction enzyme sites not specified in design file
5. **Output**: Writes complete plasmid sequence to FASTA file

## Quick Start

### Test with pUC19 Example

```bash
# Run the test script
./run_test.sh

# Or manually:
python src/main.py data/pUC19.fa data/Design_pUC19.txt Output.fa data/markers.tab
```

This will:
1. Read pUC19.fa as genomic DNA input
2. Parse Design_pUC19.txt for MCS sites and markers
3. Identify ORI from the genomic sequence
4. Build a plasmid with specified components
5. Delete EcoRI sites (not in design file)
6. Output to Output.fa

Verify EcoRI deletion:
```bash
python -c "
from src.fasta_parser import read_fasta
from src.markers_db import parse_markers_tab
from src.restriction_handler import find_restriction_sites
markers_db = parse_markers_tab('data/markers.tab')
_, seq = read_fasta('Output.fa')
sites = find_restriction_sites(seq, 'EcoRI', markers_db)
print(f'EcoRI sites: {len(sites)} (should be 0)')
"
```

## Testing

Run the test suite:

```bash
python -m pytest tests/ -v
```

Or using unittest:

```bash
python -m unittest tests.test_plasmid_maker -v
```

### Test Cases

- **FASTA Parser Tests**: Verify reading and writing FASTA files
- **Design Parser Tests**: Verify parsing of design files
- **Markers DB Tests**: Verify markers database parsing
- **ORI Finder Tests**: Verify ORI identification methods
- **Restriction Handler Tests**: Verify site finding and deletion
- **Plasmid Builder Tests**: Verify plasmid assembly
- **pUC19 Integration Test**: Full pipeline test with pUC19 data
  - Verifies EcoRI sites are deleted (not in design file)
  - Verifies plasmid construction completes successfully

## Assumptions and Design Decisions

### ORI Identification

1. **DnaA Boxes**: Primary method searches for TTATCCACA patterns and variants
   - Requires at least 3 matches within 500bp
   - Extends region to ~200-300bp (typical ORI size)

2. **AT-Rich Regions**: Fallback method finds regions with >65% AT content
   - Uses 200bp sliding window
   - Selects region with highest AT content

3. **Default Fallback**: If no patterns found, uses first 300bp
   - Conservative approach ensures ORI is always found
   - In production, would use more sophisticated methods (e.g., machine learning)

### Plasmid Construction

1. **Component Order**: ORI → Replication genes → MCS → Markers
   - Standard plasmid architecture
   - Spacers (2bp) between components

2. **Default Replication Genes**: Includes minimal sequences for plasmid maintenance
   - Simplified for this implementation
   - In production, would include full gene sequences

3. **Marker Sequences**: Uses known sequences for common markers (AmpR, lacZα)
   - Falls back gracefully if marker not found
   - Warns but continues construction

### Restriction Site Deletion

1. **Deletion Strategy**: Mutates restriction sites by changing one base
   - Changes last base of recognition sequence
   - Prevents enzyme recognition while maintaining sequence integrity

2. **Which Sites to Delete**: Deletes sites for enzymes NOT in design file
   - Design file specifies which enzymes should be in MCS
   - All other restriction sites are removed from entire plasmid

3. **Verification**: Checks that sites are successfully deleted
   - Warns if deletion cannot be verified

### Error Handling

1. **Missing Markers**: Warns but continues construction
   - Skips missing markers
   - Still produces functional plasmid

2. **Missing Restriction Sites**: Warns but continues
   - Skips enzymes not found in database
   - Uses fallback known sequences when possible

3. **Invalid Inputs**: Validates file existence and formats
   - Provides clear error messages
   - Exits gracefully with error codes

## Limitations and Future Improvements

1. **ORI Identification**: Current methods are simplified
   - Could use machine learning models trained on known ORIs
   - Could incorporate sequence conservation analysis
   - Could use organism-specific ORI databases

2. **Marker Sequences**: Limited to hardcoded sequences
   - Should use comprehensive marker sequence database
   - Could fetch from online databases (e.g., NCBI)

3. **Replication Genes**: Uses minimal placeholders
   - Should include full gene sequences
   - Should be organism-specific

4. **Restriction Site Deletion**: Simple mutation strategy
   - Could use more sophisticated site removal
   - Could preserve codon usage if in coding regions

5. **Validation**: Limited sequence validation
   - Could check for proper plasmid structure
   - Could validate against known plasmid databases

## References

- FEMS Microbiology Reviews (2013): "Broad host range plasmids" - Reference paper on plasmid replication requirements
- Standard bioinformatics methods for ORI identification
- Common plasmid architecture and design principles

## License

This project is created for educational purposes as part of a bioinformatics coding challenge.

## Author

Created for BBL434 Coding Challenge 1 - Universal Plasmid Maker
