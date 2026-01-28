"""FASTA file parser for reading genomic DNA sequences."""

from typing import List, Tuple


def read_fasta(file_path: str) -> Tuple[str, str]:
    """
    Read a FASTA file and return the header and sequence.
    
    Args:
        file_path: Path to the FASTA file
        
    Returns:
        Tuple of (header, sequence) where sequence is uppercase and without whitespace
    """
    header = ""
    sequence = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                header = line[1:]  # Remove '>' prefix
            else:
                sequence.append(line.upper())
    
    full_sequence = ''.join(sequence)
    return header, full_sequence


def write_fasta(file_path: str, header: str, sequence: str, line_length: int = 60):
    """
    Write a sequence to a FASTA file.
    
    Args:
        file_path: Path to output FASTA file
        header: Sequence header (without '>')
        sequence: DNA sequence
        line_length: Characters per line (default 60)
    """
    with open(file_path, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), line_length):
            f.write(sequence[i:i+line_length] + "\n")
