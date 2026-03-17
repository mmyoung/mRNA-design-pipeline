#!/usr/bin/env python3
"""
Step 1: Parse Input
Parse amino acid sequence from FASTA or plain text format
"""

import sys
import re
from pathlib import Path

# Standard amino acids (20 + 3 ambiguous)
VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')
VALID_AA_EXT = set('ACDEFGHIKLMNPQRSTVWYXBZJUO')  # B, Z, J, U, O are ambiguous

def parse_fasta(file_path):
    """Parse FASTA format file"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Get first word after >
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def parse_plain_text(text):
    """Parse plain amino acid sequence"""
    # Remove any whitespace and non-alpha characters
    seq = re.sub(r'[^A-Z]', '', text.upper())
    return {'input': seq}

def validate_aa_sequence(seq):
    """Validate amino acid sequence"""
    invalid = []
    for aa in seq:
        if aa not in VALID_AA_EXT:
            invalid.append(aa)
    
    if invalid:
        print(f"Warning: Invalid amino acids found: {set(invalid)}", file=sys.stderr)
        # Remove invalid characters
        seq = ''.join([aa for aa in seq if aa in VALID_AA_EXT])
    
    return seq

def main():
    if len(sys.argv) < 2:
        print("Usage: parse_input.py <input_file>", file=sys.stderr)
        sys.exit(1)
    
    input_path = Path(sys.argv[1])
    
    if not input_path.exists():
        print(f"Error: File not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    
    # Determine format and parse
    if input_path.suffix.lower() in ['.fasta', '.fa', '.fna']:
        sequences = parse_fasta(input_path)
    else:
        # Try as plain text
        content = input_path.read_text()
        sequences = parse_plain_text(content)
    
    # Validate and output
    for seq_id, seq in sequences.items():
        validated = validate_aa_sequence(seq)
        print(f">{seq_id}")
        print(validated)

if __name__ == '__main__':
    main()
