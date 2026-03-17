#!/usr/bin/env python3
"""
Step 2: Reverse Translation
Generate candidate CDS sequences from amino acid sequence
Using codon sampling strategy to avoid exponential explosion
"""

import sys
import random
import json
from pathlib import Path

# Standard genetic code ( codon -> amino acid )
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Reverse mapping: amino acid -> list of codons
AA_TO_CODONS = {}
for codon, aa in CODON_TABLE.items():
    if aa not in AA_TO_CODONS:
        AA_TO_CODONS[aa] = []
    AA_TO_CODONS[aa].append(codon)

# Remove stop codons for CDS
STOP_CODONS = ['TAA', 'TAG', 'TGA']

def reverse_translate_random(aa_seq, num_samples=100):
    """Random sampling-based reverse translation"""
    candidates = []
    
    for i in range(num_samples):
        codons = []
        for aa in aa_seq:
            if aa == '*':
                continue  # Skip stop codons in the middle
            codon_list = AA_TO_CODONS.get(aa, [])
            if codon_list:
                codons.append(random.choice(codon_list))
            else:
                # Fallback: use NNN for unknown
                codons.append('NNN')
        
        # Add stop codon at the end
        codons.append('TAA')
        cds = ''.join(codons)
        candidates.append(cds)
    
    return candidates

def reverse_translate_greedy(aa_seq, weights=None):
    """
    Greedy reverse translation using codon weights
    Prefer codons with higher weights (e.g., CAI-optimized)
    """
    if weights is None:
        # Default: use first codon (deterministic)
        weights = {codon: 1.0 for codon in CODON_TABLE}
    
    codons = []
    for aa in aa_seq:
        if aa == '*':
            continue
        codon_list = AA_TO_CODONS.get(aa, [])
        if codon_list:
            # Select codon with highest weight
            best_codon = max(codon_list, key=lambda c: weights.get(c, 0))
            codons.append(best_codon)
        else:
            codons.append('NNN')
    
    codons.append('TAA')  # Stop codon
    return [''.join(codons)]

def load_weights(weight_file=None):
    """Load codon weights from file"""
    weights = {}
    
    if weight_file and Path(weight_file).exists():
        with open(weight_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split(',')
                    if len(parts) >= 2:
                        weights[parts[0]] = float(parts[1])
    
    return weights if weights else None

def main():
    # Read amino acid sequence from stdin
    aa_seq = sys.stdin.read().strip()
    
    # Remove header if present
    if aa_seq.startswith('>'):
        aa_seq = aa_seq.split('\n', 1)[1]
    aa_seq = aa_seq.replace('\n', '').replace(' ', '')
    
    if not aa_seq:
        print("Error: Empty sequence", file=sys.stderr)
        sys.exit(1)
    
    # Get parameters
    num_samples = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    weight_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    # Load weights
    weights = load_weights(weight_file)
    
    # Generate candidates
    # Mix of random sampling and greedy optimization
    candidates = []
    
    # Add some random samples
    candidates.extend(reverse_translate_random(aa_seq, num_samples // 2))
    
    # Add greedy solution
    if weights:
        candidates.extend(reverse_translate_greedy(aa_seq, weights))
    else:
        # Random with different seeds
        for seed in range(3):
            random.seed(seed)
            candidates.extend(reverse_translate_random(aa_seq, num_samples // 4))
    
    # Remove duplicates
    candidates = list(set(candidates))
    
    # Output as FASTA
    for i, cds in enumerate(candidates):
        print(f">candidate_{i+1}")
        print(cds)

if __name__ == '__main__':
    main()
