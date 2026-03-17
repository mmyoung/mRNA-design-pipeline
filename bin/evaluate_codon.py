#!/usr/bin/env python3
"""
Step 3.1: Evaluate Codon Optimization
Calculate CAI, tAI, ENC, GC content, and other codon-related scores
"""

import sys
import json
from pathlib import Path

# Codon usage table for human (H. sapiens)
# Values are relative adaptiveness (percentage of most frequent codon)
# Source: Kazusa Codon Usage Database
HUMAN_CODON_USAGE = {
    # Phe
    'TTT': 0.45, 'TTC': 0.55,
    # Leu
    'TTA': 0.07, 'TTG': 0.13, 'CTT': 0.13, 'CTC': 0.20,
    'CTA': 0.07, 'CTG': 0.40,
    # Ile
    'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16,
    # Met (start)
    'ATG': 1.00,
    # Val
    'GTT': 0.18, 'GTC': 0.22, 'GTA': 0.11, 'GTG': 0.49,
    # Ser
    'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06,
    'AGT': 0.15, 'AGC': 0.24,
    # Pro
    'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11,
    # Thr
    'ACT': 0.24, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.12,
    # Ala
    'GCT': 0.26, 'GCC': 0.40, 'GCA': 0.23, 'GCG': 0.11,
    # Tyr
    'TAT': 0.44, 'TAC': 0.56,
    # Stop
    'TAA': 0.30, 'TAG': 0.24, 'TGA': 0.46,
    # His
    'CAT': 0.41, 'CAC': 0.59,
    # Gln
    'CAA': 0.25, 'CAG': 0.75,
    # Asn
    'AAT': 0.46, 'AAC': 0.54,
    # Lys
    'AAA': 0.42, 'AAG': 0.58,
    # Asp
    'GAT': 0.46, 'GAC': 0.54,
    # Glu
    'GAA': 0.42, 'GAG': 0.58,
    # Cys
    'TGT': 0.45, 'TGC': 0.55,
    # Trp
    'TGG': 1.00,
    # Arg
    'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21,
    'AGA': 0.20, 'AGG': 0.20,
    # Gly
    'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25,
}

def calculate_cai(cds_sequence):
    """
    Calculate Codon Adaptation Index (CAI)
    CAI = geometric mean of relative synonimous codon usage
    """
    # Remove stop codon for CAI calculation
    cds = cds_sequence
    if cds.endswith('TAA') or cds.endswith('TAG') or cds.endswith('TGA'):
        cds = cds[:-3]
    
    codons = [cds[i:i+3] for i in range(0, len(cds), 3)]
    
    # Get weights (relative adaptiveness)
    weights = [HUMAN_CODON_USAGE.get(codon, 0.01) for codon in codons]
    
    if not weights:
        return 0.0
    
    # Geometric mean
    import math
    log_sum = sum(math.log(w) for w in weights if w > 0)
    cai = math.exp(log_sum / len(weights))
    
    return cai

def calculate_tai(cds_sequence):
    """
    Calculate tRNA Adaptation Index (tAI)
    Uses tRNA gene copy numbers as weights
    """
    # Simplified tAI weights (based on tRNA gene copy numbers)
    # These are approximate values for human
    TAI_WEIGHTS = {
        'TTT': 0.44, 'TTC': 0.56, 'TTA': 0.07, 'TTG': 0.12,
        'TCT': 0.17, 'TCC': 0.21, 'TCA': 0.15, 'TCG': 0.06,
        'TAT': 0.43, 'TAC': 0.57, 'TAA': 0.30, 'TAG': 0.24,
        'TGT': 0.44, 'TGC': 0.56, 'TGA': 0.46, 'TGG': 1.00,
        'CTT': 0.12, 'CTC': 0.19, 'CTA': 0.07, 'CTG': 0.39,
        'CCT': 0.27, 'CCC': 0.33, 'CCA': 0.28, 'CCG': 0.12,
        'CAT': 0.40, 'CAC': 0.60, 'CAA': 0.24, 'CAG': 0.76,
        'CGT': 0.08, 'CGC': 0.18, 'CGA': 0.11, 'CGG': 0.20,
        'ATT': 0.35, 'ATC': 0.47, 'ATA': 0.16, 'ATG': 1.00,
        'ACT': 0.23, 'ACC': 0.36, 'ACA': 0.28, 'ACG': 0.13,
        'AAT': 0.45, 'AAC': 0.55, 'AAA': 0.41, 'AAG': 0.59,
        'AGT': 0.14, 'AGC': 0.24, 'AGA': 0.20, 'AGG': 0.20,
        'GTT': 0.17, 'GTC': 0.22, 'GTA': 0.11, 'GTG': 0.50,
        'GCT': 0.25, 'GCC': 0.40, 'GCA': 0.23, 'GCG': 0.11,
        'GAT': 0.45, 'GAC': 0.55, 'GAA': 0.41, 'GAG': 0.59,
        'GGT': 0.15, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25,
    }
    
    # Remove stop codon
    cds = cds_sequence
    if cds.endswith('TAA') or cds.endswith('TAG') or cds.endswith('TGA'):
        cds = cds[:-3]
    
    codons = [cds[i:i+3] for i in range(0, len(cds), 3)]
    weights = [TAI_WEIGHTS.get(codon, 0.01) for codon in codons]
    
    if not weights:
        return 0.0
    
    import math
    log_sum = sum(math.log(w) for w in weights if w > 0)
    tai = math.exp(log_sum / len(weights))
    
    return tai

def calculate_gc_content(cds_sequence):
    """Calculate GC content percentage"""
    cds = cds_sequence.replace('A', '').replace('T', '').replace('C', '').replace('G', '')
    if not cds_sequence:
        return 0.0
    return len(cds) / len(cds_sequence) * 100

def calculate_gc3(cds_sequence):
    """Calculate GC content at third codon position"""
    cds = cds_sequence
    if cds.endswith('TAA') or cds.endswith('TAG') or cds.endswith('TGA'):
        cds = cds[:-3]
    
    gc3 = sum(1 for i in range(2, len(cds), 3) if cds[i] in 'GC')
    total = len(cds) // 3
    
    if total == 0:
        return 0.0
    return gc3 / total * 100

def calculate_enc(cds_sequence):
    """
    Calculate Effective Number of Codons (ENC)
    Lower is better (more biased codon usage)
    """
    from collections import defaultdict
    
    cds = cds_sequence
    if cds.endswith('TAA') or cds.endswith('TAG') or cds.endswith('TGA'):
        cds = cds[:-3]
    
    # Group codons by amino acid
    aa_codons = defaultdict(list)
    for i in range(0, len(cds), 3):
        codon = cds[i:i+3]
        aa = codon_to_aa(codon)
        if aa and aa != '*':
            aa_codons[aa].append(codon)
    
    # Calculate ENC (simplified version)
    n = sum(len(codons) for codons in aa_codons.values())
    if n == 0:
        return 61  # Maximum
    
    # Simplified: use variance in codon usage
    enc = 0
    for aa, codons in aa_codons.items():
        if len(codons) > 1:
            # More uniform = higher contribution to ENC
            freq = {}
            for c in codons:
                freq[c] = freq.get(c, 0) + 1
            total = len(codons)
            variance = sum((f/total - 1/len(codons))**2 for f in freq.values()) / len(codons)
            enc += (1 + variance) * len(codons)
        else:
            enc += 1
    
    return min(enc, 61)  # Cap at 61

def codon_to_aa(codon):
    """Convert codon to amino acid"""
    CODON_MAP = {
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
    return CODON_MAP.get(codon, '')

def evaluate_codon(cds_sequence):
    """Evaluate all codon-related metrics"""
    result = {
        'id': cds_sequence.id.replace('>', ''),
        'sequence': str(cds_sequence.sequence),
    }
    
    seq = result['sequence']
    
    result['CAI'] = round(calculate_cai(seq), 4)
    result['tAI'] = round(calculate_tai(seq), 4)
    result['GC_content'] = round(calculate_gc_content(seq), 2)
    result['GC3_content'] = round(calculate_gc3(seq), 2)
    result['ENC'] = round(calculate_enc(seq), 2)
    
    return result

def main():
    # Read FASTA from stdin
    sequences = []
    current_id = None
    current_seq = []
    
    for line in sys.stdin:
        line = line.strip()
        if line.startswith('>'):
            if current_id:
                sequences.append((current_id, ''.join(current_seq)))
            current_id = line
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_id:
        sequences.append((current_id, ''.join(current_seq)))
    
    # Evaluate each sequence
    results = []
    for seq_id, seq in sequences:
        class Seq:
            pass
        s = Seq()
        s.id = seq_id
        s.sequence = seq
        result = evaluate_codon(s)
        results.append(result)
    
    # Output as JSON
    print(json.dumps(results, indent=2))

if __name__ == '__main__':
    main()
