#!/usr/bin/env python3
"""
Step 3.4: Evaluate UTR Design
Analyze 5' UTR and 3' UTR features
Kozak sequence, poly(A) signal, etc.
"""

import sys
import json
import re
from pathlib import Path

# Kozak consensus sequence (for mammals)
# GCCRCCATGG where R = purine (A/G)
KOZAK_PATTERN = r'GCC[A/G]CCATGG'

# Poly(A) signal sequences
POLYA_SIGNALS = ['AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA', 'CATAAA']

# Stronger poly(A) signals
STRONG_POLYA = ['AATAAA', 'ATTAAA']

# 5' UTR length recommendations for mRNA vaccines
IDEAL_5UTR_LENGTH = (20, 100)  # nucleotides

# 3' UTR length recommendations
IDEAL_3UTR_LENGTH = (50, 200)

def find_kozak_sequence(sequence):
    """
    Find and score Kozak sequence around start codon
    Returns position and score
    """
    # Look in first 20 nucleotides (should contain ATG)
    search_region = sequence[:25]
    
    matches = list(re.finditer(KOZAK_PATTERN, search_region, re.IGNORECASE))
    
    if matches:
        match = matches[0]
        # Score based on match quality
        score = 1.0
        # Check if ATG is at the expected position
        atg_pos = search_region.upper().find('ATG')
        if atg_pos >= 0 and abs(match.start() + 6 - atg_pos) <= 2:
            score = 1.0
        else:
            score = 0.7
        return {
            'found': True,
            'position': match.start(),
            'sequence': match.group().upper(),
            'score': score
        }
    
    # Check for ATG without strong Kozak
    atg_pos = search_region.upper().find('ATG')
    if atg_pos >= 0:
        return {
            'found': False,
            'position': atg_pos,
            'sequence': search_region[max(0, atg_pos-3):atg_pos+6],
            'score': 0.3
        }
    
    return {'found': False, 'score': 0}

def find_polya_signal(sequence):
    """
    Find poly(A) signal in sequence (near 3' end)
    """
    # Look in last 100 nucleotides
    search_region = sequence[-100:] if len(sequence) > 100 else sequence
    
    for signal in POLYA_SIGNALS:
        pos = search_region.upper().find(signal)
        if pos >= 0:
            return {
                'found': True,
                'signal': signal,
                'position': len(sequence) - len(search_region) + pos,
                'strength': 'strong' if signal in STRONG_POLYA else 'weak'
            }
    
    return {'found': False}

def calculate_gc_in_region(sequence, start, end):
    """Calculate GC content in a specific region"""
    if start >= len(sequence):
        return 0
    region = sequence[start:min(end, len(sequence))]
    if not region:
        return 0
    gc = region.count('G') + region.count('C')
    return gc / len(region) * 100

def score_utr_features(cds_sequence):
    """
    Score UTR-related features
    For this pipeline, we analyze features within/around the CDS
    """
    result = {}
    
    seq = cds_sequence.sequence if hasattr(cds_sequence, 'sequence') else str(cds_sequence)
    
    # 5' region analysis (first 50 nt)
    result['5prime_gc'] = round(calculate_gc_in_region(seq, 0, 50), 2)
    
    # Kozak sequence
    kozak = find_kozak_sequence(seq)
    result['kozak_found'] = kozak['found']
    result['kozak_score'] = kozak['score']
    if kozak.get('sequence'):
        result['kozak_sequence'] = kozak['sequence']
    
    # 3' region analysis (last 50 nt)
    result['3prime_gc'] = round(calculate_gc_in_region(seq, -50, -1), 2)
    
    # Poly(A) signal check (note: in DNA, polyA = T)
    # For mRNA, we look for AATAAA pattern in the CDS region
    # This is more relevant for 3' UTR design
    polya = find_polya_signal(seq)
    result['polya_signal'] = polya.get('signal', 'none')
    
    # Overall 5' end quality (important for translation initiation)
    # Ideal: moderate GC (~40-60%), strong Kozak
    if result['kozak_score'] >= 0.7 and 30 <= result['5prime_gc'] <= 70:
        result['5prime_quality'] = 'good'
    elif result['kozak_score'] >= 0.3:
        result['5prime_quality'] = 'moderate'
    else:
        result['5prime_quality'] = 'poor'
    
    return result

def evaluate_utr(cds_sequence):
    """Evaluate UTR-related metrics"""
    result = {
        'id': cds_sequence.id.replace('>', ''),
        'sequence': str(cds_sequence.sequence),
    }
    
    # Add UTR features
    utr_scores = score_utr_features(cds_sequence)
    result.update(utr_scores)
    
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
        result = evaluate_utr(s)
        results.append(result)
    
    # Output as JSON
    print(json.dumps(results, indent=2))

if __name__ == '__main__':
    main()
