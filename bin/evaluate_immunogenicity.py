#!/usr/bin/env python3
"""
Step 3.3: Evaluate Immunogenicity
Analyze CpG motifs, U-rich sequences, TLR recognition sites
"""

import sys
import json
import re
from pathlib import Path

# Immunogenic motifs to avoid
# Based on literature: CpG dinucleotides trigger TLR9
# U-rich sequences can trigger TLR7/8

IMMUNOGENIC_MOTIFS = {
    # CpG motifs (TLR9)
    'CpG_rare': [
        'CG',  # Any CpG - base
    ],
    # U-rich motifs (TLR7/8)
    'U_rich': [
        'UUUU',
        'UUU',
        'UGUU',
        'GUUU',
    ],
    # Known immunogenic patterns
    'danger': [
        'GGGG',  # G-rich
        'CCCCC',  # C-rich
        'AUAAU',  # AU-rich elements (ARE)
    ],
}

# Ideal ranges for mRNA vaccine design
IDEAL_RANGES = {
    'CpG_frequency': (0, 3),  # 0-3% ideal
    'GC_content': (40, 60),   # 40-60% ideal
    'U_content': (20, 35),    # 20-35% ideal
}

def calculate_cpg_frequency(sequence):
    """Calculate CpG dinucleotide frequency (%)"""
    if len(sequence) < 2:
        return 0
    
    cpg_count = sequence.count('CG')
    return cpg_count / (len(sequence) - 1) * 100

def calculate_u_frequency(sequence):
    """Calculate uracil (T in DNA) content (%)"""
    if not sequence:
        return 0
    return sequence.count('T') / len(sequence) * 100

def count_motif_matches(sequence, motifs):
    """Count occurrences of motifs"""
    total = 0
    for motif in motifs:
        total += sequence.count(motif)
    return total

def assess_tlr9_risk(cpg_freq):
    """
    Assess TLR9 activation risk based on CpG frequency
    Returns score 0-1 (lower is better)
    """
    if cpg_freq < 1:
        return 0.1  # Very low risk
    elif cpg_freq < 2:
        return 0.3  # Low
    elif cpg_freq < 4:
        return 0.5  # Moderate
    elif cpg_freq < 6:
        return 0.7  # High
    else:
        return 0.9  # Very high

def assess_tlr7_risk(u_content):
    """
    Assess TLR7/8 activation risk based on U content
    Returns score 0-1 (lower is better)
    """
    if u_content < 20:
        return 0.2
    elif u_content < 30:
        return 0.4
    elif u_content < 40:
        return 0.6
    else:
        return 0.8

def find_cpg_positions(sequence):
    """Find positions of all CpG dinucleotides"""
    positions = []
    for i in range(len(sequence) - 1):
        if sequence[i:i+2] == 'CG':
            positions.append(i)
    return positions

def find_urich_regions(sequence, min_length=4):
    """Find U-rich regions"""
    regions = []
    for motif in IMMUNOGENIC_MOTIFS['U_rich']:
        for match in re.finditer(motif, sequence):
            regions.append({
                'motif': motif,
                'start': match.start(),
                'end': match.end(),
            })
    return regions

def calculate_immunogenicity_score(sequence):
    """
    Calculate overall immunogenicity score (0-100)
    Lower is better (less immunogenic)
    """
    cpg_freq = calculate_cpg_frequency(sequence)
    u_content = calculate_u_frequency(sequence)
    
    # Risk scores
    tlr9_risk = assess_tlr9_risk(cpg_freq)
    tlr7_risk = assess_tlr7_risk(u_content)
    
    # Weighted score
    score = (tlr9_risk * 0.6 + tlr7_risk * 0.4) * 100
    
    return round(score, 2)

def evaluate_immunogenicity(cds_sequence):
    """Evaluate immunogenicity metrics"""
    result = {
        'id': cds_sequence.id.replace('>', ''),
        'sequence': str(cds_sequence.sequence),
    }
    
    seq = result['sequence']
    
    # Basic metrics
    result['CpG_frequency'] = round(calculate_cpg_frequency(seq), 2)
    result['U_content'] = round(calculate_u_frequency(seq), 2)
    
    # Risk assessments
    result['TLR9_risk'] = round(assess_tlr9_risk(result['CpG_frequency']), 3)
    result['TLR7_risk'] = round(assess_tlr7_risk(result['U_content']), 3)
    
    # Positions
    result['CpG_positions'] = find_cpg_positions(seq)
    result['U_rich_regions'] = find_urich_regions(seq)
    
    # Overall score
    result['immunogenicity_score'] = calculate_immunogenicity_score(seq)
    
    # Risk level
    if result['immunogenicity_score'] < 20:
        result['risk_level'] = 'low'
    elif result['immunogenicity_score'] < 40:
        result['risk_level'] = 'moderate'
    elif result['immunogenicity_score'] < 60:
        result['risk_level'] = 'high'
    else:
        result['risk_level'] = 'very_high'
    
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
        result = evaluate_immunogenicity(s)
        results.append(result)
    
    # Output as JSON
    print(json.dumps(results, indent=2))

if __name__ == '__main__':
    main()
