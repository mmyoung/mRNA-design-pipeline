#!/usr/bin/env python3
"""
Enhanced Combined Evaluation Script
Evaluate all metrics in one pass - including more mRNA design metrics
"""

import sys
import json
import math
import re
from collections import defaultdict

# =============================================================================
# Codon Usage Tables
# =============================================================================

# Human codon usage (relative adaptiveness)
HUMAN_CODON_USAGE = {
    'TTT': 0.45, 'TTC': 0.55, 'TTA': 0.07, 'TTG': 0.13, 'TCT': 0.18, 'TCC': 0.22,
    'TCA': 0.15, 'TCG': 0.06, 'TAT': 0.44, 'TAC': 0.56, 'TAA': 0.30, 'TAG': 0.24,
    'TGT': 0.45, 'TGC': 0.55, 'TGA': 0.46, 'TGG': 1.00, 'CTT': 0.13, 'CTC': 0.20,
    'CTA': 0.07, 'CTG': 0.40, 'CCT': 0.28, 'CCC': 0.33, 'CCA': 0.27, 'CCG': 0.11,
    'CAT': 0.41, 'CAC': 0.59, 'CAA': 0.25, 'CAG': 0.75, 'AAT': 0.46, 'AAC': 0.54,
    'AAA': 0.42, 'AAG': 0.58, 'AGT': 0.15, 'AGC': 0.24, 'AGA': 0.20, 'AGG': 0.20,
    'GTT': 0.18, 'GTC': 0.22, 'GTA': 0.11, 'GTG': 0.49, 'GCT': 0.26, 'GCC': 0.40,
    'GCA': 0.23, 'GCG': 0.11, 'GAT': 0.46, 'GAC': 0.54, 'GAA': 0.42, 'GAG': 0.58,
    'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25,
}

# tRNA Adaptation Index weights (approximate)
TAI_WEIGHTS = {
    'TTT': 0.44, 'TTC': 0.56, 'TTA': 0.07, 'TTG': 0.12, 'TCT': 0.17, 'TCC': 0.21,
    'TCA': 0.15, 'TCG': 0.06, 'TAT': 0.43, 'TAC': 0.57, 'TAA': 0.30, 'TAG': 0.24,
    'TGT': 0.44, 'TGC': 0.56, 'TGA': 0.46, 'TGG': 1.00, 'CTT': 0.12, 'CTC': 0.19,
    'CTA': 0.07, 'CTG': 0.39, 'CCT': 0.27, 'CCC': 0.33, 'CCA': 0.28, 'CCG': 0.12,
    'CAT': 0.40, 'CAC': 0.60, 'CAA': 0.24, 'CAG': 0.76, 'AAT': 0.45, 'AAC': 0.55,
    'AAA': 0.41, 'AAG': 0.59, 'AGT': 0.14, 'AGC': 0.24, 'AGA': 0.20, 'AGG': 0.20,
    'GTT': 0.17, 'GTC': 0.22, 'GTA': 0.11, 'GTG': 0.50, 'GCT': 0.25, 'GCC': 0.40,
    'GCA': 0.23, 'GCG': 0.11, 'GAT': 0.45, 'GAC': 0.55, 'GAA': 0.41, 'GAG': 0.59,
    'GGT': 0.15, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25,
}

# =============================================================================
# Metric Calculation Functions
# =============================================================================

def calculate_cai(seq):
    """Calculate Codon Adaptation Index (CAI)"""
    if seq.endswith('TAA') or seq.endswith('TAG') or seq.endswith('TGA'):
        seq = seq[:-3]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    weights = [HUMAN_CODON_USAGE.get(c, 0.01) for c in codons]
    if not weights or all(w == 0.01 for w in weights):
        return 0.0
    valid_weights = [w for w in weights if w > 0]
    if not valid_weights:
        return 0.0
    return math.exp(sum(math.log(w) for w in valid_weights) / len(valid_weights))

def calculate_tai(seq):
    """Calculate tRNA Adaptation Index (tAI)"""
    if seq.endswith('TAA') or seq.endswith('TAG') or seq.endswith('TGA'):
        seq = seq[:-3]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    weights = [TAI_WEIGHTS.get(c, 0.01) for c in codons]
    valid_weights = [w for w in weights if w > 0]
    if not valid_weights:
        return 0.0
    return math.exp(sum(math.log(w) for w in valid_weights) / len(valid_weights))

def calculate_gc(seq):
    """Calculate overall GC content"""
    if not seq:
        return 0.0
    return (seq.count('G') + seq.count('C')) / len(seq) * 100

def calculate_gc3(seq):
    """Calculate GC content at third codon position"""
    if seq.endswith('TAA') or seq.endswith('TAG') or seq.endswith('TGA'):
        seq = seq[:-3]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if not codons:
        return 0.0
    gc3 = sum(1 for c in codons if len(c) == 3 and c[2] in 'GC')
    return gc3 / len(codons) * 100 if codons else 0

def calculate_cpg(seq):
    """Calculate CpG dinucleotide frequency"""
    if len(seq) < 2:
        return 0.0
    return seq.count('CG') / (len(seq) - 1) * 100

def calculate_au_content(seq):
    """Calculate AU (A+T) content - important for mRNA stability"""
    if not seq:
        return 0.0
    return (seq.count('A') + seq.count('T')) / len(seq) * 100

def calculate_uridine_content(seq):
    """Calculate Uridine content (T in DNA)"""
    if not seq:
        return 0.0
    return seq.count('T') / len(seq) * 100

def calculate_codon_pair_score(seq):
    """
    Calculate Codon Pair Bias Score
    Positive = optimal pairs, Negative = suboptimal pairs
    Simplified version using general principles
    """
    if seq.endswith('TAA') or seq.endswith('TAG') or seq.endswith('TGA'):
        seq = seq[:-3]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if len(codons) < 2:
        return 0.0
    
    # Simplified: prefer pairs with different codon families (reduce ribosome pausing)
    # This is a heuristic - real calculation needs reference table
    score = 0
    for i in range(len(codons) - 1):
        c1, c2 = codons[i], codons[i+1]
        # Optimal: different codon families (avoid same amino acid repeats)
        if c1[:2] != c2[:2]:
            score += 0.1
        # Suboptimal: same codon (potential pausing)
        else:
            score -= 0.05
    
    return score / (len(codons) - 1) * 100 if len(codons) > 1 else 0

def calculate_enc(seq):
    """
    Effective Number of Codons (ENC)
    Lower = more biased codon usage (better for expression)
    """
    if seq.endswith('TAA') or seq.endswith('TAG') or seq.endswith('TGA'):
        seq = seq[:-3]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    
    # Group by amino acid
    aa_codons = defaultdict(list)
    for c in codons:
        if len(c) == 3:
            # Get amino acid (simplified)
            aa = codon_to_aa(c)
            if aa:
                aa_codons[aa].append(c)
    
    if not aa_codons:
        return 61.0
    
    # Calculate ENC (simplified)
    n = sum(len(codons) for codons in aa_codons.values())
    if n == 0:
        return 61.0
    
    enc = 0
    for aa, cds in aa_codons.items():
        if len(cds) > 1:
            freq = {}
            for c in cds:
                freq[c] = freq.get(c, 0) + 1
            total = len(cds)
            # Calculate variance
            variance = sum((f/total - 1/len(cds))**2 for f in freq.values()) / len(cds)
            enc += (1 + variance) * len(cds)
        else:
            enc += 1
    
    return min(enc, 61.0)

def codon_to_aa(codon):
    """Convert codon to amino acid"""
    table = {
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
    return table.get(codon, '')

def calculate_immunogenicity(seq):
    """
    Calculate immunogenicity score based on TLR recognition motifs
    Lower is better for vaccine applications
    """
    cpg = calculate_cpg(seq)
    u_content = calculate_uridine_content(seq)
    
    # TLR9 risk (CpG)
    if cpg < 1:
        tlr9 = 0.1
    elif cpg < 2:
        tlr9 = 0.3
    elif cpg < 4:
        tlr9 = 0.5
    elif cpg < 6:
        tlr9 = 0.7
    else:
        tlr9 = 0.9
    
    # TLR7/8 risk (U-rich)
    if u_content < 20:
        tlr7 = 0.2
    elif u_content < 30:
        tlr7 = 0.4
    elif u_content < 40:
        tlr7 = 0.6
    else:
        tlr7 = 0.8
    
    # Combined score
    return (tlr9 * 0.6 + tlr7 * 0.4) * 100

def find_kozak(seq):
    """Find and score Kozak sequence around start codon"""
    pattern = re.compile(r'GCC[A/G]CCATGG', re.IGNORECASE)
    match = pattern.search(seq[:25])
    if match:
        return 1.0
    # Check for ATG without perfect Kozak
    if seq[:3].upper() == 'ATG':
        return 0.3
    return 0.0

def calculate_mfe_estimate(seq):
    """
    Estimate Minimum Free Energy using a simplified model
    Real implementation would use RNAfold
    """
    if not seq:
        return 0.0
    
    # Simplified secondary structure prediction
    # Higher GC = more stable (more negative MFE)
    gc = calculate_gc(seq)
    
    # Base MFE estimate (kcal/mol per nt)
    # GC pair: ~3 kcal/mol, AU pair: ~2 kcal/mol
    gc_ratio = gc / 100
    au_ratio = (100 - gc) / 100
    
    # Average ~10 nt per stem-loop
    num_bp = len(seq) / 10
    mfe = -(gc_ratio * 3 + au_ratio * 2) * num_bp
    
    return round(mfe, 2)

def calculate_5prime_structure_score(seq):
    """
    Calculate 5'UTR/CDS initiation region structure score
    Too stable = bad (ribosome can't unwind)
    Too unstable = also bad (no regulation)
    Optimal: moderate structure
    """
    if len(seq) < 30:
        return 0.5
    
    region = seq[:30]
    gc = calculate_gc(region)
    
    # Optimal GC for ribosome initiation: 30-60%
    if 35 <= gc <= 60:
        return 1.0
    elif 25 <= gc <= 70:
        return 0.7
    else:
        return 0.4

def calculate_coding_stability(seq):
    """
    Calculate coding region stability
    Avoid extremely stable or unstable regions
    """
    gc = calculate_gc(seq)
    
    # Optimal GC for CDS: 50-60%
    if 50 <= gc <= 60:
        return 1.0
    elif 40 <= gc <= 65:
        return 0.7
    else:
        return 0.4

def calculate_translation_efficiency(seq):
    """
    Calculate translation efficiency score
    Combines CAI, tAI, ENC, and structure factors
    """
    cai = calculate_cai(seq)
    tai = calculate_tai(seq)
    enc = calculate_enc(seq)
    gc = calculate_gc(seq)
    
    # ENC normalization (lower is better)
    enc_score = 1 - (enc - 20) / (61 - 20)  # 20-61 range
    
    # GC balance
    gc_score = 1 - abs(gc - 55) / 35  # Optimal ~55%
    
    # Combined
    score = (cai * 0.3 + tai * 0.3 + enc_score * 0.2 + gc_score * 0.2)
    
    return round(score, 4)

def calculate_rscu(seq):
    """
    Calculate Relative Synonymous Codon Usage
    Returns dictionary of RSCU values per codon
    """
    if seq.endswith('TAA') or seq.endswith('TAG') or seq.endswith('TGA'):
        seq = seq[:-3]
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    
    # Group by amino acid
    aa_codons = defaultdict(list)
    for c in codons:
        if len(c) == 3:
            aa = codon_to_aa(c)
            if aa and aa != '*':
                aa_codons[aa].append(c)
    
    # Calculate RSCU
    rscu = {}
    for aa, cds in aa_codons.items():
        if len(cds) > 0:
            # Count codons
            counts = defaultdict(int)
            for c in cds:
                counts[c] += 1
            # RSCU = count / (count / num_synonymous_codons)
            num_syn = len(set(AA_TO_CODONS.get(aa, [aa])))
            for c, cnt in counts.items():
                rscu[c] = round(cnt / (len(cds) / num_syn), 3) if num_syn > 0 else 1.0
    
    return rscu

# Amino acid to codon mapping
AA_TO_CODONS = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'Y': ['TAT', 'TAC'], '*': ['TAA', 'TAG', 'TGA'], 'H': ['CAT', 'CAC'],
    'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'],
    'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'],
    'W': ['TGG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
}

# =============================================================================
# Main Evaluation Function
# =============================================================================

def evaluate_sequence(seq_id, seq):
    """Evaluate a single sequence with all metrics"""
    result = {
        'id': seq_id,
        'sequence': seq,
        # Core codon metrics
        'CAI': round(calculate_cai(seq), 4),
        'tAI': round(calculate_tai(seq), 4),
        'ENC': round(calculate_enc(seq), 2),
        # GC content
        'GC_content': round(calculate_gc(seq), 2),
        'GC3_content': round(calculate_gc3(seq), 2),
        # Immunogenicity
        'CpG_frequency': round(calculate_cpg(seq), 2),
        'AU_content': round(calculate_au_content(seq), 2),
        'Uridine_content': round(calculate_uridine_content(seq), 2),
        'Immunogenicity_score': round(calculate_immunogenicity(seq), 2),
        # Initiation
        'Kozak_score': find_kozak(seq),
        # Structure
        'MFE_estimate': calculate_mfe_estimate(seq),
        '5prime_structure_score': round(calculate_5prime_structure_score(seq), 3),
        'Coding_stability': round(calculate_coding_stability(seq), 3),
        # Efficiency
        'Translation_efficiency': calculate_translation_efficiency(seq),
        'Codon_pair_bias': round(calculate_codon_pair_score(seq), 2),
    }
    
    return result

def main():
    """Read FASTA and evaluate all sequences"""
    sequences = []
    current_id = None
    current_seq = []
    
    for line in sys.stdin:
        line = line.strip()
        if line.startswith('>'):
            if current_id:
                sequences.append((current_id, ''.join(current_seq)))
            current_id = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    
    if current_id:
        sequences.append((current_id, ''.join(current_seq)))
    
    if not sequences:
        print("Error: No sequences found", file=sys.stderr)
        sys.exit(1)
    
    # Evaluate all sequences
    results = [evaluate_sequence(seq_id, seq) for seq_id, seq in sequences]
    
    print(json.dumps(results, indent=2))

if __name__ == '__main__':
    main()
