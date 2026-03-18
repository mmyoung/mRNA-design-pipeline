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

# =============================================================================
# Human Codon Usage Tables (for optimization)
# =============================================================================

# Human codon usage frequency (COMPLETE - all 61 codons)
HUMAN_CODON_USAGE = {
    # Phe
    'TTT': 0.45, 'TTC': 0.55,
    # Leu
    'TTA': 0.07, 'TTG': 0.13, 'CTT': 0.13, 'CTC': 0.20, 'CTA': 0.07, 'CTG': 0.40,
    # Ile
    'ATT': 0.36, 'ATC': 0.48, 'ATA': 0.16,
    # Met (Start)
    'ATG': 1.00,
    # Val
    'GTT': 0.18, 'GTC': 0.22, 'GTA': 0.11, 'GTG': 0.49,
    # Ser
    'TCT': 0.18, 'TCC': 0.22, 'TCA': 0.15, 'TCG': 0.06, 'AGT': 0.15, 'AGC': 0.24,
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
    'CGT': 0.08, 'CGC': 0.19, 'CGA': 0.11, 'CGG': 0.21, 'AGA': 0.20, 'AGG': 0.20,
    # Gly
    'GGT': 0.16, 'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25,
}

# Calculate CAI weights (relative adaptiveness = freq / max_freq for each AA)
CAI_WEIGHTS = {}
for aa, codons in AA_TO_CODONS.items():
    if aa == '*':
        continue
    max_freq = max(HUMAN_CODON_USAGE.get(c, 0.01) for c in codons)
    for c in codons:
        CAI_WEIGHTS[c] = HUMAN_CODON_USAGE.get(c, 0.01) / max_freq if max_freq > 0 else 0.01

# Codons with CpG (to avoid) - these trigger TLR9
CPG_CODONS = ['CGA', 'CGG', 'TGC', 'CGT', 'CGC']

# Codons with high Uridine content (to reduce for Moderna-style)
URIDINE_CODONS = {
    'TTT': 2, 'TTC': 2, 'TTA': 1, 'TTG': 1,  # Leu
    'TCT': 1, 'TCC': 1, 'TCA': 1, 'TCG': 0,  # Ser
    'TAT': 1, 'TAC': 1,  # Tyr
    'TGT': 1, 'TGC': 1,  # Cys
    'CTT': 1, 'CTC': 1, 'CTA': 1, 'CTG': 2,  # Leu
    'CCT': 1, 'CCC': 1, 'CCA': 1, 'CCG': 0,  # Pro
    'CAT': 1, 'CAC': 1,  # His
    'CAA': 1, 'CAG': 0,  # Gln
    'ATT': 2, 'ATC': 1, 'ATA': 2,  # Ile
    'ACT': 1, 'ACC': 1, 'ACA': 1, 'ACG': 0,  # Thr
    'AAT': 1, 'AAC': 1,  # Asn
    'AAA': 2, 'AAG': 1,  # Lys
    'AGT': 1, 'AGC': 1,  # Ser
    'AGA': 1, 'AGG': 1,  # Arg
    'GTT': 1, 'GTC': 1, 'GTA': 1, 'GTG': 2,  # Val
    'GCT': 1, 'GCC': 1, 'GCA': 1, 'GCG': 0,  # Ala
    'GAT': 1, 'GAC': 1,  # Asp
    'GAA': 1, 'GAG': 0,  # Glu
    'GGT': 1, 'GGC': 1, 'GGA': 1, 'GGG': 2,  # Gly
}

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

def reverse_translate_human_optimized(aa_seq):
    """
    Human codon usage optimized reverse translation
    Uses human codon usage frequencies as weights
    """
    codons = []
    for aa in aa_seq:
        if aa == '*':
            continue
        codon_list = AA_TO_CODONS.get(aa, [])
        if codon_list:
            # Prefer codons with higher human usage frequency
            best_codon = max(codon_list, key=lambda c: HUMAN_CODON_USAGE.get(c, 0.01))
            codons.append(best_codon)
        else:
            codons.append('NNN')
    
    codons.append('TAA')  # Stop codon
    return ''.join(codons)

def reverse_translate_moderna_style(aa_seq):
    """
    Moderna-style optimization:
    1. Prefer human optimal codons
    2. Avoid CpG motifs (reduce immunogenicity)
    3. Reduce Uridine content (reduce TLR7/8 activation)
    """
    codons = []
    for aa in aa_seq:
        if aa == '*':
            continue
        codon_list = AA_TO_CODONS.get(aa, [])
        if codon_list:
            # Score each codon
            best_score = -999
            best_codon = codon_list[0]
            
            for codon in codon_list:
                score = 0
                
                # Base score: human codon usage
                score += HUMAN_CODON_USAGE.get(codon, 0.01) * 10
                
                # Penalty for CpG motifs (very high - TLR9)
                if codon in CPG_CODONS:
                    score -= 5
                
                # Penalty for high Uridine content
                uridine_count = URIDINE_CODONS.get(codon, 1)
                score -= uridine_count * 0.5
                
                if score > best_score:
                    best_score = score
                    best_codon = codon
            
            codons.append(best_codon)
        else:
            codons.append('NNN')
    
    codons.append('TAA')  # Stop codon
    return ''.join(codons)

def reverse_translate_variants(aa_seq, num_variants=10):
    """
    Generate multiple variants with different optimization strategies
    """
    variants = []
    
    # 1. Human optimized
    variants.append(('human_optimized', reverse_translate_human_optimized(aa_seq)))
    
    # 2. Moderna-style
    variants.append(('moderna_style', reverse_translate_moderna_style(aa_seq)))
    
    # 3. Random variants with different seeds
    for i in range(num_variants - 2):
        random.seed(i * 1000 + 42)
        codons = []
        for aa in aa_seq:
            if aa == '*':
                continue
            codon_list = AA_TO_CODONS.get(aa, [])
            if codon_list:
                codons.append(random.choice(codon_list))
            else:
                codons.append('NNN')
        codons.append('TAA')
        variants.append((f'random_{i+1}', ''.join(codons)))
    
    return variants

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
    # Read amino acid sequence from stdin or file
    if len(sys.argv) > 1 and sys.argv[1] != '/dev/stdin':
        # File path provided
        with open(sys.argv[1], 'r') as f:
            content = f.read()
    else:
        # Read from stdin
        content = sys.stdin.read()
    
    # Remove header if present
    if content.strip().startswith('>'):
        content = content.split('\n', 1)[1]
    aa_seq = content.replace('\n', '').replace(' ', '').strip()
    
    if not aa_seq:
        print("Error: Empty sequence", file=sys.stderr)
        sys.exit(1)
    
    # Get parameters
    num_samples = 50  # Default
    mode = 'mixed'  # Default: mixed optimization
    
    if len(sys.argv) > 2:
        try:
            num_samples = int(sys.argv[2])
        except:
            mode = sys.argv[2]
    
    # Generate candidates based on mode
    candidates = []
    
    if mode == 'human':
        # Pure human optimized
        cds = reverse_translate_human_optimized(aa_seq)
        candidates.append(('human_optimized', cds))
        
    elif mode == 'moderna':
        # Moderna-style (low CpG, low U)
        cds = reverse_translate_moderna_style(aa_seq)
        candidates.append(('moderna_style', cds))
        
    elif mode == 'top':
        # Top optimization: best CAI with iterative refinement
        variants = generate_top_candidates(aa_seq, num_samples)
        candidates.extend(variants)
        
    elif mode == 'cai':
        # Full CAI optimization
        cds = optimize_full_caI(aa_seq)
        candidates.append(('cai_full', cds))
        
    elif mode == 'balanced':
        # Balanced optimization
        cds = optimize_balanced(aa_seq)
        candidates.append(('balanced', cds))
        
    else:
        # Mixed: human optimized + moderna style + random variants
        variants = reverse_translate_variants(aa_seq, num_samples)
        candidates.extend(variants)
    
    # Remove duplicates
    unique_candidates = []
    seen = set()
    for name, seq in candidates:
        if seq not in seen:
            seen.add(seq)
            unique_candidates.append((name, seq))
    
    # Output as FASTA
    for i, (name, cds) in enumerate(unique_candidates):
        print(f">{name}")
        print(cds)

if __name__ == '__main__':
    main()

# =============================================================================
# Advanced Optimization Functions
# =============================================================================

def calculate_rscu(codon_usage):
    """Calculate Relative Synonymous Codon Usage"""
    rscu = {}
    for aa, codons in AA_TO_CODONS.items():
        if aa == '*' or aa == 'M':  # Skip stop and start
            continue
        total = sum(codon_usage.get(c, 0.01) for c in codons)
        if total > 0:
            for c in codons:
                rscu[c] = codon_usage.get(c, 0.01) * len(codons) / total
    return rscu

def calculate_cai_score(seq, codon_usage):
    """Calculate CAI for a sequence"""
    rscu = calculate_rscu(codon_usage)
    rscu_values = [rscu.get(seq[i:i+3], 0.01) for i in range(0, len(seq)-3, 3)]
    if not rscu_values:
        return 0
    # Geometric mean
    import math
    log_sum = sum(math.log(v) for v in rscu_values if v > 0)
    return math.exp(log_sum / len(rscu_values))

def optimize_caI_enhanced(aa_seq, codon_usage, iterations=5):
    """
    Enhanced CAI optimization with iterative refinement
    """
    # First pass: greedy selection based on highest codon usage
    codons = []
    for aa in aa_seq:
        if aa == '*':
            continue
        codon_list = AA_TO_CODONS.get(aa, [])
        if codon_list:
            best = max(codon_list, key=lambda c: codon_usage.get(c, 0.01))
            codons.append(best)
        else:
            codons.append('NNN')
    
    # Iterative refinement
    best_seq = ''.join(codons) + 'TAA'
    best_score = calculate_cai_score(best_seq, codon_usage)
    
    for iteration in range(iterations):
        improved = False
        for pos in range(0, len(codons), 3):  # Check every codon
            if pos >= len(codons):
                break
            aa = aa_seq[pos // 3]
            if aa == '*':
                continue
            codon_list = AA_TO_CODONS.get(aa, [])
            if not codon_list:
                continue
            
            original_codon = codons[pos // 3]
            
            for new_codon in codon_list:
                if new_codon == original_codon:
                    continue
                
                # Try replacement
                test_codons = codons[:]
                test_codons[pos // 3] = new_codon
                test_seq = ''.join(test_codons) + 'TAA'
                
                test_score = calculate_cai_score(test_seq, codon_usage)
                
                if test_score > best_score:
                    best_score = test_score
                    best_seq = test_seq
                    codons = test_codons
                    improved = True
        
        if not improved:
            break
    
    return best_seq

def optimize_full_caI(aa_seq):
    """
    Full CAI optimization - selects the best possible codon at each position
    Uses the highest frequency human codons
    """
    codons = []
    for aa in aa_seq:
        if aa == '*':
            continue
        codon_list = AA_TO_CODONS.get(aa, [])
        if codon_list:
            # Sort by human codon usage and take the best
            best = max(codon_list, key=lambda c: HUMAN_CODON_USAGE.get(c, 0))
            codons.append(best)
        else:
            codons.append('NNN')
    
    return ''.join(codons) + 'TAA'

def optimize_balanced(aa_seq):
    """
    Balanced optimization: high CAI + moderate GC + low immunogenicity
    """
    codons = []
    for aa in aa_seq:
        if aa == '*':
            continue
        codon_list = AA_TO_CODONS.get(aa, [])
        if codon_list:
            best_score = -999
            best_codon = codon_list[0]
            
            for codon in codon_list:
                score = 0
                
                # Human usage (40%)
                score += HUMAN_CODON_USAGE.get(codon, 0) * 40
                
                # Avoid CpG (25%)
                if 'CG' in codon:
                    score -= 25
                
                # Moderate GC preference (15%)
                gc_count = codon.count('G') + codon.count('C')
                score += gc_count * 5
                
                # Avoid high-U codons (20%)
                u_count = codon.count('T')
                score -= u_count * 6
                
                if score > best_score:
                    best_score = score
                    best_codon = codon
            
            codons.append(best_codon)
        else:
            codons.append('NNN')
    
    return ''.join(codons) + 'TAA'

def generate_top_candidates(aa_seq, num_candidates=20):
    """
    Generate top candidates with different optimization strategies
    """
    candidates = []
    
    # 1. Full CAI optimization
    cds = optimize_full_caI(aa_seq)
    candidates.append(('cai_optimal', cds))
    
    # 2. Balanced optimization
    cds = optimize_balanced(aa_seq)
    candidates.append(('balanced', cds))
    
    # 3-6. Iterative CAI with different iterations
    for iters in [3, 5, 10]:
        cds = optimize_caI_enhanced(aa_seq, HUMAN_CODON_USAGE, iterations=iters)
        candidates.append((f'cai_iter_{iters}', cds))
    
    # 7-10. Random variants with human preference
    for seed in range(4, 8):
        random.seed(seed * 100 + 42)
        codons = []
        for aa in aa_seq:
            if aa == '*':
                continue
            codon_list = AA_TO_CODONS.get(aa, [])
            if codon_list:
                # Weighted random: prefer high-usage codons
                weights = [HUMAN_CODON_USAGE.get(c, 0.01) ** 2 for c in codon_list]
                total = sum(weights)
                r = random.random() * total
                cumulative = 0
                selected = codon_list[0]
                for i, w in enumerate(weights):
                    cumulative += w
                    if r <= cumulative:
                        selected = codon_list[i]
                        break
                codons.append(selected)
            else:
                codons.append('NNN')
        candidates.append((f'human_weighted_{seed}', ''.join(codons) + 'TAA'))
    
    return candidates
