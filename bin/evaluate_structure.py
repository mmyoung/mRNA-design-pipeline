#!/usr/bin/env python3
"""
Step 3.2: Evaluate Secondary Structure
Calculate MFE (Minimum Free Energy) using RNAfold
"""

import sys
import json
import subprocess
import tempfile
import os
from pathlib import Path

def run_rnafold(sequence):
    """
    Run RNAfold to predict secondary structure
    Returns MFE (kcal/mol)
    """
    try:
        # Use RNAfold from ViennaRNA
        result = subprocess.run(
            ['RNAfold', '--noPS'],
            input=sequence,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        output = result.stdout
        lines = output.strip().split('\n')
        
        if len(lines) >= 2:
            # Second line contains structure and energy
            structure_line = lines[1]
            # Extract MFE from parenthesis (e.g., "((-...))  -5.30")
            parts = structure_line.rsplit('(', 1)
            if len(parts) == 2:
                mfe_str = parts[1].strip()
                try:
                    mfe = float(mfe_str)
                    return mfe
                except ValueError:
                    pass
        
        return None
        
    except FileNotFoundError:
        print("Warning: RNAfold not found. Using fallback calculation.", file=sys.stderr)
        return calculate_gc_structure(sequence)
    except Exception as e:
        print(f"Warning: RNAfold failed: {e}", file=sys.stderr)
        return calculate_gc_structure(sequence)

def calculate_gc_structure(sequence):
    """
    Fallback: Calculate structure score based on GC content
    and local stability estimation
    """
    # Simple heuristic: higher GC = more stable
    gc_count = sequence.count('G') + sequence.count('C')
    total = len(sequence)
    
    if total == 0:
        return 0
    
    gc_ratio = gc_count / total
    
    # Estimate MFE (rough approximation)
    # GC base pairs contribute ~3 kcal/mol
    # AT base pairs contribute ~2 kcal/mol
    estimated_mfe = -(gc_ratio * 3 + (1 - gc_ratio) * 2) * (total / 10)
    
    return round(estimated_mfe, 2)

def calculate_local_structure(sequence, window=30):
    """
    Calculate local structure stability in sliding windows
    Focus on 5' UTR and translation initiation region
    """
    # First ~50 nucleotides are critical for translation
    # Look at structure in this region
    
    region = sequence[:min(50, len(sequence))]
    gc_count = region.count('G') + region.count('C')
    
    if len(region) == 0:
        return 0
    
    # Higher GC in 5' region can inhibit ribosome binding
    # But some structure is needed for regulation
    gc_ratio = gc_count / len(region)
    
    # Score: moderate GC (~40-60%) is optimal
    if 0.35 <= gc_ratio <= 0.65:
        return 1.0
    elif 0.25 <= gc_ratio <= 0.75:
        return 0.7
    else:
        return 0.4

def calculate_cAI(sequence):
    """
    Calculate codon-level adaptiveness index
    Similar to CAI but at finer resolution
    """
    # Simplified version
    return 0.5  # Placeholder

def evaluate_structure(cds_sequence):
    """Evaluate secondary structure metrics"""
    result = {
        'id': cds_sequence.id.replace('>', ''),
        'sequence': str(cds_sequence.sequence),
    }
    
    seq = result['sequence']
    
    # Calculate MFE using RNAfold
    result['MFE'] = run_rnafold(seq)
    
    # Local structure at 5' end
    result['local_structure_score'] = calculate_local_structure(seq)
    
    # Overall GC content (related to stability)
    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100 if seq else 0
    result['structure_stability'] = 'high' if gc > 60 else 'moderate' if gc > 40 else 'low'
    
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
        result = evaluate_structure(s)
        results.append(result)
    
    # Output as JSON
    print(json.dumps(results, indent=2))

if __name__ == '__main__':
    main()
