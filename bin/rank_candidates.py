#!/usr/bin/env python3
"""
Step 4: Rank Candidates
Combine all scores and rank by weighted composite score
"""

import sys
import json
from pathlib import Path

# Default weights (can be configured)
DEFAULT_WEIGHTS = {
    'CAI': 0.15,           # Codon Adaptation Index
    'tAI': 0.10,           # tRNA Adaptation Index
    'GC_content': 0.10,    # GC content (balanced)
    'MFE': 0.15,           # Secondary structure (more negative = better)
    'immunogenicity_score': 0.20,  # Lower is better
    'kozak_score': 0.15,   # Kozak sequence strength
    'GC3_content': 0.10,  # GC at 3rd position
    'ENC': 0.05,           # Effective Number of Codons (lower = better)
}

def normalize_score(value, min_val, max_val, reverse=False):
    """
    Normalize score to 0-100 range
    reverse=True when lower is better
    """
    if max_val == min_val:
        return 50
    
    normalized = (value - min_val) / (max_val - min_val) * 100
    
    if reverse:
        normalized = 100 - normalized
    
    return max(0, min(100, normalized))

def calculate_composite_score(scores, weights):
    """
    Calculate weighted composite score
    """
    total_score = 0
    total_weight = 0
    
    for metric, weight in weights.items():
        if metric in scores and scores[metric] is not None:
            value = scores[metric]
            
            # Handle metrics where lower is better
            reverse = metric in ['immunogenicity_score', 'ENC']
            
            # Get typical ranges for normalization
            ranges = {
                'CAI': (0.5, 1.0),
                'tAI': (0.3, 0.9),
                'GC_content': (30, 70),
                'MFE': (-100, -5),
                'immunogenicity_score': (0, 100),
                'kozak_score': (0, 1),
                'GC3_content': (20, 80),
                'ENC': (20, 61),
            }
            
            min_v, max_v = ranges.get(metric, (0, 100))
            normalized = normalize_score(value, min_v, max_v, reverse)
            
            total_score += normalized * weight
            total_weight += weight
    
    if total_weight == 0:
        return 0
    
    return round(total_score / total_weight, 2)

def rank_candidates(score_files):
    """
    Read all score files, combine, and rank candidates
    """
    # Read from stdin (file paths as JSON)
    input_data = sys.stdin.read().strip()
    
    if not input_data:
        print("Error: No input data", file=sys.stderr)
        sys.exit(1)
    
    # Parse input - assume it's a list of score objects
    try:
        all_scores = json.loads(input_data)
    except json.JSONDecodeError:
        # Try reading as file path
        with open(input_data, 'r') as f:
            all_scores = json.load(f)
    
    # Handle single or multiple score sets
    if not isinstance(all_scores, list):
        all_scores = [all_scores]
    
    # Merge scores by candidate ID
    merged = {}
    for score_set in all_scores:
        if isinstance(score_set, list):
            for item in score_set:
                cand_id = item.get('id', 'unknown')
                if cand_id not in merged:
                    merged[cand_id] = {}
                merged[cand_id].update(item)
        else:
            cand_id = score_set.get('id', 'unknown')
            merged[cand_id] = score_set
    
    # Calculate composite scores
    results = []
    for cand_id, scores in merged.items():
        # Use default or custom weights
        weights = DEFAULT_WEIGHTS.copy()
        
        # Remove non-numeric keys for scoring
        clean_scores = {k: v for k, v in scores.items() 
                        if isinstance(v, (int, float))}
        
        composite = calculate_composite_score(clean_scores, weights)
        
        scores['composite_score'] = composite
        scores['rank'] = 0  # Will be set after sorting
        results.append(scores)
    
    # Sort by composite score (descending)
    results.sort(key=lambda x: x.get('composite_score', 0), reverse=True)
    
    # Assign ranks
    for i, result in enumerate(results):
        result['rank'] = i + 1
    
    return results

def main():
    # Get file paths from command line arguments
    if len(sys.argv) < 2:
        print("Usage: rank_candidates.py <codon_scores.json> [<structure_scores.json> <immuno_scores.json> <utr_scores.json>]", file=sys.stderr)
        sys.exit(1)
    
    # Read all score files
    all_scores = []
    
    for arg in sys.argv[1:]:
        if arg.endswith('.json') and Path(arg).exists():
            with open(arg, 'r') as f:
                data = json.load(f)
                if isinstance(data, list):
                    all_scores.extend(data)
                else:
                    all_scores.append(data)
    
    if not all_scores:
        print("Error: No score data found", file=sys.stderr)
        sys.exit(1)
    
    # Rank
    ranked = rank_candidates(all_scores)
    
    # Output top candidates
    output = {
        'total_candidates': len(ranked),
        'weights_used': DEFAULT_WEIGHTS,
        'top_10': ranked[:10] if len(ranked) > 10 else ranked,
        'all_candidates': ranked
    }
    
    print(json.dumps(output, indent=2))

if __name__ == '__main__':
    main()
