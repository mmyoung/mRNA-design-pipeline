#!/usr/bin/env python3
"""
Export ranked sequences to FASTA format
"""

import sys
import json

def main():
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'r') as f:
            data = json.load(f)
    else:
        data = json.load(sys.stdin)
    
    # Sort by rank if available, otherwise by score
    if data and 'rank' in data[0]:
        sorted_data = sorted(data, key=lambda x: x.get('rank', 999))
    elif data and 'Composite_score' in data[0]:
        sorted_data = sorted(data, key=lambda x: x.get('Composite_score', 0), reverse=True)
    else:
        sorted_data = data
    
    for i, entry in enumerate(sorted_data):
        seq_id = entry.get('id', f"sequence_{i+1}")
        seq = entry.get('sequence', '')
        score = entry.get('composite_score', entry.get('score', 'N/A'))
        rank = entry.get('rank', i+1)
        
        print(f">rank{rank}_score{score}_{seq_id}")
        # Print sequence in 80-character lines
        for j in range(0, len(seq), 80):
            print(seq[j:j+80])

if __name__ == '__main__':
    main()
