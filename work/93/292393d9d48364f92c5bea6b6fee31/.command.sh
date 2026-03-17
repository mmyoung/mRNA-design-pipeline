#!/bin/bash -ue
python3 /Users/lin/.openclaw/workspace/mRNA_design_pipeline/bin/rank_candidates.py codon_scores.json structure_scores.json immuno_scores.json utr_scores.json > ranked.json
