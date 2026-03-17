#!/bin/bash -ue
python3 /Users/lin/.openclaw/workspace/mRNA_design_pipeline/bin/evaluate_immunogenicity.py < candidates.fasta > immuno_scores.json
