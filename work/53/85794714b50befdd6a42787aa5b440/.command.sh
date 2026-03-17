#!/bin/bash -ue
python3 /Users/lin/.openclaw/workspace/mRNA_design_pipeline/bin/reverse_translate.py < parsed_seq.fasta 100 > candidates.fasta
