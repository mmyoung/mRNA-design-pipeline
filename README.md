# 🧬 mRNA Design Pipeline

A Nextflow-based pipeline for designing and evaluating mRNA vaccine candidates.

## Overview

This pipeline takes an amino acid sequence as input and generates ranked mRNA candidates based on multiple efficiency metrics:

- **Codon Optimization**: CAI, tAI, ENC, GC content
- **Secondary Structure**: MFE prediction via RNAfold
- **Immunogenicity**: CpG frequency, TLR risk assessment
- **UTR Features**: Kozak sequence, poly(A) signals

## Installation

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (v21+)
- Python 3.10+
- ViennaRNA (for RNAfold)
- Conda (recommended for dependency management)

### Setup

```bash
# Clone or download this pipeline
cd mRNA_design_pipeline

# Make scripts executable
chmod +x bin/*.py

# Test with sample data
echo "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH" > test_input.fasta
```

## Usage

### Basic Usage

```bash
nextflow run main.nf --input test_input.fasta
```

### With Custom Parameters

```bash
nextflow run main.nf \
  --input test_input.fasta \
  --sample_seqs 200 \
  --output_dir results
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | Required | Input amino acid sequence (FASTA) |
| `--output_dir` | results | Output directory |
| `--sample_seqs` | 100 | Number of candidate sequences to generate |
| `--weight_cai` | 0.20 | CAI weight in scoring |
| `--weight_tai` | 0.15 | tAI weight in scoring |
| `--weight_gc` | 0.15 | GC content weight |
| `--weight_mfe` | 0.15 | MFE weight |
| `--weight_cpg` | 0.15 | Immunogenicity weight |
| `--weight_kozak` | 0.10 | Kozak score weight |
| `--weight_polyA` | 0.10 | Poly(A) signal weight |

## Output

The pipeline generates:

- `results/report.html` - Interactive HTML report with visualizations
- `results/ranked_candidates.json` - All candidates with scores
- `results/top_candidates.fasta` - Top candidate sequences

## Pipeline Architecture

```
Input (AA sequence)
    │
    ▼
┌─────────────────┐
│  Parse Input   │  → Validate AA sequence
└─────────────────┘
    │
    ▼
┌─────────────────┐
│Reverse Translate│  → Generate candidate CDS
└─────────────────┘
    │
    ├────────────────┬────────────────┐
    ▼                ▼                ▼
┌──────────┐   ┌────────────┐   ┌──────────────┐
│ Codon   │   │ Structure  │   │Immunogenicity│
│ Scoring │   │ Prediction │   │  Analysis    │
└──────────┘   └────────────┘   └──────────────┘
    │                │                │
    └────────────────┼────────────────┘
                     ▼
            ┌─────────────────┐
            │  Rank Candidates │  → Composite score
            └─────────────────┘
                     │
                     ▼
            ┌─────────────────┐
            │ Generate Report │  → HTML output
            └─────────────────┘
```

## Metrics Explained

### Codon Optimization

- **CAI** (Codon Adaptation Index): 0-1, higher is better
- **tAI** (tRNA Adaptation Index): 0-1, higher is better
- **ENC** (Effective Number of Codons): 20-61, lower is better
- **GC3**: GC content at 3rd codon position

### Structure

- **MFE** (Minimum Free Energy): kcal/mol, more negative = more stable

### Immunogenicity

- **CpG Frequency**: %, lower is better
- **TLR Risk**: 0-1, lower is better
- **Immunogenicity Score**: 0-100, lower is better

## Development

### Adding New Metrics

1. Create new evaluation script in `bin/`
2. Add to `main.nf` workflow
3. Update scoring weights in `bin/rank_candidates.py`
4. Update report template

### Running Tests

```bash
# Test individual modules
python bin/parse_input.py test_input.fasta
python bin/evaluate_codon.py < candidates.fasta
```

## License

MIT License

## Citation

If you use this pipeline in your research, please cite:

> mRNA Design Pipeline - A computational tool for mRNA vaccine design
