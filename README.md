# 🧬 mRNA Design Pipeline

<p align="center">
  <a href="https://nextflow.io"><img src="https://img.shields.io/badge/Nextflow-25.10+-grey.svg?style=flat&logo=nextflow&logoColor=white"></a>
  <a href="https://www.python.org"><img src="https://img.shields.io/badge/Python-3.10+-blue.svg"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-green.svg"></a>
</p>

A Nextflow-based computational pipeline for designing and optimizing mRNA vaccine candidates. The pipeline evaluates multiple metrics that affect mRNA expression efficiency and generates ranked candidates with comprehensive HTML reports.

## 📖 Overview

This pipeline takes an **amino acid sequence** as input and generates **optimized mRNA candidates** based on multiple efficiency metrics:

```
Input (AA sequence) → Reverse Translation → Multi-dimensional Evaluation → Ranked Output
```

## ✨ Features

### Multiple Optimization Modes
- **Human Optimized** - Uses human codon usage preferences
- **Moderna Style** - Low CpG + Low Uridine (reduced immunogenicity)
- **CAI Optimal** - Maximum codon adaptation index (CAI ≈ 1.0)
- **Balanced** - Trade-off between CAI and immunogenicity
- **Top** - Generate multiple variants and pick the best

### Multi-dimensional Evaluation
- **Codon Optimization**: CAI, tAI, ENC, Codon Pair Bias
- **GC Content**: Overall GC%, GC3%
- **Immunogenicity**: CpG frequency, Uridine content, TLR risk assessment
- **Structure**: MFE estimation, 5' end structure score
- **Translation Efficiency**: Composite scoring
- **Initiation**: Kozak sequence analysis

### Verified with Real Vaccine Proteins
- ✅ SARS-CoV-2 Spike protein (mRNA-1273 / Spikevax)
- ✅ RSV F protein (mRNA-1345)

## 🚀 Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) v21.0+ (optional, can run Python scripts directly)
- Python 3.10+

### Installation

```bash
# Clone the repository
git clone https://github.com/mmyoung/mRNA-design-pipeline.git
cd mRNA-design-pipeline

# Install Python dependencies (if needed)
pip install biopython
```

### Run the Pipeline

```bash
# Activate environment
conda activate mrna-pipeline

# Run with test data (Spike protein)
nextflow run main.nf

# Run with custom input
nextflow run main.nf --input your_sequence.fasta

# Or use Python directly for specific modes:
python bin/parse_input.py test_spike.fasta | python bin/reverse_translate.py /dev/stdin top > candidates.fasta
python bin/evaluate_all.py < candidates.fasta > scores.json
python bin/rank_and_report.py scores.json > report.html
```

## 📋 Optimization Modes

| Mode | Command | Description |
|------|---------|-------------|
| `top` | `--mode top` | Generate multiple variants, pick best (default) |
| `human` | `--mode human` | Pure human codon usage optimization |
| `moderna` | `--mode moderna` | Low CpG + Low Uridine (reduces immune response) |
| `cai` | `--mode cai` | Maximum CAI optimization (CAI ≈ 1.0) |
| `balanced` | `--mode balanced` | Trade-off between CAI and immunogenicity |

### Example Usage

```bash
# Generate CAI-optimal sequence
python bin/parse_input.py your_protein.fasta | python bin/reverse_translate.py /dev/stdin cai > optimized.fasta

# Generate Moderna-style sequence (low immunogenicity)
python bin/parse_input.py your_protein.fasta | python bin/reverse_translate.py /dev/stdin moderna > moderna_style.fasta
```

## 📋 Input Format

The pipeline accepts FASTA format:

```fasta
>spike_protein|SARS-CoV-2
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFS
NVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIV
```

## 📊 Output

### Generated Files
- `results/all_scores.json` - Detailed metrics for all candidates
- `results/ranked_sequences.fasta` - Sequences ranked by score
- `results/report.html` - Interactive HTML report

### Report Contents
- Summary statistics
- Score distribution charts
- Top candidates table with all metrics
- Radar chart comparison

## 📈 Evaluation Metrics

### Codon Optimization (Higher is Better)
| Metric | Description | Target | Achieved |
|--------|-------------|--------|----------|
| CAI | Codon Adaptation Index | > 0.8 | **~1.0** ✓ |
| tAI | tRNA Adaptation Index | > 0.7 | ~0.99 |
| ENC | Effective Number of Codons | 20-35 | 61 |

### Immunogenicity (Lower is Better)
| Metric | Description | Target | Achieved |
|--------|-------------|--------|----------|
| CpG Frequency | CpG dinucleotide % | < 3% | ~6% |
| Uridine Content | T content % | 15-25% | **~12-15%** |

### Structure
| Metric | Description | Target |
|--------|-------------|--------|
| GC Content | Overall GC % | 50-65% |
| MFE | Minimum Free Energy | -15 ~ -30 |

## 🔬 Test Results

### SARS-CoV-2 Spike Protein (1273 aa)

| Metric | Before | After Optimization |
|--------|--------|-------------------|
| CAI | 0.18 | **1.00** |
| Uridine | 26.6% | **14.6%** |
| GC Content | 46% | 64% |

### RSV F Protein (1336 aa)

| Metric | Before | After Optimization |
|--------|--------|-------------------|
| CAI | ~0.18 | **0.99** |
| Uridine | ~26% | **12.1%** |
| GC Content | ~46% | 65% |

## 🏗️ Pipeline Architecture

```
┌─────────────────┐
│  Parse Input   │  → Validate AA sequence
└─────────────────┘
          │
          ▼
┌─────────────────┐
│Reverse Translate│  → Generate optimized CDS (multiple modes)
└─────────────────┘
          │
          ▼
┌─────────────────┐
│  Evaluate All   │  → CAI, tAI, GC%, CpG, Uridine, etc.
└─────────────────┘
          │
          ▼
┌─────────────────┐
│ Rank & Report   │  → Composite score + HTML + FASTA
└─────────────────┘
```

## 📁 Project Structure

```
mRNA-design-pipeline/
├── main.nf                    # Nextflow main workflow
├── nextflow.config             # Configuration
├── environment.yml            # Conda environment
├── bin/
│   ├── parse_input.py         # Input parsing (FASTA/plain text)
│   ├── reverse_translate.py   # Generate candidates with multiple modes
│   ├── evaluate_all.py        # Calculate all metrics
│   ├── rank_and_report.py    # Generate HTML report
│   └── export_fasta.py       # Export ranked sequences
├── test_spike.fasta          # SARS-CoV-2 Spike protein (1273 aa)
├── test_rsv.fasta            # RSV F protein (1336 aa)
└── results/                  # Output directory
    ├── all_scores.json
    ├── ranked_sequences.fasta
    └── report.html
```

## 📚 References

The scoring weights are based on publicly available information from:

- Moderna Therapeutics mRNA platform publications
- Nature Biotechnology (2021) - mRNA design optimization
- Nat Rev Drug Discov (2020) - mRNA vaccine design
- WHO/EMA - Approved mRNA vaccine sequences

## 🔧 Development

### Adding New Metrics

1. Add metric calculation in `bin/evaluate_all.py`
2. Update weights in `bin/rank_and_report.py`
3. Update HTML template for visualization

### Human Codon Usage Table

The pipeline uses a complete human codon usage table with 61 codons:

```python
# Example weights (relative adaptiveness)
HUMAN_CODON_USAGE = {
    'TTC': 1.00,  # Phe (optimal)
    'CTG': 1.00,  # Leu (optimal)
    'GTG': 1.00,  # Val (optimal)
    # ... all 61 codons
}
```

## ⚠️ Limitations

- MFE is estimated from GC content (not full RNAfold calculation)
- Limited to human codon usage currently
- Weights are based on literature; experimental validation recommended

## 🤝 Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## 📜 License

MIT License - see [LICENSE](LICENSE) for details.

## 👤 Author

- Lin Yang - [younglin113@163.com]

---

<p align="center">Built with ❤️ for mRNA vaccine research</p>
