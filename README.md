# 🧬 mRNA Design Pipeline

<p align="center">
  <a href="https://nextflow.io"><img src="https://img.shields.io/badge/Nextflow-25.10+-grey.svg?style=flat&logo=nextflow&logoColor=white"></a>
  <a href="https://www.python.org"><img src="https://img.shields.io/badge/Python-3.10+-blue.svg"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-green.svg"></a>
</p>

A Nextflow-based computational pipeline for designing and optimizing mRNA vaccine candidates. The pipeline evaluates multiple metrics that affect mRNA expression efficiency and generates ranked candidates with comprehensive HTML reports.

## 📖 Overview

This pipeline takes an **amino acid sequence** as input and generates **ranked mRNA candidates** based on multiple efficiency metrics:

```
Input (AA sequence) → Reverse Translation → Multi-dimensional Evaluation → Ranked Output
```

## ✨ Features

### Multi-dimensional Evaluation
- **Codon Optimization**: CAI, tAI, ENC, Codon Pair Bias
- **GC Content**: Overall GC%, GC3%
- **Immunogenicity**: CpG frequency, Uridine content, TLR risk assessment
- **Structure**: MFE estimation, 5' end structure score
- **Translation Efficiency**: Composite scoring
- **Initiation**: Kozak sequence analysis

### Output
- Interactive HTML report with visualizations
- Detailed metrics for each candidate
- Top candidates ranked by composite score

## 🚀 Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) v21.0+
- [Conda](https://docs.conda.io/) (recommended) or Python 3.10+
- Java 11+

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/mRNA-design-pipeline.git
cd mRNA-design-pipeline

# Create conda environment
conda env create -f environment.yml
conda activate mrna-pipeline
```

### Run the Pipeline

```bash
# Activate environment
conda activate mrna-pipeline

# Run with test data
nextflow run main.nf

# Run with custom input
nextflow run main.nf --input your_sequence.fasta

# Generate more candidates
nextflow run main.nf --input your_sequence.fasta --sample_seqs 100
```

## 📋 Input Format

The pipeline accepts FASTA format:

```fasta
>my_sequence
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH
```

Or plain amino acid sequence:

```bash
echo "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH" > input.txt
nextflow run main.nf --input input.txt
```

## 📊 Output

### Generated Files
- `results/report.html` - Interactive HTML report
- `work/` - Intermediate files (Nextflow cache)

### Report Contents
- Summary statistics
- Score distribution charts
- Top candidates table with all metrics
- Radar chart comparison

## ⚙️ Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | `test_input.fasta` | Input amino acid sequence (FASTA) |
| `--sample_seqs` | 50 | Number of candidate sequences to generate |

## 📈 Evaluation Metrics

### Codon Optimization
| Metric | Description | Target |
|--------|-------------|--------|
| CAI | Codon Adaptation Index | > 0.8 |
| tAI | tRNA Adaptation Index | > 0.7 |
| ENC | Effective Number of Codons | 20-35 |

### Immunogenicity (Lower is Better)
| Metric | Description | Target |
|--------|-------------|--------|
| CpG Frequency | CpG dinucleotide % | < 3% |
| Uridine Content | T content % | 20-30% |

### Structure
| Metric | Description | Target |
|--------|-------------|--------|
| GC Content | Overall GC % | 50-60% |
| MFE | Minimum Free Energy | -15 ~ -30 |
| 5' Structure | Initiation region accessibility | High |

## 🏗️ Pipeline Architecture

```
┌─────────────────┐
│  Parse Input   │  → Validate AA sequence
└─────────────────┘
    │
    ▼
┌─────────────────┐
│Reverse Translate│  → Generate candidate CDS
└─────────────────┘
    │
    ▼
┌─────────────────┐
│  Evaluate All   │  → CAI, tAI, GC%, CpG, etc.
└─────────────────┘
    │
    ▼
┌─────────────────┐
│ Rank & Report   │  → Composite score + HTML
└─────────────────┘
```

## 🔧 Development

### Project Structure
```
mRNA-design-pipeline/
├── main.nf                 # Nextflow main workflow
├── nextflow.config         # Configuration
├── environment.yml        # Conda environment
├── bin/
│   ├── parse_input.py         # Input parsing
│   ├── reverse_translate.py   # Generate candidates
│   ├── evaluate_all.py        # Calculate all metrics
│   └── rank_and_report.py     # Generate HTML report
├── templates/             # HTML templates
├── data/                  # Reference data (codon tables)
└── test_input.fasta       # Test data
```

### Adding New Metrics

1. Add metric calculation in `bin/evaluate_all.py`
2. Update weights in `bin/rank_and_report.py`
3. Update HTML template for visualization

## 📚 References

The scoring weights are based on publicly available information from:

- Moderna Therapeutics mRNA platform publications
- Nature Biotechnology (2021) - mRNA design optimization
- Nat Rev Drug Discov (2020) - mRNA vaccine design

## ⚠️ Limitations

- MFE is estimated from GC content (not full RNAfold calculation)
- Limited to human codon usage currently
- Weights are based on literature; experimental validation recommended

## 🤝 Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## � License

MIT License - see [LICENSE](LICENSE) for details.

## 👤 Author

- Your Name - [your.email@example.com]

---

<p align="center">Built with ❤️ for mRNA vaccine research</p>
