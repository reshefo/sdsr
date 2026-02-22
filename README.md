# SDSR - Spectral Divide-and-Conquer for Species Tree Reconstruction

A scalable phylogenetic reconstruction algorithm that uses spectral clustering and divide-and-conquer strategies to build species trees from gene alignments.

## About

SDSR (Spectral Divide-and-conquer approach for Species tree Reconstruction) is designed for efficient reconstruction of large-scale phylogenetic trees. The method partitions species into smaller subgroups using spectral clustering, reconstructs subtrees independently, and merges them into a complete species tree.

## Installation

### Prerequisites
- Python 3.7 or higher
- Linux or WSL

### Setup

1. Clone the repository:
```bash
git clone https://github.com/reshefo/sdsr.git
cd sdsr
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Quick Test

To verify the installation is working, simply run:
```bash
python sdsr_runner.py
```

This will run SDSR with the default test alignments included in the repository.

The reconstructed tree will be saved as `sdsr_tree.newick`.

### Running with Your Data

- Place your alignment files in the `alignments/` directory
- Alignments should be in a supported format (FASTA, PHYLIP, etc.)
- Open `sdsr_runner.py` and configure the paramerts as needed.
- Run the reconstruction:
```bash
python sdsr_runner.py
```
- The reconstructed tree will be saved as `sdsr_tree.newick`



### Troubleshooting
If the binaries don't exacutable:

```bash
# For RAxML
chmod +x spectraltree/spectraltree_libs/raxml_bins/raxmlHPC-SSE3-linux

# For ASTRAL
chmod +x spectraltree/spectraltree_libs/astral_bin/astral.4.10.12.jar
```