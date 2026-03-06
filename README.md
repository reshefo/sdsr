# SDSR - Spectral Divide-and-Conquer for Species Tree Reconstruction

A scalable phylogenetic reconstruction algorithm that uses spectral clustering and divide-and-conquer strategies to build species trees from gene alignments.

## About

SDSR (Spectral Divide-and-conquer approach for Species tree Reconstruction) is designed for efficient reconstruction of large-scale phylogenetic trees. The method partitions species into smaller subgroups using spectral clustering, reconstructs subtrees independently, and merges them into a complete species tree.

## Installation

### Prerequisites
- Python 3.7 or higher (recommended: 3.11 or 3.12)
- Linux or WSL

> Note: Some dependencies (e.g., `ete3`) rely on the `cgi` module which was removed in Python 3.13+.  
> Therefore Python 3.13+ is currently not supported.

### Setup

1. Clone the repository:
```bash
git clone https://github.com/reshefo/sdsr.git
cd sdsr
```

2. Create venv in your home directory (recommended)
```bash
python3 -m venv ~/sdsr_venv
source ~/sdsr_venv/bin/activate
```

3. Install Python dependencies:
```bash
pip3 install -r requirements.txt
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

- Place your alignment files in the `alignments/` directory.
  - Default format: multiple MSAs in a single PHYLIP file.
  - To use multiple separate PHYLIP files, set the `format_type` parameter to `"phylip_multiple_genes_from_different_files"`.
- Open `sdsr_runner.py` and configure the other parameters as needed.
- Run the reconstruction:
```bash
python sdsr_runner.py
```
- The reconstructed tree will be saved as `sdsr_tree.newick`


### Troubleshooting
If the binaries aren't executable:

```bash
# For RAxML
chmod +x spectraltree/spectraltree_libs/raxml_bins/raxmlHPC-SSE3-linux

# For ASTRAL
chmod +x spectraltree/spectraltree_libs/astral_bin/astral4
```