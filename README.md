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

2. Optional: create and activate virtual environment
```bash
python3 -m venv sdsr_venv
source sdsr_venv/bin/activate
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
chmod +x spectraltree/libs/raxmlHPC_bin/raxmlHPC-SSE3-linux

chmod +x spectraltree/libs/astral_bin/*
```

If ASTRAL fails with a `GLIBC` version error (e.g., `GLIBC_2.XX not found`), your system's glibc version may be too old. This can happen on older Linux distributions (e.g., Ubuntu 18.04 or earlier). Upgrading your OS or glibc should resolve the issue.


## Citation

If you use SDSR in your research, please cite:

```bibtex
@article{reshef2026sdsr,
  title={SDSR: A Spectral Divide-and-Conquer Approach for Species Tree Reconstruction},
  author={Reshef, Ortal and Glassman, Ofer and Zuk, Or and Aizenbud, Yariv and Nadler, Boaz and Jaffe, Ariel},
  year={2026},
  eprint={2603.10215},
  archivePrefix={arXiv},
  primaryClass={q-bio}
}
```

For more details, see the full paper: [https://arxiv.org/abs/2603.10215](https://arxiv.org/abs/2603.10215)