# Foldseek Pocket Miner

<p align="center">
  <img src="puppy-dog-eyes-sad.gif" width="150">
  <br>
  <em>Please ⭐ star this repo if it helps you!</em>
</p>

A Python tool that leverages Foldseek to discover similar protein structures, extract binding pockets with their ligands, and generate comprehensive PyMOL sessions for structural analysis.

## Features

- **Structure Search**: Use Foldseek to search PDB100 database for structurally similar proteins
- **Automatic Download**: Fetch protein structures from RCSB PDB
- **Structure Alignment**: Align similar structures to the query protein
- **Ligand Extraction**: Extract ligands from aligned structures to identify binding pockets
- **Conservation Analysis**: Calculate and visualize sequence conservation across hits
- **Pharmacophore Generation**: Generate ligand-based pharmacophore models
- **PyMOL Session**: Create a comprehensive `.pse` file with:
  1. Input structure
  2. Sequence conservation coloring
  3. Aligned ligands forming a pharmacophore

## Installation

### Prerequisites

1. **Foldseek**: Install the Foldseek command-line tool
```bash
# macOS (Apple Silicon)
conda install -c conda-forge -c bioconda foldseek

# or download binary
wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz
tar xzf foldseek-osx-universal.tar.gz
export PATH="$(pwd)/foldseek/bin:$PATH"
```

2. **PyMOL**: Install PyMOL (open-source or commercial)
```bash
# macOS with Homebrew
brew install pymol

# or conda
conda install -c conda-forge pymol-open-source
```

### Install the package

```bash
# Clone the repository
git clone https://github.com/JianJinglin/foldseek-pocket-miner.git
cd foldseek-pocket-miner

# Install in development mode
pip install -e .

# Or install with dev dependencies
pip install -e ".[dev]"
```

## Usage

### Command Line Interface

```bash
# Search by PDB ID
foldseek-pocket-miner --pdb 1ABC --chain A --output results/

# Search by sequence (will use AlphaFold or ESMFold for structure prediction)
foldseek-pocket-miner --sequence "MVLSPADKTN..." --output results/

# Search with a local structure file
foldseek-pocket-miner --structure my_protein.pdb --output results/

# Customize search parameters
foldseek-pocket-miner --pdb 1ABC \
    --chain A \
    --max-hits 100 \
    --evalue 1e-5 \
    --min-seq-id 0.3 \
    --output results/
```

### Python API

```python
from foldseek_pocket_miner import PocketMiner

# Initialize the miner
miner = PocketMiner(output_dir="results/")

# Run the full pipeline
session_path = miner.run(
    pdb_id="1ABC",
    chain="A",
    max_hits=50,
    evalue_threshold=1e-3
)

print(f"PyMOL session saved to: {session_path}")
```

## Output

The tool generates:
- `results/structures/` - Downloaded PDB files
- `results/aligned/` - Aligned structure files
- `results/ligands/` - Extracted ligand files
- `results/conservation.csv` - Conservation scores per residue
- `results/pocket_analysis.pse` - PyMOL session file

## Dependencies

- Python >= 3.9
- Foldseek (command-line tool)
- PyMOL
- BioPython
- NumPy
- Pandas
- Requests
- tqdm

## License

MIT License

## Citation

If you use this tool, please cite:

**Foldseek:**
```bibtex
@article{van2024fast,
  title={Fast and accurate protein structure search with Foldseek},
  author={van Kempen, Michel and Kim, Stephanie S and Tumescheit, Charlotte and Mirdita, Milot and Lee, Jeongjae and Gilchrist, Cameron LM and S{\"o}ding, Johannes and Steinegger, Martin},
  journal={Nature Biotechnology},
  volume={42},
  number={2},
  pages={243--246},
  year={2024},
  publisher={Nature Publishing Group US New York}
}
```

**Foldseek Pocket Miner:**
```bibtex
@software{jian2024foldseek_pocket_miner,
  author = {Jian, Jinglin},
  title = {Foldseek Pocket Miner: A Tool for Protein Binding Pocket Discovery},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/JianJinglin/foldseek-pocket-miner}
}
```

## Acknowledgments

- RCSB PDB for structure data
- Foldseek developers

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=JianJinglin/foldseek-pocket-miner&type=Date)](https://star-history.com/#JianJinglin/foldseek-pocket-miner&Date)
