#!/usr/bin/env python3
"""
Example usage of FoldseekPocketMiner.

This script demonstrates how to use the PocketMiner API to:
1. Search for similar structures using Foldseek
2. Download and align structures
3. Extract ligands and calculate conservation
4. Generate a PyMOL session
"""

from pathlib import Path
from foldseek_pocket_miner import PocketMiner

# Define output directory
OUTPUT_DIR = Path("./example_output")

# Example 1: Search by PDB ID
def example_pdb_search():
    """Search for similar structures to a known PDB entry."""
    print("Example 1: Search by PDB ID")
    print("-" * 40)

    # Initialize miner
    miner = PocketMiner(
        output_dir=str(OUTPUT_DIR / "pdb_search"),
        database="pdb100"
    )

    # Run the pipeline
    result = miner.run(
        pdb_id="3HTB",  # Human CDK2 kinase
        chain="A",
        max_hits=50,
        evalue_threshold=1e-3,
        only_with_ligands=True,
        session_name="cdk2_pocket"
    )

    print(f"Found {len(result.foldseek_hits)} similar structures")
    print(f"Extracted {len(result.ligands)} ligands")
    print(f"Session saved to: {result.session_path}")

    return result


# Example 2: Search with a local structure file
def example_structure_file(structure_path: str):
    """Search using a local structure file."""
    print("\nExample 2: Search with local structure")
    print("-" * 40)

    miner = PocketMiner(
        output_dir=str(OUTPUT_DIR / "structure_search"),
        database="pdb100"
    )

    result = miner.run(
        structure_path=structure_path,
        max_hits=100,
        only_with_ligands=True,
        session_name="my_protein_pocket"
    )

    return result


# Example 3: Quick search with minimal parameters
def example_quick_search():
    """Quick search with default parameters."""
    print("\nExample 3: Quick search")
    print("-" * 40)

    miner = PocketMiner(output_dir=str(OUTPUT_DIR / "quick"))

    # Use the quick run method
    session_path = miner.run_quick(
        pdb_id="1ATP",  # cAMP-dependent protein kinase
        chain="E",
        max_hits=30
    )

    print(f"Session saved to: {session_path}")
    return session_path


# Example 4: Using individual components
def example_components():
    """Demonstrate using individual components."""
    print("\nExample 4: Using individual components")
    print("-" * 40)

    from foldseek_pocket_miner import (
        FoldseekSearcher,
        StructureDownloader,
        LigandExtractor,
        ConservationAnalyzer
    )

    # Search
    searcher = FoldseekSearcher(database="pdb100")
    # hits = searcher.search("my_structure.pdb", max_hits=50)

    # Download
    downloader = StructureDownloader(output_dir="./structures")
    # path = downloader.download_structure("1ABC", chain="A")

    # Extract ligands
    extractor = LigandExtractor(output_dir="./ligands")
    # ligands = extractor.extract_ligands("structure.pdb")

    # Analyze conservation
    analyzer = ConservationAnalyzer(output_dir="./conservation")
    # conservation = analyzer.calculate_conservation(query, aligned_list)

    print("Components initialized successfully")


if __name__ == "__main__":
    # Run Example 1 (PDB search)
    try:
        result = example_pdb_search()
    except Exception as e:
        print(f"Example 1 failed: {e}")
        print("Make sure foldseek and pymol are installed")

    # Quick search example
    try:
        example_quick_search()
    except Exception as e:
        print(f"Example 3 failed: {e}")

    # Components example (always works)
    example_components()

    print("\n" + "=" * 40)
    print("Examples completed!")
    print("=" * 40)
