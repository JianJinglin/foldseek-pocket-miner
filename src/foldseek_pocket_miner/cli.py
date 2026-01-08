"""
Command-line interface for FoldseekPocketMiner.
"""

import argparse
import logging
import sys
from pathlib import Path

from .pocket_miner import PocketMiner


def setup_logging(verbose: bool = False) -> None:
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def main():
    """Main entry point for CLI."""
    parser = argparse.ArgumentParser(
        description="Foldseek Pocket Miner - Find similar structures and extract binding pockets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Search by PDB ID
  foldseek-pocket-miner --pdb 1ABC --chain A --output results/

  # Search with local structure file
  foldseek-pocket-miner --structure my_protein.pdb --output results/

  # Customize search parameters
  foldseek-pocket-miner --pdb 1ABC --chain A --max-hits 100 --evalue 1e-5 --output results/
        """
    )

    # Input options (mutually exclusive)
    input_group = parser.add_argument_group('Input options (choose one)')
    input_group.add_argument(
        '--pdb', '-p',
        type=str,
        help='PDB ID of the query structure (e.g., 1ABC)'
    )
    input_group.add_argument(
        '--structure', '-s',
        type=str,
        help='Path to query structure file (PDB or mmCIF)'
    )
    input_group.add_argument(
        '--sequence', '-q',
        type=str,
        help='Amino acid sequence (will use structure prediction)'
    )

    # Chain option
    parser.add_argument(
        '--chain', '-c',
        type=str,
        default='A',
        help='Chain ID in query structure (default: A)'
    )

    # Output options
    parser.add_argument(
        '--output', '-o',
        type=str,
        default='./pocket_miner_output',
        help='Output directory (default: ./pocket_miner_output)'
    )
    parser.add_argument(
        '--name', '-n',
        type=str,
        default='pocket_analysis',
        help='Name for output session file (default: pocket_analysis)'
    )

    # Search parameters
    search_group = parser.add_argument_group('Search parameters')
    search_group.add_argument(
        '--max-hits', '-m',
        type=int,
        default=100,
        help='Maximum number of Foldseek hits (default: 100)'
    )
    search_group.add_argument(
        '--evalue', '-e',
        type=float,
        default=1e-3,
        help='E-value threshold for Foldseek (default: 1e-3)'
    )
    search_group.add_argument(
        '--min-seq-id',
        type=float,
        default=0.0,
        help='Minimum sequence identity (0-1, default: 0.0)'
    )
    search_group.add_argument(
        '--database', '-d',
        type=str,
        default='pdb100',
        choices=['pdb100', 'afdb', 'afdb50', 'cath50', 'mgnify_esm30'],
        help='Foldseek database to search (default: pdb100)'
    )

    # Tool paths
    tools_group = parser.add_argument_group('Tool paths')
    tools_group.add_argument(
        '--foldseek-path',
        type=str,
        default='foldseek',
        help='Path to foldseek executable (default: foldseek)'
    )
    tools_group.add_argument(
        '--pymol-path',
        type=str,
        default='pymol',
        help='Path to PyMOL executable (default: pymol)'
    )

    # Additional options
    parser.add_argument(
        '--no-ligand-filter',
        action='store_true',
        help='Include structures without ligands'
    )
    parser.add_argument(
        '--threads', '-t',
        type=int,
        default=4,
        help='Number of threads for parallel operations (default: 4)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 0.1.0'
    )

    args = parser.parse_args()

    # Validate input
    if not any([args.pdb, args.structure, args.sequence]):
        parser.error("Must provide one of: --pdb, --structure, or --sequence")

    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    # Initialize PocketMiner
    logger.info("Initializing FoldseekPocketMiner")
    miner = PocketMiner(
        output_dir=args.output,
        foldseek_path=args.foldseek_path,
        pymol_path=args.pymol_path,
        database=args.database,
        n_threads=args.threads
    )

    try:
        # Run pipeline
        result = miner.run(
            pdb_id=args.pdb,
            chain=args.chain,
            structure_path=args.structure,
            sequence=args.sequence,
            max_hits=args.max_hits,
            evalue_threshold=args.evalue,
            min_seq_id=args.min_seq_id,
            only_with_ligands=not args.no_ligand_filter,
            session_name=args.name
        )

        # Print summary
        print("\n" + "=" * 60)
        print("POCKET MINER RESULTS")
        print("=" * 60)
        print(f"Query structure: {result.query_structure}")
        print(f"Foldseek hits: {len(result.foldseek_hits)}")
        print(f"Downloaded structures: {len([v for v in result.downloaded_structures.values() if v])}")
        print(f"Aligned structures: {len(result.aligned_structures)}")
        print(f"Extracted ligands: {len(result.ligands)}")
        print(f"Conservation scores: {len(result.conservation)} residues")
        print(f"\nOutput directory: {result.output_dir}")
        print(f"PyMOL session: {result.session_path}")
        print("=" * 60)

        if result.ligands:
            print("\nExtracted ligands:")
            for lig in result.ligands[:10]:  # Show first 10
                print(f"  - {lig.resname} from {lig.source_pdb} chain {lig.chain}")
            if len(result.ligands) > 10:
                print(f"  ... and {len(result.ligands) - 10} more")

        print(f"\nOpen the session in PyMOL:")
        print(f"  pymol {result.session_path}")

        return 0

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
