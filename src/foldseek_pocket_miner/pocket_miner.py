"""
Main PocketMiner class - Orchestrates the full pipeline.
"""

import os
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any
from dataclasses import dataclass, field

from .foldseek_search import FoldseekSearcher, FoldseekHit
from .structure_downloader import StructureDownloader
from .structure_aligner import StructureAligner
from .ligand_extractor import LigandExtractor, Ligand
from .conservation import ConservationAnalyzer, ResidueConservation
from .pymol_session import PyMOLSessionGenerator, SessionConfig

logger = logging.getLogger(__name__)


@dataclass
class PocketMinerResult:
    """Results from PocketMiner pipeline."""
    query_structure: str
    foldseek_hits: List[FoldseekHit]
    downloaded_structures: Dict[str, str]
    aligned_structures: List[str]
    ligands: List[Ligand]
    conservation: List[ResidueConservation]
    session_path: str
    output_dir: str


class PocketMiner:
    """
    Main class for pocket mining workflow.

    Pipeline:
    1. Search for similar structures using Foldseek
    2. Download matching structures from RCSB
    3. Align structures to query
    4. Extract ligands from aligned structures
    5. Calculate sequence conservation
    6. Generate PyMOL session
    """

    def __init__(
        self,
        output_dir: str,
        foldseek_path: str = "foldseek",
        pymol_path: str = "pymol",
        database: str = "pdb100",
        n_threads: int = 4
    ):
        """
        Initialize PocketMiner.

        Args:
            output_dir: Base output directory
            foldseek_path: Path to foldseek executable
            pymol_path: Path to PyMOL executable
            database: Foldseek database to search
            n_threads: Number of threads for parallel operations
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        self.structures_dir = self.output_dir / "structures"
        self.aligned_dir = self.output_dir / "aligned"
        self.ligands_dir = self.output_dir / "ligands"

        for d in [self.structures_dir, self.aligned_dir, self.ligands_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # Initialize components
        self.searcher = FoldseekSearcher(
            foldseek_path=foldseek_path,
            database=database
        )
        self.downloader = StructureDownloader(
            output_dir=str(self.structures_dir),
            max_workers=n_threads
        )
        self.aligner = StructureAligner(
            output_dir=str(self.aligned_dir),
            pymol_path=pymol_path
        )
        self.ligand_extractor = LigandExtractor(
            output_dir=str(self.ligands_dir)
        )
        self.conservation_analyzer = ConservationAnalyzer(
            output_dir=str(self.output_dir)
        )
        self.session_generator = PyMOLSessionGenerator(
            output_dir=str(self.output_dir),
            pymol_path=pymol_path
        )

        self.n_threads = n_threads

    def run(
        self,
        pdb_id: Optional[str] = None,
        chain: Optional[str] = None,
        structure_path: Optional[str] = None,
        sequence: Optional[str] = None,
        max_hits: int = 100,
        evalue_threshold: float = 1e-3,
        min_seq_id: float = 0.0,
        only_with_ligands: bool = True,
        session_name: str = "pocket_analysis"
    ) -> PocketMinerResult:
        """
        Run the full pocket mining pipeline.

        Args:
            pdb_id: PDB ID of query structure (provide this or structure_path)
            chain: Chain ID in query structure
            structure_path: Path to query structure file
            sequence: Amino acid sequence (for structure prediction)
            max_hits: Maximum Foldseek hits
            evalue_threshold: E-value threshold for Foldseek
            min_seq_id: Minimum sequence identity
            only_with_ligands: Only keep hits that have ligands
            session_name: Name for output session

        Returns:
            PocketMinerResult with all outputs
        """
        logger.info("Starting PocketMiner pipeline")

        # Step 1: Get query structure
        logger.info("Step 1: Preparing query structure")
        query_path = self._prepare_query(pdb_id, chain, structure_path, sequence)

        # Step 2: Search for similar structures
        logger.info("Step 2: Searching for similar structures with Foldseek")
        hits = self.searcher.search(
            structure_path=query_path,
            max_hits=max_hits,
            evalue=evalue_threshold,
            min_seq_id=min_seq_id
        )
        logger.info(f"Found {len(hits)} Foldseek hits")

        if not hits:
            logger.warning("No hits found. Cannot proceed.")
            return PocketMinerResult(
                query_structure=query_path,
                foldseek_hits=[],
                downloaded_structures={},
                aligned_structures=[],
                ligands=[],
                conservation=[],
                session_path="",
                output_dir=str(self.output_dir)
            )

        # Step 3: Download structures
        logger.info("Step 3: Downloading structures from RCSB")
        pdb_ids = [hit.pdb_id for hit in hits]
        chains = [hit.chain for hit in hits]
        downloaded = self.downloader.download_batch(pdb_ids, chains)

        # Filter successful downloads
        valid_paths = {k: v for k, v in downloaded.items() if v is not None}
        logger.info(f"Successfully downloaded {len(valid_paths)} structures")

        # Step 4: Filter by ligand presence (optional)
        if only_with_ligands:
            logger.info("Step 4: Filtering structures with ligands")
            structures_with_ligands = self.downloader.filter_by_ligand(
                list(valid_paths.keys())
            )
            # Update valid paths
            ligand_pdb_ids = {item[0] for item in structures_with_ligands}
            valid_paths = {k: v for k, v in valid_paths.items()
                          if k.split('_')[0] in ligand_pdb_ids}
            logger.info(f"Kept {len(valid_paths)} structures with ligands")

        # Step 5: Align structures
        logger.info("Step 5: Aligning structures to query")
        target_paths = list(valid_paths.values())
        alignment_results = self.aligner.align_batch(
            query_path=query_path,
            target_paths=target_paths,
            query_chain=chain
        )
        aligned_paths = [r.aligned_path for r in alignment_results]
        logger.info(f"Successfully aligned {len(aligned_paths)} structures")

        # Step 6: Extract ligands
        logger.info("Step 6: Extracting ligands from aligned structures")
        all_ligands = []
        for aligned_path in aligned_paths:
            ligands = self.ligand_extractor.extract_ligands(aligned_path)
            for lig in ligands:
                lig_path = self.ligand_extractor.save_ligand(lig)
                all_ligands.append(lig)
        logger.info(f"Extracted {len(all_ligands)} ligands")

        # Step 7: Calculate conservation
        logger.info("Step 7: Calculating sequence conservation")
        conservation = self.conservation_analyzer.calculate_conservation(
            query_path=query_path,
            aligned_paths=aligned_paths,
            query_chain=chain
        )
        self.conservation_analyzer.save_conservation(conservation)
        logger.info(f"Calculated conservation for {len(conservation)} residues")

        # Step 8: Generate PyMOL session
        logger.info("Step 8: Generating PyMOL session")
        ligand_files = [str(self.ligands_dir / f"{lig.source_pdb}_{lig.resname}_{lig.chain}{lig.resseq}.pdb")
                       for lig in all_ligands]
        ligand_files = [f for f in ligand_files if Path(f).exists()]

        conservation_scores = self.conservation_analyzer.get_conservation_bfactors(conservation)

        session_path = self.session_generator.generate_session(
            query_structure=query_path,
            ligand_files=ligand_files,
            conservation_scores=conservation_scores,
            output_name=session_name,
            aligned_structures=aligned_paths[:5]  # Include top 5 aligned structures
        )

        logger.info(f"Pipeline complete! Session saved to: {session_path}")

        return PocketMinerResult(
            query_structure=query_path,
            foldseek_hits=hits,
            downloaded_structures=downloaded,
            aligned_structures=aligned_paths,
            ligands=all_ligands,
            conservation=conservation,
            session_path=session_path,
            output_dir=str(self.output_dir)
        )

    def _prepare_query(
        self,
        pdb_id: Optional[str],
        chain: Optional[str],
        structure_path: Optional[str],
        sequence: Optional[str]
    ) -> str:
        """
        Prepare the query structure.

        Args:
            pdb_id: PDB ID to download
            chain: Chain ID
            structure_path: Existing structure file
            sequence: Amino acid sequence

        Returns:
            Path to query structure file
        """
        if structure_path:
            # Use existing file
            if not Path(structure_path).exists():
                raise FileNotFoundError(f"Structure file not found: {structure_path}")
            return structure_path

        elif pdb_id:
            # Download from RCSB
            query_downloader = StructureDownloader(
                output_dir=str(self.output_dir / "query")
            )
            path = query_downloader.download_structure(pdb_id, chain)
            if not path:
                raise RuntimeError(f"Failed to download query structure: {pdb_id}")
            return path

        elif sequence:
            # TODO: Implement structure prediction using ESMFold or AlphaFold
            raise NotImplementedError(
                "Structure prediction from sequence not yet implemented. "
                "Please provide a PDB ID or structure file."
            )

        else:
            raise ValueError(
                "Must provide either pdb_id, structure_path, or sequence"
            )

    def run_quick(
        self,
        pdb_id: str,
        chain: str = "A",
        max_hits: int = 50
    ) -> str:
        """
        Quick run with default parameters.

        Args:
            pdb_id: PDB ID of query
            chain: Chain ID
            max_hits: Maximum hits

        Returns:
            Path to generated session file
        """
        result = self.run(
            pdb_id=pdb_id,
            chain=chain,
            max_hits=max_hits,
            only_with_ligands=True,
            session_name=f"{pdb_id}_{chain}_pocket"
        )
        return result.session_path
