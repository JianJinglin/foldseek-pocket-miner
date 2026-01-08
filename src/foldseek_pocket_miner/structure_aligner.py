"""
Structure alignment module - Align protein structures using various methods.
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class AlignmentResult:
    """Result of structure alignment."""
    query_path: str
    target_path: str
    aligned_path: str
    rmsd: float
    aligned_residues: int
    seq_identity: float
    transformation_matrix: Optional[np.ndarray] = None


class StructureAligner:
    """Align protein structures using various methods."""

    def __init__(
        self,
        output_dir: str,
        method: str = "pymol",
        pymol_path: str = "pymol"
    ):
        """
        Initialize the structure aligner.

        Args:
            output_dir: Directory to save aligned structures
            method: Alignment method ('pymol', 'foldseek', 'tmalign')
            pymol_path: Path to PyMOL executable
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.method = method
        self.pymol_path = pymol_path

    def align(
        self,
        query_path: str,
        target_path: str,
        query_chain: Optional[str] = None,
        target_chain: Optional[str] = None
    ) -> Optional[AlignmentResult]:
        """
        Align target structure to query structure.

        Args:
            query_path: Path to query structure
            target_path: Path to target structure to align
            query_chain: Chain ID in query structure
            target_chain: Chain ID in target structure

        Returns:
            AlignmentResult object, or None if alignment failed
        """
        if self.method == "pymol":
            return self._align_pymol(query_path, target_path, query_chain, target_chain)
        elif self.method == "foldseek":
            return self._align_foldseek(query_path, target_path)
        elif self.method == "tmalign":
            return self._align_tmalign(query_path, target_path)
        else:
            raise ValueError(f"Unknown alignment method: {self.method}")

    def _align_pymol(
        self,
        query_path: str,
        target_path: str,
        query_chain: Optional[str] = None,
        target_chain: Optional[str] = None
    ) -> Optional[AlignmentResult]:
        """Align using PyMOL's super command."""
        query_name = Path(query_path).stem
        target_name = Path(target_path).stem
        output_path = self.output_dir / f"{target_name}_aligned.pdb"

        # Build PyMOL script
        query_sel = f"{query_name} and chain {query_chain}" if query_chain else query_name
        target_sel = f"{target_name} and chain {target_chain}" if target_chain else target_name

        script = f"""
import pymol
from pymol import cmd

# Load structures
cmd.load("{query_path}", "{query_name}")
cmd.load("{target_path}", "{target_name}")

# Perform alignment using super (structure-based)
result = cmd.super("{target_sel}", "{query_sel}")

# result is (RMSD, aligned_atoms, cycles, RMSD_initial, aligned_atoms_initial, ...)
rmsd = result[0]
aligned_atoms = result[1]

# Save aligned structure
cmd.save("{output_path}", "{target_name}")

# Print results for parsing
print(f"ALIGNMENT_RESULT: rmsd={{rmsd}}, atoms={{aligned_atoms}}")

cmd.quit()
"""

        script_path = self.output_dir / f"align_{target_name}.py"
        with open(script_path, 'w') as f:
            f.write(script)

        try:
            result = subprocess.run(
                [self.pymol_path, "-cq", str(script_path)],
                capture_output=True,
                text=True,
                timeout=120
            )

            # Parse results from output
            rmsd = 0.0
            aligned_atoms = 0

            for line in result.stdout.split('\n'):
                if 'ALIGNMENT_RESULT:' in line:
                    parts = line.split('ALIGNMENT_RESULT:')[1].strip()
                    for part in parts.split(','):
                        key, val = part.strip().split('=')
                        if key == 'rmsd':
                            rmsd = float(val)
                        elif key == 'atoms':
                            aligned_atoms = int(val)

            if output_path.exists():
                return AlignmentResult(
                    query_path=query_path,
                    target_path=target_path,
                    aligned_path=str(output_path),
                    rmsd=rmsd,
                    aligned_residues=aligned_atoms // 4,  # Approximate residues from atoms
                    seq_identity=0.0  # Not available from super
                )
            else:
                logger.warning(f"Alignment failed: output file not created")
                return None

        except subprocess.TimeoutExpired:
            logger.warning(f"Alignment timed out for {target_path}")
            return None
        except Exception as e:
            logger.warning(f"Alignment failed: {e}")
            return None
        finally:
            if script_path.exists():
                script_path.unlink()

    def _align_foldseek(
        self,
        query_path: str,
        target_path: str
    ) -> Optional[AlignmentResult]:
        """Align using Foldseek structurealign."""
        target_name = Path(target_path).stem
        output_path = self.output_dir / f"{target_name}_aligned.pdb"

        # Foldseek can also do structural alignment
        cmd = [
            "foldseek",
            "structurealign",
            query_path,
            target_path,
            str(output_path),
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            if output_path.exists():
                return AlignmentResult(
                    query_path=query_path,
                    target_path=target_path,
                    aligned_path=str(output_path),
                    rmsd=0.0,  # Parse from output if available
                    aligned_residues=0,
                    seq_identity=0.0
                )
            return None

        except Exception as e:
            logger.warning(f"Foldseek alignment failed: {e}")
            return None

    def _align_tmalign(
        self,
        query_path: str,
        target_path: str
    ) -> Optional[AlignmentResult]:
        """Align using TM-align."""
        target_name = Path(target_path).stem
        output_path = self.output_dir / f"{target_name}_aligned.pdb"

        cmd = [
            "TMalign",
            target_path,
            query_path,
            "-o", str(output_path)
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            # Parse TM-align output
            rmsd = 0.0
            tm_score = 0.0
            aligned_residues = 0

            for line in result.stdout.split('\n'):
                if 'RMSD' in line and 'aligned residues' in line:
                    # Parse: RMSD=  2.05, TM-score= 0.71234, aligned residues=  145
                    parts = line.split(',')
                    for part in parts:
                        if 'RMSD=' in part:
                            rmsd = float(part.split('=')[1].strip())
                        elif 'aligned residues=' in part:
                            aligned_residues = int(part.split('=')[1].strip())

            if output_path.exists():
                return AlignmentResult(
                    query_path=query_path,
                    target_path=target_path,
                    aligned_path=str(output_path),
                    rmsd=rmsd,
                    aligned_residues=aligned_residues,
                    seq_identity=0.0
                )
            return None

        except FileNotFoundError:
            logger.warning("TMalign not found. Please install TM-align.")
            return None
        except Exception as e:
            logger.warning(f"TM-align failed: {e}")
            return None

    def align_batch(
        self,
        query_path: str,
        target_paths: List[str],
        query_chain: Optional[str] = None,
        target_chains: Optional[List[str]] = None,
        show_progress: bool = True
    ) -> List[AlignmentResult]:
        """
        Align multiple structures to a query.

        Args:
            query_path: Path to query structure
            target_paths: List of paths to target structures
            query_chain: Chain ID in query
            target_chains: List of chain IDs in targets
            show_progress: Whether to show progress bar

        Returns:
            List of successful AlignmentResult objects
        """
        from tqdm import tqdm

        if target_chains is None:
            target_chains = [None] * len(target_paths)

        results = []
        iterator = zip(target_paths, target_chains)

        if show_progress:
            iterator = tqdm(list(iterator), desc="Aligning structures")

        for target_path, target_chain in iterator:
            result = self.align(query_path, target_path, query_chain, target_chain)
            if result:
                results.append(result)

        logger.info(f"Successfully aligned {len(results)}/{len(target_paths)} structures")
        return results


def calculate_transformation_matrix(
    coords1: np.ndarray,
    coords2: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Calculate optimal rotation and translation to superimpose coords2 onto coords1.

    Uses Kabsch algorithm.

    Args:
        coords1: Reference coordinates (N x 3)
        coords2: Coordinates to transform (N x 3)

    Returns:
        Tuple of (rotation matrix, translation vector, RMSD)
    """
    # Center coordinates
    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)

    centered1 = coords1 - centroid1
    centered2 = coords2 - centroid2

    # Compute covariance matrix
    H = centered2.T @ centered1

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Rotation matrix
    R = Vt.T @ U.T

    # Handle reflection case
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # Translation
    t = centroid1 - R @ centroid2

    # Calculate RMSD
    transformed = (R @ coords2.T).T + t
    rmsd = np.sqrt(np.mean(np.sum((coords1 - transformed) ** 2, axis=1)))

    return R, t, rmsd
