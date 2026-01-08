"""
Conservation analysis module - Calculate sequence conservation from structural alignments.
"""

import os
import logging
from pathlib import Path
from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass
from collections import Counter

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Standard amino acid codes
AA_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# Amino acid groupings for conservation scoring
AA_GROUPS = {
    'A': 'aliphatic', 'V': 'aliphatic', 'L': 'aliphatic', 'I': 'aliphatic', 'M': 'aliphatic',
    'F': 'aromatic', 'W': 'aromatic', 'Y': 'aromatic',
    'K': 'positive', 'R': 'positive', 'H': 'positive',
    'D': 'negative', 'E': 'negative',
    'S': 'polar', 'T': 'polar', 'N': 'polar', 'Q': 'polar',
    'C': 'special', 'G': 'special', 'P': 'special',
}


@dataclass
class ResidueConservation:
    """Conservation data for a single residue."""
    resseq: int
    resname: str
    chain: str
    conservation_score: float  # 0-1, higher = more conserved
    frequency: Dict[str, float]  # Amino acid frequencies at this position
    num_sequences: int


class ConservationAnalyzer:
    """Calculate sequence conservation from aligned structures."""

    def __init__(self, output_dir: Optional[str] = None):
        """
        Initialize the conservation analyzer.

        Args:
            output_dir: Directory to save conservation data
        """
        self.output_dir = Path(output_dir) if output_dir else None
        if self.output_dir:
            self.output_dir.mkdir(parents=True, exist_ok=True)

    def extract_sequence(
        self,
        pdb_path: str,
        chain: Optional[str] = None
    ) -> Dict[int, str]:
        """
        Extract sequence from a PDB file.

        Args:
            pdb_path: Path to PDB file
            chain: Optional chain to extract

        Returns:
            Dictionary mapping residue number to amino acid one-letter code
        """
        sequence = {}

        with open(pdb_path, 'r') as f:
            seen_residues = set()

            for line in f:
                if line.startswith('ATOM'):
                    res_chain = line[21:22].strip()
                    if chain and res_chain != chain:
                        continue

                    resseq = int(line[22:26].strip())
                    resname = line[17:20].strip()

                    if (res_chain, resseq) not in seen_residues:
                        seen_residues.add((res_chain, resseq))
                        aa = AA_CODES.get(resname, 'X')
                        sequence[resseq] = aa

        return sequence

    def calculate_conservation(
        self,
        query_path: str,
        aligned_paths: List[str],
        query_chain: Optional[str] = None,
        method: str = "shannon"
    ) -> List[ResidueConservation]:
        """
        Calculate conservation scores for query structure.

        Args:
            query_path: Path to query structure
            aligned_paths: List of paths to aligned structures
            query_chain: Chain in query structure
            method: Conservation scoring method ('shannon', 'identity', 'similarity')

        Returns:
            List of ResidueConservation objects
        """
        # Extract query sequence
        query_seq = self.extract_sequence(query_path, query_chain)

        if not query_seq:
            logger.warning(f"No sequence found in query: {query_path}")
            return []

        # Extract aligned sequences
        aligned_seqs = []
        for path in aligned_paths:
            seq = self.extract_sequence(path)
            if seq:
                aligned_seqs.append(seq)

        if not aligned_seqs:
            logger.warning("No aligned sequences found")
            return []

        # Calculate conservation per position
        conservation = []

        for resseq, query_aa in sorted(query_seq.items()):
            # Collect amino acids at this position from aligned structures
            aas_at_position = [query_aa]  # Include query

            for seq in aligned_seqs:
                if resseq in seq:
                    aas_at_position.append(seq[resseq])

            # Calculate frequency
            counts = Counter(aas_at_position)
            total = len(aas_at_position)
            frequency = {aa: count / total for aa, count in counts.items()}

            # Calculate conservation score
            if method == "shannon":
                score = self._shannon_entropy_score(frequency)
            elif method == "identity":
                score = self._identity_score(frequency, query_aa)
            elif method == "similarity":
                score = self._similarity_score(frequency, query_aa)
            else:
                score = self._shannon_entropy_score(frequency)

            conservation.append(ResidueConservation(
                resseq=resseq,
                resname=query_aa,
                chain=query_chain or 'A',
                conservation_score=score,
                frequency=frequency,
                num_sequences=total
            ))

        return conservation

    def _shannon_entropy_score(self, frequency: Dict[str, float]) -> float:
        """
        Calculate conservation score based on Shannon entropy.

        Lower entropy = higher conservation.

        Args:
            frequency: Amino acid frequencies

        Returns:
            Conservation score (0-1, higher = more conserved)
        """
        entropy = 0.0
        for freq in frequency.values():
            if freq > 0:
                entropy -= freq * np.log2(freq)

        # Normalize by maximum possible entropy (log2(20) for 20 amino acids)
        max_entropy = np.log2(20)
        normalized_entropy = entropy / max_entropy

        # Convert to conservation score (1 - normalized entropy)
        return 1.0 - normalized_entropy

    def _identity_score(self, frequency: Dict[str, float], query_aa: str) -> float:
        """
        Calculate simple identity-based conservation score.

        Args:
            frequency: Amino acid frequencies
            query_aa: Query amino acid

        Returns:
            Conservation score (0-1)
        """
        return frequency.get(query_aa, 0.0)

    def _similarity_score(self, frequency: Dict[str, float], query_aa: str) -> float:
        """
        Calculate similarity-based conservation score.

        Considers amino acid properties.

        Args:
            frequency: Amino acid frequencies
            query_aa: Query amino acid

        Returns:
            Conservation score (0-1)
        """
        query_group = AA_GROUPS.get(query_aa, 'other')

        score = 0.0
        for aa, freq in frequency.items():
            aa_group = AA_GROUPS.get(aa, 'other')
            if aa == query_aa:
                score += freq  # Full score for identity
            elif aa_group == query_group:
                score += freq * 0.7  # Partial score for similar
            else:
                score += freq * 0.2  # Low score for different

        return min(1.0, score)

    def to_dataframe(
        self,
        conservation: List[ResidueConservation]
    ) -> pd.DataFrame:
        """Convert conservation data to DataFrame."""
        data = []
        for res in conservation:
            data.append({
                'resseq': res.resseq,
                'resname': res.resname,
                'chain': res.chain,
                'conservation_score': res.conservation_score,
                'num_sequences': res.num_sequences,
            })
        return pd.DataFrame(data)

    def save_conservation(
        self,
        conservation: List[ResidueConservation],
        filename: str = "conservation.csv"
    ) -> str:
        """Save conservation data to CSV."""
        if not self.output_dir:
            raise ValueError("No output directory specified")

        df = self.to_dataframe(conservation)
        output_path = self.output_dir / filename
        df.to_csv(output_path, index=False)

        logger.info(f"Saved conservation data to {output_path}")
        return str(output_path)

    def get_conservation_bfactors(
        self,
        conservation: List[ResidueConservation]
    ) -> Dict[int, float]:
        """
        Convert conservation scores to B-factor values for PyMOL coloring.

        Args:
            conservation: List of conservation data

        Returns:
            Dictionary mapping residue number to B-factor value
        """
        # Scale scores to 0-100 range for B-factor coloring
        return {res.resseq: res.conservation_score * 100 for res in conservation}


def calculate_pocket_conservation(
    conservation: List[ResidueConservation],
    pocket_residues: List[int]
) -> Tuple[float, float]:
    """
    Calculate average conservation for pocket residues.

    Args:
        conservation: Full conservation data
        pocket_residues: List of residue numbers in the pocket

    Returns:
        Tuple of (average conservation, standard deviation)
    """
    pocket_scores = []

    for res in conservation:
        if res.resseq in pocket_residues:
            pocket_scores.append(res.conservation_score)

    if not pocket_scores:
        return 0.0, 0.0

    return np.mean(pocket_scores), np.std(pocket_scores)
