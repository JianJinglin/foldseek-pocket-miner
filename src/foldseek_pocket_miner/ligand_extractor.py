"""
Ligand extractor module - Extract ligands from PDB structures and analyze pharmacophores.
"""

import os
import logging
from pathlib import Path
from typing import List, Optional, Dict, Set, Tuple
from dataclasses import dataclass, field
from collections import defaultdict

import numpy as np

logger = logging.getLogger(__name__)

# Common molecules to exclude (not drug-like)
EXCLUDE_RESIDUES = {
    # Water
    'HOH', 'WAT', 'DOD', 'H2O',
    # Ions
    'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'CO', 'NI', 'CD',
    'BR', 'IOD', 'F', 'SO4', 'PO4', 'NO3',
    # Common buffer components
    'ACT', 'GOL', 'EDO', 'PEG', 'DMS', 'MPD', 'TRS', 'MES', 'EPE',
    'BME', 'DTT', 'CIT', 'TAR', 'FMT', 'ACE', 'NH4',
    # Lipids
    'PLM', 'OLA', 'MYR',
    # Modified amino acids (usually covalent)
    'MSE', 'SEP', 'TPO', 'PTR',
}


@dataclass
class Atom:
    """Represents an atom in a structure."""
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    bfactor: float = 0.0
    element: str = ""

    @property
    def coords(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z])


@dataclass
class Ligand:
    """Represents a ligand molecule."""
    resname: str
    chain: str
    resseq: int
    atoms: List[Atom] = field(default_factory=list)
    source_pdb: str = ""
    source_chain: str = ""

    @property
    def center(self) -> np.ndarray:
        """Calculate center of mass."""
        if not self.atoms:
            return np.zeros(3)
        coords = np.array([[a.x, a.y, a.z] for a in self.atoms])
        return np.mean(coords, axis=0)

    @property
    def num_atoms(self) -> int:
        return len(self.atoms)

    def to_pdb_string(self) -> str:
        """Convert ligand to PDB format string."""
        lines = []
        for atom in self.atoms:
            line = (
                f"HETATM{atom.serial:5d} {atom.name:<4s}{atom.altloc:1s}"
                f"{self.resname:>3s} {atom.chain:1s}{self.resseq:4d}    "
                f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
                f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}          "
                f"{atom.element:>2s}"
            )
            lines.append(line)
        return '\n'.join(lines)


@dataclass
class PharmacophoreFeature:
    """Represents a pharmacophore feature point."""
    feature_type: str  # 'donor', 'acceptor', 'aromatic', 'hydrophobic', 'positive', 'negative'
    position: np.ndarray
    radius: float = 1.5
    direction: Optional[np.ndarray] = None


class LigandExtractor:
    """Extract ligands from protein structures."""

    def __init__(
        self,
        output_dir: str,
        exclude_residues: Optional[Set[str]] = None,
        min_atoms: int = 5
    ):
        """
        Initialize the ligand extractor.

        Args:
            output_dir: Directory to save extracted ligands
            exclude_residues: Set of residue names to exclude
            min_atoms: Minimum number of atoms for a valid ligand
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.exclude_residues = exclude_residues or EXCLUDE_RESIDUES
        self.min_atoms = min_atoms

    def extract_ligands(
        self,
        pdb_path: str,
        chain: Optional[str] = None
    ) -> List[Ligand]:
        """
        Extract ligands from a PDB file.

        Args:
            pdb_path: Path to PDB file
            chain: Optional chain to focus on (ligands near this chain)

        Returns:
            List of Ligand objects
        """
        pdb_path = Path(pdb_path)
        if not pdb_path.exists():
            logger.warning(f"PDB file not found: {pdb_path}")
            return []

        # Parse PDB file
        ligands_dict: Dict[Tuple[str, str, int], List[Atom]] = defaultdict(list)

        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    try:
                        atom = self._parse_hetatm(line)
                        resname = atom.resname.strip()

                        # Skip excluded residues
                        if resname in self.exclude_residues:
                            continue

                        key = (resname, atom.chain, atom.resseq)
                        ligands_dict[key].append(atom)

                    except Exception as e:
                        logger.debug(f"Failed to parse HETATM line: {e}")
                        continue

        # Create Ligand objects
        ligands = []
        for (resname, lig_chain, resseq), atoms in ligands_dict.items():
            if len(atoms) >= self.min_atoms:
                ligand = Ligand(
                    resname=resname,
                    chain=lig_chain,
                    resseq=resseq,
                    atoms=atoms,
                    source_pdb=pdb_path.stem,
                    source_chain=lig_chain
                )
                ligands.append(ligand)

        logger.debug(f"Extracted {len(ligands)} ligands from {pdb_path}")
        return ligands

    def _parse_hetatm(self, line: str) -> Atom:
        """Parse a HETATM line from PDB file."""
        return Atom(
            serial=int(line[6:11].strip()),
            name=line[12:16].strip(),
            altloc=line[16:17].strip(),
            resname=line[17:20].strip(),
            chain=line[21:22].strip(),
            resseq=int(line[22:26].strip()),
            x=float(line[30:38].strip()),
            y=float(line[38:46].strip()),
            z=float(line[46:54].strip()),
            occupancy=float(line[54:60].strip()) if len(line) > 54 else 1.0,
            bfactor=float(line[60:66].strip()) if len(line) > 60 else 0.0,
            element=line[76:78].strip() if len(line) > 76 else ""
        )

    def save_ligand(
        self,
        ligand: Ligand,
        suffix: str = ""
    ) -> str:
        """
        Save a ligand to a PDB file.

        Args:
            ligand: Ligand object to save
            suffix: Optional suffix for filename

        Returns:
            Path to saved file
        """
        filename = f"{ligand.source_pdb}_{ligand.resname}_{ligand.chain}{ligand.resseq}{suffix}.pdb"
        output_path = self.output_dir / filename

        with open(output_path, 'w') as f:
            f.write(f"REMARK   Ligand {ligand.resname} from {ligand.source_pdb}\n")
            f.write(ligand.to_pdb_string())
            f.write("\nEND\n")

        return str(output_path)

    def extract_from_batch(
        self,
        pdb_paths: List[str],
        save: bool = True
    ) -> Dict[str, List[Ligand]]:
        """
        Extract ligands from multiple PDB files.

        Args:
            pdb_paths: List of PDB file paths
            save: Whether to save extracted ligands

        Returns:
            Dictionary mapping PDB path to list of ligands
        """
        from tqdm import tqdm

        results = {}

        for pdb_path in tqdm(pdb_paths, desc="Extracting ligands"):
            ligands = self.extract_ligands(pdb_path)

            if save:
                for ligand in ligands:
                    self.save_ligand(ligand)

            results[pdb_path] = ligands

        total_ligands = sum(len(v) for v in results.values())
        logger.info(f"Extracted {total_ligands} ligands from {len(pdb_paths)} structures")

        return results

    def cluster_ligands(
        self,
        ligands: List[Ligand],
        distance_threshold: float = 3.0
    ) -> List[List[Ligand]]:
        """
        Cluster ligands based on spatial proximity.

        Args:
            ligands: List of ligands to cluster
            distance_threshold: Maximum distance between cluster centers

        Returns:
            List of ligand clusters
        """
        if not ligands:
            return []

        # Calculate centers
        centers = np.array([lig.center for lig in ligands])

        # Simple hierarchical clustering
        from scipy.cluster.hierarchy import fcluster, linkage

        if len(ligands) == 1:
            return [[ligands[0]]]

        Z = linkage(centers, method='average')
        clusters_idx = fcluster(Z, t=distance_threshold, criterion='distance')

        # Group ligands by cluster
        clusters: Dict[int, List[Ligand]] = defaultdict(list)
        for lig, cluster_id in zip(ligands, clusters_idx):
            clusters[cluster_id].append(lig)

        return list(clusters.values())


class PharmacophoreGenerator:
    """Generate pharmacophore features from ligands."""

    # Pharmacophore feature definitions based on atom types
    DONOR_ATOMS = {'N', 'O'}  # with H attached
    ACCEPTOR_ATOMS = {'N', 'O', 'S'}  # with lone pairs
    AROMATIC_RESIDUES = {'PHE', 'TYR', 'TRP', 'HIS'}
    POSITIVE_ATOMS = {'N'}  # in specific contexts
    NEGATIVE_ATOMS = {'O'}  # in COO-

    def __init__(self):
        pass

    def generate_from_ligand(self, ligand: Ligand) -> List[PharmacophoreFeature]:
        """
        Generate pharmacophore features from a ligand.

        Args:
            ligand: Ligand to analyze

        Returns:
            List of PharmacophoreFeature objects
        """
        features = []

        for atom in ligand.atoms:
            element = atom.element.strip().upper() if atom.element else atom.name[0].upper()

            # Simple heuristics for feature assignment
            if element in ['N', 'O']:
                # Could be donor or acceptor
                features.append(PharmacophoreFeature(
                    feature_type='hbond',
                    position=atom.coords,
                    radius=1.5
                ))

            elif element == 'C':
                # Check if aromatic (simplified)
                if 'AR' in atom.name.upper() or any(c.isdigit() for c in atom.name):
                    pass  # Could add aromatic features

        # Add center as hydrophobic feature if molecule is large enough
        if len(ligand.atoms) > 6:
            features.append(PharmacophoreFeature(
                feature_type='hydrophobic',
                position=ligand.center,
                radius=2.0
            ))

        return features

    def generate_from_cluster(
        self,
        ligands: List[Ligand],
        feature_types: Optional[List[str]] = None
    ) -> List[PharmacophoreFeature]:
        """
        Generate consensus pharmacophore from a cluster of ligands.

        Args:
            ligands: Cluster of aligned ligands
            feature_types: Optional list of feature types to include

        Returns:
            List of consensus pharmacophore features
        """
        all_features = []

        for ligand in ligands:
            features = self.generate_from_ligand(ligand)
            all_features.extend(features)

        # Cluster features by position
        if not all_features:
            return []

        positions = np.array([f.position for f in all_features])

        # Simple clustering of feature positions
        from scipy.cluster.hierarchy import fcluster, linkage

        if len(positions) == 1:
            return all_features

        Z = linkage(positions, method='average')
        clusters = fcluster(Z, t=2.0, criterion='distance')

        # Create consensus features
        consensus = []
        cluster_features: Dict[int, List[PharmacophoreFeature]] = defaultdict(list)

        for feat, cluster_id in zip(all_features, clusters):
            cluster_features[cluster_id].append(feat)

        for cluster_id, feats in cluster_features.items():
            # Use centroid as consensus position
            positions = np.array([f.position for f in feats])
            centroid = np.mean(positions, axis=0)

            # Use most common feature type
            types = [f.feature_type for f in feats]
            most_common = max(set(types), key=types.count)

            consensus.append(PharmacophoreFeature(
                feature_type=most_common,
                position=centroid,
                radius=1.5
            ))

        return consensus
