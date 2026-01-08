"""
Structure downloader module - Download protein structures from RCSB PDB.
"""

import os
import logging
import gzip
from pathlib import Path
from typing import List, Optional, Dict, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)


class StructureDownloader:
    """Download protein structures from RCSB PDB."""

    RCSB_BASE_URL = "https://files.rcsb.org/download"
    RCSB_API_URL = "https://data.rcsb.org/rest/v1/core"

    def __init__(
        self,
        output_dir: str,
        format: str = "pdb",
        max_workers: int = 4
    ):
        """
        Initialize the structure downloader.

        Args:
            output_dir: Directory to save downloaded structures
            format: Structure format ('pdb' or 'cif')
            max_workers: Maximum number of parallel downloads
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.format = format.lower()
        self.max_workers = max_workers

        if self.format not in ['pdb', 'cif', 'mmcif']:
            raise ValueError(f"Unsupported format: {format}. Use 'pdb' or 'cif'.")

    def download_structure(
        self,
        pdb_id: str,
        chain: Optional[str] = None,
        overwrite: bool = False
    ) -> Optional[str]:
        """
        Download a single structure from RCSB PDB.

        Args:
            pdb_id: PDB ID (4 characters)
            chain: Optional chain ID (not used for download, but for naming)
            overwrite: Whether to overwrite existing files

        Returns:
            Path to downloaded file, or None if download failed
        """
        pdb_id = pdb_id.upper()

        # Determine output filename
        suffix = "_" + chain if chain else ""
        if self.format == 'pdb':
            filename = f"{pdb_id}{suffix}.pdb"
            url = f"{self.RCSB_BASE_URL}/{pdb_id}.pdb"
        else:
            filename = f"{pdb_id}{suffix}.cif"
            url = f"{self.RCSB_BASE_URL}/{pdb_id}.cif"

        output_path = self.output_dir / filename

        # Check if file already exists
        if output_path.exists() and not overwrite:
            logger.debug(f"Structure already exists: {output_path}")
            return str(output_path)

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            with open(output_path, 'w') as f:
                f.write(response.text)

            logger.debug(f"Downloaded: {pdb_id} -> {output_path}")
            return str(output_path)

        except requests.exceptions.HTTPError as e:
            # Try compressed version
            try:
                gz_url = url + ".gz"
                response = requests.get(gz_url, timeout=30)
                response.raise_for_status()

                content = gzip.decompress(response.content).decode('utf-8')
                with open(output_path, 'w') as f:
                    f.write(content)

                logger.debug(f"Downloaded (gzip): {pdb_id} -> {output_path}")
                return str(output_path)

            except Exception:
                logger.warning(f"Failed to download {pdb_id}: {e}")
                return None

        except Exception as e:
            logger.warning(f"Failed to download {pdb_id}: {e}")
            return None

    def download_batch(
        self,
        pdb_ids: List[str],
        chains: Optional[List[str]] = None,
        overwrite: bool = False,
        show_progress: bool = True
    ) -> Dict[str, Optional[str]]:
        """
        Download multiple structures in parallel.

        Args:
            pdb_ids: List of PDB IDs
            chains: Optional list of chain IDs (parallel to pdb_ids)
            overwrite: Whether to overwrite existing files
            show_progress: Whether to show progress bar

        Returns:
            Dictionary mapping PDB ID to downloaded file path (None if failed)
        """
        if chains is None:
            chains = [None] * len(pdb_ids)

        if len(chains) != len(pdb_ids):
            raise ValueError("Number of chains must match number of PDB IDs")

        results = {}

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(
                    self.download_structure, pdb_id, chain, overwrite
                ): (pdb_id, chain)
                for pdb_id, chain in zip(pdb_ids, chains)
            }

            iterator = as_completed(futures)
            if show_progress:
                iterator = tqdm(iterator, total=len(futures), desc="Downloading structures")

            for future in iterator:
                pdb_id, chain = futures[future]
                try:
                    path = future.result()
                    key = f"{pdb_id}_{chain}" if chain else pdb_id
                    results[key] = path
                except Exception as e:
                    logger.error(f"Error downloading {pdb_id}: {e}")
                    results[pdb_id] = None

        # Log summary
        success = sum(1 for v in results.values() if v is not None)
        logger.info(f"Downloaded {success}/{len(pdb_ids)} structures")

        return results

    def get_structure_info(self, pdb_id: str) -> Optional[Dict]:
        """
        Get metadata about a structure from RCSB.

        Args:
            pdb_id: PDB ID

        Returns:
            Dictionary with structure metadata
        """
        url = f"{self.RCSB_API_URL}/entry/{pdb_id}"

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.warning(f"Failed to get info for {pdb_id}: {e}")
            return None

    def get_ligands(self, pdb_id: str) -> List[Dict]:
        """
        Get list of ligands in a structure.

        Args:
            pdb_id: PDB ID

        Returns:
            List of ligand information dictionaries
        """
        url = f"https://data.rcsb.org/graphql"

        query = """
        query ($pdb_id: String!) {
            entry(entry_id: $pdb_id) {
                nonpolymer_entities {
                    nonpolymer_comp {
                        chem_comp {
                            id
                            name
                            formula
                            formula_weight
                        }
                    }
                    rcsb_nonpolymer_entity {
                        pdbx_description
                    }
                }
            }
        }
        """

        try:
            response = requests.post(
                url,
                json={"query": query, "variables": {"pdb_id": pdb_id}},
                timeout=30
            )
            response.raise_for_status()
            data = response.json()

            ligands = []
            entities = data.get('data', {}).get('entry', {}).get('nonpolymer_entities', [])

            for entity in entities:
                if entity and entity.get('nonpolymer_comp'):
                    comp = entity['nonpolymer_comp'].get('chem_comp', {})
                    ligands.append({
                        'id': comp.get('id'),
                        'name': comp.get('name'),
                        'formula': comp.get('formula'),
                        'formula_weight': comp.get('formula_weight'),
                    })

            return ligands

        except Exception as e:
            logger.warning(f"Failed to get ligands for {pdb_id}: {e}")
            return []

    def filter_by_ligand(
        self,
        pdb_ids: List[str],
        exclude_common: bool = True
    ) -> List[Tuple[str, List[Dict]]]:
        """
        Filter structures that contain ligands.

        Args:
            pdb_ids: List of PDB IDs to check
            exclude_common: Exclude common molecules (water, ions, etc.)

        Returns:
            List of (pdb_id, ligands) tuples for structures with ligands
        """
        # Common molecules to exclude
        common_ids = {
            'HOH', 'WAT', 'DOD',  # Water
            'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'CO',  # Ions
            'SO4', 'PO4', 'ACT', 'GOL', 'EDO', 'PEG', 'DMS',  # Buffers/solvents
        }

        results = []

        for pdb_id in tqdm(pdb_ids, desc="Checking ligands"):
            ligands = self.get_ligands(pdb_id)

            if exclude_common:
                ligands = [l for l in ligands if l.get('id') not in common_ids]

            if ligands:
                results.append((pdb_id, ligands))

        logger.info(f"Found {len(results)}/{len(pdb_ids)} structures with ligands")
        return results
