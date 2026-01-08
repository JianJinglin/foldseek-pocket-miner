"""
Foldseek search module - Interface to Foldseek command line tool for structure similarity search.
"""

import subprocess
import tempfile
import os
import logging
from pathlib import Path
from typing import Optional, List, Dict, Any
from dataclasses import dataclass

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class FoldseekHit:
    """Represents a single Foldseek search hit."""
    query: str
    target: str
    pdb_id: str
    chain: str
    evalue: float
    score: float
    seq_identity: float
    alignment_length: int
    query_start: int
    query_end: int
    target_start: int
    target_end: int

    @classmethod
    def from_row(cls, row: pd.Series) -> "FoldseekHit":
        """Create a FoldseekHit from a DataFrame row."""
        # Parse target to extract PDB ID and chain
        target = row['target']
        # Target format is usually like "1abc_A" or "pdb1abc_A"
        if '_' in target:
            parts = target.split('_')
            pdb_id = parts[0].replace('pdb', '').upper()
            chain = parts[1] if len(parts) > 1 else 'A'
        else:
            pdb_id = target[:4].upper()
            chain = target[4:] if len(target) > 4 else 'A'

        return cls(
            query=row.get('query', ''),
            target=target,
            pdb_id=pdb_id,
            chain=chain,
            evalue=float(row.get('evalue', 0)),
            score=float(row.get('score', 0)),
            seq_identity=float(row.get('fident', 0)),
            alignment_length=int(row.get('alnlen', 0)),
            query_start=int(row.get('qstart', 0)),
            query_end=int(row.get('qend', 0)),
            target_start=int(row.get('tstart', 0)),
            target_end=int(row.get('tend', 0)),
        )


class FoldseekSearcher:
    """Interface to Foldseek for structure similarity search."""

    # Foldseek output format columns
    OUTPUT_COLUMNS = [
        'query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen',
        'qstart', 'qend', 'tstart', 'tend', 'evalue', 'score'
    ]

    def __init__(
        self,
        foldseek_path: str = "foldseek",
        database: str = "pdb100",
        temp_dir: Optional[str] = None
    ):
        """
        Initialize the Foldseek searcher.

        Args:
            foldseek_path: Path to foldseek executable
            database: Database to search (pdb100, afdb, etc.)
            temp_dir: Temporary directory for intermediate files
        """
        self.foldseek_path = foldseek_path
        self.database = database
        self.temp_dir = temp_dir or tempfile.gettempdir()

        # Verify foldseek is available
        self._verify_foldseek()

    def _verify_foldseek(self) -> None:
        """Verify that foldseek is installed and accessible."""
        try:
            result = subprocess.run(
                [self.foldseek_path, "--version"],
                capture_output=True,
                text=True
            )
            logger.info(f"Foldseek version: {result.stdout.strip()}")
        except FileNotFoundError:
            raise RuntimeError(
                f"Foldseek not found at '{self.foldseek_path}'. "
                "Please install foldseek or provide the correct path."
            )

    def search(
        self,
        structure_path: str,
        max_hits: int = 100,
        evalue: float = 1e-3,
        min_seq_id: float = 0.0,
        coverage: float = 0.0,
        threads: int = 4
    ) -> List[FoldseekHit]:
        """
        Search for similar structures using Foldseek.

        Args:
            structure_path: Path to the query structure (PDB or mmCIF)
            max_hits: Maximum number of hits to return
            evalue: E-value threshold
            min_seq_id: Minimum sequence identity (0-1)
            coverage: Minimum coverage (0-1)
            threads: Number of threads to use

        Returns:
            List of FoldseekHit objects
        """
        structure_path = Path(structure_path)
        if not structure_path.exists():
            raise FileNotFoundError(f"Structure file not found: {structure_path}")

        # Create temporary output file
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.tsv', delete=False, dir=self.temp_dir
        ) as tmp_out:
            output_path = tmp_out.name

        try:
            # Build foldseek command
            # Using easy-search which handles database download automatically
            cmd = [
                self.foldseek_path,
                "easy-search",
                str(structure_path),
                self.database,  # Will use remote database or local if available
                output_path,
                os.path.join(self.temp_dir, "foldseek_tmp"),
                "--max-seqs", str(max_hits),
                "-e", str(evalue),
                "--min-seq-id", str(min_seq_id),
                "-c", str(coverage),
                "--threads", str(threads),
                "--format-output", ",".join(self.OUTPUT_COLUMNS),
            ]

            logger.info(f"Running Foldseek: {' '.join(cmd)}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                logger.error(f"Foldseek stderr: {result.stderr}")
                raise RuntimeError(f"Foldseek search failed: {result.stderr}")

            # Parse results
            hits = self._parse_results(output_path)
            logger.info(f"Found {len(hits)} hits")

            return hits

        finally:
            # Cleanup
            if os.path.exists(output_path):
                os.remove(output_path)

    def _parse_results(self, output_path: str) -> List[FoldseekHit]:
        """Parse Foldseek output file."""
        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            return []

        df = pd.read_csv(
            output_path,
            sep='\t',
            names=self.OUTPUT_COLUMNS,
            header=None
        )

        hits = []
        for _, row in df.iterrows():
            try:
                hit = FoldseekHit.from_row(row)
                hits.append(hit)
            except Exception as e:
                logger.warning(f"Failed to parse hit: {e}")
                continue

        return hits

    def search_web(
        self,
        structure_path: str,
        max_hits: int = 100,
        database: str = "pdb100"
    ) -> List[FoldseekHit]:
        """
        Search using Foldseek web server API (alternative when local database not available).

        Args:
            structure_path: Path to the query structure
            max_hits: Maximum number of hits
            database: Database to search

        Returns:
            List of FoldseekHit objects
        """
        import requests

        # Read structure file
        with open(structure_path, 'r') as f:
            structure_content = f.read()

        # Submit to Foldseek web server
        url = "https://search.foldseek.com/api/ticket"

        files = {
            'q': ('query.pdb', structure_content),
        }
        data = {
            'database[]': database,
            'mode': '3diaa',
        }

        response = requests.post(url, files=files, data=data)
        response.raise_for_status()

        ticket = response.json()
        ticket_id = ticket['id']

        # Poll for results
        import time
        result_url = f"https://search.foldseek.com/api/ticket/{ticket_id}"

        while True:
            response = requests.get(result_url)
            result = response.json()

            if result['status'] == 'COMPLETE':
                break
            elif result['status'] == 'ERROR':
                raise RuntimeError("Foldseek web search failed")

            time.sleep(2)

        # Parse results
        hits = []
        for entry in result.get('results', [])[:max_hits]:
            # Parse web results into FoldseekHit format
            # This may need adjustment based on actual API response format
            pass

        return hits


def download_pdb100_database(output_dir: str, foldseek_path: str = "foldseek") -> str:
    """
    Download the PDB100 database for local searches.

    Args:
        output_dir: Directory to store the database
        foldseek_path: Path to foldseek executable

    Returns:
        Path to the downloaded database
    """
    db_path = os.path.join(output_dir, "pdb100")

    cmd = [
        foldseek_path,
        "databases",
        "PDB",
        db_path,
        os.path.join(output_dir, "tmp"),
    ]

    logger.info(f"Downloading PDB100 database to {db_path}")
    subprocess.run(cmd, check=True)

    return db_path
