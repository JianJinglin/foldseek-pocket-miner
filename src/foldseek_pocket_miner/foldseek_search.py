"""
Foldseek search module - Interface to Foldseek web API for structure similarity search.
"""

import subprocess
import tempfile
import os
import logging
import time
from pathlib import Path
from typing import Optional, List, Dict, Any
from dataclasses import dataclass

import requests
import urllib3
import pandas as pd

# Suppress SSL warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

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
        target = row['target']
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

    @classmethod
    def from_web_result(cls, result: Dict, query_name: str = "query") -> "FoldseekHit":
        """Create a FoldseekHit from Foldseek web API result."""
        target = result.get('target', '')

        # Parse target ID - formats vary:
        # New format: "3hpq-assembly2.cif.gz_B Crystal structure description..."
        # Old format: "pdb100_1abc_A" or "1abc_A"
        pdb_id = ''
        chain = 'A'

        # Try new format first: PDB ID is before first '-' or before first '_'
        if '-' in target:
            # Format: "3hpq-assembly2.cif.gz_B description"
            pdb_part = target.split('-')[0]
            pdb_id = pdb_part.upper()
            # Chain is after the last underscore in the filename part (before space)
            filename_part = target.split()[0] if ' ' in target else target
            if '_' in filename_part:
                chain = filename_part.split('_')[-1]
        elif '_' in target:
            # Old format: "pdb100_1abc_A" or "1abc_A"
            parts = target.split('_')
            if len(parts) >= 2:
                # Skip database prefix if present
                if parts[0] in ['pdb100', 'afdb', 'afdb50']:
                    pdb_id = parts[1].upper()
                    chain = parts[2] if len(parts) > 2 else 'A'
                else:
                    pdb_id = parts[0].upper()
                    chain = parts[1] if len(parts) > 1 else 'A'
        else:
            pdb_id = target[:4].upper()

        # Log the parsed values for debugging
        logger.debug(f"Parsed target '{target[:50]}...' -> pdb_id='{pdb_id}', chain='{chain}'")

        return cls(
            query=query_name,
            target=target,
            pdb_id=pdb_id,
            chain=chain,
            evalue=float(result.get('evalue', result.get('eval', 0))),
            score=float(result.get('score', result.get('bits', 0))),
            seq_identity=float(result.get('seqId', result.get('pident', 0))) / 100.0,
            alignment_length=int(result.get('alnLength', result.get('alnlen', 0))),
            query_start=int(result.get('qStartPos', result.get('qstart', 0))),
            query_end=int(result.get('qEndPos', result.get('qend', 0))),
            target_start=int(result.get('dbStartPos', result.get('tstart', 0))),
            target_end=int(result.get('dbEndPos', result.get('tend', 0))),
        )


class FoldseekSearcher:
    """Interface to Foldseek for structure similarity search."""

    # Foldseek Web API endpoints
    WEB_API_BASE = "https://search.foldseek.com/api"

    # Foldseek output format columns for local search
    OUTPUT_COLUMNS = [
        'query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen',
        'qstart', 'qend', 'tstart', 'tend', 'evalue', 'score'
    ]

    def __init__(
        self,
        foldseek_path: str = "foldseek",
        database: str = "pdb100",
        temp_dir: Optional[str] = None,
        use_web_api: bool = True
    ):
        """
        Initialize the Foldseek searcher.

        Args:
            foldseek_path: Path to foldseek executable
            database: Database to search (pdb100, afdb, etc.)
            temp_dir: Temporary directory for intermediate files
            use_web_api: Use web API instead of local foldseek (recommended)
        """
        self.foldseek_path = foldseek_path
        self.database = database
        self.temp_dir = temp_dir or tempfile.gettempdir()
        self.use_web_api = use_web_api

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

        # Use web API (recommended - no local database needed)
        if self.use_web_api:
            return self.search_web(
                structure_path=str(structure_path),
                max_hits=max_hits,
                database=self.database,
                evalue=evalue
            )

        # Local search (requires database download)
        return self._search_local(
            structure_path=str(structure_path),
            max_hits=max_hits,
            evalue=evalue,
            min_seq_id=min_seq_id,
            coverage=coverage,
            threads=threads
        )

    def _search_local(
        self,
        structure_path: str,
        max_hits: int,
        evalue: float,
        min_seq_id: float,
        coverage: float,
        threads: int
    ) -> List[FoldseekHit]:
        """Run local Foldseek search."""
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.tsv', delete=False, dir=self.temp_dir
        ) as tmp_out:
            output_path = tmp_out.name

        try:
            cmd = [
                self.foldseek_path,
                "easy-search",
                structure_path,
                self.database,
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

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                logger.error(f"Foldseek stderr: {result.stderr}")
                raise RuntimeError(f"Foldseek search failed: {result.stderr}")

            hits = self._parse_results(output_path)
            logger.info(f"Found {len(hits)} hits")
            return hits

        finally:
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
        database: str = "pdb100",
        evalue: float = 1e-3,
        timeout: int = 300
    ) -> List[FoldseekHit]:
        """
        Search using Foldseek web server API.

        Args:
            structure_path: Path to the query structure
            max_hits: Maximum number of hits
            database: Database to search (pdb100, afdb50, afdb-swissprot, etc.)
            evalue: E-value threshold
            timeout: Maximum time to wait for results (seconds)

        Returns:
            List of FoldseekHit objects
        """
        logger.info(f"Searching Foldseek web API with database: {database}")

        # Read structure file
        with open(structure_path, 'r') as f:
            structure_content = f.read()

        # Map database names to Foldseek API database IDs
        db_mapping = {
            'pdb100': 'pdb100',
            'pdb': 'pdb100',
            'afdb': 'afdb50',
            'afdb50': 'afdb50',
            'afdb-swissprot': 'afdb-swissprot',
            'afdb-proteome': 'afdb-proteome',
        }
        api_database = db_mapping.get(database.lower(), 'pdb100')

        # Submit job to Foldseek web server
        submit_url = f"{self.WEB_API_BASE}/ticket"

        files = {
            'q': ('query.pdb', structure_content, 'application/octet-stream'),
        }
        data = {
            'database[]': api_database,
            'mode': '3diaa',  # 3Di+AA mode for better sensitivity
        }

        try:
            response = requests.post(
                submit_url,
                files=files,
                data=data,
                timeout=60,
                verify=False
            )
            response.raise_for_status()
            ticket = response.json()
            ticket_id = ticket['id']
            logger.info(f"Foldseek job submitted, ticket ID: {ticket_id}")

        except Exception as e:
            logger.error(f"Failed to submit Foldseek job: {e}")
            raise RuntimeError(f"Failed to submit Foldseek search: {e}")

        # Poll for results
        status_url = f"{self.WEB_API_BASE}/ticket/{ticket_id}"
        start_time = time.time()

        while True:
            if time.time() - start_time > timeout:
                raise RuntimeError(f"Foldseek search timed out after {timeout}s")

            try:
                response = requests.get(status_url, timeout=30, verify=False)
                result = response.json()
                status = result.get('status', '')

                if status == 'COMPLETE':
                    logger.info("Foldseek search completed")
                    break
                elif status == 'ERROR':
                    error_msg = result.get('error', 'Unknown error')
                    raise RuntimeError(f"Foldseek search failed: {error_msg}")
                elif status in ['PENDING', 'RUNNING']:
                    logger.debug(f"Foldseek status: {status}")
                    time.sleep(2)
                else:
                    logger.debug(f"Foldseek status: {status}")
                    time.sleep(2)

            except requests.exceptions.RequestException as e:
                logger.warning(f"Error polling status: {e}")
                time.sleep(2)

        # Fetch and parse results
        try:
            result_url = f"{self.WEB_API_BASE}/result/{ticket_id}/0"
            response = requests.get(result_url, timeout=60, verify=False)
            response.raise_for_status()
            results_data = response.json()

        except Exception as e:
            logger.error(f"Failed to fetch results: {e}")
            raise RuntimeError(f"Failed to fetch Foldseek results: {e}")

        # Parse results
        hits = []

        # Log the response structure for debugging
        logger.info(f"Response keys: {list(results_data.keys())}")

        # Handle different API response formats
        # New format: results -> [{'db': ..., 'alignments': [[hit1, hit2, ...]], ...}]
        # Old format: alignments -> [hit1, hit2, ...]
        results = results_data.get('results', [])
        if results:
            logger.info(f"Found {len(results)} result sets")
            # Extract alignments from each result set
            all_alignments = []
            for result_set in results:
                if isinstance(result_set, dict):
                    db_alignments = result_set.get('alignments', [])
                    if db_alignments:
                        all_alignments.extend(db_alignments)
            alignments = all_alignments
            logger.info(f"Total alignment groups: {len(alignments)}")
        else:
            # Fallback to old format
            alignments = results_data.get('alignments', [])

        for alignment in alignments:
            # Log first alignment structure for debugging
            if alignments.index(alignment) == 0:
                if isinstance(alignment, list) and len(alignment) > 0:
                    logger.info(f"First hit data keys: {list(alignment[0].keys()) if isinstance(alignment[0], dict) else type(alignment[0])}")
                    logger.info(f"First hit data sample: {str(alignment[0])[:500]}")
                elif isinstance(alignment, dict):
                    logger.info(f"First alignment keys: {list(alignment.keys())}")
                    logger.info(f"First alignment data sample: {str(alignment)[:500]}")

            # Each alignment can have multiple hits
            if isinstance(alignment, list):
                for hit_data in alignment[:max_hits]:
                    try:
                        hit = FoldseekHit.from_web_result(hit_data)
                        if hit.evalue <= evalue:
                            hits.append(hit)
                    except Exception as e:
                        logger.debug(f"Failed to parse hit: {e}")
                        continue
            elif isinstance(alignment, dict):
                try:
                    hit = FoldseekHit.from_web_result(alignment)
                    if hit.evalue <= evalue:
                        hits.append(hit)
                except Exception as e:
                    logger.debug(f"Failed to parse hit: {e}")
                    continue

        # Limit results
        hits = hits[:max_hits]

        # Remove duplicates by PDB ID
        seen = set()
        unique_hits = []
        for hit in hits:
            key = f"{hit.pdb_id}_{hit.chain}"
            if key not in seen:
                seen.add(key)
                unique_hits.append(hit)

        logger.info(f"Found {len(unique_hits)} unique hits from Foldseek web API")
        return unique_hits


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
