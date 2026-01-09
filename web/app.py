"""
Flask web application for Foldseek Pocket Miner.
"""

import os
import json
import uuid
import logging
import threading
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional, List, Tuple, Set
from collections import Counter, defaultdict
import numpy as np

from flask import Flask, render_template, request, jsonify, send_file, send_from_directory
from flask_cors import CORS
import requests

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent))

from foldseek_pocket_miner.structure_downloader import StructureDownloader
from foldseek_pocket_miner.foldseek_search import FoldseekSearcher
from foldseek_pocket_miner.ligand_extractor import LigandExtractor, EXCLUDE_RESIDUES

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Try to import Kamaji for ligand classification
KAMAJI_AVAILABLE = False
try:
    from openbabel import openbabel as ob
    from openbabel import pybel
    from kamaji import Kamaji
    from kamaji.pdbparser import MultiStatePDB
    KAMAJI_AVAILABLE = True
    logger.info("Kamaji loaded successfully for ligand classification")
except ImportError as e:
    logger.warning(f"Kamaji not available (missing dependencies: {e}). Using fallback classification.")

# Initialize Flask app
app = Flask(__name__)
CORS(app)

# Configuration
UPLOAD_FOLDER = Path(__file__).parent / "uploads"
RESULTS_FOLDER = Path(__file__).parent / "results"
UPLOAD_FOLDER.mkdir(exist_ok=True)
RESULTS_FOLDER.mkdir(exist_ok=True)

app.config['UPLOAD_FOLDER'] = str(UPLOAD_FOLDER)
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50MB max

# Job storage
jobs: Dict[str, Dict[str, Any]] = {}

# Amino acid codes
AA_CODES = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# Ligand category definitions - aligned with Kamaji types
# Kamaji types: protein, dna, rna, cofactor, glyco, ligand, water, metal, metal_catalytic, salt, modifier, additive
LIGAND_CATEGORIES = {
    'ligand': {
        'name': 'Small Molecules',
        'description': 'Drug-like compounds and generic ligands',
        'kamaji_type': 'ligand',
    },
    'cofactor': {
        'name': 'Cofactors',
        'description': 'ATP, NAD, FAD, Heme and other cofactors',
        'kamaji_type': 'cofactor',
    },
    'glyco': {
        'name': 'Carbohydrates',
        'description': 'Sugars and glycosylation groups',
        'kamaji_type': 'glyco',
    },
    'metal_catalytic': {
        'name': 'Catalytic Metals',
        'description': 'Mg, Mn, Fe, Co, Ni, Cu, Zn',
        'kamaji_type': 'metal_catalytic',
    },
    'metal': {
        'name': 'Metals',
        'description': 'Non-catalytic metal ions',
        'kamaji_type': 'metal',
    },
    'modifier': {
        'name': 'Modifiers',
        'description': 'Post-translational modifications',
        'kamaji_type': 'modifier',
    },
    'additive': {
        'name': 'Additives',
        'description': 'Crystallization additives (GOL, PEG, etc.)',
        'kamaji_type': 'additive',
    },
    'salt': {
        'name': 'Salts',
        'description': 'Na, Cl, Ca, K ions',
        'kamaji_type': 'salt',
    },
    # Fallback categories for when Kamaji is not available
    'small_molecule': {
        'name': 'Small Molecules',
        'description': 'Drug-like compounds (fallback)',
    },
    'nucleotide': {
        'name': 'Nucleotides',
        'description': 'ATP, GTP, NAD, FAD (fallback)',
    },
    'saccharide': {
        'name': 'Saccharides',
        'description': 'Sugars and carbohydrates (fallback)',
    },
    'metal_complex': {
        'name': 'Metal Complexes',
        'description': 'Metal-containing cofactors (fallback)',
    },
}

# Known nucleotide-like ligands
NUCLEOTIDE_LIGANDS = {
    'ATP', 'ADP', 'AMP', 'ANP', 'ACP', 'AGS',
    'GTP', 'GDP', 'GMP', 'GNP', 'GSP',
    'CTP', 'CDP', 'CMP',
    'UTP', 'UDP', 'UMP',
    'NAD', 'NAP', 'NDP', 'NAI',
    'FAD', 'FMN',
    'SAM', 'SAH',
    'COA', 'ACO',
}

# Known saccharide ligands
SACCHARIDE_LIGANDS = {
    'GAL', 'GLC', 'MAN', 'FUC', 'SIA', 'NAG', 'NDG', 'BMA', 'BGC',
    'GLA', 'A2G', 'RAM', 'XYS', 'RIB', 'FRU',
}

# Known metal-containing ligands
METAL_LIGANDS = {
    'HEM', 'HEC', 'HEA', 'HEB', 'HEG',
    'CLF', 'CLA', 'BCL', 'BCB', 'BPH', 'BPB',
    'FES', 'F3S', 'SF4', 'FEO',
}


def categorize_ligand(resname: str) -> str:
    """
    Categorize a ligand based on its residue name.

    Returns category key from LIGAND_CATEGORIES.
    """
    resname = resname.upper().strip()

    if resname in NUCLEOTIDE_LIGANDS:
        return 'nucleotide'
    elif resname in SACCHARIDE_LIGANDS:
        return 'saccharide'
    elif resname in METAL_LIGANDS:
        return 'metal_complex'
    else:
        return 'small_molecule'


def cluster_ligands_by_similarity(ligands: List[Dict]) -> Dict[str, List[Dict]]:
    """
    Cluster ligands by chemical category (fallback without Kamaji).

    Groups ligands into categories based on their chemical nature.
    """
    clusters = defaultdict(list)

    for lig in ligands:
        category = categorize_ligand(lig.get('name', ''))
        lig['category'] = category
        lig['category_name'] = LIGAND_CATEGORIES.get(category, {}).get('name', 'Other')
        clusters[category].append(lig)

    return dict(clusters)


def classify_pdb_with_kamaji(pdb_content: str) -> Dict[str, List[Dict]]:
    """
    Use Kamaji to classify all residues in a PDB structure.

    Returns a dictionary mapping category names to lists of residue info.
    Categories: ligand, cofactor, glyco, metal_catalytic, metal, modifier, additive, salt
    """
    if not KAMAJI_AVAILABLE:
        logger.debug("Kamaji not available, skipping classification")
        return {}

    try:
        # Write PDB to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_path = f.name

        try:
            # Parse PDB with Kamaji's MultiStatePDB
            with open(temp_path, 'r') as f:
                rawpdb = f.readlines()

            mmol = MultiStatePDB(rawpdb)
            if mmol.model_count == 0:
                return {}

            mmol.set_state(model=0, alt_mode='a')
            mol = mmol.get_structure()
            pdb_info = mmol.pdb_info

            # Run Kamaji classification
            kamaji = Kamaji()
            kamaji.set_molecule(mol, pdb_info=pdb_info)
            profile = kamaji.profile

            # Extract classified residues
            results = {}
            relevant_types = ['ligand', 'cofactor', 'glyco', 'metal_catalytic',
                             'metal', 'modifier', 'additive', 'salt']

            types_found = profile.get_types() or []

            for res_type in relevant_types:
                if res_type not in types_found:
                    continue

                residues = profile.get_type(res_type)
                if not residues:
                    continue

                results[res_type] = []
                for res_id, res_data in residues.items():
                    res_info = {
                        'name': res_data.get('name', ''),
                        'chain': res_data.get('chain', ''),
                        'num': res_data.get('num', 0),
                        'type': res_type,
                        'category_name': LIGAND_CATEGORIES.get(res_type, {}).get('name', res_type),
                    }
                    # Add extra info from Kamaji
                    if 'info' in res_data:
                        info = res_data['info']
                        if isinstance(info, dict):
                            res_info['kamaji_info'] = {
                                'class': info.get('class', ''),
                                'size': info.get('size', ''),
                                'mw': info.get('mw', 0),
                            }
                    results[res_type].append(res_info)

            return results

        finally:
            os.unlink(temp_path)

    except Exception as e:
        logger.warning(f"Kamaji classification failed: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return {}


def cluster_ligands_with_kamaji(ligands: List[Dict], pdb_content: str) -> Dict[str, List[Dict]]:
    """
    Cluster ligands using Kamaji classification.

    Falls back to simple categorization if Kamaji is not available.
    """
    if not KAMAJI_AVAILABLE:
        return cluster_ligands_by_similarity(ligands)

    try:
        # Get Kamaji classification for the PDB
        kamaji_results = classify_pdb_with_kamaji(pdb_content)

        if not kamaji_results:
            # Fallback to simple categorization
            return cluster_ligands_by_similarity(ligands)

        # Build a lookup from residue name to Kamaji category
        resname_to_category = {}
        for category, residues in kamaji_results.items():
            for res in residues:
                key = res.get('name', '')
                if key:
                    resname_to_category[key] = category

        # Categorize ligands based on Kamaji results
        clusters = defaultdict(list)
        for lig in ligands:
            resname = lig.get('name', '')
            # Look up category from Kamaji results
            category = resname_to_category.get(resname, 'ligand')
            lig['category'] = category
            lig['category_name'] = LIGAND_CATEGORIES.get(category, {}).get('name', 'Other')
            clusters[category].append(lig)

        return dict(clusters)

    except Exception as e:
        logger.warning(f"Kamaji clustering failed: {e}, using fallback")
        return cluster_ligands_by_similarity(ligands)


def fetch_taxonomy_for_pdb(pdb_id: str) -> Dict[str, Any]:
    """
    Fetch taxonomy/species information from RCSB for a PDB entry.

    Returns dict with taxonomy_id, scientific_name, common_name, lineage.
    """
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
        response = requests.get(url, timeout=10, verify=False)
        response.raise_for_status()
        data = response.json()

        # Extract source organism info
        polymer_entities = data.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids', [])

        # Get taxonomy from first polymer entity
        if polymer_entities:
            entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.lower()}/{polymer_entities[0]}"
            entity_response = requests.get(entity_url, timeout=10, verify=False)
            if entity_response.ok:
                entity_data = entity_response.json()
                source_organisms = entity_data.get('rcsb_entity_source_organism', [])

                if source_organisms:
                    org = source_organisms[0]
                    return {
                        'taxonomy_id': org.get('ncbi_taxonomy_id'),
                        'scientific_name': org.get('ncbi_scientific_name', ''),
                        'common_name': org.get('common_name', ''),
                        'source_type': org.get('source_type', ''),
                    }

        return {}
    except Exception as e:
        logger.debug(f"Failed to fetch taxonomy for {pdb_id}: {e}")
        return {}


def fetch_taxonomy_batch(pdb_ids: List[str]) -> Dict[str, Dict]:
    """
    Fetch taxonomy information for multiple PDB entries using GraphQL.

    Returns dict mapping pdb_id to taxonomy info.
    """
    if not pdb_ids:
        return {}

    # Use GraphQL for batch query
    graphql_url = "https://data.rcsb.org/graphql"

    # Build query for multiple entries
    query = """
    query getTaxonomy($ids: [String!]!) {
        entries(entry_ids: $ids) {
            rcsb_id
            polymer_entities {
                rcsb_entity_source_organism {
                    ncbi_taxonomy_id
                    ncbi_scientific_name
                    common_name
                    source_type
                }
            }
        }
    }
    """

    results = {}

    try:
        # Query in batches of 50
        batch_size = 50
        for i in range(0, len(pdb_ids), batch_size):
            batch = pdb_ids[i:i + batch_size]

            response = requests.post(
                graphql_url,
                json={'query': query, 'variables': {'ids': batch}},
                timeout=30,
                verify=False
            )

            if response.ok:
                data = response.json()
                entries = data.get('data', {}).get('entries', [])

                for entry in entries:
                    if not entry:
                        continue
                    pdb_id = entry.get('rcsb_id', '').upper()

                    # Get first source organism from first polymer entity
                    polymer_entities = entry.get('polymer_entities', [])
                    if polymer_entities:
                        for pe in polymer_entities:
                            organisms = pe.get('rcsb_entity_source_organism', [])
                            if organisms:
                                org = organisms[0]
                                results[pdb_id] = {
                                    'taxonomy_id': org.get('ncbi_taxonomy_id'),
                                    'scientific_name': org.get('ncbi_scientific_name', ''),
                                    'common_name': org.get('common_name', ''),
                                    'source_type': org.get('source_type', ''),
                                }
                                break

    except Exception as e:
        logger.warning(f"GraphQL taxonomy fetch failed: {e}")
        # Fallback to individual queries
        for pdb_id in pdb_ids[:20]:  # Limit fallback
            tax = fetch_taxonomy_for_pdb(pdb_id)
            if tax:
                results[pdb_id.upper()] = tax

    return results


def group_hits_by_species(hits: List[Dict], taxonomy_data: Dict[str, Dict]) -> Dict[str, List[Dict]]:
    """
    Group foldseek hits by species/organism.

    Returns dict mapping species name to list of hits.
    """
    groups = defaultdict(list)

    for hit in hits:
        pdb_id = hit.get('pdb_id', '').upper()
        tax = taxonomy_data.get(pdb_id, {})
        species = tax.get('scientific_name') or 'Unknown'
        hit['species'] = species
        hit['taxonomy_id'] = tax.get('taxonomy_id')
        hit['common_name'] = tax.get('common_name') or ''
        groups[species].append(hit)

    return dict(groups)


def calculate_conservation_by_species(
    query_seq: Dict[int, str],
    aligned_seqs: List[Dict[int, str]],
    species_filter: Optional[Set[str]] = None,
    exclude_species: Optional[Set[str]] = None
) -> Dict[int, float]:
    """
    Calculate conservation scores with species filtering.

    Args:
        query_seq: Query sequence as dict of resseq -> aa
        aligned_seqs: List of aligned sequences with metadata
        species_filter: Only include these species (if provided)
        exclude_species: Exclude these species (e.g., {'Homo sapiens'})

    Returns dict mapping residue number to conservation score (0-1).
    """
    conservation = {}

    for resseq, query_aa in query_seq.items():
        aas = [query_aa]

        for seq_data in aligned_seqs:
            # Check species filter
            if isinstance(seq_data, dict):
                species = seq_data.get('species', '')
                seq = seq_data.get('sequence', {})

                if species_filter and species not in species_filter:
                    continue
                if exclude_species and species in exclude_species:
                    continue

                if resseq in seq:
                    aas.append(seq[resseq])
            else:
                # Legacy format: just sequence dict
                if resseq in seq_data:
                    aas.append(seq_data[resseq])

        if len(aas) <= 1:
            conservation[resseq] = 0.5
            continue

        counts = Counter(aas)
        total = len(aas)

        entropy = 0.0
        for count in counts.values():
            freq = count / total
            if freq > 0:
                entropy -= freq * np.log2(freq)

        max_entropy = np.log2(min(20, total))
        if max_entropy > 0:
            normalized = entropy / max_entropy
        else:
            normalized = 0

        conservation[resseq] = 1.0 - normalized

    return conservation


class JobStatus:
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


def extract_sequence_from_pdb(pdb_content: str, chain: str = None) -> Dict[int, str]:
    """Extract amino acid sequence from PDB content."""
    sequence = {}
    seen = set()

    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            res_chain = line[21:22].strip()
            if chain and res_chain != chain:
                continue

            resseq = int(line[22:26].strip())
            resname = line[17:20].strip()

            if (res_chain, resseq) not in seen:
                seen.add((res_chain, resseq))
                aa = AA_CODES.get(resname, 'X')
                sequence[resseq] = aa

    return sequence


def calculate_conservation(query_seq: Dict[int, str], aligned_seqs: List[Dict[int, str]]) -> Dict[int, float]:
    """
    Calculate conservation score for each residue.

    Returns dict mapping residue number to conservation score (0-1).
    Higher score = more conserved.
    """
    conservation = {}

    for resseq, query_aa in query_seq.items():
        # Collect amino acids at this position
        aas = [query_aa]
        for seq in aligned_seqs:
            if resseq in seq:
                aas.append(seq[resseq])

        # Calculate Shannon entropy-based conservation
        if len(aas) <= 1:
            conservation[resseq] = 0.5  # Unknown
            continue

        # Count frequencies
        counts = Counter(aas)
        total = len(aas)

        # Shannon entropy
        entropy = 0.0
        for count in counts.values():
            freq = count / total
            if freq > 0:
                entropy -= freq * np.log2(freq)

        # Normalize (max entropy = log2(20) for 20 amino acids)
        max_entropy = np.log2(min(20, total))
        if max_entropy > 0:
            normalized = entropy / max_entropy
        else:
            normalized = 0

        # Convert to conservation (1 - entropy)
        conservation[resseq] = 1.0 - normalized

    return conservation


def extract_ligands_from_pdb(pdb_content: str, source_name: str = "") -> List[Dict]:
    """Extract ligand atoms from PDB content."""
    ligands = {}

    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            resname = line[17:20].strip()
            if resname in EXCLUDE_RESIDUES:
                continue

            chain = line[21:22].strip()
            resseq = int(line[22:26].strip())
            key = (resname, chain, resseq)

            if key not in ligands:
                ligands[key] = {
                    'name': resname,
                    'chain': chain,
                    'resseq': resseq,
                    'source': source_name,
                    'atoms': []
                }

            # Extract atom info
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            ligands[key]['atoms'].append({'x': x, 'y': y, 'z': z, 'line': line})

    # Filter ligands with enough atoms
    result = []
    for lig in ligands.values():
        if len(lig['atoms']) >= 5:
            result.append(lig)

    return result


def align_structure_pymol(query_pdb: str, target_pdb: str, pymol_path: str = "pymol") -> Tuple[str, float]:
    """
    Align target structure to query using PyMOL.
    Returns aligned PDB content and RMSD.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        query_path = os.path.join(tmpdir, "query.pdb")
        target_path = os.path.join(tmpdir, "target.pdb")
        aligned_path = os.path.join(tmpdir, "aligned.pdb")

        with open(query_path, 'w') as f:
            f.write(query_pdb)
        with open(target_path, 'w') as f:
            f.write(target_pdb)

        script = f'''
from pymol import cmd
cmd.load("{query_path}", "query")
cmd.load("{target_path}", "target")
result = cmd.align("target", "query")
cmd.save("{aligned_path}", "target")
print(f"RMSD: {{result[0]}}")
cmd.quit()
'''
        script_path = os.path.join(tmpdir, "align.py")
        with open(script_path, 'w') as f:
            f.write(script)

        try:
            result = subprocess.run(
                [pymol_path, "-cq", script_path],
                capture_output=True,
                text=True,
                timeout=60
            )

            rmsd = 0.0
            for line in result.stdout.split('\n'):
                if 'RMSD:' in line:
                    rmsd = float(line.split(':')[1].strip())

            if os.path.exists(aligned_path):
                with open(aligned_path, 'r') as f:
                    return f.read(), rmsd
        except Exception as e:
            logger.warning(f"PyMOL alignment failed: {e}")

        return target_pdb, 0.0


def run_pipeline(job_id: str, pdb_id: str = None, chain: str = "A",
                structure_path: str = None, max_hits: int = 50):
    """Run the pocket mining pipeline with conservation and pharmacophore analysis."""
    try:
        jobs[job_id]["status"] = JobStatus.RUNNING
        jobs[job_id]["progress"] = 0
        jobs[job_id]["current_step"] = "Initializing..."

        output_dir = RESULTS_FOLDER / job_id
        output_dir.mkdir(exist_ok=True)

        # Step 1: Download query structure
        jobs[job_id]["current_step"] = "Downloading query structure..."
        jobs[job_id]["progress"] = 5

        if pdb_id and not structure_path:
            downloader = StructureDownloader(output_dir=str(output_dir / "query"))
            structure_path = downloader.download_structure(pdb_id, chain)
            if not structure_path:
                raise Exception(f"Failed to download {pdb_id}")

        with open(structure_path, 'r') as f:
            query_pdb = f.read()

        jobs[job_id]["query_structure"] = structure_path
        jobs[job_id]["query_pdb"] = query_pdb
        jobs[job_id]["progress"] = 10

        # Extract query sequence
        query_seq = extract_sequence_from_pdb(query_pdb, chain)
        logger.info(f"Query has {len(query_seq)} residues")

        # Step 2: Foldseek search
        jobs[job_id]["current_step"] = "Searching similar structures (Foldseek)..."
        jobs[job_id]["progress"] = 15

        searcher = FoldseekSearcher(database="pdb100", use_web_api=True)
        hits = searcher.search(
            structure_path=structure_path,
            max_hits=max_hits,
            evalue=1e-3
        )

        # Step 2.5: Fetch taxonomy for all hits
        jobs[job_id]["current_step"] = "Fetching taxonomy information..."
        pdb_ids = [h.pdb_id for h in hits]
        taxonomy_data = fetch_taxonomy_batch(pdb_ids)
        logger.info(f"Fetched taxonomy for {len(taxonomy_data)} structures")

        hits_list = [
            {"pdb_id": h.pdb_id, "chain": h.chain, "evalue": h.evalue,
             "score": h.score, "seq_identity": h.seq_identity}
            for h in hits
        ]

        # Add taxonomy info to hits
        for hit in hits_list:
            pdb_id = hit['pdb_id'].upper()
            tax = taxonomy_data.get(pdb_id, {})
            hit['species'] = tax.get('scientific_name') or 'Unknown'
            hit['taxonomy_id'] = tax.get('taxonomy_id')
            hit['common_name'] = tax.get('common_name') or ''

        jobs[job_id]["hits"] = hits_list
        jobs[job_id]["taxonomy_data"] = taxonomy_data

        # Group hits by species
        species_groups = group_hits_by_species(hits_list.copy(), taxonomy_data)
        jobs[job_id]["species_groups"] = {k: len(v) for k, v in species_groups.items()}
        jobs[job_id]["progress"] = 32

        if not hits:
            jobs[job_id]["current_step"] = "No similar structures found"
            jobs[job_id]["status"] = JobStatus.COMPLETED
            return

        # Step 3: Download structures
        jobs[job_id]["current_step"] = f"Downloading {len(hits)} structures..."
        jobs[job_id]["progress"] = 35

        struct_downloader = StructureDownloader(output_dir=str(output_dir / "structures"))
        pdb_ids = [h.pdb_id for h in hits]
        chains_list = [h.chain for h in hits]
        downloaded = struct_downloader.download_batch(pdb_ids, chains_list, show_progress=False)

        valid_paths = {k: v for k, v in downloaded.items() if v}
        jobs[job_id]["downloaded_count"] = len(valid_paths)
        logger.info(f"Downloaded {len(valid_paths)} structures")
        jobs[job_id]["progress"] = 45

        # Step 4: Align structures and calculate conservation
        jobs[job_id]["current_step"] = "Aligning structures & calculating conservation..."
        jobs[job_id]["progress"] = 50

        aligned_seqs = []
        aligned_seqs_with_species = []  # For species-filtered conservation
        aligned_pdbs = []
        all_ligands = []

        # Check if PyMOL is available
        pymol_available = False
        try:
            result = subprocess.run(["pymol", "-c", "-Q"], capture_output=True, timeout=5)
            pymol_available = True
        except:
            # Try common paths
            for path in ["/opt/homebrew/bin/pymol", "/usr/local/bin/pymol",
                        "/Applications/PyMOL.app/Contents/MacOS/PyMOL"]:
                if os.path.exists(path):
                    pymol_available = True
                    break

        processed = 0
        for key, path in list(valid_paths.items())[:30]:  # Limit to 30 structures
            if not path:
                continue

            try:
                with open(path, 'r') as f:
                    target_pdb = f.read()

                # Align if PyMOL available
                if pymol_available:
                    aligned_pdb, rmsd = align_structure_pymol(query_pdb, target_pdb)
                else:
                    aligned_pdb = target_pdb  # Use unaligned

                # Extract PDB ID from key
                source_pdb_id = key.split('_')[0]

                # Extract sequence for conservation
                target_chain = key.split('_')[1] if '_' in key else None
                seq = extract_sequence_from_pdb(aligned_pdb, target_chain)
                if seq:
                    aligned_seqs.append(seq)
                    # Store with species info for filtered analysis
                    species = taxonomy_data.get(source_pdb_id.upper(), {}).get('scientific_name', 'Unknown')
                    aligned_seqs_with_species.append({
                        'sequence': seq,
                        'species': species,
                        'pdb_id': source_pdb_id
                    })

                # Extract ligands
                ligands = extract_ligands_from_pdb(aligned_pdb, source_pdb_id)

                # Run Kamaji classification on this structure
                kamaji_categories = {}
                if KAMAJI_AVAILABLE:
                    try:
                        kamaji_results = classify_pdb_with_kamaji(aligned_pdb)
                        for category, residues in kamaji_results.items():
                            for res in residues:
                                resname = res.get('name', '')
                                if resname:
                                    kamaji_categories[resname] = {
                                        'category': category,
                                        'category_name': LIGAND_CATEGORIES.get(category, {}).get('name', category),
                                        'info': res.get('kamaji_info', {})
                                    }
                    except Exception as e:
                        logger.debug(f"Kamaji classification failed for {source_pdb_id}: {e}")

                for lig in ligands:
                    # Build ligand PDB string
                    lig_pdb = f"REMARK Ligand {lig['name']} from {lig['source']}\n"
                    for atom in lig['atoms']:
                        lig_pdb += atom['line'] + '\n'
                    lig_pdb += "END\n"

                    # Get Kamaji category if available
                    kamaji_cat = kamaji_categories.get(lig['name'], {})
                    category = kamaji_cat.get('category', categorize_ligand(lig['name']))
                    category_name = kamaji_cat.get('category_name', LIGAND_CATEGORIES.get(category, {}).get('name', 'Other'))

                    all_ligands.append({
                        'name': lig['name'],
                        'source': lig['source'],
                        'chain': lig['chain'],
                        'atoms': len(lig['atoms']),
                        'pdb': lig_pdb,
                        'category': category,
                        'category_name': category_name,
                        'kamaji_info': kamaji_cat.get('info', {})
                    })

                aligned_pdbs.append({
                    'pdb_id': source_pdb_id,
                    'pdb': aligned_pdb
                })

                processed += 1
                jobs[job_id]["progress"] = 50 + int((processed / min(30, len(valid_paths))) * 30)

            except Exception as e:
                logger.warning(f"Failed to process {key}: {e}")
                continue

        logger.info(f"Processed {processed} structures, found {len(all_ligands)} ligands")

        # Step 5: Calculate conservation
        jobs[job_id]["current_step"] = "Computing sequence conservation..."
        jobs[job_id]["progress"] = 82

        conservation = calculate_conservation(query_seq, aligned_seqs)

        # Convert to list for JSON
        conservation_data = [
            {"resseq": resseq, "score": score, "aa": query_seq.get(resseq, 'X')}
            for resseq, score in sorted(conservation.items())
        ]
        jobs[job_id]["conservation"] = conservation_data
        logger.info(f"Calculated conservation for {len(conservation)} residues")

        # Step 6: Prepare visualization data
        jobs[job_id]["current_step"] = "Preparing visualization..."
        jobs[job_id]["progress"] = 90

        # Deduplicate ligands by name
        seen_ligands = set()
        unique_ligands = []
        for lig in all_ligands:
            key = f"{lig['name']}_{lig['source']}"
            if key not in seen_ligands and len(unique_ligands) < 30:
                seen_ligands.add(key)
                unique_ligands.append(lig)

        # Categorize ligands (using Kamaji classification if available)
        ligands_info = [
            {
                "name": l["name"],
                "source": l["source"],
                "chain": l["chain"],
                "atoms": l["atoms"],
                "category": l.get("category", "ligand"),
                "category_name": l.get("category_name", "Small Molecules"),
                "kamaji_info": l.get("kamaji_info", {})
            }
            for l in unique_ligands
        ]

        # Cluster ligands by their pre-assigned Kamaji categories
        ligand_clusters = defaultdict(list)
        for lig in ligands_info:
            category = lig.get("category", "ligand")
            ligand_clusters[category].append(lig)

        jobs[job_id]["ligand_clusters"] = {
            cat: [{"name": lig["name"], "source": lig["source"], "category_name": lig.get("category_name", ""),
                   "kamaji_info": lig.get("kamaji_info", {})}
                  for lig in ligs]
            for cat, ligs in ligand_clusters.items()
        }
        jobs[job_id]["kamaji_available"] = KAMAJI_AVAILABLE

        jobs[job_id]["ligand_pdbs"] = unique_ligands
        jobs[job_id]["ligands"] = ligands_info
        jobs[job_id]["aligned_count"] = len(aligned_seqs)
        jobs[job_id]["aligned_seqs_with_species"] = aligned_seqs_with_species

        jobs[job_id]["progress"] = 100
        jobs[job_id]["current_step"] = "Complete!"
        jobs[job_id]["status"] = JobStatus.COMPLETED
        jobs[job_id]["completed_at"] = datetime.now().isoformat()

        logger.info(f"Job {job_id} completed: {len(hits)} hits, {len(unique_ligands)} ligands, {len(conservation)} conserved residues")

    except Exception as e:
        logger.error(f"Job {job_id} failed: {e}")
        import traceback
        jobs[job_id]["status"] = JobStatus.FAILED
        jobs[job_id]["error"] = str(e)
        jobs[job_id]["traceback"] = traceback.format_exc()


@app.route('/')
def index():
    """Render main page."""
    return render_template('index.html')


@app.route('/api/submit', methods=['POST'])
def submit_job():
    """Submit a new pocket mining job."""
    data = request.json

    pdb_id = data.get('pdb_id', '').strip().upper()
    chain = data.get('chain', 'A').strip().upper()
    max_hits = int(data.get('max_hits', 50))

    if not pdb_id:
        return jsonify({"error": "PDB ID is required"}), 400

    if len(pdb_id) != 4:
        return jsonify({"error": "Invalid PDB ID format"}), 400

    job_id = str(uuid.uuid4())[:8]
    jobs[job_id] = {
        "id": job_id,
        "pdb_id": pdb_id,
        "chain": chain,
        "max_hits": max_hits,
        "status": JobStatus.PENDING,
        "progress": 0,
        "current_step": "Queued",
        "created_at": datetime.now().isoformat(),
        "hits": [],
        "ligands": [],
        "conservation": [],
    }

    thread = threading.Thread(
        target=run_pipeline,
        args=(job_id, pdb_id, chain, None, max_hits)
    )
    thread.daemon = True
    thread.start()

    return jsonify({"job_id": job_id, "message": "Job submitted successfully"})


@app.route('/api/status/<job_id>')
def get_status(job_id):
    """Get job status."""
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404

    job = jobs[job_id]
    return jsonify({
        "id": job["id"],
        "status": job["status"],
        "progress": job.get("progress", 0),
        "current_step": job.get("current_step", ""),
        "error": job.get("error"),
    })


@app.route('/api/results/<job_id>')
def get_results(job_id):
    """Get job results."""
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404

    job = jobs[job_id]

    if job["status"] != JobStatus.COMPLETED:
        return jsonify({"error": "Job not completed yet"}), 400

    return jsonify({
        "id": job["id"],
        "pdb_id": job["pdb_id"],
        "chain": job["chain"],
        "status": job["status"],
        "hits": job.get("hits", []),
        "ligands": job.get("ligands", []),
        "ligand_pdbs": job.get("ligand_pdbs", []),
        "ligand_clusters": job.get("ligand_clusters", {}),
        "conservation": job.get("conservation", []),
        "downloaded_count": job.get("downloaded_count", 0),
        "aligned_count": job.get("aligned_count", 0),
        "query_pdb": job.get("query_pdb", ""),
        "species_groups": job.get("species_groups", {}),
        "taxonomy_data": job.get("taxonomy_data", {}),
        "kamaji_available": job.get("kamaji_available", False),
    })


@app.route('/api/conservation/<job_id>', methods=['POST'])
def get_filtered_conservation(job_id):
    """
    Get conservation analysis filtered by species.

    Request body:
    {
        "include_species": ["Cytomegalovirus", "HIV-1"],  // optional
        "exclude_species": ["Homo sapiens"],  // optional
    }
    """
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404

    job = jobs[job_id]
    data = request.json or {}

    include_species = set(data.get('include_species', []))
    exclude_species = set(data.get('exclude_species', []))

    # Get stored aligned sequences with species info
    aligned_seqs_with_species = job.get("aligned_seqs_with_species", [])

    if not aligned_seqs_with_species:
        # Fallback to regular conservation
        return jsonify({
            "conservation": job.get("conservation", []),
            "filtered": False
        })

    # Get query sequence
    query_pdb = job.get("query_pdb", "")
    chain = job.get("chain", "A")
    query_seq = extract_sequence_from_pdb(query_pdb, chain)

    # Calculate filtered conservation
    conservation = calculate_conservation_by_species(
        query_seq,
        aligned_seqs_with_species,
        species_filter=include_species if include_species else None,
        exclude_species=exclude_species if exclude_species else None
    )

    conservation_data = [
        {"resseq": resseq, "score": score, "aa": query_seq.get(resseq, 'X')}
        for resseq, score in sorted(conservation.items())
    ]

    return jsonify({
        "conservation": conservation_data,
        "filtered": True,
        "include_species": list(include_species),
        "exclude_species": list(exclude_species)
    })


@app.route('/api/pdb/<pdb_id>')
def fetch_pdb(pdb_id):
    """Fetch PDB structure from RCSB."""
    import urllib3
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        response = requests.get(url, timeout=30, verify=False)
        response.raise_for_status()
        return response.text, 200, {'Content-Type': 'text/plain'}
    except Exception as e:
        return jsonify({"error": str(e)}), 404


@app.route('/static/<path:filename>')
def serve_static(filename):
    """Serve static files."""
    return send_from_directory('static', filename)


def set_bfactors_for_conservation(pdb_content: str, conservation: Dict[int, float], chain: str = None) -> str:
    """
    Set B-factors in PDB to conservation scores for PyMOL visualization.
    Conservation scores (0-1) are scaled to 0-100 for B-factor column.
    """
    lines = []
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM'):
            res_chain = line[21:22].strip()
            if chain and res_chain and res_chain != chain:
                lines.append(line)
                continue

            try:
                resseq = int(line[22:26].strip())
                score = conservation.get(resseq, 0.5)
                # Scale to 0-100 for B-factor
                bfactor = score * 100.0
                # Replace B-factor column (60-66)
                new_line = line[:60] + f"{bfactor:6.2f}" + line[66:]
                lines.append(new_line)
            except:
                lines.append(line)
        else:
            lines.append(line)

    return '\n'.join(lines)


@app.route('/api/download_pse/<job_id>')
def download_pse(job_id):
    """Generate and download PyMOL session file."""
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404

    job = jobs[job_id]

    if job["status"] != JobStatus.COMPLETED:
        return jsonify({"error": "Job not completed yet"}), 400

    # Create temporary directory for PyMOL session generation
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_id = job["pdb_id"]
        chain = job.get("chain", "A")

        # Get conservation data as dict
        conservation = {c["resseq"]: c["score"] for c in job.get("conservation", [])}

        # Write query structure with B-factors set to conservation scores
        query_pdb = job.get("query_pdb", "")
        if conservation:
            query_pdb = set_bfactors_for_conservation(query_pdb, conservation, chain)

        query_path = os.path.join(tmpdir, f"{pdb_id}_query.pdb")
        with open(query_path, 'w') as f:
            f.write(query_pdb)

        # Write ligands
        ligand_paths = []
        for i, lig in enumerate(job.get("ligand_pdbs", [])):
            lig_path = os.path.join(tmpdir, f"ligand_{i}_{lig['name']}_{lig['source']}.pdb")
            with open(lig_path, 'w') as f:
                f.write(lig.get('pdb', ''))
            ligand_paths.append((lig_path, lig['name'], lig['source']))

        # Create PyMOL script
        pse_path = os.path.join(tmpdir, f"{pdb_id}_pocket_miner.pse")

        script = f'''
from pymol import cmd, util

# Load query structure
cmd.load("{query_path}", "query_{pdb_id}")

# Set up conservation coloring (B-factor spectrum: blue=low, white=mid, red=high)
cmd.show("cartoon", "query_{pdb_id}")
cmd.spectrum("b", "blue_white_red", "query_{pdb_id} and name CA", minimum=0, maximum=100)

# Color cartoon by B-factor (conservation)
cmd.set("cartoon_color", "spectrum", "query_{pdb_id}")
util.cbaw("query_{pdb_id}")  # Color by atom type for sticks
cmd.spectrum("b", "blue_white_red", "query_{pdb_id}", minimum=0, maximum=100)

# Create groups
cmd.group("protein", "query_{pdb_id}")
'''

        # Add ligands
        if ligand_paths:
            script += '\n# Load ligands\n'
            for lig_path, lig_name, lig_source in ligand_paths:
                obj_name = f"lig_{lig_name}_{lig_source}".replace('-', '_')
                script += f'cmd.load("{lig_path}", "{obj_name}")\n'
                script += f'cmd.show("sticks", "{obj_name}")\n'
                script += f'util.cbag("{obj_name}")  # Color by atom type\n'

            script += '\n# Group all ligands\n'
            ligand_names = [f"lig_{l[1]}_{l[2]}".replace('-', '_') for l in ligand_paths]
            script += f'cmd.group("ligands", "{" ".join(ligand_names)}")\n'

        # Finalize visualization
        script += f'''
# Set view and style
cmd.bg_color("white")
cmd.set("ray_shadow", 0)
cmd.set("depth_cue", 0)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_flat_sheets", 1)
cmd.set("stick_radius", 0.15)
cmd.set("sphere_scale", 0.2)

# Create surface (hidden by default)
cmd.create("protein_surface", "query_{pdb_id}")
cmd.show("surface", "protein_surface")
cmd.set("surface_color", "white", "protein_surface")
cmd.set("transparency", 0.7, "protein_surface")
cmd.disable("protein_surface")

# Center view
cmd.zoom("query_{pdb_id}")

# Save session
cmd.save("{pse_path}")
cmd.quit()
'''

        script_path = os.path.join(tmpdir, "create_session.py")
        with open(script_path, 'w') as f:
            f.write(script)

        # Run PyMOL to create session
        pymol_paths = ["pymol", "/opt/homebrew/bin/pymol", "/usr/local/bin/pymol",
                       "/Applications/PyMOL.app/Contents/MacOS/PyMOL"]

        pymol_cmd = None
        for path in pymol_paths:
            try:
                result = subprocess.run([path, "--version"], capture_output=True, timeout=5)
                pymol_cmd = path
                break
            except:
                continue

        if not pymol_cmd:
            return jsonify({"error": "PyMOL not available on server"}), 500

        try:
            result = subprocess.run(
                [pymol_cmd, "-cq", script_path],
                capture_output=True,
                text=True,
                timeout=120,
                cwd=tmpdir
            )

            if result.returncode != 0:
                logger.error(f"PyMOL error: {result.stderr}")
                return jsonify({"error": f"PyMOL failed: {result.stderr}"}), 500

            if not os.path.exists(pse_path):
                return jsonify({"error": "Failed to create PyMOL session"}), 500

            # Copy to results folder for serving
            output_path = RESULTS_FOLDER / job_id / f"{pdb_id}_pocket_miner.pse"
            import shutil
            shutil.copy(pse_path, output_path)

            return send_file(
                output_path,
                as_attachment=True,
                download_name=f"{pdb_id}_pocket_miner.pse",
                mimetype='application/octet-stream'
            )

        except subprocess.TimeoutExpired:
            return jsonify({"error": "PyMOL session creation timed out"}), 500
        except Exception as e:
            logger.error(f"Failed to create PyMOL session: {e}")
            return jsonify({"error": str(e)}), 500


if __name__ == '__main__':
    import os
    debug = os.environ.get('FLASK_DEBUG', 'true').lower() == 'true'
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=debug)
