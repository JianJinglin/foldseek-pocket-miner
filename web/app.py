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


# Cache for ligand metadata
_ligand_metadata_cache: Dict[str, Dict] = {}


def fetch_ligand_metadata(ligand_id: str) -> Dict[str, Any]:
    """
    Fetch ligand metadata from RCSB API.

    Returns molecular weight, formula, description etc.
    """
    ligand_id = ligand_id.upper().strip()

    # Check cache
    if ligand_id in _ligand_metadata_cache:
        return _ligand_metadata_cache[ligand_id]

    try:
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
        response = requests.get(url, timeout=5, verify=False)

        if response.ok:
            data = response.json()

            # Extract relevant info
            chem_comp = data.get('chem_comp', {})
            descriptor = data.get('rcsb_chem_comp_descriptor', {})

            metadata = {
                'id': ligand_id,
                'name': chem_comp.get('name', ''),
                'formula': chem_comp.get('formula', ''),
                'molecular_weight': chem_comp.get('formula_weight', 0),
                'type': chem_comp.get('type', ''),
                'description': chem_comp.get('name', ''),  # Full chemical name
                'smiles': descriptor.get('SMILES', '') if descriptor else '',
            }

            # Cache the result
            _ligand_metadata_cache[ligand_id] = metadata
            return metadata

    except Exception as e:
        logger.debug(f"Failed to fetch ligand metadata for {ligand_id}: {e}")

    # Return empty dict on failure
    return {'id': ligand_id}


def fetch_ligand_metadata_batch(ligand_ids: List[str]) -> Dict[str, Dict]:
    """
    Fetch metadata for multiple ligands.

    Returns dict mapping ligand ID to metadata.
    """
    results = {}
    uncached = []

    # Check cache first
    for lid in ligand_ids:
        lid_upper = lid.upper().strip()
        if lid_upper in _ligand_metadata_cache:
            results[lid_upper] = _ligand_metadata_cache[lid_upper]
        else:
            uncached.append(lid_upper)

    # Fetch uncached (up to 10 to avoid too many requests)
    for lid in uncached[:10]:
        metadata = fetch_ligand_metadata(lid)
        results[lid] = metadata

    return results


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


def classify_organism_type(scientific_name: str, taxonomy_id: Optional[int], lineage: List[str] = None) -> str:
    """
    Classify organism into high-level categories: viral, mammal, human, bacteria, other.
    """
    if not scientific_name:
        return 'unknown'

    name_lower = scientific_name.lower()

    # Check for human
    if 'homo sapiens' in name_lower or taxonomy_id == 9606:
        return 'human'

    # Check for viruses (common patterns)
    viral_keywords = ['virus', 'viral', 'phage', 'viridae', 'virales',
                      'herpes', 'hiv', 'cmv', 'cytomegalovirus', 'influenza',
                      'coronavirus', 'sars', 'hepatitis', 'retrovirus']
    if any(kw in name_lower for kw in viral_keywords):
        return 'viral'

    # Check for mammals
    mammal_keywords = ['mus musculus', 'rattus', 'mouse', 'rat', 'bovine',
                       'bos taurus', 'sus scrofa', 'pig', 'rabbit', 'oryctolagus',
                       'canis', 'dog', 'felis', 'cat', 'primate', 'macaca',
                       'monkey', 'chimpanzee', 'gorilla']
    if any(kw in name_lower for kw in mammal_keywords):
        return 'mammal'

    # Check for bacteria
    bacteria_keywords = ['escherichia', 'bacillus', 'streptococcus', 'staphylococcus',
                         'salmonella', 'pseudomonas', 'mycobacterium', 'clostridium']
    if any(kw in name_lower for kw in bacteria_keywords):
        return 'bacteria'

    # Check for fungi/yeast
    fungi_keywords = ['saccharomyces', 'candida', 'aspergillus', 'yeast', 'fungi']
    if any(kw in name_lower for kw in fungi_keywords):
        return 'fungi'

    return 'other'


def fetch_taxonomy_batch(pdb_ids: List[str]) -> Dict[str, Dict]:
    """
    Fetch taxonomy and citation information for multiple PDB entries using GraphQL.

    Returns dict mapping pdb_id to taxonomy info including:
    - taxonomy_id, scientific_name, common_name, source_type
    - organism_type (human, mammal, viral, bacteria, fungi, other)
    - citation info (doi, pubmed_id, title, year, authors)
    """
    if not pdb_ids:
        return {}

    # Use GraphQL for batch query
    graphql_url = "https://data.rcsb.org/graphql"

    # Build query for multiple entries - include citation info
    query = """
    query getTaxonomyAndCitation($ids: [String!]!) {
        entries(entry_ids: $ids) {
            rcsb_id
            struct {
                title
            }
            rcsb_primary_citation {
                pdbx_database_id_DOI
                pdbx_database_id_PubMed
                title
                year
                rcsb_authors
                journal_abbrev
            }
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

                    result = {
                        'taxonomy_id': None,
                        'scientific_name': '',
                        'common_name': '',
                        'source_type': '',
                        'organism_type': 'unknown',
                        'structure_title': '',
                        'citation': None,
                    }

                    # Get structure title
                    struct = entry.get('struct', {})
                    if struct:
                        result['structure_title'] = struct.get('title', '')

                    # Get citation info
                    citation = entry.get('rcsb_primary_citation', {})
                    if citation:
                        doi = citation.get('pdbx_database_id_DOI', '')
                        pubmed_id = citation.get('pdbx_database_id_PubMed', '')
                        result['citation'] = {
                            'doi': doi,
                            'pubmed_id': pubmed_id,
                            'title': citation.get('title', ''),
                            'year': citation.get('year'),
                            'authors': citation.get('rcsb_authors', []),
                            'journal': citation.get('journal_abbrev', ''),
                            'doi_url': f"https://doi.org/{doi}" if doi else '',
                            'pubmed_url': f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/" if pubmed_id else '',
                        }

                    # Get first source organism from first polymer entity
                    polymer_entities = entry.get('polymer_entities', [])
                    if polymer_entities:
                        for pe in polymer_entities:
                            organisms = pe.get('rcsb_entity_source_organism', [])
                            if organisms:
                                org = organisms[0]
                                result['taxonomy_id'] = org.get('ncbi_taxonomy_id')
                                result['scientific_name'] = org.get('ncbi_scientific_name', '')
                                result['common_name'] = org.get('common_name', '')
                                result['source_type'] = org.get('source_type', '')
                                # Classify organism type
                                result['organism_type'] = classify_organism_type(
                                    result['scientific_name'],
                                    result['taxonomy_id']
                                )
                                break

                    results[pdb_id] = result

    except Exception as e:
        logger.warning(f"GraphQL taxonomy fetch failed: {e}")
        # Fallback to individual queries
        for pdb_id in pdb_ids[:20]:  # Limit fallback
            tax = fetch_taxonomy_for_pdb(pdb_id)
            if tax:
                tax['organism_type'] = classify_organism_type(
                    tax.get('scientific_name', ''),
                    tax.get('taxonomy_id')
                )
                tax['citation'] = None
                results[pdb_id.upper()] = tax

    return results


def fetch_citation_chain(doi: str = None, pubmed_id: str = None, max_depth: int = 3) -> Dict:
    """
    Fetch citation chain using Semantic Scholar API.
    Returns papers that cite this paper (citing) and papers this paper cites (cited).

    Args:
        doi: DOI of the paper
        pubmed_id: PubMed ID of the paper
        max_depth: Maximum depth of citation chain to fetch

    Returns:
        Dict with 'paper', 'citing' (papers that cite this), 'cited' (papers this cites)
    """
    if not doi and not pubmed_id:
        return {'error': 'No DOI or PubMed ID provided'}

    # Build query identifier for Semantic Scholar
    paper_id = None
    if doi:
        paper_id = f"DOI:{doi}"
    elif pubmed_id:
        paper_id = f"PMID:{pubmed_id}"

    result = {
        'paper': None,
        'citing': [],  # Papers that cite this paper
        'cited': [],   # Papers that this paper cites
        'chain': []    # Citation chain visualization
    }

    try:
        # Fetch paper details
        paper_url = f"https://api.semanticscholar.org/graph/v1/paper/{paper_id}"
        params = {
            'fields': 'paperId,title,year,authors,venue,citationCount,citations.paperId,citations.title,citations.year,citations.authors,references.paperId,references.title,references.year,references.authors'
        }

        resp = requests.get(paper_url, params=params, timeout=10)
        if resp.status_code != 200:
            return {'error': f'Paper not found in Semantic Scholar: {resp.status_code}'}

        paper_data = resp.json()

        # Main paper info
        result['paper'] = {
            'id': paper_data.get('paperId', ''),
            'title': paper_data.get('title', ''),
            'year': paper_data.get('year'),
            'authors': [a.get('name', '') for a in paper_data.get('authors', [])[:3]],
            'venue': paper_data.get('venue', ''),
            'citation_count': paper_data.get('citationCount', 0)
        }

        # Papers that cite this paper (incoming citations)
        citations = paper_data.get('citations', []) or []
        for cit in citations[:10]:  # Limit to top 10
            if cit:
                result['citing'].append({
                    'id': cit.get('paperId', ''),
                    'title': cit.get('title', ''),
                    'year': cit.get('year'),
                    'authors': [a.get('name', '') for a in (cit.get('authors', []) or [])[:2]]
                })

        # Papers this paper cites (outgoing references)
        references = paper_data.get('references', []) or []
        for ref in references[:10]:  # Limit to top 10
            if ref:
                result['cited'].append({
                    'id': ref.get('paperId', ''),
                    'title': ref.get('title', ''),
                    'year': ref.get('year'),
                    'authors': [a.get('name', '') for a in (ref.get('authors', []) or [])[:2]]
                })

        # Build simple chain visualization
        # Format: [citing papers] -> [this paper] -> [cited papers]
        chain = []

        # Add citing papers (newer papers that cite this one)
        for cit in result['citing'][:3]:
            chain.append({
                'title': cit['title'][:50] + '...' if len(cit.get('title', '')) > 50 else cit.get('title', ''),
                'year': cit.get('year'),
                'type': 'citing'
            })

        # Add the main paper
        chain.append({
            'title': result['paper']['title'][:50] + '...' if len(result['paper'].get('title', '')) > 50 else result['paper'].get('title', ''),
            'year': result['paper'].get('year'),
            'type': 'main'
        })

        # Add cited papers (older papers that this one cites)
        for ref in result['cited'][:3]:
            chain.append({
                'title': ref['title'][:50] + '...' if len(ref.get('title', '')) > 50 else ref.get('title', ''),
                'year': ref.get('year'),
                'type': 'cited'
            })

        result['chain'] = chain

    except requests.exceptions.Timeout:
        return {'error': 'Timeout fetching citation data'}
    except Exception as e:
        logger.warning(f"Citation chain fetch failed: {e}")
        return {'error': str(e)}

    return result


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


def extract_ligands_from_pdb(pdb_content: str, source_name: str = "", target_chain: str = None) -> List[Dict]:
    """
    Extract ligand atoms from PDB content.

    Args:
        pdb_content: PDB file content
        source_name: Name/ID of the source structure
        target_chain: If specified, only extract ligands from this chain
    """
    ligands = {}

    for line in pdb_content.split('\n'):
        if line.startswith('HETATM'):
            resname = line[17:20].strip()
            if resname in EXCLUDE_RESIDUES:
                continue

            chain = line[21:22].strip()

            # Filter by chain if specified
            if target_chain and chain != target_chain:
                continue

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


def align_structure_pymol(query_pdb: str, target_pdb: str, pymol_path: str = "pymol", target_chain: str = None) -> Tuple[str, float]:
    """
    Align target structure to query using PyMOL.
    Returns aligned PDB content and RMSD.

    Args:
        query_pdb: Query PDB content
        target_pdb: Target PDB content
        pymol_path: Path to PyMOL executable
        target_chain: If specified, only align and save this chain from target
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        query_path = os.path.join(tmpdir, "query.pdb")
        target_path = os.path.join(tmpdir, "target.pdb")
        aligned_path = os.path.join(tmpdir, "aligned.pdb")

        with open(query_path, 'w') as f:
            f.write(query_pdb)
        with open(target_path, 'w') as f:
            f.write(target_pdb)

        # Build selection string for target chain
        target_sel = f"target and chain {target_chain}" if target_chain else "target"

        script = f'''
from pymol import cmd
cmd.load("{query_path}", "query")
cmd.load("{target_path}", "target")
# Align only the specified chain (or whole structure if no chain specified)
result = cmd.align("{target_sel}", "query")
# Save only the aligned chain
cmd.save("{aligned_path}", "{target_sel}")
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

        # Add taxonomy, citation, and organism_type info to hits
        for hit in hits_list:
            pdb_id = hit['pdb_id'].upper()
            tax = taxonomy_data.get(pdb_id, {})
            hit['species'] = tax.get('scientific_name') or 'Unknown'
            hit['taxonomy_id'] = tax.get('taxonomy_id')
            hit['common_name'] = tax.get('common_name') or ''
            hit['organism_type'] = tax.get('organism_type', 'unknown')
            hit['structure_title'] = tax.get('structure_title', '')
            hit['citation'] = tax.get('citation')

        jobs[job_id]["hits"] = hits_list
        jobs[job_id]["taxonomy_data"] = taxonomy_data

        # Group hits by species
        species_groups = group_hits_by_species(hits_list.copy(), taxonomy_data)
        jobs[job_id]["species_groups"] = {k: len(v) for k, v in species_groups.items()}

        # Group hits by organism type (human, mammal, viral, bacteria, etc.)
        organism_type_groups = defaultdict(int)
        for hit in hits_list:
            org_type = hit.get('organism_type', 'unknown')
            organism_type_groups[org_type] += 1
        jobs[job_id]["organism_type_groups"] = dict(organism_type_groups)
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

                # Extract PDB ID and chain from key (format: PDBID_CHAIN, e.g., 5CNN_A)
                source_pdb_id = key.split('_')[0]
                target_chain = key.split('_')[1] if '_' in key else None

                # Align if PyMOL available - only align the specific chain
                if pymol_available:
                    aligned_pdb, rmsd = align_structure_pymol(query_pdb, target_pdb, target_chain=target_chain)
                else:
                    aligned_pdb = target_pdb  # Use unaligned

                # Extract sequence for conservation
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

                # Extract ligands - only from the matched chain to avoid ligands from other chains
                logger.info(f"Extracting ligands from {source_pdb_id}, chain={target_chain}, key={key}")
                ligands = extract_ligands_from_pdb(aligned_pdb, source_pdb_id, target_chain)
                logger.info(f"Found {len(ligands)} ligands from chain {target_chain}")

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
                    'chain': target_chain,
                    'key': key,  # e.g., "5CNN_A"
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

        # Fetch metadata from RCSB for unique ligand IDs
        unique_ligand_names = list(set(lig["name"] for lig in ligands_info))
        rcsb_metadata = fetch_ligand_metadata_batch(unique_ligand_names)

        def get_ligand_display_info(lig):
            """Get display info for a ligand from Kamaji or RCSB."""
            kamaji_info = lig.get("kamaji_info", {})
            rcsb_info = rcsb_metadata.get(lig["name"].upper(), {})

            # Prefer Kamaji MW if available, else use RCSB
            mw = kamaji_info.get("mw", 0) if kamaji_info else 0
            if not mw:
                mw = rcsb_info.get("molecular_weight", 0)

            # Get description from Kamaji class or RCSB name
            description = ""
            if kamaji_info and kamaji_info.get("class"):
                description = kamaji_info["class"]
            elif rcsb_info.get("description"):
                description = rcsb_info["description"]

            # Get formula from RCSB
            formula = rcsb_info.get("formula", "")

            return {
                "name": lig["name"],
                "source": lig["source"],
                "category_name": lig.get("category_name", ""),
                "kamaji_info": kamaji_info,
                "molecular_weight": mw,
                "formula": formula,
                "description": description,
            }

        jobs[job_id]["ligand_clusters"] = {
            cat: [get_ligand_display_info(lig) for lig in ligs]
            for cat, ligs in ligand_clusters.items()
        }
        jobs[job_id]["kamaji_available"] = KAMAJI_AVAILABLE

        jobs[job_id]["ligand_pdbs"] = unique_ligands
        jobs[job_id]["ligands"] = ligands_info
        jobs[job_id]["aligned_count"] = len(aligned_seqs)
        jobs[job_id]["aligned_seqs_with_species"] = aligned_seqs_with_species
        jobs[job_id]["aligned_pdbs"] = aligned_pdbs  # Store aligned PDB data for download

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
        "organism_type_groups": job.get("organism_type_groups", {}),
        "taxonomy_data": job.get("taxonomy_data", {}),
        "kamaji_available": job.get("kamaji_available", False),
    })


@app.route('/api/conservation/<job_id>', methods=['POST'])
def get_filtered_conservation(job_id):
    """
    Get conservation analysis filtered by species or selected hits.

    Request body:
    {
        "include_species": ["Cytomegalovirus", "HIV-1"],  // optional
        "exclude_species": ["Homo sapiens"],  // optional
        "selected_hits": [{"pdb_id": "1ABC", "chain": "A"}, ...]  // optional - specific hits to include
    }
    """
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404

    job = jobs[job_id]
    data = request.json or {}

    include_species = set(data.get('include_species', []))
    exclude_species = set(data.get('exclude_species', []))
    selected_hits = data.get('selected_hits', [])  # List of {pdb_id, chain} dicts

    # Get stored aligned sequences with species info
    aligned_seqs_with_species = job.get("aligned_seqs_with_species", [])

    if not aligned_seqs_with_species:
        # Fallback to regular conservation
        return jsonify({
            "conservation": job.get("conservation", []),
            "filtered": False
        })

    # Filter by selected hits if provided
    if selected_hits:
        # Build set of selected hit keys
        selected_keys = set(f"{h['pdb_id'].upper()}_{h['chain'].upper()}" for h in selected_hits)

        # Filter aligned sequences to only include selected hits
        aligned_seqs_with_species = [
            seq_info for seq_info in aligned_seqs_with_species
            if f"{seq_info.get('pdb_id', '').upper()}_{seq_info.get('chain', '').upper()}" in selected_keys
        ]

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
        "exclude_species": list(exclude_species),
        "selected_hits_count": len(selected_hits) if selected_hits else len(aligned_seqs_with_species)
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

        # Write ligands with Kamaji category info
        ligand_paths = []

        # Create lookup for ligand category info
        ligands_info = job.get("ligands", [])
        ligand_category_lookup = {}
        for lig_info in ligands_info:
            key = f"{lig_info.get('name', '')}_{lig_info.get('source', '')}"
            ligand_category_lookup[key] = {
                'category': lig_info.get('category', 'ligand'),
                'category_name': lig_info.get('category_name', 'Small Molecules')
            }

        for i, lig in enumerate(job.get("ligand_pdbs", [])):
            lig_path = os.path.join(tmpdir, f"ligand_{i}_{lig['name']}_{lig['source']}.pdb")
            with open(lig_path, 'w') as f:
                f.write(lig.get('pdb', ''))

            # Get Kamaji category for this ligand
            lig_key = f"{lig['name']}_{lig.get('source', '')}"
            cat_info = ligand_category_lookup.get(lig_key, {})
            category = cat_info.get('category', 'ligand')
            category_name = cat_info.get('category_name', 'Small Molecules')

            ligand_paths.append((lig_path, lig['name'], lig['source'], category, category_name))

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

        # Add ligands grouped by Kamaji category
        if ligand_paths:
            script += '\n# Load ligands\n'

            # Group ligands by Kamaji category
            ligands_by_category = {}
            for lig_path, lig_name, lig_source, category, category_name in ligand_paths:
                obj_name = f"lig_{lig_name}_{lig_source}".replace('-', '_')
                script += f'cmd.load("{lig_path}", "{obj_name}")\n'
                script += f'cmd.show("sticks", "{obj_name}")\n'
                script += f'util.cbag("{obj_name}")  # Color by atom type\n'

                # Group by Kamaji category
                cat_key = category if category else 'ligand'
                if cat_key not in ligands_by_category:
                    ligands_by_category[cat_key] = []
                ligands_by_category[cat_key].append(obj_name)

            # Create Kamaji category subgroups for ligands
            script += '\n# Group ligands by Kamaji category\n'
            lig_cat_group_names = []
            for category, lig_names in ligands_by_category.items():
                group_name = f"lig_{category}"
                lig_cat_group_names.append(group_name)
                script += f'cmd.group("{group_name}", "{" ".join(lig_names)}")\n'

            # Create main ligands group
            script += f'cmd.group("ligands", "{" ".join(lig_cat_group_names)}")\n'

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

# Create surface colored by conservation (blue=low, white=mid, red=high)
cmd.create("protein_surface", "query_{pdb_id}")
cmd.show("surface", "protein_surface")
cmd.spectrum("b", "blue_white_red", "protein_surface", minimum=0, maximum=100)
cmd.set("transparency", 0.4, "protein_surface")
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


@app.route('/api/citation_chain', methods=['GET'])
def get_citation_chain():
    """Get citation chain for a paper by DOI or PubMed ID."""
    doi = request.args.get('doi', '')
    pubmed_id = request.args.get('pubmed_id', '')

    if not doi and not pubmed_id:
        return jsonify({"error": "Please provide doi or pubmed_id parameter"}), 400

    result = fetch_citation_chain(doi=doi, pubmed_id=pubmed_id)

    if 'error' in result:
        return jsonify(result), 404

    return jsonify(result)


@app.route('/api/download_selected/<job_id>', methods=['POST'])
def download_selected(job_id):
    """Generate and download PyMOL session file with selected structures and ligands."""
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404

    job = jobs[job_id]

    if job["status"] != JobStatus.COMPLETED:
        return jsonify({"error": "Job not completed yet"}), 400

    data = request.get_json() or {}
    selected_hits = data.get('selected_hits', [])
    selected_ligand_hits = data.get('selected_ligand_hits', [])

    if not selected_hits:
        return jsonify({"error": "No structures selected"}), 400

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

        # Get selected hits data
        all_hits = job.get("hits", [])
        taxonomy_data = job.get("taxonomy_data", {})
        aligned_pdbs = job.get("aligned_pdbs", [])  # Get pre-aligned structures

        # Create lookup for aligned PDBs by key (pdb_id_chain, e.g., "5CNN_A")
        aligned_pdb_lookup = {}
        for apdb in aligned_pdbs:
            # Use key if available, otherwise construct from pdb_id and chain
            key = apdb.get('key', f"{apdb['pdb_id']}_{apdb.get('chain', '')}").lower()
            aligned_pdb_lookup[key] = apdb['pdb']

        # Extract PDB IDs from "pdb_chain" format for ligand matching
        selected_ligand_pdb_ids = set(h.split('_')[0].lower() for h in selected_ligand_hits)

        # Write selected aligned structures - use pre-aligned PDBs from job data
        hit_paths = []
        for hit in all_hits:
            hit_pdb_id = hit.get('pdb_id', '').lower()
            hit_chain = hit.get('chain', '')
            hit_key = f"{hit_pdb_id}_{hit_chain}".lower()

            # Check if this specific pdb_chain combination is selected
            if hit_key in [h.lower() for h in selected_hits]:
                # Use pre-aligned structure from job data
                hit_pdb = aligned_pdb_lookup.get(hit_key)

                if hit_pdb:
                    hit_path = os.path.join(tmpdir, f"hit_{hit_key}.pdb")
                    with open(hit_path, 'w') as f:
                        f.write(hit_pdb)
                    # Get taxonomy info
                    tax_info = taxonomy_data.get(hit_pdb_id.upper(), {})
                    species = tax_info.get('scientific_name', 'Unknown')
                    organism_type = tax_info.get('organism_type', 'unknown')
                    hit_paths.append((hit_path, hit_key, hit_pdb_id, species, organism_type))
                else:
                    logger.warning(f"Aligned structure not found for {hit_key}")

        # Write ligands from selected ligand hits only, grouped by Kamaji category
        ligand_paths = []

        # Create lookup for ligand category info
        ligands_info = job.get("ligands", [])
        ligand_category_lookup = {}
        for lig_info in ligands_info:
            key = f"{lig_info.get('name', '')}_{lig_info.get('source', '')}"
            ligand_category_lookup[key] = {
                'category': lig_info.get('category', 'ligand'),
                'category_name': lig_info.get('category_name', 'Small Molecules')
            }

        for i, lig in enumerate(job.get("ligand_pdbs", [])):
            lig_source = lig.get('source', '').lower()
            # Only include ligands from selected structures (match by PDB ID)
            if lig_source in selected_ligand_pdb_ids:
                lig_path = os.path.join(tmpdir, f"ligand_{i}_{lig['name']}_{lig_source}.pdb")
                with open(lig_path, 'w') as f:
                    f.write(lig.get('pdb', ''))

                # Get Kamaji category for this ligand
                lig_key = f"{lig['name']}_{lig.get('source', '')}"
                cat_info = ligand_category_lookup.get(lig_key, {})
                category = cat_info.get('category', 'ligand')
                category_name = cat_info.get('category_name', 'Small Molecules')

                ligand_paths.append((lig_path, lig['name'], lig_source, category, category_name))

        # Create PyMOL script
        pse_path = os.path.join(tmpdir, f"{pdb_id}_selected.pse")

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

        # Add aligned structures grouped by organism type
        if hit_paths:
            script += '\n# Load aligned structures\n'

            # Group structures by organism type
            structures_by_org = {}
            for hit_path, hit_id, hit_pdb_id, species, org_type in hit_paths:
                obj_name = f"hit_{hit_id}".replace('-', '_')
                script += f'cmd.load("{hit_path}", "{obj_name}")\n'
                script += f'cmd.show("cartoon", "{obj_name}")\n'
                script += f'cmd.color("gray70", "{obj_name}")\n'

                # Group by organism type
                org_key = org_type if org_type else 'unknown'
                if org_key not in structures_by_org:
                    structures_by_org[org_key] = []
                structures_by_org[org_key].append(obj_name)

            # Create organism type subgroups
            script += '\n# Group structures by organism type\n'
            org_group_names = []
            for org_type, struct_names in structures_by_org.items():
                group_name = f"struct_{org_type}"
                org_group_names.append(group_name)
                script += f'cmd.group("{group_name}", "{" ".join(struct_names)}")\n'

            # Create main aligned_structures group
            script += f'cmd.group("aligned_structures", "{" ".join(org_group_names)}")\n'
            script += 'cmd.disable("aligned_structures")  # Hidden by default\n'

        # Add ligands grouped by Kamaji category
        if ligand_paths:
            script += '\n# Load ligands\n'

            # Group ligands by Kamaji category
            ligands_by_category = {}
            for lig_path, lig_name, lig_source, category, category_name in ligand_paths:
                obj_name = f"lig_{lig_name}_{lig_source}".replace('-', '_')
                script += f'cmd.load("{lig_path}", "{obj_name}")\n'
                script += f'cmd.show("sticks", "{obj_name}")\n'
                script += f'util.cbag("{obj_name}")  # Color by atom type\n'

                # Group by Kamaji category
                cat_key = category if category else 'ligand'
                if cat_key not in ligands_by_category:
                    ligands_by_category[cat_key] = []
                ligands_by_category[cat_key].append(obj_name)

            # Create Kamaji category subgroups for ligands
            script += '\n# Group ligands by Kamaji category\n'
            lig_cat_group_names = []
            for category, lig_names in ligands_by_category.items():
                group_name = f"lig_{category}"
                lig_cat_group_names.append(group_name)
                script += f'cmd.group("{group_name}", "{" ".join(lig_names)}")\n'

            # Create main ligands group
            script += f'cmd.group("ligands", "{" ".join(lig_cat_group_names)}")\n'

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

# Create surface colored by conservation (blue=low, white=mid, red=high)
cmd.create("protein_surface", "query_{pdb_id}")
cmd.show("surface", "protein_surface")
cmd.spectrum("b", "blue_white_red", "protein_surface", minimum=0, maximum=100)
cmd.set("transparency", 0.4, "protein_surface")
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
            output_path = RESULTS_FOLDER / job_id / f"{pdb_id}_selected.pse"
            import shutil
            shutil.copy(pse_path, output_path)

            return send_file(
                output_path,
                as_attachment=True,
                download_name=f"{pdb_id}_selected.pse",
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
