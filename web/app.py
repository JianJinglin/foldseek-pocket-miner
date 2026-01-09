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
from typing import Dict, Any, Optional, List, Tuple
from collections import Counter
import numpy as np

from flask import Flask, render_template, request, jsonify, send_file, send_from_directory
from flask_cors import CORS

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from foldseek_pocket_miner.structure_downloader import StructureDownloader
from foldseek_pocket_miner.foldseek_search import FoldseekSearcher
from foldseek_pocket_miner.ligand_extractor import LigandExtractor, EXCLUDE_RESIDUES

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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

        jobs[job_id]["hits"] = [
            {"pdb_id": h.pdb_id, "chain": h.chain, "evalue": h.evalue,
             "score": h.score, "seq_identity": h.seq_identity}
            for h in hits
        ]
        jobs[job_id]["progress"] = 30

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

                # Extract sequence for conservation
                target_chain = key.split('_')[1] if '_' in key else None
                seq = extract_sequence_from_pdb(aligned_pdb, target_chain)
                if seq:
                    aligned_seqs.append(seq)

                # Extract ligands
                source_pdb_id = key.split('_')[0]
                ligands = extract_ligands_from_pdb(aligned_pdb, source_pdb_id)

                for lig in ligands:
                    # Build ligand PDB string
                    lig_pdb = f"REMARK Ligand {lig['name']} from {lig['source']}\n"
                    for atom in lig['atoms']:
                        lig_pdb += atom['line'] + '\n'
                    lig_pdb += "END\n"

                    all_ligands.append({
                        'name': lig['name'],
                        'source': lig['source'],
                        'chain': lig['chain'],
                        'atoms': len(lig['atoms']),
                        'pdb': lig_pdb
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

        jobs[job_id]["ligand_pdbs"] = unique_ligands
        jobs[job_id]["ligands"] = [
            {"name": l["name"], "source": l["source"], "chain": l["chain"], "atoms": l["atoms"]}
            for l in unique_ligands
        ]
        jobs[job_id]["aligned_count"] = len(aligned_seqs)

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
        "conservation": job.get("conservation", []),
        "downloaded_count": job.get("downloaded_count", 0),
        "aligned_count": job.get("aligned_count", 0),
        "query_pdb": job.get("query_pdb", ""),
    })


@app.route('/api/pdb/<pdb_id>')
def fetch_pdb(pdb_id):
    """Fetch PDB structure from RCSB."""
    import requests
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
