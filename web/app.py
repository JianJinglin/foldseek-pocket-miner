"""
Flask web application for Foldseek Pocket Miner.
"""

import os
import json
import uuid
import logging
import threading
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional

from flask import Flask, render_template, request, jsonify, send_file, send_from_directory
from flask_cors import CORS

# Add parent directory to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from foldseek_pocket_miner import PocketMiner
from foldseek_pocket_miner.structure_downloader import StructureDownloader
from foldseek_pocket_miner.foldseek_search import FoldseekSearcher
from foldseek_pocket_miner.ligand_extractor import LigandExtractor

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

# Job storage (in production, use Redis or database)
jobs: Dict[str, Dict[str, Any]] = {}


class JobStatus:
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


def run_pipeline(job_id: str, pdb_id: str = None, chain: str = "A",
                structure_path: str = None, max_hits: int = 50):
    """Run the pocket mining pipeline in a background thread."""
    try:
        jobs[job_id]["status"] = JobStatus.RUNNING
        jobs[job_id]["progress"] = 0
        jobs[job_id]["current_step"] = "Initializing..."

        output_dir = RESULTS_FOLDER / job_id
        output_dir.mkdir(exist_ok=True)

        # Step 1: Prepare query structure
        jobs[job_id]["current_step"] = "Preparing query structure..."
        jobs[job_id]["progress"] = 5

        if pdb_id and not structure_path:
            downloader = StructureDownloader(output_dir=str(output_dir / "query"))
            structure_path = downloader.download_structure(pdb_id, chain)
            if not structure_path:
                raise Exception(f"Failed to download {pdb_id}")

        jobs[job_id]["query_structure"] = structure_path
        jobs[job_id]["progress"] = 10

        # Step 2: Run Foldseek search
        jobs[job_id]["current_step"] = "Searching for similar structures with Foldseek..."
        jobs[job_id]["progress"] = 15

        searcher = FoldseekSearcher(
            database="pdb100",
            use_web_api=True  # Use Foldseek web API (no local database needed)
        )

        hits = searcher.search(
            structure_path=structure_path,
            max_hits=max_hits,
            evalue=1e-3
        )

        jobs[job_id]["hits"] = [
            {
                "pdb_id": h.pdb_id,
                "chain": h.chain,
                "evalue": h.evalue,
                "score": h.score,
                "seq_identity": h.seq_identity
            }
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

        struct_downloader = StructureDownloader(
            output_dir=str(output_dir / "structures")
        )

        pdb_ids = [h.pdb_id for h in hits]
        chains = [h.chain for h in hits]
        downloaded = struct_downloader.download_batch(pdb_ids, chains, show_progress=False)

        valid_paths = {k: v for k, v in downloaded.items() if v}
        jobs[job_id]["downloaded_count"] = len(valid_paths)
        jobs[job_id]["progress"] = 50

        # Step 4: Filter by ligands
        jobs[job_id]["current_step"] = "Filtering structures with ligands..."
        jobs[job_id]["progress"] = 55

        structures_with_ligands = struct_downloader.filter_by_ligand(list(valid_paths.keys()))
        ligand_pdb_ids = {item[0] for item in structures_with_ligands}

        # Store ligand info
        jobs[job_id]["ligand_info"] = [
            {"pdb_id": pdb_id, "ligands": [l["id"] for l in ligs]}
            for pdb_id, ligs in structures_with_ligands
        ]
        jobs[job_id]["progress"] = 65

        # Step 5: Extract ligands from downloaded structures
        jobs[job_id]["current_step"] = "Extracting ligands..."
        jobs[job_id]["progress"] = 70

        extractor = LigandExtractor(output_dir=str(output_dir / "ligands"))
        all_ligands = []

        for key, path in valid_paths.items():
            if path and key.split('_')[0] in ligand_pdb_ids:
                ligands = extractor.extract_ligands(path)
                for lig in ligands:
                    lig_path = extractor.save_ligand(lig)
                    all_ligands.append({
                        "name": lig.resname,
                        "source": lig.source_pdb,
                        "chain": lig.chain,
                        "atoms": lig.num_atoms,
                        "path": lig_path
                    })

        jobs[job_id]["ligands"] = all_ligands
        jobs[job_id]["progress"] = 85

        # Step 6: Prepare visualization data
        jobs[job_id]["current_step"] = "Preparing visualization..."
        jobs[job_id]["progress"] = 90

        # Read query structure for visualization
        with open(structure_path, 'r') as f:
            jobs[job_id]["query_pdb"] = f.read()

        # Read ligand structures
        ligand_pdbs = []
        for lig_info in all_ligands[:20]:  # Limit to 20 ligands for visualization
            if os.path.exists(lig_info["path"]):
                with open(lig_info["path"], 'r') as f:
                    ligand_pdbs.append({
                        "name": lig_info["name"],
                        "source": lig_info["source"],
                        "pdb": f.read()
                    })
        jobs[job_id]["ligand_pdbs"] = ligand_pdbs

        jobs[job_id]["progress"] = 100
        jobs[job_id]["current_step"] = "Complete!"
        jobs[job_id]["status"] = JobStatus.COMPLETED
        jobs[job_id]["completed_at"] = datetime.now().isoformat()

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

    # Create job
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
    }

    # Start background thread
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
        "ligand_info": job.get("ligand_info", []),
        "downloaded_count": job.get("downloaded_count", 0),
        "query_pdb": job.get("query_pdb", ""),
        "ligand_pdbs": job.get("ligand_pdbs", []),
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


@app.route('/api/ligand/<job_id>/<int:ligand_idx>')
def get_ligand_pdb(job_id, ligand_idx):
    """Get ligand PDB content."""
    if job_id not in jobs:
        return jsonify({"error": "Job not found"}), 404

    job = jobs[job_id]
    ligands = job.get("ligands", [])

    if ligand_idx >= len(ligands):
        return jsonify({"error": "Ligand not found"}), 404

    ligand_path = ligands[ligand_idx].get("path")
    if ligand_path and os.path.exists(ligand_path):
        with open(ligand_path, 'r') as f:
            return f.read(), 200, {'Content-Type': 'text/plain'}

    return jsonify({"error": "Ligand file not found"}), 404


@app.route('/static/<path:filename>')
def serve_static(filename):
    """Serve static files."""
    return send_from_directory('static', filename)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
