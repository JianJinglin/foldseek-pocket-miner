"""
PyMOL session generator module - Create comprehensive PyMOL sessions with conservation and ligands.
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import List, Optional, Dict
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class SessionConfig:
    """Configuration for PyMOL session generation."""
    show_surface: bool = True
    surface_transparency: float = 0.7
    conservation_spectrum: str = "blue_white_red"  # Low to high conservation
    ligand_representation: str = "sticks"
    ligand_color_by: str = "source"  # 'source', 'element', 'pharmacophore'
    show_pocket_residues: bool = True
    pocket_cutoff: float = 5.0  # Angstroms from ligand


class PyMOLSessionGenerator:
    """Generate PyMOL sessions with conservation coloring and ligands."""

    def __init__(
        self,
        output_dir: str,
        pymol_path: str = "pymol",
        config: Optional[SessionConfig] = None
    ):
        """
        Initialize the session generator.

        Args:
            output_dir: Directory to save sessions
            pymol_path: Path to PyMOL executable
            config: Session configuration
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.pymol_path = pymol_path
        self.config = config or SessionConfig()

    def generate_session(
        self,
        query_structure: str,
        ligand_files: List[str],
        conservation_scores: Optional[Dict[int, float]] = None,
        output_name: str = "pocket_analysis",
        aligned_structures: Optional[List[str]] = None
    ) -> str:
        """
        Generate a PyMOL session with all components.

        Args:
            query_structure: Path to query structure
            ligand_files: List of paths to ligand files
            conservation_scores: Dictionary mapping residue number to conservation score (0-1)
            output_name: Name for output session file
            aligned_structures: Optional list of aligned structures to include

        Returns:
            Path to generated .pse file
        """
        output_path = self.output_dir / f"{output_name}.pse"
        script_path = self.output_dir / f"{output_name}_script.py"

        # Generate PyMOL script
        script = self._generate_script(
            query_structure=query_structure,
            ligand_files=ligand_files,
            conservation_scores=conservation_scores,
            output_path=str(output_path),
            aligned_structures=aligned_structures
        )

        # Write script
        with open(script_path, 'w') as f:
            f.write(script)

        # Run PyMOL
        try:
            result = subprocess.run(
                [self.pymol_path, "-cq", str(script_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode != 0:
                logger.warning(f"PyMOL warnings: {result.stderr}")

            if output_path.exists():
                logger.info(f"Generated PyMOL session: {output_path}")
                return str(output_path)
            else:
                raise RuntimeError(f"PyMOL session not created: {result.stderr}")

        except subprocess.TimeoutExpired:
            raise RuntimeError("PyMOL session generation timed out")
        finally:
            # Optionally clean up script
            # script_path.unlink()
            pass

    def _generate_script(
        self,
        query_structure: str,
        ligand_files: List[str],
        conservation_scores: Optional[Dict[int, float]],
        output_path: str,
        aligned_structures: Optional[List[str]] = None
    ) -> str:
        """Generate PyMOL Python script."""
        query_name = Path(query_structure).stem

        script_lines = [
            "from pymol import cmd, stored",
            "import os",
            "",
            "# Initialize PyMOL",
            "cmd.reinitialize()",
            "",
            f"# Load query structure",
            f'cmd.load("{query_structure}", "query")',
            "",
        ]

        # Load aligned structures if provided
        if aligned_structures:
            script_lines.append("# Load aligned structures")
            for i, struct_path in enumerate(aligned_structures[:10]):  # Limit to 10
                struct_name = f"aligned_{i}"
                script_lines.append(f'cmd.load("{struct_path}", "{struct_name}")')
            script_lines.append("")

        # Load ligands
        script_lines.append("# Load ligands")
        ligand_names = []
        for i, lig_path in enumerate(ligand_files):
            lig_name = f"ligand_{i}"
            ligand_names.append(lig_name)
            script_lines.append(f'cmd.load("{lig_path}", "{lig_name}")')
        script_lines.append("")

        # Color scheme for ligands
        colors = [
            "salmon", "lime", "cyan", "yellow", "orange",
            "pink", "palegreen", "lightblue", "wheat", "violet"
        ]

        # Set representations
        script_lines.extend([
            "# Set representations",
            'cmd.show("cartoon", "query")',
            'cmd.color("gray80", "query")',
            "",
        ])

        # Color ligands
        script_lines.append("# Color and display ligands")
        for i, lig_name in enumerate(ligand_names):
            color = colors[i % len(colors)]
            script_lines.extend([
                f'cmd.show("sticks", "{lig_name}")',
                f'cmd.color("{color}", "{lig_name}")',
                f'cmd.set("stick_radius", 0.15, "{lig_name}")',
            ])
        script_lines.append("")

        # Apply conservation coloring if provided
        if conservation_scores:
            script_lines.extend([
                "# Apply conservation coloring",
                "stored.conservation = {}",
            ])

            # Add conservation scores
            for resseq, score in conservation_scores.items():
                # Scale to 0-100 for b-factor
                bfactor = score * 100
                script_lines.append(f"stored.conservation[{resseq}] = {bfactor}")

            script_lines.extend([
                "",
                "# Set b-factors based on conservation",
                'cmd.alter("query and name CA", "b=stored.conservation.get(resi, 50)")',
                "",
                "# Color by conservation (spectrum)",
                'cmd.spectrum("b", "blue_white_red", "query", minimum=0, maximum=100)',
                "",
            ])

        # Find and show pocket residues
        if ligand_names:
            script_lines.extend([
                "# Create pocket selection",
                'cmd.select("pocket", f"query within {self.config.pocket_cutoff} of (' + ' or '.join(ligand_names) + ')")',
                'cmd.show("sticks", "pocket")',
                "",
            ])

        # Surface settings
        if self.config.show_surface:
            script_lines.extend([
                "# Show surface",
                'cmd.show("surface", "query")',
                f'cmd.set("transparency", {self.config.surface_transparency}, "query")',
                'cmd.set("surface_color", "white", "query")',
                "",
            ])

        # Center and orient view
        script_lines.extend([
            "# Center view on ligands",
            f'cmd.center("{" or ".join(ligand_names)}")',
            'cmd.zoom("all", buffer=5)',
            "",
        ])

        # Add labels and annotations
        script_lines.extend([
            "# Add legend/notes",
            'cmd.set("label_size", 20)',
            'cmd.set("label_color", "black")',
            "",
        ])

        # Rendering settings
        script_lines.extend([
            "# Rendering settings",
            'cmd.set("ray_shadows", "off")',
            'cmd.set("antialias", 2)',
            'cmd.set("depth_cue", 0)',
            'cmd.bg_color("white")',
            "",
        ])

        # Save session
        script_lines.extend([
            "# Save session",
            f'cmd.save("{output_path}")',
            "",
            "# Quit",
            "cmd.quit()",
        ])

        return "\n".join(script_lines)

    def generate_quick_session(
        self,
        structures: List[str],
        output_name: str = "quick_view"
    ) -> str:
        """
        Generate a quick session with multiple structures.

        Args:
            structures: List of structure files
            output_name: Output session name

        Returns:
            Path to session file
        """
        output_path = self.output_dir / f"{output_name}.pse"
        script_path = self.output_dir / f"{output_name}_script.py"

        script_lines = [
            "from pymol import cmd",
            "",
            "cmd.reinitialize()",
            "",
        ]

        colors = [
            "marine", "green", "orange", "yellow", "cyan",
            "salmon", "pink", "teal", "purple", "lime"
        ]

        for i, struct_path in enumerate(structures):
            name = Path(struct_path).stem
            color = colors[i % len(colors)]
            script_lines.extend([
                f'cmd.load("{struct_path}", "{name}")',
                f'cmd.show("cartoon", "{name}")',
                f'cmd.color("{color}", "{name}")',
            ])

        script_lines.extend([
            "",
            'cmd.zoom("all")',
            'cmd.bg_color("white")',
            f'cmd.save("{output_path}")',
            "cmd.quit()",
        ])

        with open(script_path, 'w') as f:
            f.write("\n".join(script_lines))

        subprocess.run(
            [self.pymol_path, "-cq", str(script_path)],
            capture_output=True,
            timeout=120
        )

        return str(output_path)


def create_conservation_pymol_script(
    pdb_path: str,
    conservation_scores: Dict[int, float],
    output_script: str
) -> str:
    """
    Create a standalone PyMOL script for conservation coloring.

    Args:
        pdb_path: Path to PDB file
        conservation_scores: Dictionary of residue number to conservation score
        output_script: Path to output script file

    Returns:
        Path to script file
    """
    lines = [
        "from pymol import cmd, stored",
        "",
        f'cmd.load("{pdb_path}", "protein")',
        "",
        "# Conservation scores (0-100 scale)",
        "stored.conservation = {",
    ]

    for resseq, score in sorted(conservation_scores.items()):
        lines.append(f"    {resseq}: {score * 100:.1f},")

    lines.extend([
        "}",
        "",
        '# Apply conservation to B-factors',
        'cmd.alter("protein", "b=stored.conservation.get(int(resi), 50)")',
        "",
        '# Color by conservation',
        'cmd.spectrum("b", "blue_white_red", "protein", minimum=0, maximum=100)',
        'cmd.show("cartoon", "protein")',
        "",
        '# Create legend note',
        'print("Conservation coloring: blue=variable, white=moderate, red=conserved")',
    ])

    with open(output_script, 'w') as f:
        f.write("\n".join(lines))

    return output_script
