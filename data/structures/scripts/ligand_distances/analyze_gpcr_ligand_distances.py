#!/usr/bin/env python3
"""
Analyze GPCR structures to calculate distances between receptor residues and bound ligands.
"""

import os
import numpy as np
from collections import defaultdict
import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
CIF_DIR = REPO_ROOT / "structures" / "raw" / "experimental"

try:
    from Bio.PDB import MMCIFParser, Selection
    from Bio.PDB.Polypeptide import is_aa
except ImportError:
    print("BioPython not found. Installing...")
    import subprocess
    subprocess.check_call(["pip", "install", "biopython"])
    from Bio.PDB import MMCIFParser, Selection
    from Bio.PDB.Polypeptide import is_aa


def get_ligand_common_name(pdb_code):
    """
    Placeholder for ligand common name - returns empty string for manual entry.

    Parameters:
    -----------
    pdb_code : str
        The 3-letter PDB ligand code

    Returns:
    --------
    str : Empty string (to be filled manually)
    """
    return ""


def get_residue_atoms(residue):
    """Get all atoms from a residue."""
    return list(residue.get_atoms())


def calculate_min_distance(residue, ligand):
    """
    Calculate minimum distance between any atom in a residue and any atom in a ligand.

    Parameters:
    -----------
    residue : Bio.PDB.Residue
        The receptor residue
    ligand : Bio.PDB.Residue
        The ligand residue

    Returns:
    --------
    float : Minimum distance in Angstroms
    """
    res_atoms = get_residue_atoms(residue)
    lig_atoms = get_residue_atoms(ligand)

    if not res_atoms or not lig_atoms:
        return None

    min_dist = float('inf')
    for res_atom in res_atoms:
        for lig_atom in lig_atoms:
            dist = res_atom - lig_atom  # BioPython calculates Euclidean distance
            if dist < min_dist:
                min_dist = dist

    return min_dist


def calculate_min_sidechain_distance(residue, ligand):
    """
    Calculate minimum distance between sidechain atoms in a residue and any atom in a ligand.
    Excludes backbone atoms (N, CA, C, O) except for Glycine, which has no sidechain.
    For Glycine, uses all atoms (backbone only).

    Parameters:
    -----------
    residue : Bio.PDB.Residue
        The receptor residue
    ligand : Bio.PDB.Residue
        The ligand residue

    Returns:
    --------
    float : Minimum distance in Angstroms
    """
    res_name = residue.get_resname()

    # For Glycine, use all atoms (it has no sidechain)
    if res_name == 'GLY':
        res_atoms = get_residue_atoms(residue)
    else:
        # Backbone atom names to exclude
        backbone_atoms = {'N', 'CA', 'C', 'O'}
        # Get sidechain atoms (everything except backbone)
        res_atoms = [atom for atom in residue.get_atoms() if atom.get_name() not in backbone_atoms]

    lig_atoms = get_residue_atoms(ligand)

    if not res_atoms or not lig_atoms:
        return None

    min_dist = float('inf')
    for res_atom in res_atoms:
        for lig_atom in lig_atoms:
            dist = res_atom - lig_atom
            if dist < min_dist:
                min_dist = dist

    return min_dist


def identify_chains(structure):
    """
    Identify receptor, ligand, and G-protein chains in the structure.

    Returns:
    --------
    dict with keys: 'receptor_chains', 'ligand_residues', 'gprotein_chains'
    """
    model = structure[0]

    receptor_chains = []
    gprotein_chains = []
    ligand_residues = []

    for chain in model:
        chain_id = chain.get_id()
        residues = list(chain)

        if not residues:
            continue

        # Count amino acids vs hetero atoms
        aa_count = sum(1 for r in residues if is_aa(r, standard=True))
        hetero_count = sum(1 for r in residues if r.id[0].startswith('H_'))

        # If mostly amino acids and long enough, it's likely receptor or G-protein
        if aa_count > 50:
            # Typically receptor is chain A or R, but we'll categorize by length
            # GPCR receptors are usually ~300-400 residues
            # G-proteins have multiple chains (alpha ~350, beta ~340, gamma ~70)
            if aa_count > 200 and aa_count < 500:
                receptor_chains.append(chain_id)
            else:
                gprotein_chains.append(chain_id)

        # Identify ligands (HETATM residues that are not water)
        for residue in residues:
            res_id = residue.id
            res_name = residue.get_resname()

            # Ligands are HETATM records, not water or common ions
            if res_id[0].startswith('H_') and res_name not in ['HOH', 'WAT', 'NA', 'CL', 'MG', 'K', 'CA', 'ZN']:
                # Further filter: ligands usually have 10-100 atoms
                atom_count = len(list(residue.get_atoms()))
                if 5 < atom_count < 200:
                    ligand_residues.append((chain_id, residue))

    return {
        'receptor_chains': receptor_chains,
        'ligand_residues': ligand_residues,
        'gprotein_chains': gprotein_chains
    }


def determine_species(structure, receptor_chains):
    """
    Determine if structure is mouse or human based on conserved ASP position.
    Mouse has ASP147, Human has ASP149.

    Returns:
    --------
    str: 'mouse', 'human', or 'unknown'
    int: offset to convert to human numbering (mouse = +2, human = 0)
    """
    model = structure[0]
    has_asp147 = False
    has_asp149 = False

    for chain_id in receptor_chains:
        chain = model[chain_id]
        for residue in chain:
            if residue.get_resname() == 'ASP':
                res_num = residue.id[1]
                if res_num == 147:
                    has_asp147 = True
                elif res_num == 149:
                    has_asp149 = True

    if has_asp147 and not has_asp149:
        return 'mouse', 2
    elif has_asp149 and not has_asp147:
        return 'human', 0
    else:
        return 'unknown', 0


def analyze_structure(cif_file, output_dir):
    """
    Analyze a single GPCR structure file.

    Parameters:
    -----------
    cif_file : str
        Path to the mmCIF file
    output_dir : str
        Directory to save output CSV files
    """
    parser = MMCIFParser(QUIET=True)
    structure_name = Path(cif_file).stem

    print(f"\nAnalyzing {structure_name}...")

    try:
        structure = parser.get_structure(structure_name, cif_file)
    except Exception as e:
        print(f"  Error parsing {structure_name}: {e}")
        return None

    # Identify components
    components = identify_chains(structure)

    receptor_chains = components['receptor_chains']
    ligand_residues = components['ligand_residues']
    gprotein_chains = components['gprotein_chains']

    print(f"  Receptor chains: {receptor_chains}")
    print(f"  G-protein chains: {gprotein_chains}")
    print(f"  Ligands found: {len(ligand_residues)}")

    if not receptor_chains:
        print(f"  WARNING: No receptor chain identified in {structure_name}")
        return None

    if not ligand_residues:
        print(f"  WARNING: No ligand found in {structure_name}")
        return None

    # Determine species
    species, human_offset = determine_species(structure, receptor_chains)
    print(f"  Species: {species.upper()} (offset to human: +{human_offset})")

    # Calculate distances from ALL receptor chains to ligands
    model = structure[0]
    results = []

    for lig_chain_id, ligand in ligand_residues:
        lig_name = ligand.get_resname()
        lig_id = ligand.id[1]

        print(f"  Calculating distances for ligand: {lig_name} (chain {lig_chain_id}, residue {lig_id})")

        # Check all receptor chains for this ligand
        for receptor_chain_id in receptor_chains:
            receptor_chain = model[receptor_chain_id]

            for residue in receptor_chain:
                if not is_aa(residue, standard=True):
                    continue

                res_num = residue.id[1]
                res_name = residue.get_resname()

                # Calculate both all-atom and sidechain distances
                min_dist = calculate_min_distance(residue, ligand)
                min_sidechain_dist = calculate_min_sidechain_distance(residue, ligand)

                if min_dist is not None:
                    results.append({
                        'pdb_id': structure_name,
                        'species': species,
                        'receptor_chain': receptor_chain_id,
                        'residue_number': res_num,
                        'residue_number_human': res_num + human_offset,
                        'residue_name': res_name,
                        'ligand_name': lig_name,
                        'ligand_common_name': get_ligand_common_name(lig_name),
                        'ligand_chain': lig_chain_id,
                        'ligand_resid': lig_id,
                        'distance_angstrom': round(min_dist, 3),
                        'distance_sidechain_angstrom': round(min_sidechain_dist, 3) if min_sidechain_dist is not None else None,
                        'gprotein_chains': ','.join(gprotein_chains)
                    })

    if results:
        # Save individual CSV
        output_file = os.path.join(output_dir, f"{structure_name}_distances.csv")
        with open(output_file, 'w', newline='') as f:
            fieldnames = ['pdb_id', 'species', 'receptor_chain', 'residue_number',
                         'residue_number_human', 'residue_name', 'ligand_name',
                         'ligand_common_name', 'ligand_chain', 'ligand_resid',
                         'distance_angstrom', 'distance_sidechain_angstrom',
                         'gprotein_chains']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        print(f"  Saved {len(results)} distance measurements to {output_file}")

    return results


def generate_summary_report(all_results, output_dir):
    """
    Generate a summary report of all structures analyzed.
    """
    if not all_results:
        print("No results to summarize.")
        return

    # Combine all results
    combined_file = os.path.join(output_dir, "all_structures_combined.csv")
    with open(combined_file, 'w', newline='') as f:
        fieldnames = ['pdb_id', 'species', 'receptor_chain', 'residue_number',
                     'residue_number_human', 'residue_name', 'ligand_name',
                     'ligand_common_name', 'ligand_chain', 'ligand_resid',
                     'distance_angstrom', 'distance_sidechain_angstrom',
                     'gprotein_chains']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for results in all_results:
            writer.writerows(results)

    print(f"\nCombined results saved to {combined_file}")

    # Generate summary statistics
    summary_file = os.path.join(output_dir, "summary_report.txt")
    with open(summary_file, 'w') as f:
        f.write("GPCR Ligand Distance Analysis - Summary Report\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Total structures analyzed: {len(all_results)}\n")

        total_measurements = sum(len(r) for r in all_results)
        f.write(f"Total distance measurements: {total_measurements}\n\n")

        # Per-structure summary
        f.write("Per-Structure Summary:\n")
        f.write("-" * 80 + "\n")

        # Count species
        mouse_count = sum(1 for r in all_results if r and r[0]['species'] == 'mouse')
        human_count = sum(1 for r in all_results if r and r[0]['species'] == 'human')

        f.write(f"Mouse structures: {mouse_count}\n")
        f.write(f"Human structures: {human_count}\n\n")

        for results in all_results:
            if not results:
                continue

            pdb_id = results[0]['pdb_id']
            species = results[0]['species']
            ligand_name = results[0]['ligand_name']
            gproteins = results[0]['gprotein_chains']

            distances = [r['distance_angstrom'] for r in results]
            min_dist = min(distances)
            max_dist = max(distances)
            mean_dist = np.mean(distances)

            # Find residues within 5 Angstroms (binding site)
            binding_site = [r for r in results if r['distance_angstrom'] <= 5.0]

            f.write(f"\n{pdb_id}:\n")
            f.write(f"  Species: {species.upper()}\n")
            f.write(f"  Ligand: {ligand_name}\n")
            f.write(f"  G-protein chains: {gproteins}\n")
            f.write(f"  Distance range: {min_dist:.2f} - {max_dist:.2f} Å\n")
            f.write(f"  Mean distance: {mean_dist:.2f} Å\n")
            f.write(f"  Residues in binding site (≤5Å): {len(binding_site)}\n")

            if binding_site:
                f.write(f"  Binding site residues (human numbering): ")
                binding_res = sorted([(r['residue_name'], r['residue_number_human'], r['distance_angstrom'])
                                     for r in binding_site], key=lambda x: x[2])
                f.write(", ".join([f"{r[0]}{r[1]}({r[2]:.1f}Å)" for r in binding_res[:10]]))
                if len(binding_site) > 10:
                    f.write(f" ... and {len(binding_site)-10} more")
                f.write("\n")

    print(f"Summary report saved to {summary_file}")


def main():
    """Main analysis pipeline."""
    # Create output directory
    output_dir = str(REPO_ROOT / "structures" / "processed" / "ligand_distances_all")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Output directory: {output_dir}")

    # Find all CIF files
    cif_files = sorted([f for f in os.listdir(str(CIF_DIR)) if f.endswith('.cif')])

    print(f"Found {len(cif_files)} CIF files to analyze\n")

    all_results = []

    for cif_file in cif_files:
        cif_path = os.path.join(str(CIF_DIR), cif_file)
        results = analyze_structure(cif_path, output_dir)
        if results:
            all_results.append(results)

    # Generate summary report
    generate_summary_report(all_results, output_dir)

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print(f"Results saved in: {output_dir}")


if __name__ == "__main__":
    main()
