#!/usr/bin/env python3
"""
Identify all residue contacts within 5 Angstroms of ligands for each MOR structure.

This script processes all CIF structure files and identifies residues that are within
5Å of any ligand atom, creating comprehensive contact maps for each structure.
"""

import os
import sys
import numpy as np
import pandas as pd
from collections import defaultdict
from pathlib import Path

try:
    from Bio.PDB import MMCIFParser
    from Bio.PDB.Polypeptide import is_aa
except ImportError:
    print("BioPython not found. Installing...")
    import subprocess
    subprocess.check_call(["pip", "install", "biopython"])
    from Bio.PDB import MMCIFParser
    from Bio.PDB.Polypeptide import is_aa


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
    float : Minimum distance in Angstroms, or None if calculation fails
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

    return min_dist if min_dist != float('inf') else None


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
    float : Minimum distance in Angstroms, or None if calculation fails
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

    return min_dist if min_dist != float('inf') else None


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
                # Further filter: ligands usually have 5-200 atoms
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


def load_gpcrdb_annotations(annotations_file):
    """
    Load GPCRdb annotations to add structural context to contacts.
    
    Parameters:
    -----------
    annotations_file : str
        Path to GPCRdb_OPRM1_table.csv
        
    Returns:
    --------
    dict : Mapping of human position to GPCRdb annotation
    """
    if not os.path.exists(annotations_file):
        print(f"Warning: GPCRdb annotations file not found: {annotations_file}")
        return {}
    
    try:
        df = pd.read_csv(annotations_file)
        # Create mapping: human position -> (GPCRdb, SSE, wt_aa)
        annotations = {}
        for _, row in df.iterrows():
            pos = row['pos']
            gpcrdb = row.get('GPCRdb', '')
            sse = row.get('SSE', '')
            wt_aa = row.get('wt_aa', '')
            annotations[pos] = {
                'GPCRdb': gpcrdb,
                'SSE': sse,
                'wt_aa': wt_aa
            }
        print(f"Loaded GPCRdb annotations for {len(annotations)} positions")
        return annotations
    except Exception as e:
        print(f"Warning: Could not load GPCRdb annotations: {e}")
        return {}


def analyze_structure_contacts(cif_file, output_dir, gpcrdb_annotations, contact_cutoff=5.0):
    """
    Analyze a single GPCR structure file and identify contacts within cutoff distance.

    Parameters:
    -----------
    cif_file : str
        Path to the mmCIF file
    output_dir : str
        Directory to save output CSV files
    gpcrdb_annotations : dict
        GPCRdb annotation mapping
    contact_cutoff : float
        Distance cutoff in Angstroms (default: 5.0)

    Returns:
    --------
    list : List of contact dictionaries
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

    # Calculate contacts from ALL receptor chains to ligands
    model = structure[0]
    contacts = []

    for lig_chain_id, ligand in ligand_residues:
        lig_name = ligand.get_resname()
        lig_id = ligand.id[1]

        print(f"  Identifying contacts for ligand: {lig_name} (chain {lig_chain_id}, residue {lig_id})")

        # Check all receptor chains for this ligand
        for receptor_chain_id in receptor_chains:
            receptor_chain = model[receptor_chain_id]

            for residue in receptor_chain:
                if not is_aa(residue, standard=True):
                    continue

                res_num = residue.id[1]
                res_num_human = res_num + human_offset
                res_name = residue.get_resname()

                # Calculate both all-atom and sidechain distances
                min_dist = calculate_min_distance(residue, ligand)
                min_sidechain_dist = calculate_min_sidechain_distance(residue, ligand)

                # Only include contacts within cutoff distance
                if min_dist is not None and min_dist <= contact_cutoff:
                    # Get GPCRdb annotations if available
                    annotation = gpcrdb_annotations.get(res_num_human, {})
                    
                    contact = {
                        'pdb_id': structure_name,
                        'species': species,
                        'receptor_chain': receptor_chain_id,
                        'residue_number': res_num,
                        'residue_number_human': res_num_human,
                        'residue_name': res_name,
                        'GPCRdb': annotation.get('GPCRdb', ''),
                        'SSE': annotation.get('SSE', ''),
                        'ligand_name': lig_name,
                        'ligand_chain': lig_chain_id,
                        'ligand_resid': lig_id,
                        'distance_angstrom': round(min_dist, 3),
                        'distance_sidechain_angstrom': round(min_sidechain_dist, 3) if min_sidechain_dist is not None else None,
                        'gprotein_chains': ','.join(gprotein_chains) if gprotein_chains else ''
                    }
                    contacts.append(contact)

    if contacts:
        # Sort by distance
        contacts.sort(key=lambda x: x['distance_angstrom'])
        
        # Save individual CSV
        output_file = os.path.join(output_dir, f"{structure_name}_contacts_5A.csv")
        df_contacts = pd.DataFrame(contacts)
        df_contacts.to_csv(output_file, index=False)

        print(f"  Found {len(contacts)} contacts within {contact_cutoff}Å")
        print(f"  Saved to {output_file}")

    return contacts


def generate_summary_report(all_contacts, output_dir, contact_cutoff=5.0):
    """
    Generate a summary report of all contacts identified.
    """
    if not all_contacts:
        print("No contacts to summarize.")
        return

    # Combine all contacts
    combined_file = os.path.join(output_dir, f"all_structures_contacts_{contact_cutoff}A.csv")
    df_all = pd.DataFrame(all_contacts)
    df_all.to_csv(combined_file, index=False)
    print(f"\nCombined contacts saved to {combined_file}")

    # Generate summary statistics
    summary_file = os.path.join(output_dir, f"contacts_summary_{contact_cutoff}A.txt")
    with open(summary_file, 'w') as f:
        f.write("MOR Ligand Contact Analysis - Summary Report\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Contact cutoff distance: {contact_cutoff} Angstroms\n\n")

        f.write(f"Total structures analyzed: {len(set(c['pdb_id'] for c in all_contacts))}\n")
        f.write(f"Total contacts identified: {len(all_contacts)}\n\n")

        # Per-structure summary
        f.write("Per-Structure Summary:\n")
        f.write("-" * 80 + "\n")

        # Group by structure
        by_structure = defaultdict(list)
        for contact in all_contacts:
            by_structure[contact['pdb_id']].append(contact)

        # Count species
        mouse_count = sum(1 for contacts in by_structure.values() 
                         if contacts and contacts[0]['species'] == 'mouse')
        human_count = sum(1 for contacts in by_structure.values()
                         if contacts and contacts[0]['species'] == 'human')

        f.write(f"Mouse structures: {mouse_count}\n")
        f.write(f"Human structures: {human_count}\n\n")

        for pdb_id in sorted(by_structure.keys()):
            contacts = by_structure[pdb_id]
            if not contacts:
                continue

            species = contacts[0]['species']
            ligand_name = contacts[0]['ligand_name']
            gproteins = contacts[0]['gprotein_chains']

            distances = [c['distance_angstrom'] for c in contacts]
            min_dist = min(distances)
            max_dist = max(distances)
            mean_dist = np.mean(distances)

            # Get unique positions
            unique_positions = sorted(set(c['residue_number_human'] for c in contacts))

            f.write(f"\n{pdb_id}:\n")
            f.write(f"  Species: {species.upper()}\n")
            f.write(f"  Ligand: {ligand_name}\n")
            f.write(f"  G-protein chains: {gproteins if gproteins else 'None'}\n")
            f.write(f"  Total contacts: {len(contacts)}\n")
            f.write(f"  Unique positions: {len(unique_positions)}\n")
            f.write(f"  Distance range: {min_dist:.2f} - {max_dist:.2f} Å\n")
            f.write(f"  Mean distance: {mean_dist:.2f} Å\n")
            
            # List top 10 closest contacts
            f.write(f"  Closest contacts (human numbering):\n")
            for i, contact in enumerate(contacts[:10], 1):
                pos = contact['residue_number_human']
                res = contact['residue_name']
                dist = contact['distance_angstrom']
                sse = contact.get('SSE', 'N/A')
                f.write(f"    {i:2d}. {res}{pos} ({sse}): {dist:.2f} Å\n")
            if len(contacts) > 10:
                f.write(f"    ... and {len(contacts)-10} more contacts\n")

        # Position frequency analysis
        f.write("\n\nPosition Frequency Analysis:\n")
        f.write("-" * 80 + "\n")
        f.write("Residues that contact ligands in multiple structures:\n\n")

        # Count how many structures each position appears in
        position_counts = defaultdict(set)
        for contact in all_contacts:
            pos = contact['residue_number_human']
            pdb = contact['pdb_id']
            position_counts[pos].add(pdb)

        # Sort by frequency
        position_freq = [(pos, len(structures), sorted(structures)) 
                         for pos, structures in position_counts.items()]
        position_freq.sort(key=lambda x: x[1], reverse=True)

        f.write(f"Total unique positions with contacts: {len(position_freq)}\n\n")
        f.write("Top 20 most frequently contacted positions:\n")
        for i, (pos, count, structures) in enumerate(position_freq[:20], 1):
            # Get residue name from first contact
            sample_contact = next(c for c in all_contacts if c['residue_number_human'] == pos)
            res_name = sample_contact['residue_name']
            sse = sample_contact.get('SSE', 'N/A')
            gpcrdb = sample_contact.get('GPCRdb', '')
            f.write(f"  {i:2d}. {res_name}{pos} ({sse}, {gpcrdb}): {count} structures\n")

    print(f"Summary report saved to {summary_file}")


def main():
    """Main analysis pipeline."""
    # Get script directory and set paths
    script_dir = Path(__file__).parent

    # Path to GPCRdb annotations (relative to repo root)
    repo_root = Path(__file__).resolve().parents[3]
    annotations_file = repo_root / "data" / "external" / "annotations" / "GPCRdb_OPRM1_table.csv"

    # Create output directory
    output_dir = repo_root / "data" / "processed" / "structural" / "ligand_contacts"
    output_dir.mkdir(exist_ok=True)

    print(f"Output directory: {output_dir}")
    print(f"GPCRdb annotations: {annotations_file}")

    # Load GPCRdb annotations
    gpcrdb_annotations = load_gpcrdb_annotations(str(annotations_file))

    # Find all CIF files
    cif_files = sorted([f for f in (repo_root / "data" / "raw" / "structural" / "pdb").glob("*.cif")])
    
    if not cif_files:
        print(f"ERROR: No CIF files found in {repo_root / 'data' / 'raw' / 'structural' / 'pdb'}")
        print("Please ensure CIF files are present in the structural/pdb directory")
        sys.exit(1)

    print(f"\nFound {len(cif_files)} CIF files to analyze")

    all_contacts = []
    contact_cutoff = 5.0

    for cif_file in cif_files:
        contacts = analyze_structure_contacts(
            str(cif_file), 
            str(output_dir), 
            gpcrdb_annotations,
            contact_cutoff=contact_cutoff
        )
        if contacts:
            all_contacts.extend(contacts)

    # Generate summary report
    if all_contacts:
        generate_summary_report(all_contacts, str(output_dir), contact_cutoff=contact_cutoff)

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print(f"Results saved in: {output_dir}")
    print(f"\nTotal contacts identified: {len(all_contacts)}")
    print(f"Unique structures: {len(set(c['pdb_id'] for c in all_contacts))}")
    print(f"Unique positions: {len(set(c['residue_number_human'] for c in all_contacts))}")


if __name__ == "__main__":
    main()
