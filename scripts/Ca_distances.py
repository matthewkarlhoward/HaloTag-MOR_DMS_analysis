from Bio.PDB import PDBParser

def calculate_distance(pdb_file1, pdb_file2):
    parser = PDBParser(QUIET=True)

    # Parse PDB structures
    structure1 = parser.get_structure('inactive', 4dkl)
    structure2 = parser.get_structure('active', 5c1m)

    distances = []

    # Iterate through chains, residues, and atoms
    for model1, model2 in zip(structure1, structure2):
        for chain1, chain2 in zip(model1, model2):
            for residue1, residue2 in zip(chain1, chain2):
                if 'CA' in residue1 and 'CA' in residue2:
                    atom1 = residue1['CA']
                    atom2 = residue2['CA']

                    # Calculate Euclidean distance between Cα atoms
                    distance = atom1 - atom2
                    distances.append(distance)

    return distances

# Example usage
pdb_file1 = '4dkl.pdb'
pdb_file2 = '5c1m.pdb'

distances = calculate_distance(pdb_file1, pdb_file2)

print("Distances between Cα atoms:")
for distance in distances:
    print(distance)
