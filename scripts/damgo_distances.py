import csv
import numpy as np

# Assuming "mor" is the object name for the first structure and "P" is the chain name for the second structure
object_name_mor = "chain R"
chain_name_p = "chain P"
cutoff_distance = 100.0  # Set your desired cutoff distance in angstroms

# Get the coordinates of all atoms in "mor"
coords_mor = np.array([atom.coord for atom in cmd.get_model(f"({object_name_mor})").atom])

# Get the coordinates of all atoms in chain P
coords_chain_p = np.array([atom.coord for atom in cmd.get_model(f"({chain_name_p})").atom])

# Calculate distances efficiently using NumPy broadcasting
distances = np.linalg.norm(coords_mor[:, np.newaxis, :] - coords_chain_p, axis=2)

# Find the minimum distance for each residue in "mor"
min_distances = np.min(distances, axis=1)

# Get the corresponding PyMOL residue numbers
residue_numbers_mor = [atom.resi for atom in cmd.get_model(f"({object_name_mor})").atom]

# Create a dictionary to store the valid shortest distances for each residue in "mor"
damage_distances = {str(resi): {"pos": str(resi), "damgo_distance": min_distance}
                    for resi, min_distance in zip(residue_numbers_mor, min_distances)
                    if min_distance < cutoff_distance}

# Save the results to a CSV file
csv_filename = "damgo_distances.csv"
with open(csv_filename, mode='w', newline='') as csv_file:
    fieldnames = ["pos", "damgo_distance"]
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

    # Write the header
    writer.writeheader()

    # Write the data
    for result in damage_distances.values():
        writer.writerow(result)

print(f"Results saved to {csv_filename}")
