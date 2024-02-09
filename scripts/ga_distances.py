import csv
import numpy as np

# Assuming "mor" and "alpha" are the object names
object_name_mor = "chain R"
object_name_alpha = "chain A"
cutoff_distance = 100.0  # Set your desired cutoff distance in angstroms

# Get the coordinates of CA atoms in "mor"
coords_mor = np.array([atom.coord for atom in cmd.get_model(f"({object_name_mor} and name CA)").atom])

# Get the coordinates of CA atoms in "alpha"
coords_alpha = np.array([atom.coord for atom in cmd.get_model(f"({object_name_alpha} and name CA)").atom])

# Calculate distances efficiently using NumPy broadcasting
distances = np.linalg.norm(coords_mor[:, np.newaxis, :] - coords_alpha, axis=2)

# Find the minimum distance for each residue
min_distances = np.min(distances, axis=1)

# Identify residues with distances below the cutoff
valid_residues = np.where(min_distances < cutoff_distance)[0]

# Create a dictionary to store the valid shortest distances for each residue in "mor"
shortest_distances = {str(resi): {"pos": str(resi), "Ga_distance": min_distances[resi]} for resi in valid_residues}

# Save the results to a CSV file
csv_filename = "ga_distances.csv"
with open(csv_filename, mode='w', newline='') as csv_file:
    fieldnames = ["pos", "Ga_distance"]
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

    # Write the header
    writer.writeheader()

    # Write the data
    for result in shortest_distances.values():
        writer.writerow(result)

print(f"Results saved to {csv_filename}")
