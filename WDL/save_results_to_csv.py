import sys
import csv

# Get the disease-gene associations from command-line arguments
disease_gene_associations = sys.argv[1:]

# Output CSV file name
output_file = "disease_gene_associations.csv"

# Function to save results to a CSV file
def save_results_to_csv(disease_gene_associations, output_file):
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Disease-Gene Associations"])
        for association in disease_gene_associations:
            writer.writerow([association])

# Save results to CSV
save_results_to_csv(disease_gene_associations, output_file)

# Confirm the results were saved
print(f"Results saved to {output_file}")
