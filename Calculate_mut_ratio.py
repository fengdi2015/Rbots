# Calculate the ratio of S and NS for each Sample ID and Gene Symbol

# Count occurrences of S and NS per Sample ID and Gene Symbol
mutation_counts = df_updated.pivot_table(index=["Sample ID", "Gene Symbol"], 
                                         columns="Mutation Type", 
                                         aggfunc="size", fill_value=0)

# Ensure both S and NS columns exist
if "S" not in mutation_counts.columns:
    mutation_counts["S"] = 0
if "NS" not in mutation_counts.columns:
    mutation_counts["NS"] = 0

# Calculate ratio of S to NS
mutation_counts["S/NS Ratio"] = mutation_counts["S"] / (mutation_counts["NS"] + 1e-6)  # Avoid division by zero

# Reset index for display
mutation_ratios = mutation_counts.reset_index()

# Display DataFrame
tools.display_dataframe_to_user(name="Mutation Type Ratios", dataframe=mutation_ratios)

# Re-import necessary libraries after execution state reset
import pandas as pd
import random
import ace_tools as tools

# Define gene symbols
genes = ["TP53", "BRCA1", "EGFR", "KRAS", "PTEN"]

# Create mock data
updated_data = []
sample_id = 1

# Generate multiple mutations per sample
for gene in genes:
    num_samples = random.randint(3, 6)  # Each gene has 3-6 mutations
    for _ in range(num_samples):
        sample_label = f"Sample_{sample_id}"
        mutation_types = random.choices(["S", "NS"], k=random.randint(1, 3))  # Each sample can have multiple mutation types
        for mutation in mutation_types:
            updated_data.append([sample_label, gene, mutation])
        sample_id += 1

# Create DataFrame
df_updated = pd.DataFrame(updated_data, columns=["Sample ID", "Gene Symbol", "Mutation Type"])

# Count occurrences of S and NS per Sample ID and Gene Symbol
mutation_counts = df_updated.pivot_table(index=["Sample ID", "Gene Symbol"], 
                                         columns="Mutation Type", 
                                         aggfunc="size", fill_value=0)

# Ensure both S and NS columns exist
if "S" not in mutation_counts.columns:
    mutation_counts["S"] = 0
if "NS" not in mutation_counts.columns:
    mutation_counts["NS"] = 0

# Calculate ratio of S to NS
mutation_counts["S/NS Ratio"] = mutation_counts["S"] / (mutation_counts["NS"] + 1e-6)  # Avoid division by zero

# Reset index for display
mutation_ratios = mutation_counts.reset_index()

# Display DataFrame
tools.display_dataframe_to_user(name="Mutation Type Ratios", dataframe=mutation_ratios)
