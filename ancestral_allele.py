# this file takes in the intersected bed files containing the hal snps data and the snpeff annotated variants
import os
import pandas as pd
import sys

## Input file
input_file = sys.argv[1]

# Load the dataset

df = pd.read_csv(input_file, sep=",")  # Assuming tab-separated file, change if needed

# Define function to infer the ancestral allele
def infer_ancestral_allele(row):
    """
    Infers the ancestral allele using weasel, polecat, and steppe polecat data.
    - If weasel, polecat, and steppe polecat have the same allele, that is inferred as ancestral.
    - If weasel has a different allele than the others, assume weasel carries a unique mutation.
    """
    alleles = {row["Ferret"], row["Steppe"], row["Weasel"]}  # Unique alleles in outgroups
    alleles.discard("")  # Remove missing data if any

    if len(alleles) == 1:
        return list(alleles)[0]  # Single allele found, treat as ancestral
    return None  # Unclear ancestral state

# Apply function to infer ancestral alleles
df["Ancestral_Allele"] = df.apply(infer_ancestral_allele, axis=1)

# Identify derived mutations in ferrets
df["Derived_In_Polecat"] = df.apply(lambda row: row["Polecat"] != row["Ancestral_Allele"] if row["Ancestral_Allele"] else "Unclear", axis=1)

# Save the processed data
output_file = "filtered_variants_with_ancestral_info.csv"
df.to_csv(output_file, sep="\t", index=False)

print(f"Processed data saved to {output_file}")
