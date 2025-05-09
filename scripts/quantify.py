import pandas as pd
import os
from snakemake.shell import shell

# Aggregate all count files into one matrix
count_files = snakemake.input
output_file = snakemake.output[0]

# Read first file to get gene IDs
first_file = count_files[0]
count_data = pd.read_csv(first_file, sep="\t", comment="#", index_col=0)
combined_counts = count_data.iloc[:, -1].to_frame()
combined_counts.columns = [os.path.basename(first_file).replace(".counts.txt", "")]

# Add remaining samples
for f in count_files[1:]:
    sample_name = os.path.basename(f).replace(".counts.txt", "")
    count_data = pd.read_csv(f, sep="\t", comment="#", index_col=0)
    combined_counts[sample_name] = count_data.iloc[:, -1]

# Save combined counts
combined_counts.to_csv(output_file, sep="\t")
