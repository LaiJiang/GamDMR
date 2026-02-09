#collect thhe misssing taks IDs from stage 1
#output: 8_col_stage1.txt containg all IDs need to be run for stage 2 with more hours and stringent cpg selections.
import os
import re

# Directory containing the .RData files
DATA_DIR_stage1 = "/home/laj773/scratch/UQAC/meth/results/somnibus/"

DATA_DIR_stage2 = "/home/laj773/scratch/UQAC/meth/results/somnibus_stage2/"

# Pattern to extract four-digit numbers
pattern = re.compile(r"_(\d{4})_")

# Collect all numbers present in filenames (list may contain duplicates)
stage1_runs = []
for fname in os.listdir(DATA_DIR_stage1):
    if fname.endswith("_somnibus.RData") or fname.endswith("_somnbibus.RData"):
        match = pattern.search(fname)
        if match:
            stage1_runs.append(int(match.group(1)))

stage2_runs = []
for fname in os.listdir(DATA_DIR_stage2):
    if fname.endswith("_somnibus.RData") or fname.endswith("_somnbibus.RData"):
        match = pattern.search(fname)
        if match:
            stage2_runs.append(int(match.group(1)))



# Ensure uniqueness
present = set(stage1_runs + stage2_runs)

# Define full range of task IDs
all_ids = set(range(1, 1738))  # 1 through 1737 inclusive


# Compute missing IDs
missing = sorted(all_ids - present)


# Output file path
output_file = "/home/laj773/scratch/UQAC/meth/results/9_regional/10_col_stage2.txt"


# Write missing IDs to a text file
with open(output_file, "w") as out:
    for task_id in missing:
        out.write(f"{task_id}\n")

print(f"Found {len(missing)} missing IDs after stage 1 + stage 2 somnibus. Written to {output_file}.")
