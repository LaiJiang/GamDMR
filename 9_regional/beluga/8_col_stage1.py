#collect thhe misssing taks IDs from stage 1
#output: 8_col_stage1.txt containg all IDs need to be run for stage 2 with more hours and stringent cpg selections.
import os
import re

# Directory containing the .RData files
DATA_DIR = "/home/laj773/scratch/UQAC/meth/results/somnibus/"

# Pattern to extract four-digit numbers
pattern = re.compile(r"_(\d{4})_")

# Collect all numbers present in filenames (list may contain duplicates)
present_list = []
for fname in os.listdir(DATA_DIR):
    if fname.endswith("_somnibus.RData") or fname.endswith("_somnbibus.RData"):
        match = pattern.search(fname)
        if match:
            present_list.append(int(match.group(1)))


# Ensure uniqueness
present = set(present_list)

# Define full range of task IDs
all_ids = set(range(1, 1738))  # 1 through 1737 inclusive

# Compute missing IDs
missing = sorted(all_ids - present)

# Write missing IDs to a text file, one per line
with open("/home/laj773/scratch/UQAC/meth/results/9_regional/8_col_stage1.txt", "w") as out:
    for task_id in missing:
        out.write(f"{task_id}\n")

print(f"Found {len(missing)} missing IDs. Written to missing_ids.txt.")
