#!/bin/bash

# Set the directory where the chunk_XXXX.csv files are located
CHUNK_DIR="/home/laj773/scratch/UQAC/meth/data/meth_split/"

# Output file
OUTFILE="/home/laj773/scratch/UQAC/meth/scr/10_sanity/results/first_row_chunk.txt"



# Empty the output file if it already exists
> "$OUTFILE"

# Loop from 1 to 3000
for i in $(seq -w 1 1800); do
    FILE="$CHUNK_DIR/chunk_${i}.csv"
    if [ -f "$FILE" ]; then
        awk 'NR==2 {print $1 "\t" $2}' "$FILE" >> "$OUTFILE"
    else
        echo "Warning: $FILE does not exist." >&2
    fi
done
