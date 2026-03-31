from pathlib import Path
import pandas as pd

INPUT_FILE = Path("data/peek_methylation_sva.csv")
OUTPUT_FILE = Path("data/peek_methylation_sva_var_filtered.tsv")
VAR_THRESHOLD = 0.035
CHUNK_SIZE = 10_000
METH_SUFFIX = "_meth"


def filter_rows_by_variance(
    input_file: Path,
    output_file: Path,
    var_threshold: float = VAR_THRESHOLD,
    chunk_size: int = CHUNK_SIZE,
) -> None:
    first_chunk = True

    for chunk in pd.read_csv(input_file, sep="\t", chunksize=chunk_size):
        meth_cols = [col for col in chunk.columns if col.endswith(METH_SUFFIX)]
        if not meth_cols:
            raise ValueError(f"No columns ending with '{METH_SUFFIX}' were found.")

        filtered_chunk = chunk.loc[
            chunk[meth_cols].var(axis=1, skipna=True) >= var_threshold
        ]

        filtered_chunk.to_csv(
            output_file,
            sep="\t",
            index=False,
            mode="w" if first_chunk else "a",
            header=first_chunk,
        )
        first_chunk = False

    print(f"Filtered data saved to: {output_file}")


if __name__ == "__main__":
    filter_rows_by_variance(INPUT_FILE, OUTPUT_FILE)
