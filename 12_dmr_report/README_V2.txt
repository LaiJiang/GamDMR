CHUNK-BUNDLED PARALLEL DMR CHARACTERISTICS PIPELINE — V2
=========================================================

Why V2?
-------
Rorqual counts every job-array task against the submitted-job limit.  The
prepared DMR manifest contains 1733 source methylation chunks.  Therefore:

    --array=1-1733%100

still represents 1733 submitted tasks and is rejected with:

    AssocMaxSubmitJobLimit

The "%100" only limits the number RUNNING simultaneously; it does not reduce
the number of submitted/pending array tasks.

V2 bundles several source chunks into each SLURM array task.

Default:
    CHUNKS_PER_TASK=10

For 1733 chunks:
    ceil(1733 / 10) = 174 array tasks

Each array task processes its 10 chunks sequentially.  Each source chunk is
still read exactly once, and one partial RDS is still written per source chunk.
The final summarizer is unchanged.

Files
-----
2_prepare_DMR_chunk_manifest.R       unchanged
2_submit_DMR_pipeline.sh             UPDATED
3_extract_DMR_statistics_array.R     UPDATED
3_extract_DMR_statistics_array.sh    UPDATED
4_summarize_DMR_statistics.R         unchanged
4_summarize_DMR_statistics.sh        unchanged

Recommended run
---------------
Replace the V1 scripts with these V2 versions, then:

    cd ~/scratch/UQAC/meth/scr/15_revision/9_dmr_report

    chmod +x 2_submit_DMR_pipeline.sh
    chmod +x 3_extract_DMR_statistics_array.sh
    chmod +x 4_summarize_DMR_statistics.sh

    ./2_submit_DMR_pipeline.sh

With the current 1733-chunk manifest the launcher should print approximately:

    Source methylation chunks:    1733
    Chunks per array task:        10
    Total submitted array tasks:  174
    Maximum concurrently running: 100

Alternative settings
--------------------
More conservative job count:

    CHUNKS_PER_TASK=20 ./2_submit_DMR_pipeline.sh

This gives about 87 array tasks.

Lower filesystem concurrency:

    CHUNKS_PER_TASK=10 MAX_CONCURRENT=50 ./2_submit_DMR_pipeline.sh

Both:

    CHUNKS_PER_TASK=20 MAX_CONCURRENT=50 ./2_submit_DMR_pipeline.sh

Resources
---------
Each bundled array task:
    1 CPU
    6 GB RAM
    4 hours

Chunks are processed sequentially, so bundling increases wall time per array
task but does not multiply memory use.

Final collector:
    1 CPU
    16 GB RAM
    1 hour

Monitoring
----------
    squeue -u $USER

After completion:
    cat ~/scratch/UQAC/meth/results/15_revision/9_dmr_report/DMR_statistics/DMR_characteristics_report.txt
