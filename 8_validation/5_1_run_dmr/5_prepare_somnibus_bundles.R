#!/usr/bin/env Rscript

#!/usr/bin/env Rscript

user_lib <- path.expand("~/scratch/R/library")
.libPaths(c(user_lib, .libPaths()))

cat("R library paths:\n")
print(.libPaths())

suppressPackageStartupMessages({
  library(data.table)
  library(SOMNiBUS)
})

cat(
  "SOMNiBUS version:",
  as.character(packageVersion("SOMNiBUS")),
  "\n"
)

# --------------------------- configuration -----------------------------------
PATH_wk <- path.expand("~/scratch/UQAC/meth/")

region_file_path <- file.path(
  PATH_wk, "results/15_revision/5_cv/region_somnibus_M12.csv"
)

pheno_rds_path <- file.path(
  PATH_wk, "scr/6_whole/pheno_file.rds"
)

chunk_dir <- file.path(PATH_wk, "data/meth_split")
output_root <- file.path(
  PATH_wk, "results/15_revision/5_cv/2_somnibus/prep"
)
bundle_dir <- file.path(output_root, "bundles")
manifest_dir <- file.path(output_root, "manifests")

N_jobs <- as.integer(Sys.getenv("N_JOBS", "900"))
gap_bp <- as.integer(Sys.getenv("SOMNIBUS_GAP_BP", "250"))
min_cpgs <- as.integer(Sys.getenv("SOMNIBUS_MIN_CPGS", "51"))
max_cpgs <- as.integer(Sys.getenv("SOMNIBUS_MAX_CPGS", "2000"))
overwrite <- identical(Sys.getenv("OVERWRITE", "0"), "1")

# Match the previous SOMNiBUS analysis.
selected_covariates <- c("AA_only")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Pass SLURM_ARRAY_TASK_ID.")
job_id <- suppressWarnings(as.integer(args[1]))
if (is.na(job_id) || job_id < 1L || job_id > N_jobs) {
  stop("Invalid SLURM_ARRAY_TASK_ID.")
}

dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

bundle_file <- file.path(
  bundle_dir, sprintf("somnibus_bundle_job_%04d.rds", job_id)
)
manifest_file <- file.path(
  manifest_dir, sprintf("somnibus_manifest_job_%04d.tsv", job_id)
)

if (!overwrite && (file.exists(bundle_file) || file.exists(manifest_file))) {
  stop("Output already exists. Set OVERWRITE=1 to replace it.")
}

normalize_id <- function(x) {
  x <- gsub("\\.", "-", as.character(x))
  sub("^X", "", x)
}

# Whole chunks remain together, so a chunk is read by one preparation job only.
assign_chunks <- function(region_dt, n_jobs) {
  load_dt <- region_dt[
    ,
    .(
      load = sum(pmax(fifelse(is.na(n_cpgs), 1L, n_cpgs), 1L)),
      n_parent_regions = .N
    ),
    by = data_chunk_id
  ]
  setorder(load_dt, -load, data_chunk_id)

  job_load <- numeric(n_jobs)
  load_dt[, prep_job_id := NA_integer_]

  for (i in seq_len(nrow(load_dt))) {
    j <- which.min(job_load)
    load_dt$prep_job_id[i] <- j
    job_load[j] <- job_load[j] + load_dt$load[i]
  }

  setorder(load_dt, data_chunk_id)
  load_dt
}

# Convert all selected CpGs from one chunk in one data.table::melt call.
convert_chunk_long <- function(wide_dt, pheno_dt) {
  meth_cols <- grep("_meth$", names(wide_dt), value = TRUE)
  if (length(meth_cols) == 0L) stop("No *_meth columns.")

  raw_ids <- sub("_meth$", "", meth_cols)
  tot_cols <- paste0(raw_ids, "_tot")
  missing_tot <- setdiff(tot_cols, names(wide_dt))
  if (length(missing_tot) > 0L) {
    stop("Missing matching *_tot columns: ",
         paste(head(missing_tot, 10L), collapse = ", "))
  }

  keep_cols <- c("chr", "start", meth_cols, tot_cols)
  wide_dt <- wide_dt[, ..keep_cols]

  long_dt <- melt(
    wide_dt,
    id.vars = c("chr", "start"),
    measure.vars = list(meth_cols, tot_cols),
    variable.name = "sample_index",
    value.name = c("meth_prop", "Total_Counts"),
    variable.factor = FALSE,
    na.rm = FALSE
  )

  long_dt[, ID := normalize_id(raw_ids[as.integer(sample_index)])]
  long_dt[, sample_index := NULL]
  long_dt[, meth_prop := suppressWarnings(as.numeric(meth_prop))]
  long_dt[, Total_Counts := suppressWarnings(as.numeric(Total_Counts))]

  long_dt <- long_dt[
    is.finite(meth_prop) &
      is.finite(Total_Counts) &
      Total_Counts > 0 &
      meth_prop >= 0 &
      meth_prop <= 1
  ]

  long_dt[
    ,
    `:=`(
      Meth_Counts = as.numeric(round(meth_prop * Total_Counts)),
      Position = as.integer(start)
    )
  ]

  long_dt <- long_dt[
    Meth_Counts >= 0 & Meth_Counts <= Total_Counts
  ]

  long_dt[, c("meth_prop", "start") := NULL]

  merge(long_dt, pheno_dt, by = "ID", all = FALSE, sort = FALSE)
}

manifest_row <- function(
  status, parent, chunk_id, message = NA_character_,
  bundle_index = NA_integer_, omnibus_region_id = NA_character_,
  child_index = NA_integer_, region_start = NA_integer_,
  region_end = NA_integer_, n_cpgs_out = NA_integer_,
  n_rows = NA_integer_, n_samples = NA_integer_, n_fids = NA_integer_
) {
  data.table(
    prep_job_id = job_id,
    bundle_file = bundle_file,
    bundle_index = bundle_index,
    omnibus_region_id = omnibus_region_id,
    data_chunk_id = as.integer(chunk_id),
    parent_region_id = as.integer(parent$region_id),
    chr = as.integer(parent$chr),
    parent_start = as.integer(parent$region_start),
    parent_end = as.integer(parent$region_end),
    parent_n_cpgs = as.integer(parent$n_cpgs),
    child_index = child_index,
    region_start = region_start,
    region_end = region_end,
    n_cpgs = n_cpgs_out,
    n_rows = n_rows,
    n_samples = n_samples,
    n_fids = n_fids,
    status = status,
    message = message
  )
}

# --------------------------- inputs ------------------------------------------
region_dt <- fread(region_file_path)
required_region <- c(
  "data_chunk_id", "chr", "region_start", "region_end",
  "n_cpgs", "region_id"
)
missing_region <- setdiff(required_region, names(region_dt))
if (length(missing_region) > 0L) {
  stop("Region file is missing: ", paste(missing_region, collapse = ", "))
}

region_dt <- region_dt[
  complete.cases(region_dt[, ..required_region])
]
region_dt[
  ,
  (required_region) := lapply(.SD, as.integer),
  .SDcols = required_region
]
if (anyDuplicated(region_dt$region_id)) stop("region_id must be unique.")
if (any(region_dt$region_start > region_dt$region_end)) {
  stop("At least one region has start > end.")
}
setorder(region_dt, data_chunk_id, region_start, region_end, region_id)

chunk_assignment <- assign_chunks(region_dt, N_jobs)
assigned_chunks <- chunk_assignment[
  prep_job_id == job_id, data_chunk_id
]

if (job_id == 1L) {
  fwrite(
    chunk_assignment,
    file.path(output_root, "somnibus_chunk_job_assignment.tsv"),
    sep = "\t"
  )
}

pheno_dt <- as.data.table(
  copy(readRDS(pheno_rds_path))
)

###########
required_pheno <- c("ID", "FID", selected_covariates)
missing_pheno <- setdiff(required_pheno, names(pheno_dt))
if (length(missing_pheno) > 0L) {
  stop("pheno_file is missing: ", paste(missing_pheno, collapse = ", "))
}

pheno_dt <- pheno_dt[, ..required_pheno]
pheno_dt[, ID := normalize_id(ID)]
pheno_dt[, FID := as.character(FID)]
pheno_dt[, AA_only := suppressWarnings(as.integer(as.character(AA_only)))]

if (any(!is.na(pheno_dt$AA_only) & !pheno_dt$AA_only %in% c(0L, 1L))) {
  stop("AA_only must contain only 0, 1, or NA.")
}
if (anyDuplicated(pheno_dt$ID)) stop("Duplicated IDs in pheno_file.")

sample_map <- unique(pheno_dt[, .(ID, FID, AA_only)])

# --------------------------- process chunks ----------------------------------
bundle_regions <- list()
manifest_list <- list()

for (chunk_id in assigned_chunks) {
  cat("\nChunk", chunk_id, "\n")
  parents <- region_dt[data_chunk_id == chunk_id]
  chunk_file <- file.path(chunk_dir, sprintf("chunk_%04d.csv", chunk_id))

  if (!file.exists(chunk_file)) {
    for (i in seq_len(nrow(parents))) {
      manifest_list[[length(manifest_list) + 1L]] <- manifest_row(
        "MISSING_CHUNK_FILE", parents[i], chunk_id, chunk_file
      )
    }
    next
  }

  chunk_dt <- tryCatch(fread(chunk_file), error = function(e) e)
  if (inherits(chunk_dt, "error")) {
    for (i in seq_len(nrow(parents))) {
      manifest_list[[length(manifest_list) + 1L]] <- manifest_row(
        "CHUNK_READ_ERROR", parents[i], chunk_id,
        conditionMessage(chunk_dt)
      )
    }
    next
  }

  if (!all(c("chr", "start") %in% names(chunk_dt))) {
    stop("Chunk ", chunk_id, " lacks chr/start.")
  }
  chunk_dt[, chr := as.integer(chr)]
  chunk_dt[, start := as.integer(start)]

  # Select the union of parent-region CpGs before converting the chunk.
  keep <- rep(FALSE, nrow(chunk_dt))
  for (i in seq_len(nrow(parents))) {
    p <- parents[i]
    keep <- keep | (
      chunk_dt$chr == p$chr &
        chunk_dt$start >= p$region_start &
        chunk_dt$start <= p$region_end
    )
  }
  selected_wide <- chunk_dt[keep]
  rm(chunk_dt, keep)

  if (nrow(selected_wide) == 0L) {
    for (i in seq_len(nrow(parents))) {
      manifest_list[[length(manifest_list) + 1L]] <- manifest_row(
        "NO_MATCHING_CPG_ROWS", parents[i], chunk_id
      )
    }
    next
  }

  chunk_long <- tryCatch(
    convert_chunk_long(selected_wide, pheno_dt),
    error = function(e) e
  )
  rm(selected_wide)

  if (inherits(chunk_long, "error")) {
    for (i in seq_len(nrow(parents))) {
      manifest_list[[length(manifest_list) + 1L]] <- manifest_row(
        "LONG_FORMAT_ERROR", parents[i], chunk_id,
        conditionMessage(chunk_long)
      )
    }
    next
  }

  for (i in seq_len(nrow(parents))) {
    p <- parents[i]
    cat("  parent region", p$region_id, "\n")

    parent_long <- chunk_long[
      chr == p$chr &
        Position >= p$region_start &
        Position <= p$region_end
    ]

    if (nrow(parent_long) == 0L) {
      manifest_list[[length(manifest_list) + 1L]] <- manifest_row(
        "NO_COMPLETE_DATA", p, chunk_id
      )
      next
    }

    parent_n_cpgs_observed <- uniqueN(parent_long$Position)
    parent_n_samples <- uniqueN(parent_long$ID)
    parent_n_fids <- uniqueN(parent_long$FID)

    if (parent_n_cpgs_observed < min_cpgs) {
      manifest_list[[length(manifest_list) + 1L]] <- manifest_row(
        "BELOW_MIN_CPGS", p, chunk_id,
        paste0("Observed CpGs = ", parent_n_cpgs_observed,
               "; min.cpgs = ", min_cpgs),
        n_cpgs_out = parent_n_cpgs_observed,
        n_rows = nrow(parent_long),
        n_samples = parent_n_samples,
        n_fids = parent_n_fids
      )
      next
    }

    # FID is saved separately in sample_map, not passed to SOMNiBUS.
    model_cols <- c(
      "Meth_Counts", "Total_Counts", "Position", "ID",
      selected_covariates
    )
    parent_model <- parent_long[, ..model_cols]
    setorder(parent_model, Position, ID)

    children <- tryCatch(
      splitDataByRegion(
        parent_model,
        gap = gap_bp,
        min.cpgs = min_cpgs,
        max.cpgs = max_cpgs,
        verbose = FALSE
      ),
      error = function(e) e
    )

    if (inherits(children, "error")) {
      manifest_list[[length(manifest_list) + 1L]] <- manifest_row(
        "SPLIT_ERROR", p, chunk_id, conditionMessage(children),
        n_cpgs_out = parent_n_cpgs_observed,
        n_rows = nrow(parent_long),
        n_samples = parent_n_samples,
        n_fids = parent_n_fids
      )
      next
    }

    if (length(children) == 0L) {
      manifest_list[[length(manifest_list) + 1L]] <- manifest_row(
        "NO_SOMNIBUS_REGION", p, chunk_id,
        "splitDataByRegion() returned no region.",
        n_cpgs_out = parent_n_cpgs_observed,
        n_rows = nrow(parent_long),
        n_samples = parent_n_samples,
        n_fids = parent_n_fids
      )
      next
    }

    for (child_index in seq_along(children)) {
      child <- as.data.table(children[[child_index]])
      setorder(child, Position, ID)

      omnibus_region_id <- paste0(
        "P", p$region_id, "_S", sprintf("%03d", child_index)
      )
      bundle_regions[[length(bundle_regions) + 1L]] <- as.data.frame(child)
      names(bundle_regions)[length(bundle_regions)] <- omnibus_region_id
      bundle_index <- length(bundle_regions)

      child_ids <- unique(child$ID)
      child_n_fids <- uniqueN(sample_map[ID %in% child_ids, FID])

      manifest_list[[length(manifest_list) + 1L]] <- manifest_row(
        "READY", p, chunk_id,
        bundle_index = bundle_index,
        omnibus_region_id = omnibus_region_id,
        child_index = child_index,
        region_start = min(child$Position),
        region_end = max(child$Position),
        n_cpgs_out = uniqueN(child$Position),
        n_rows = nrow(child),
        n_samples = uniqueN(child$ID),
        n_fids = child_n_fids
      )
    }
  }

  rm(chunk_long)
  invisible(gc(verbose = FALSE))
}

if (length(manifest_list) > 0L) {
  manifest_dt <- rbindlist(manifest_list, use.names = TRUE, fill = TRUE)
} else {
  manifest_dt <- data.table(
    prep_job_id = integer(), bundle_file = character(),
    bundle_index = integer(), omnibus_region_id = character(),
    data_chunk_id = integer(), parent_region_id = integer(),
    chr = integer(), parent_start = integer(), parent_end = integer(),
    parent_n_cpgs = integer(), child_index = integer(),
    region_start = integer(), region_end = integer(),
    n_cpgs = integer(), n_rows = integer(), n_samples = integer(),
    n_fids = integer(), status = character(), message = character()
  )
}

setorder(
  manifest_dt, data_chunk_id, parent_region_id, child_index,
  na.last = TRUE
)

bundle <- list(
  metadata = list(
    prep_job_id = job_id,
    n_jobs = N_jobs,
    source_region_file = region_file_path,
    selected_covariates = selected_covariates,
    split_settings = list(
      gap = gap_bp, min.cpgs = min_cpgs, max.cpgs = max_cpgs
    ),
    assigned_chunks = assigned_chunks,
    created = format(Sys.time(), tz = "America/Toronto", usetz = TRUE)
  ),
  sample_map = as.data.frame(sample_map),
  manifest_ready = as.data.frame(manifest_dt[status == "READY"]),
  regions = bundle_regions
)

saveRDS(bundle, bundle_file, compress = FALSE)
fwrite(manifest_dt, manifest_file, sep = "\t", quote = FALSE, na = "NA")

cat("\nFinished job", job_id, "\n")
cat("Bundle:", bundle_file, "\n")
cat("Manifest:", manifest_file, "\n")
cat("READY regions:", length(bundle_regions), "\n")
if (nrow(manifest_dt) > 0L) print(manifest_dt[, .N, by = status])
