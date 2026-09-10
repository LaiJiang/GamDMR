#!/usr/bin/env Rscript

# =============================================================================
# 4_summarize_DMR_statistics.R
#
# Collect all chunk-level partial sufficient statistics and create:
#
# 1) Supplementary_Table_DMR_characteristics.csv
#       one row per DMR across GAM-DMR, SOMNiBUS, BSmooth, DMRcate
#
# 2) DMR_characteristics_summary_by_method.csv
#       method-level mean/median/IQR/range summaries
#
# 3) DMR_characteristics_report.txt
#       ready-to-paste values/sentences for response letter/manuscript
#
# Region methylation definition:
#   First average raw 0-1 CpG methylation proportions within each DMR for each
#   subject, then average those subject-level regional means within AA/control.
#   Therefore each subject contributes equally to a group-level regional mean.
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

BASE <- path.expand("~/scratch/UQAC/meth")

OUT_ROOT <- if (length(args) >= 1L) args[1L] else file.path(
    BASE,
    "results/15_revision/9_dmr_report/DMR_statistics"
)

WORK_DIR <- file.path(OUT_ROOT, "work")
PARTIAL_DIR <- file.path(OUT_ROOT, "partials")

DMR_RDS <- file.path(WORK_DIR, "dmr_table.rds")
MANIFEST_FILE <- file.path(WORK_DIR, "chunk_manifest.csv")
SUBJECT_FILE <- file.path(WORK_DIR, "analysis_subjects.csv")

for (f in c(DMR_RDS, MANIFEST_FILE, SUBJECT_FILE)) {
    if (!file.exists(f)) stop("Required prepared file missing: ", f)
}

dmr <- readRDS(DMR_RDS)
manifest <- fread(MANIFEST_FILE)
subjects <- fread(SUBJECT_FILE)

n_dmr <- nrow(dmr)
n_subjects <- nrow(subjects)

cat("============================================================\n")
cat("Collecting DMR methylation statistics\n")
cat("============================================================\n")
cat("DMRs:          ", n_dmr, "\n", sep = "")
cat("Subjects:      ", n_subjects, "\n", sep = "")
cat("Source chunks: ", nrow(manifest), "\n\n", sep = "")

# -------------------------------------------------------------------------
# Require every current manifest partial
# -------------------------------------------------------------------------

expected_partial <- file.path(
    PARTIAL_DIR,
    sprintf("chunk_%04d.rds", manifest$data_chunk_id)
)

missing <- expected_partial[!file.exists(expected_partial)]

if (length(missing)) {
    stop(
        length(missing),
        " expected partial files are missing. First few:\n",
        paste(head(missing, 20L), collapse = "\n")
    )
}

# -------------------------------------------------------------------------
# Global sufficient statistics
# -------------------------------------------------------------------------

sum_beta <- matrix(
    0,
    nrow = n_dmr,
    ncol = n_subjects
)

n_beta <- matrix(
    0L,
    nrow = n_dmr,
    ncol = n_subjects
)

# Unique CpG coordinates for every DMR.
seen_cpg <- vector("list", n_dmr)

for (ii in seq_len(nrow(manifest))) {

    chunk_id <- manifest$data_chunk_id[ii]
    f <- expected_partial[ii]

    p <- readRDS(f)

    if (!identical(as.character(p$subject_ids), as.character(subjects$ID))) {
        stop("Subject order mismatch in partial: ", f)
    }

    if (!identical(as.integer(p$AA_only), as.integer(subjects$AA_only))) {
        stop("AA status mismatch in partial: ", f)
    }

    if (as.integer(p$data_chunk_id) != as.integer(chunk_id)) {
        stop("Chunk ID mismatch in partial: ", f)
    }

    for (entry in p$entries) {

        rid <- as.integer(entry$DMR_row_id)

        if (is.na(rid) || rid < 1L || rid > n_dmr) {
            stop("Invalid DMR_row_id in partial: ", f)
        }

        # The original split methylation files are expected to partition CpG
        # rows. Detect any unexpected duplicate CpG coordinate across chunks
        # rather than silently double-counting it.
        if (length(entry$cpg_keys)) {

            old <- seen_cpg[[rid]]

            if (length(old)) {
                dup <- intersect(old, entry$cpg_keys)

                if (length(dup)) {
                    stop(
                        "Duplicate CpG coordinate found for DMR_row_id ",
                        rid,
                        " across source chunks. Example: ",
                        dup[1L],
                        ". Stopping to prevent double counting."
                    )
                }
            }

            seen_cpg[[rid]] <- c(old, entry$cpg_keys)
        }

        sum_beta[rid, ] <- sum_beta[rid, ] + as.numeric(entry$sum_beta)
        n_beta[rid, ] <- n_beta[rid, ] + as.integer(entry$n_beta)
    }

    if (ii == 1L || ii %% 50L == 0L || ii == nrow(manifest)) {
        cat(
            sprintf(
                "Collected %d/%d chunks; %s\n",
                ii,
                nrow(manifest),
                format(Sys.time())
            )
        )
    }

    rm(p)
}

# -------------------------------------------------------------------------
# Subject-level regional methylation
# -------------------------------------------------------------------------

subject_region_mean <- sum_beta / n_beta
subject_region_mean[n_beta == 0L] <- NA_real_

case_idx <- which(subjects$AA_only == 1L)
ctrl_idx <- which(subjects$AA_only == 0L)

safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    mean(x)
}

safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    median(x)
}

safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(NA_real_)
    sd(x)
}

safe_q <- function(x, p) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    as.numeric(quantile(x, p, type = 7, names = FALSE))
}

n_cpg <- lengths(seen_cpg)

mean_all <- median_all <- rep(NA_real_, n_dmr)
mean_AA <- median_AA <- sd_AA <- rep(NA_real_, n_dmr)
mean_ctrl <- median_ctrl <- sd_ctrl <- rep(NA_real_, n_dmr)
n_all_obs <- n_AA_obs <- n_ctrl_obs <- integer(n_dmr)

for (rid in seq_len(n_dmr)) {

    z <- subject_region_mean[rid, ]
    z_AA <- z[case_idx]
    z_ctrl <- z[ctrl_idx]

    mean_all[rid] <- safe_mean(z)
    median_all[rid] <- safe_median(z)

    mean_AA[rid] <- safe_mean(z_AA)
    median_AA[rid] <- safe_median(z_AA)
    sd_AA[rid] <- safe_sd(z_AA)

    mean_ctrl[rid] <- safe_mean(z_ctrl)
    median_ctrl[rid] <- safe_median(z_ctrl)
    sd_ctrl[rid] <- safe_sd(z_ctrl)

    n_all_obs[rid] <- sum(is.finite(z))
    n_AA_obs[rid] <- sum(is.finite(z_AA))
    n_ctrl_obs[rid] <- sum(is.finite(z_ctrl))
}

# -------------------------------------------------------------------------
# Per-DMR supplementary table
# -------------------------------------------------------------------------

supp <- copy(dmr)

supp[, `:=`(
    n_CpGs = n_cpg,

    mean_methylation_all = mean_all,
    median_methylation_all = median_all,

    mean_methylation_AA = mean_AA,
    median_methylation_AA = median_AA,
    SD_methylation_AA = sd_AA,

    mean_methylation_control = mean_ctrl,
    median_methylation_control = median_ctrl,
    SD_methylation_control = sd_ctrl,

    methylation_difference_AA_minus_control = mean_AA - mean_ctrl,
    absolute_methylation_difference = abs(mean_AA - mean_ctrl),

    n_subjects_with_methylation = n_all_obs,
    n_AA_with_methylation = n_AA_obs,
    n_control_with_methylation = n_ctrl_obs
)]

setcolorder(
    supp,
    c(
        "DMR_row_id",
        "DMR_method_id",
        "method",
        "DMR_chr",
        "DMR_region_start",
        "DMR_region_end",
        "DMR_width_bp",
        "n_CpGs",
        "mean_methylation_all",
        "mean_methylation_AA",
        "mean_methylation_control",
        "methylation_difference_AA_minus_control",
        "absolute_methylation_difference",
        "median_methylation_all",
        "median_methylation_AA",
        "median_methylation_control",
        "SD_methylation_AA",
        "SD_methylation_control",
        "n_subjects_with_methylation",
        "n_AA_with_methylation",
        "n_control_with_methylation",
        "n_source_chunks"
    )
)

if ("DMR_method_index" %in% names(supp)) {
    supp[, DMR_method_index := NULL]
}

supp_file <- file.path(
    OUT_ROOT,
    "Supplementary_Table_DMR_characteristics.csv"
)

fwrite(supp, supp_file)

# QC file for any DMR with zero CpGs.
zero_file <- file.path(OUT_ROOT, "DMRs_with_zero_extracted_CpGs.csv")
fwrite(supp[n_CpGs == 0L], zero_file)

# -------------------------------------------------------------------------
# Method-level summary
# -------------------------------------------------------------------------

method_order <- c("GAM-DMR", "SOMNiBUS", "BSmooth", "DMRcate")

method_summary <- supp[, {

    cpg <- as.numeric(n_CpGs)
    width <- as.numeric(DMR_width_bp)
    aa <- mean_methylation_AA
    ctrl <- mean_methylation_control
    delta <- methylation_difference_AA_minus_control
    abs_delta <- absolute_methylation_difference

    list(
        N_DMRs = .N,
        N_DMRs_with_CpGs = sum(cpg > 0, na.rm = TRUE),
        N_DMRs_with_zero_CpGs = sum(cpg == 0 | is.na(cpg)),

        CpGs_mean = safe_mean(cpg),
        CpGs_SD = safe_sd(cpg),
        CpGs_median = safe_median(cpg),
        CpGs_Q1 = safe_q(cpg, 0.25),
        CpGs_Q3 = safe_q(cpg, 0.75),
        CpGs_IQR = safe_q(cpg, 0.75) - safe_q(cpg, 0.25),
        CpGs_min = if (all(is.na(cpg))) NA_real_ else min(cpg, na.rm = TRUE),
        CpGs_max = if (all(is.na(cpg))) NA_real_ else max(cpg, na.rm = TRUE),

        Width_bp_mean = safe_mean(width),
        Width_bp_SD = safe_sd(width),
        Width_bp_median = safe_median(width),
        Width_bp_Q1 = safe_q(width, 0.25),
        Width_bp_Q3 = safe_q(width, 0.75),
        Width_bp_IQR = safe_q(width, 0.75) - safe_q(width, 0.25),
        Width_bp_min = if (all(is.na(width))) NA_real_ else min(width, na.rm = TRUE),
        Width_bp_max = if (all(is.na(width))) NA_real_ else max(width, na.rm = TRUE),

        Mean_regional_methylation_AA = safe_mean(aa),
        Median_regional_methylation_AA = safe_median(aa),

        Mean_regional_methylation_control = safe_mean(ctrl),
        Median_regional_methylation_control = safe_median(ctrl),

        Mean_AA_minus_control_difference = safe_mean(delta),
        Median_AA_minus_control_difference = safe_median(delta),

        Mean_absolute_methylation_difference = safe_mean(abs_delta),
        Median_absolute_methylation_difference = safe_median(abs_delta)
    )

}, by = method]

method_summary[, order_tmp := match(method, method_order)]
setorder(method_summary, order_tmp)
method_summary[, order_tmp := NULL]

summary_file <- file.path(
    OUT_ROOT,
    "DMR_characteristics_summary_by_method.csv"
)

fwrite(method_summary, summary_file)

# -------------------------------------------------------------------------
# Ready-to-paste text report
# -------------------------------------------------------------------------

fmt <- function(x, digits = 2L) {
    ifelse(
        is.na(x),
        "NA",
        formatC(x, format = "f", digits = digits, big.mark = ",")
    )
}

fmt0 <- function(x) fmt(x, 0L)

report_file <- file.path(
    OUT_ROOT,
    "DMR_characteristics_report.txt"
)

con <- file(report_file, open = "wt")

writeLines(
    c(
        "DMR CHARACTERISTICS REPORT",
        "==========================",
        "",
        paste0(
            "Analytic subjects: ",
            n_subjects,
            " (AA=",
            length(case_idx),
            ", controls=",
            length(ctrl_idx),
            ")"
        ),
        "",
        "Definition:",
        paste0(
            "For each DMR and subject, raw methylation proportions were averaged ",
            "over observed CpGs in that interval. AA/control regional methylation ",
            "was then calculated as the mean of subject-level regional means."
        ),
        "",
        "METHOD-LEVEL VALUES",
        "-------------------"
    ),
    con
)

for (m in method_order) {

    z <- method_summary[method == m]
    if (!nrow(z)) next

    writeLines(
        c(
            "",
            paste0(m, ":"),
            paste0("  N DMRs = ", fmt0(z$N_DMRs)),
            paste0(
                "  CpGs/region: mean ",
                fmt(z$CpGs_mean, 2),
                "; median ",
                fmt(z$CpGs_median, 1),
                " (IQR ",
                fmt(z$CpGs_Q1, 1),
                "--",
                fmt(z$CpGs_Q3, 1),
                "; range ",
                fmt0(z$CpGs_min),
                "--",
                fmt0(z$CpGs_max),
                ")"
            ),
            paste0(
                "  Width (bp): mean ",
                fmt(z$Width_bp_mean, 1),
                "; median ",
                fmt(z$Width_bp_median, 1),
                " (IQR ",
                fmt(z$Width_bp_Q1, 1),
                "--",
                fmt(z$Width_bp_Q3, 1),
                "; range ",
                fmt0(z$Width_bp_min),
                "--",
                fmt0(z$Width_bp_max),
                ")"
            ),
            paste0(
                "  Mean regional methylation: AA=",
                fmt(z$Mean_regional_methylation_AA, 4),
                "; controls=",
                fmt(z$Mean_regional_methylation_control, 4)
            ),
            paste0(
                "  Median absolute AA-control methylation difference = ",
                fmt(z$Median_absolute_methylation_difference, 4)
            ),
            paste0(
                "  DMRs with zero extracted CpGs = ",
                fmt0(z$N_DMRs_with_zero_CpGs)
            )
        ),
        con
    )
}

gam <- method_summary[method == "GAM-DMR"]

if (nrow(gam) == 1L) {

    sentence <- paste0(
        "The ",
        fmt0(gam$N_DMRs),
        " GAM-DMRs contained a mean of ",
        fmt(gam$CpGs_mean, 2),
        " CpGs and a median of ",
        fmt(gam$CpGs_median, 1),
        " CpGs (IQR ",
        fmt(gam$CpGs_Q1, 1),
        "--",
        fmt(gam$CpGs_Q3, 1),
        "; range ",
        fmt0(gam$CpGs_min),
        "--",
        fmt0(gam$CpGs_max),
        ") per region and had a mean genomic width of ",
        fmt(gam$Width_bp_mean, 1),
        " bp and a median width of ",
        fmt(gam$Width_bp_median, 1),
        " bp (IQR ",
        fmt(gam$Width_bp_Q1, 1),
        "--",
        fmt(gam$Width_bp_Q3, 1),
        "; range ",
        fmt0(gam$Width_bp_min),
        "--",
        fmt0(gam$Width_bp_max),
        "). Mean regional methylation on the original proportion scale was ",
        fmt(gam$Mean_regional_methylation_AA, 4),
        " in AA cases and ",
        fmt(gam$Mean_regional_methylation_control, 4),
        " in controls, with a median absolute between-group methylation difference of ",
        fmt(gam$Median_absolute_methylation_difference, 4),
        "."
    )

    writeLines(
        c(
            "",
            "",
            "READY-TO-PASTE GAM-DMR RESULTS SENTENCE",
            "---------------------------------------",
            sentence
        ),
        con
    )
}

writeLines(
    c(
        "",
        "",
        "OUTPUT FILES",
        "------------",
        paste0("Supplementary per-DMR table: ", supp_file),
        paste0("Method summary CSV: ", summary_file),
        paste0("Ready-to-paste report: ", report_file),
        paste0("Zero-CpG QC table: ", zero_file),
        "",
        "Scale note:",
        paste0(
            "Descriptive methylation values are raw 0-1 methylation proportions. ",
            "GAM-DMR inference remains based on arcsine-square-root transformed ",
            "methylation."
        )
    ),
    con
)

close(con)

# Cache the completed result for future table formatting without rescanning.
saveRDS(
    list(
        per_DMR = supp,
        by_method = method_summary,
        subjects = subjects
    ),
    file.path(OUT_ROOT, "DMR_characteristics_results.rds"),
    compress = "xz"
)

cat("\n============================================================\n")
cat("Completed successfully\n")
cat("============================================================\n")
cat("Supplementary table:\n  ", supp_file, "\n", sep = "")
cat("Method summary:\n  ", summary_file, "\n", sep = "")
cat("Report:\n  ", report_file, "\n", sep = "")
cat("Zero-CpG QC:\n  ", zero_file, "\n\n", sep = "")
print(method_summary)
