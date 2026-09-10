#!/usr/bin/env python3

"""
4_summarize_epigen.py

Combine the outputs of 2_epigen_dmr_overlap.py across four DMR methods and
generate supplementary tables.

Expected folder layout under --root:

  GAM_DMR/
  SOMNiBUS/
  BSmooth/
  DMRcate/

Each folder should contain:
  EPIGEN_DMR_meQTL_overlap_summary.csv
  EPIGEN_DMR_meQTL_associations.tsv.gz

Outputs:
  Supplementary_Table_EPIGEN_meQTL_overlap_by_method.csv
  Supplementary_Table_EPIGEN_meQTL_overlap_by_DMR.csv
  EPIGEN_meQTL_overlap_summary.txt
"""

import argparse
import csv
import gzip
import math
import os
from statistics import median


METHODS = [
    ("GAM-DMR",  "GAM_DMR"),
    ("SOMNiBUS", "SOMNiBUS"),
    ("BSmooth",  "BSmooth"),
    ("DMRcate",  "DMRcate"),
]


def yes(x):
    return str(x).strip().lower() == "yes"


def to_int(x, default=0):
    try:
        return int(float(x))
    except Exception:
        return default


def pct(n, d):
    return 100.0 * n / d if d else math.nan


def fmt_pct(x):
    return "" if not math.isfinite(x) else f"{x:.2f}"


def read_summary(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def scan_associations(path):
    """
    Count globally unique CpGs/SNPs/LD clumps for a METHOD.

    This avoids summing per-DMR unique counts, which could double-count the
    same CpG or SNP if two DMR intervals overlap.
    """
    all_cpg = set()
    cis_cpg = set()
    trans_cpg = set()

    all_snp = set()
    cis_snp = set()
    trans_snp = set()

    ld_clumps = set()
    n_rows = 0

    if not os.path.exists(path):
        raise FileNotFoundError(path)

    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")

        for row in r:
            n_rows += 1

            cpg = row.get("CpG", "").strip()
            snp = row.get("SNP", "").strip()
            typ = row.get("type", "").strip().lower()
            ld = row.get("LD_clump", "").strip()

            if cpg:
                all_cpg.add(cpg)
                if typ == "cis":
                    cis_cpg.add(cpg)
                elif typ == "trans":
                    trans_cpg.add(cpg)

            if snp:
                all_snp.add(snp)
                if typ == "cis":
                    cis_snp.add(snp)
                elif typ == "trans":
                    trans_snp.add(snp)

            if ld:
                ld_clumps.add(ld)

    return {
        "N_association_DMR_overlap_rows": n_rows,
        "N_unique_meQTL_target_CpGs": len(all_cpg),
        "N_unique_cis_meQTL_target_CpGs": len(cis_cpg),
        "N_unique_trans_meQTL_target_CpGs": len(trans_cpg),
        "N_unique_meQTL_SNPs": len(all_snp),
        "N_unique_cis_meQTL_SNPs": len(cis_snp),
        "N_unique_trans_meQTL_SNPs": len(trans_snp),
        "N_unique_LD_clumps": len(ld_clumps),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root",
        required=True,
        help="Root EPIGEN output directory containing the four method folders",
    )
    args = ap.parse_args()

    root = os.path.abspath(args.root)

    method_rows = []
    combined_dmr_rows = []

    for method, slug in METHODS:

        method_dir = os.path.join(root, slug)

        summary_file = os.path.join(
            method_dir,
            "EPIGEN_DMR_meQTL_overlap_summary.csv",
        )

        assoc_file = os.path.join(
            method_dir,
            "EPIGEN_DMR_meQTL_associations.tsv.gz",
        )

        if not os.path.exists(summary_file):
            raise FileNotFoundError(
                f"Missing summary file for {method}: {summary_file}"
            )

        rows = read_summary(summary_file)

        n_dmr = len(rows)
        n_any = sum(yes(r.get("EPIGEN_meQTL_overlap")) for r in rows)
        n_cis = sum(yes(r.get("EPIGEN_cis_meQTL_overlap")) for r in rows)
        n_trans = sum(yes(r.get("EPIGEN_trans_meQTL_overlap")) for r in rows)

        positive_cpg_counts = [
            to_int(r.get("n_unique_meQTL_CpGs", 0))
            for r in rows
            if yes(r.get("EPIGEN_meQTL_overlap"))
        ]

        positive_cis_counts = [
            to_int(r.get("n_cis_meQTL_CpGs", 0))
            for r in rows
            if yes(r.get("EPIGEN_cis_meQTL_overlap"))
        ]

        assoc_stats = scan_associations(assoc_file)

        method_row = {
            "Method": method,
            "N_DMRs": n_dmr,

            "N_DMRs_with_any_EPIGEN_meQTL": n_any,
            "Percent_DMRs_with_any_EPIGEN_meQTL": fmt_pct(pct(n_any, n_dmr)),

            "N_DMRs_with_cis_meQTL": n_cis,
            "Percent_DMRs_with_cis_meQTL": fmt_pct(pct(n_cis, n_dmr)),

            "N_DMRs_with_trans_meQTL": n_trans,
            "Percent_DMRs_with_trans_meQTL": fmt_pct(pct(n_trans, n_dmr)),

            "Median_unique_meQTL_CpGs_per_positive_DMR":
                "" if not positive_cpg_counts else f"{median(positive_cpg_counts):.1f}",

            "Median_cis_meQTL_CpGs_per_cis_positive_DMR":
                "" if not positive_cis_counts else f"{median(positive_cis_counts):.1f}",
        }

        method_row.update(assoc_stats)
        method_rows.append(method_row)

        # Combined DMR-level table.
        for r in rows:
            rr = {"Method": method}
            rr.update(r)
            combined_dmr_rows.append(rr)

    # ------------------------------------------------------------------
    # Method-level supplementary comparison table
    # ------------------------------------------------------------------

    method_out = os.path.join(
        root,
        "Supplementary_Table_EPIGEN_meQTL_overlap_by_method.csv",
    )

    with open(method_out, "w", newline="", encoding="utf-8") as f:
        fields = list(method_rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(method_rows)

    # ------------------------------------------------------------------
    # Combined DMR-level supplementary table
    # ------------------------------------------------------------------

    dmr_out = os.path.join(
        root,
        "Supplementary_Table_EPIGEN_meQTL_overlap_by_DMR.csv",
    )

    # Preserve the union of all fields in stable order.
    fields = ["Method"]
    seen = {"Method"}

    for r in combined_dmr_rows:
        for k in r.keys():
            if k not in seen:
                fields.append(k)
                seen.add(k)

    with open(dmr_out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(combined_dmr_rows)

    # ------------------------------------------------------------------
    # Human-readable summary
    # ------------------------------------------------------------------

    txt_out = os.path.join(root, "EPIGEN_meQTL_overlap_summary.txt")

    with open(txt_out, "w", encoding="utf-8") as f:
        f.write("EPIGEN significant meQTL target-CpG overlap with DMRs\n")
        f.write("=====================================================\n\n")

        for r in method_rows:
            f.write(f"{r['Method']}\n")
            f.write(f"  DMRs: {r['N_DMRs']}\n")
            f.write(
                "  Any meQTL overlap: "
                f"{r['N_DMRs_with_any_EPIGEN_meQTL']} "
                f"({r['Percent_DMRs_with_any_EPIGEN_meQTL']}%)\n"
            )
            f.write(
                "  Cis-meQTL overlap: "
                f"{r['N_DMRs_with_cis_meQTL']} "
                f"({r['Percent_DMRs_with_cis_meQTL']}%)\n"
            )
            f.write(
                "  Trans-meQTL overlap: "
                f"{r['N_DMRs_with_trans_meQTL']} "
                f"({r['Percent_DMRs_with_trans_meQTL']}%)\n"
            )
            f.write(
                "  Unique EPIGEN meQTL target CpGs: "
                f"{r['N_unique_meQTL_target_CpGs']}\n"
            )
            f.write("\n")

        f.write(
            "Interpretation note: meQTL_full.txt.gz contains significant EPIGEN "
            "meQTL associations. A DMR classified as 'No' therefore means that "
            "no significant EPIGEN meQTL target CpG in this file fell within "
            "that DMR; it does not establish absence of genetic regulation.\n"
        )

    # ------------------------------------------------------------------
    # Console output
    # ------------------------------------------------------------------

    print("\nEPIGEN meQTL overlap summary")
    print("=" * 78)

    for r in method_rows:
        print(
            f"{r['Method']:10s} | "
            f"N={r['N_DMRs']:5d} | "
            f"any={r['N_DMRs_with_any_EPIGEN_meQTL']:5d} "
            f"({r['Percent_DMRs_with_any_EPIGEN_meQTL']:>6s}%) | "
            f"cis={r['N_DMRs_with_cis_meQTL']:5d} "
            f"({r['Percent_DMRs_with_cis_meQTL']:>6s}%) | "
            f"trans={r['N_DMRs_with_trans_meQTL']:5d} "
            f"({r['Percent_DMRs_with_trans_meQTL']:>6s}%)"
        )

    print("\nCreated:")
    print(" ", method_out)
    print(" ", dmr_out)
    print(" ", txt_out)


if __name__ == "__main__":
    main()
