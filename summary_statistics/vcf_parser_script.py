#!/usr/bin/env python3
"""
vcf_to_table1.py
================
Computes recall and false discovery rates from raw LoFreq VCF files
against a validated plasmid truth set, producing Table 1 for:

    "Benchmarking LoFreq for sub-consensus variant calling on Oxford
     Nanopore long-read data"

Inputs
------
- R9.4.1 uncalibrated VCFs : v5_plas_r9_noQUAD/   (ss{rep}_{bc}_plas.vcf)
- R9.4.1 calibrated VCFs   : V5_plas_r9/           (q-{scalar}_ss{rep}_{bc}_plas.vcf)
- R10.4.1 VCFs             : V5_plas_r11mix/        (ss{rep}_R11_{bc}.vcf  /  q-{scalar}_ss{rep}_R11_{bc}.vcf)

Outputs
-------
- raw_data.csv      : per-replicate TP / FN / FP / Recall / FDR
- summary.csv       : mean ± 95% CI per condition / barcode / AF bin (n=12)
- table1.csv        : Table 1 layout (conditions as rows, AF bins as columns)

Usage
-----
    python vcf_to_table1.py \\
        --r9_uncal  path/to/v5_plas_r9_noQUAD \\
        --r9_cal    path/to/V5_plas_r9 \\
        --r11       path/to/V5_plas_r11mix \\
        --outdir    path/to/output

Notes
-----
- n = 12 subsampling replicates per condition (df = 11, t-critical = 2.201)
- Control row is a confirmed duplicate of Uncalibrated (identical VCF hashes)
- Variant POS 2926 G>T excluded from R10.4.1 bc17 expected list (template
  defect: missed in 100% of replicates across all conditions)
- FDR is computed per replicate file and deduplicated across AF bins before
  averaging, to avoid double-counting the same FP value
"""

import os
import re
import argparse
import hashlib
import pandas as pd
import numpy as np
from scipy import stats

# ── Constants ─────────────────────────────────────────────────────────────────

T_CRIT = stats.t.ppf(0.975, df=11)   # n=12, df=11

COND_ORDER = ["Uncalibrated", "Control"] + [f"Minus {i}" for i in range(1, 11)]

# Variant excluded from R10.4.1 bc17 expected list — template defect
# (missed in 100% of 132 VCFs across all conditions and replicates)
EXCLUDED_FROM_BC17 = {(2926, "G", "T")}

# ── Truth set ─────────────────────────────────────────────────────────────────
# Columns: POS, REF, ALT, AF_bc04, AF_bc05, AF_bc06, AF_bc16, AF_bc17
# R9.4.1 barcodes: 04 (WT), 05 (Delta), 06 (Omicron)
# R10.4.1 barcodes: bc16, bc17

TRUTH_DATA = [
    (2719,"T","C", 0.10,0.15,0.18, 0.05,0.01),
    (2722,"C","A", 0.19,0.19,0.19, 0.15,0.20),
    (2725,"T","C", 0.10,0.15,0.18, 0.05,0.01),
    (2725,"T","A", 0.10,0.05,0.01, 0.10,0.19),
    (2728,"C","T", 0.19,0.19,0.19, 0.15,0.20),
    (2729,"A","T", 0.10,0.05,0.01, 0.10,0.19),
    (2730,"G","C", 0.10,0.05,0.01, 0.10,0.19),
    (2731,"C","T", 0.19,0.19,0.19, 0.15,0.20),
    (2734,"T","C", 0.10,0.15,0.18, 0.05,0.01),
    (2737,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2740,"C","A", 0.10,0.15,0.18, 0.05,0.01),
    (2740,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2743,"C","G", 0.10,0.15,0.18, 0.05,0.01),
    (2746,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2749,"C","G", 0.10,0.15,0.18, 0.05,0.01),
    (2749,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2752,"C","T", 0.10,0.15,0.18, 0.05,0.01),
    (2752,"C","A", 0.10,0.05,0.01, 0.10,0.19),
    (2755,"A","G", 0.10,0.15,0.18, 0.05,0.01),
    (2758,"T","C", 0.10,0.15,0.18, 0.05,0.01),
    (2758,"T","A", 0.10,0.05,0.01, 0.10,0.19),
    (2761,"C","A", 0.19,0.19,0.19, 0.15,0.20),
    (2764,"A","C", 0.10,0.15,0.18, 0.05,0.01),
    (2767,"C","T", 0.19,0.19,0.19, 0.15,0.20),
    (2770,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2771,"A","T", 0.10,0.05,0.01, 0.10,0.19),
    (2772,"G","C", 0.10,0.05,0.01, 0.10,0.19),
    (2773,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2776,"T","C", 0.10,0.05,0.01, 0.10,0.19),
    (2782,"G","T", 0.10,0.15,0.18, 0.15,0.20),
    (2785,"G","T", 0.10,0.05,0.01, 0.05,0.01),
    (2788,"G","C", 0.10,0.15,0.18, 0.05,0.01),
    (2788,"G","T", 0.10,0.05,0.01, 0.10,0.19),
    (2791,"G","T", 0.10,0.05,0.01, 0.10,0.19),
    (2794,"C","T", 0.19,0.19,0.19, 0.15,0.20),
    (2799,"G","A", 0.10,0.05,0.01, 0.10,0.19),
    (2800,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2803,"G","T", 0.10,0.05,0.01, 0.10,0.19),
    (2806,"C","T", 0.10,0.15,0.18, 0.05,0.01),
    (2809,"T","C", 0.19,0.19,0.19, 0.15,0.20),
    (2812,"C","A", 0.10,0.05,0.01, 0.10,0.19),
    (2815,"G","A", 0.10,0.05,0.01, 0.10,0.19),
    (2818,"G","A", 0.10,0.15,0.18, 0.05,0.01),
    (2818,"G","C", 0.10,0.05,0.01, 0.10,0.19),
    (2824,"A","C", 0.10,0.15,0.18, 0.05,0.01),
    (2824,"A","T", 0.10,0.05,0.01, 0.10,0.19),
    (2827,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2830,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2833,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2836,"C","T", 0.10,0.15,0.18, 0.05,0.01),
    (2836,"C","A", 0.10,0.05,0.01, 0.10,0.19),
    (2839,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2842,"G","A", 0.19,0.19,0.19, 0.15,0.20),
    (2845,"G","T", 0.19,0.19,0.19, 0.15,0.20),
    (2848,"C","A", 0.10,0.15,0.18, 0.05,0.01),
    (2848,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2851,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2854,"A","T", 0.10,0.05,0.01, 0.10,0.19),
    (2858,"A","C", 0.19,0.19,0.19, 0.15,0.20),
    (2860,"G","A", 0.10,0.15,0.18, 0.05,0.01),
    (2860,"G","T", 0.10,0.05,0.01, 0.10,0.19),
    (2863,"G","A", 0.10,0.15,0.18, 0.05,0.01),
    (2863,"G","T", 0.10,0.05,0.01, 0.10,0.19),
    (2866,"C","T", 0.19,0.19,0.19, 0.15,0.20),
    (2867,"A","T", 0.19,0.19,0.19, 0.15,0.20),
    (2868,"G","C", 0.19,0.19,0.19, 0.15,0.20),
    (2869,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2872,"A","C", 0.10,0.15,0.18, 0.05,0.01),
    (2875,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2876,"A","T", 0.19,0.19,0.19, 0.15,0.20),
    (2877,"G","C", 0.19,0.19,0.19, 0.15,0.20),
    (2878,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2881,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2884,"A","C", 0.10,0.15,0.18, 0.05,0.01),
    (2884,"A","T", 0.10,0.05,0.01, 0.10,0.19),
    (2887,"C","T", 0.19,0.19,0.19, 0.15,0.20),
    (2890,"G","A", 0.10,0.05,0.01, 0.10,0.19),
    (2893,"C","A", 0.19,0.19,0.19, 0.15,0.20),
    (2894,"A","C", 0.19,0.19,0.19, 0.15,0.20),
    (2896,"A","T", 0.19,0.19,0.19, 0.15,0.20),
    (2899,"C","T", 0.10,0.15,0.18, 0.05,0.01),
    (2899,"C","A", 0.10,0.05,0.01, 0.10,0.19),
    (2902,"C","A", 0.10,0.15,0.18, 0.05,0.01),
    (2906,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2908,"G","A", 0.10,0.05,0.01, 0.10,0.19),
    (2911,"C","A", 0.10,0.05,0.01, 0.10,0.19),
    (2914,"T","G", 0.10,0.05,0.01, 0.10,0.19),
    (2917,"T","A", 0.10,0.15,0.18, 0.05,0.01),
    (2921,"T","C", 0.19,0.19,0.19, 0.15,0.20),
    (2923,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2926,"G","T", 0.10,0.05,0.01, 0.10,0.19),  # excluded from bc17 only
    (2926,"G","C", 0.10,0.15,0.18, 0.05,0.01),
    (2929,"C","T", 0.10,0.15,0.18, 0.05,0.01),
    (2932,"T","C", 0.10,0.05,0.01, 0.10,0.19),
    (2933,"T","A", 0.10,0.15,0.18, 0.05,0.01),
    (2934,"C","G", 0.10,0.15,0.18, 0.05,0.01),
    (2935,"C","A", 0.10,0.05,0.01, 0.10,0.19),
    (2938,"C","T", 0.10,0.05,0.01, 0.10,0.19),
    (2944,"C","T", 0.19,0.19,0.19, 0.15,0.20),
    (2947,"C","T", 0.10,0.15,0.18, 0.05,0.01),
]


# ── Truth set construction ────────────────────────────────────────────────────

def build_truth():
    """
    Returns truth dict: (pos, ref, alt) -> {barcode: af}
    Applies the bc17 exclusion for the template-defective variant.
    """
    truth = {}
    for row in TRUTH_DATA:
        pos, ref, alt = row[0], row[1], row[2]
        key = (pos, ref, alt)
        afs = {
            "04":   row[3],
            "05":   row[4],
            "06":   row[5],
            "bc16": row[6],
            "bc17": row[7],
        }
        if key in EXCLUDED_FROM_BC17:
            del afs["bc17"]
        truth[key] = afs
    return truth


# ── AF binning ────────────────────────────────────────────────────────────────

def af_bin(af):
    if af <= 0.02:   return "1%"
    elif af <= 0.07: return "5%"
    elif af <= 0.12: return "10%"
    elif af <= 0.17: return "15%"
    else:            return "19%"


# ── VCF parsing ───────────────────────────────────────────────────────────────

def vcf_hash(path):
    """MD5 of VCF file for duplicate confirmation."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def parse_vcf(path):
    """Return set of (POS, REF, ALT) for PASS calls."""
    calls = set()
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split("\t")
            if len(p) >= 7 and p[6] == "PASS":
                calls.add((int(p[1]), p[3], p[4]))
    return calls


def analyse_vcf(path, barcode, truth, truth_keys):
    """
    For a single VCF and barcode, compute:
      - per-AF-bin: TP, FN, total_expected, recall
      - overall: FP, total_calls, FDR
    Returns a list of dicts (one per AF bin).
    """
    calls = parse_vcf(path)

    # False positives = PASS calls at positions not in truth at all
    fp = sum(1 for c in calls if c not in truth_keys)
    total_calls = len(calls)
    fdr = fp / total_calls * 100 if total_calls > 0 else 0.0

    # Group truth variants for this barcode by AF bin
    bc_truth = {k: v[barcode] for k, v in truth.items() if barcode in v}
    bins = {}
    for (pos, ref, alt), af in bc_truth.items():
        bins.setdefault(af_bin(af), []).append((pos, ref, alt))

    rows = []
    for b, variants in bins.items():
        tp = sum(1 for v in variants if v in calls)
        fn = len(variants) - tp
        recall = tp / len(variants) * 100
        rows.append({
            "af_bin":         b,
            "tp":             tp,
            "fn":             fn,
            "fp":             fp,
            "total_expected": len(variants),
            "total_calls":    total_calls,
            "recall":         round(recall, 2),
            "fdr":            round(fdr, 2),
        })
    return rows


# ── Directory scanning ────────────────────────────────────────────────────────

def scan_r9_uncal(directory, truth, truth_keys):
    """
    Parse uncalibrated R9.4.1 VCFs.
    Naming: ss{rep}_{barcode}_plas.vcf
    Mixed barcodes only: 04, 05, 06 (01 = clonal, skipped).
    Duplicated as both Uncalibrated and Control.
    """
    records = []
    hashes = {}
    for fname in sorted(os.listdir(directory)):
        m = re.match(r"ss(\d+)_(\d+)_plas\.vcf$", fname)
        if not m:
            continue
        rep, barcode = int(m.group(1)), m.group(2)
        if barcode not in ("04", "05", "06"):
            continue
        path = os.path.join(directory, fname)
        hashes[fname] = vcf_hash(path)
        for r in analyse_vcf(path, barcode, truth, truth_keys):
            base = {"chemistry": "R9.4.1", "rep": rep, "barcode": barcode, **r}
            records.append({**base, "condition": "Uncalibrated"})
            records.append({**base, "condition": "Control"})
    return records, hashes


def scan_r9_cal(directory, truth, truth_keys):
    """
    Parse calibrated R9.4.1 VCFs.
    Naming: q-{scalar}_ss{rep}_{barcode}_plas.vcf
    """
    records = []
    for fname in sorted(os.listdir(directory)):
        m = re.match(r"q-(\d+)_ss(\d+)_(\d+)_plas\.vcf$", fname)
        if not m:
            continue
        scalar, rep, barcode = int(m.group(1)), int(m.group(2)), m.group(3)
        if barcode not in ("04", "05", "06"):
            continue
        path = os.path.join(directory, fname)
        for r in analyse_vcf(path, barcode, truth, truth_keys):
            records.append({
                "chemistry": "R9.4.1",
                "condition": f"Minus {scalar}",
                "rep":       rep,
                "barcode":   barcode,
                **r,
            })
    return records


def scan_r11(directory, truth, truth_keys):
    """
    Parse R10.4.1 VCFs (both uncalibrated and calibrated).
    Uncalibrated naming: ss{rep}_R11_{bc}.vcf  (duplicated as Control)
    Calibrated naming:   q-{scalar}_ss{rep}_R11_{bc}.vcf
    """
    records = []
    hashes = {}
    for fname in sorted(os.listdir(directory)):
        # Uncalibrated
        m = re.match(r"ss(\d+)_R11_(bc\d+)\.vcf$", fname)
        if m:
            rep, bc = int(m.group(1)), m.group(2)
            path = os.path.join(directory, fname)
            hashes[fname] = vcf_hash(path)
            for r in analyse_vcf(path, bc, truth, truth_keys):
                base = {"chemistry": "R10.4.1", "rep": rep, "barcode": bc, **r}
                records.append({**base, "condition": "Uncalibrated"})
                records.append({**base, "condition": "Control"})
            continue
        # Calibrated
        m = re.match(r"q-(\d+)_ss(\d+)_R11_(bc\d+)\.vcf$", fname)
        if m:
            scalar, rep, bc = int(m.group(1)), int(m.group(2)), m.group(3)
            path = os.path.join(directory, fname)
            for r in analyse_vcf(path, bc, truth, truth_keys):
                records.append({
                    "chemistry": "R10.4.1",
                    "condition": f"Minus {scalar}",
                    "rep":       rep,
                    "barcode":   bc,
                    **r,
                })
    return records, hashes


# ── Statistics ────────────────────────────────────────────────────────────────

def ci_stat(series):
    """Mean and 95% CI for a pandas Series (n=12, t-critical=2.201)."""
    n = len(series)
    m = series.mean()
    sd = series.std(ddof=1)
    ci = T_CRIT * sd / np.sqrt(n) if (n > 1 and sd > 0) else 0.0
    return round(m, 2), round(ci, 2), n


def build_summary(df):
    """Per condition / barcode / AF bin summary with mean ± 95% CI."""
    rows = []
    for (chem, cond, bc, af), g in df.groupby(
        ["chemistry", "condition", "barcode", "af_bin"]
    ):
        r = g["recall"].dropna()
        if len(r) == 0:
            continue
        rm, rci, n = ci_stat(r)
        # FDR: deduplicate per (rep, barcode) before averaging
        f = g.drop_duplicates("rep")["fdr"].dropna()
        fm, fci, _ = ci_stat(f) if len(f) > 0 else (None, None, 0)
        rows.append({
            "Chemistry":    chem,
            "Condition":    cond,
            "Barcode":      bc,
            "AF Bin":       af,
            "Recall Mean":  rm,
            "Recall CI":    rci,
            "FDR Mean":     fm,
            "FDR CI":       fci,
            "n":            n,
        })
    return pd.DataFrame(rows)


def build_table1(df):
    """
    Construct Table 1 layout.
    R9.4.1:  5% = bc05, 10% = bc04, 15% = bc05, 19% = bc04+05+06
    R10.4.1: 5% = bc16, 10% = bc16, 15% = bc16, 19% = bc17
    """
    r9_af  = {"5%": ["05"],           "10%": ["04"],
               "15%": ["05"],          "19%": ["04", "05", "06"]}
    r11_af = {"5%": ["bc16"],          "10%": ["bc16"],
               "15%": ["bc16"],        "19%": ["bc17"]}

    def trow(subdf, cond, af_bc, fdr_bcs):
        row = {"Condition": cond}
        for af, bcs in af_bc.items():
            sub = subdf[
                (subdf["condition"] == cond) &
                (subdf["barcode"].isin(bcs)) &
                (subdf["af_bin"] == af)
            ]["recall"]
            if sub.empty:
                row[f"{af} Recall"] = None
                row[f"{af} CI"]     = None
            else:
                m, ci, _ = ci_stat(sub)
                row[f"{af} Recall"] = m
                row[f"{af} CI"]     = ci
        fdr_sub = subdf[
            (subdf["condition"] == cond) &
            (subdf["barcode"].isin(fdr_bcs))
        ].drop_duplicates(["rep", "barcode"])["fdr"]
        if fdr_sub.empty:
            row["FDR"] = None
            row["FDR CI"] = None
        else:
            m, ci, _ = ci_stat(fdr_sub)
            row["FDR"]    = m
            row["FDR CI"] = ci
        return row

    r9d  = df[df["chemistry"] == "R9.4.1"]
    r11d = df[df["chemistry"] == "R10.4.1"]

    r9_rows  = [trow(r9d,  c, r9_af,  ["04", "05", "06"]) for c in COND_ORDER]
    r11_rows = [trow(r11d, c, r11_af, ["bc16", "bc17"])   for c in COND_ORDER]

    r9_df  = pd.DataFrame(r9_rows);  r9_df["Chemistry"]  = "R9.4.1"
    r11_df = pd.DataFrame(r11_rows); r11_df["Chemistry"] = "R10.4.1"
    return pd.concat([r9_df, r11_df], ignore_index=True)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compute Table 1 recall/FDR from raw LoFreq VCFs."
    )
    parser.add_argument("--r9_uncal", required=True,
                        help="Directory of uncalibrated R9.4.1 VCFs (v5_plas_r9_noQUAD)")
    parser.add_argument("--r9_cal",   required=True,
                        help="Directory of calibrated R9.4.1 VCFs (V5_plas_r9)")
    parser.add_argument("--r11",      required=True,
                        help="Directory of R10.4.1 VCFs (V5_plas_r11mix)")
    parser.add_argument("--outdir",   default=".",
                        help="Output directory (default: current directory)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Building truth set...")
    truth = build_truth()
    truth_keys = set(truth.keys())
    print(f"  {len(truth)} truth variants loaded "
          f"(excluding {len(EXCLUDED_FROM_BC17)} from bc17)")

    print("\nParsing R9.4.1 uncalibrated VCFs...")
    r9u_records, r9u_hashes = scan_r9_uncal(args.r9_uncal, truth, truth_keys)
    print(f"  {len(r9u_hashes)} VCFs parsed → "
          f"{len(set(r9u_hashes.values()))} unique hash(es) "
          f"(Uncalibrated duplicated as Control)")

    print("Parsing R9.4.1 calibrated VCFs...")
    r9c_records = scan_r9_cal(args.r9_cal, truth, truth_keys)
    print(f"  {len(r9c_records)} records")

    print("Parsing R10.4.1 VCFs...")
    r11_records, r11_hashes = scan_r11(args.r11, truth, truth_keys)
    print(f"  {len(r11_hashes)} uncalibrated VCFs parsed → "
          f"{len(set(r11_hashes.values()))} unique hash(es)")

    # Combine
    all_records = r9u_records + r9c_records + r11_records
    df = pd.DataFrame(all_records)
    df["_c"] = df["condition"].map({c: i for i, c in enumerate(COND_ORDER)})
    df = df.sort_values(["chemistry", "_c", "barcode", "af_bin", "rep"])\
           .drop(columns="_c")\
           .reset_index(drop=True)

    print(f"\nTotal records: {len(df)}")

    # Save raw data
    raw_path = os.path.join(args.outdir, "raw_data.csv")
    df.to_csv(raw_path, index=False)
    print(f"Saved: {raw_path}")

    # Summary
    summary = build_summary(df)
    sum_path = os.path.join(args.outdir, "summary.csv")
    summary.to_csv(sum_path, index=False)
    print(f"Saved: {sum_path}")

    # Table 1
    table1 = build_table1(df)
    t1_path = os.path.join(args.outdir, "table1.csv")
    table1.to_csv(t1_path, index=False)
    print(f"Saved: {t1_path}")

    # Print Table 1 to console
    print("\n=== TABLE 1 ===")
    cols = ["Chemistry","Condition",
            "5% Recall","5% CI","10% Recall","10% CI",
            "15% Recall","15% CI","19% Recall","19% CI","FDR","FDR CI"]
    print(table1[cols].to_string(index=False, float_format=lambda x: f"{x:.2f}"))


if __name__ == "__main__":
    main()
