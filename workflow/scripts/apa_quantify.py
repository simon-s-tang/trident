"""
3' UTR isoform quantification with apatools, from per-sample kallisto output
to the table of proximal/distal isoform pairs. No statistical test is run.

Steps follow apatools' bulkapa.py up to create_table_of_pairs: merge kallisto
counts, correct short proximal isoforms, drop length outliers, keep reliably
expressed multi-isoform genes, bulk correction, per-sample thresholds with a
group-wise expression filter, optional quantile normalisation, pairing.
"""
import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
from bulkanalysis.filtering import filter_samples_based_on_thresholds, gmm_threshold
from bulkanalysis.normalization import quantile_normalize
from bulkanalysis.plot import plot_gene_expression_distribution
from gtfparse import read_gtf as gtfparse_read_gtf
from kallistomanager.loader import merge_kallisto_counts

import apatools.filtering
from apatools.apanalysis import create_table_of_pairs
from apatools.filtering import (
    add_length_info,
    bulk_correction,
    check_proximal_abundances_in_raw_isoforms_counts,
    correct_counts_for_small_proximals,
    remove_length_outliers,
    stats_raw,
)


# Run-time fixes to apatools (github.com/idiap/apatools, commit 359bda1):
# 1. gtfparse >= 2 returns polars tables; apatools.filtering expects pandas.
def read_gtf(path):
    gtf = gtfparse_read_gtf(path)
    return gtf if isinstance(gtf, pd.DataFrame) else gtf.to_pandas()


apatools.filtering.read_gtf = read_gtf


# 2. apatools takes an isoform's gene as id[:-2]; for IDs like
#    GENE.TRANSCRIPT.SOURCE the gene is the part before the first dot.
def isoform_gene(isoform_id):
    return isoform_id.split(".")[0]


def get_multi_utrs_isoforms_from_reliably_expressed_genes(
    df_counts: pd.DataFrame, df_counts_genes: pd.DataFrame
):
    reliably_expressed_genes = list(df_counts_genes.index)

    df_counts["gene"] = [isoform_gene(x) for x in df_counts.index]
    df_counts_iso = df_counts[df_counts["gene"].isin(reliably_expressed_genes)]

    multi_utr_genes = (
        df_counts_iso.groupby("gene")
        .count()[df_counts_iso.groupby("gene").count()[df_counts_iso.columns[0]] > 1]
        .index
    )
    df_counts_iso = df_counts_iso[df_counts_iso["gene"].isin(multi_utr_genes)]

    print(
        f"After filtering isoforms of genes that are not reliably expressed, \
          we have {len(df_counts_iso)} reliably expressed isoforms."
    )

    return df_counts_iso


# Command line
def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--kallisto", nargs="+", required=True,
                        help="abundance.tsv per sample, each in a folder named after the sample")
    parser.add_argument("--genes", required=True,
                        help="reliably expressed genes (scripts/apa_gene_list.py)")
    parser.add_argument("--samples", required=True, help="samples.tsv")
    parser.add_argument("--group-column", required=True,
                        help="column of samples.tsv that defines the groups")
    parser.add_argument("--gtf-focus", required=True,
                        help="GTF of the quantified 3' UTR windows (the kallisto index)")
    parser.add_argument("--gtf-no-focus", required=True,
                        help="GTF of the full-length 3' UTR isoforms")
    parser.add_argument("--focus-nt", type=int, required=True)
    parser.add_argument("--minimal-required-length", type=int, required=True)
    parser.add_argument("--quantile-threshold", type=float, required=True)
    parser.add_argument("--background-quantile", choices=["true", "false"], required=True)
    parser.add_argument("--pct-in-treatment", type=float, required=True)
    parser.add_argument("--quantile-normalization", choices=["true", "false"], required=True)
    parser.add_argument("--figures-extension", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    figures = os.path.join(args.outdir, "Figures")
    os.makedirs(figures, exist_ok=True)

    # apatools reports through print(); send it to message.log
    log = open(os.path.join(args.outdir, "message.log"), "w")
    sys.stdout = log

    # Groups from the sample sheet
    sheet = pd.read_csv(args.samples, sep="\t", dtype=str)
    if args.group_column not in sheet.columns:
        sys.exit(f"{args.samples} has no column {args.group_column}")
    treatments = {
        group: list(rows["sample_name"])
        for group, rows in sheet.groupby(args.group_column, sort=False)
    }
    all_samples = [s for samples in treatments.values() for s in samples]
    print(f"Groups: {treatments}")

    # 1. Load kallisto counts and the gene list
    merge_kallisto_counts(save_dir=args.outdir, data_list=args.kallisto)
    df_counts = pd.read_csv(os.path.join(args.outdir, "counts_matrix.csv"), index_col=0)
    df_filtered_genes = pd.read_csv(args.genes, index_col=0)

    df_filtered_genes = df_filtered_genes[all_samples]
    df_counts = df_counts[all_samples]

    print(f"Before filtering, we have {len(df_counts)} isoforms.")
    reps = df_counts.columns

    # 2. Correct isoforms annotated shorter than the quantified window
    df_counts = correct_counts_for_small_proximals(
        df_counts=df_counts,
        reps=reps,
        gtf_file_quantification=args.gtf_focus,
        focus_nt=args.focus_nt,
        minimal_required_length=args.minimal_required_length,
    )

    # 3. Drop genes with outlier lengths
    df_counts = remove_length_outliers(df_counts[reps])

    stats_raw(
        df_counts,
        save_file=os.path.join(figures, f"reads_per_sample.{args.figures_extension}"),
    )

    # 4. Keep reliably expressed genes with more than one isoform
    df_counts_iso = get_multi_utrs_isoforms_from_reliably_expressed_genes(
        df_counts, df_filtered_genes
    )
    print(
        f"Before filtering, we have {len(df_counts_iso)} isoforms and "
        f"{len(df_counts_iso.groupby('gene').sum())} genes."
    )

    # 5. Isoform lengths and proximal abundance check
    df_counts_iso = add_length_info(df_counts_iso, path_gtf_no_length_focus=args.gtf_no_focus)
    _, df_counts_iso = check_proximal_abundances_in_raw_isoforms_counts(
        df_counts_iso=df_counts_iso, path_to_results=args.outdir, reps=reps
    )

    # 6. Bulk correction
    df_counts_iso_corr = bulk_correction(
        df_counts_iso,
        columns_to_correct=list(reps),
        replace_col=True,
        gene_id="gene",
        length_id="length",
        clip=True,
    )

    # 7. Per-sample thresholds, then the group-wise expression filter
    df_counts_iso_log = np.log2(
        df_counts_iso_corr.drop(columns=["gene", "length", "strand"], axis=1).set_index(
            "transcript_id"
        )
        + 1
    )
    thresholds = {}
    for col in df_counts_iso_log.columns:
        thresholds[col] = gmm_threshold(
            data=np.array(df_counts_iso_log[df_counts_iso_log[col] != 0][col]),
            quantile_threshold=args.quantile_threshold,
            background_quantile=args.background_quantile == "true",
            density=True,
            show_plot=False,
            data_name=col,
            saving_folder=figures,
            extension=args.figures_extension,
            return_int=False,
        )
    print(thresholds)

    df_counts_iso_corr_filtered = filter_samples_based_on_thresholds(
        df_genes=df_counts_iso_corr,
        treatments=treatments,
        thresholds=thresholds,
        pct_in_treatment=args.pct_in_treatment,
        save_file=os.path.join(
            figures, f"isoform_exp_after_filtering.{args.figures_extension}"
        ),
    )
    print(f"{df_counts_iso_corr_filtered.shape[0]} isoforms left after filtering.")
    df_counts_iso_corr_filtered.set_index("transcript_id", inplace=True)

    # 8. Optional quantile normalisation
    if args.quantile_normalization == "true":
        df_counts_iso_corr_normalized = quantile_normalize(df_counts_iso_corr_filtered[reps])
        plot_gene_expression_distribution(
            df_counts_iso_corr_normalized,
            title="Isoform expression after filtering and normalization",
            save_file=os.path.join(
                figures,
                f"isoform_expression_after_filtering_quantile.{args.figures_extension}",
            ),
        )
        df_counts_iso_corr_normalized["gene"] = [
            isoform_gene(x) for x in df_counts_iso_corr_normalized.index
        ]
        df_counts_iso_corr_normalized = add_length_info(
            df_counts_iso_corr_normalized, path_gtf_no_length_focus=args.gtf_no_focus
        )
    else:
        df_counts_iso_corr_normalized = df_counts_iso_corr_filtered.copy()

    # 9. Proximal/distal pairs (create_table_of_pairs writes table_of_pairs.csv)
    table_of_pairs = create_table_of_pairs(
        df_counts_iso_corr_normalized.reset_index(),
        treatments=treatments,
        replicates=False,
        path_to_results=args.outdir,
    )
    print(f"We analyze {len(table_of_pairs)} pairs of isoforms.")

    df_counts_iso_corr_normalized.to_csv(
        os.path.join(args.outdir, "df_counts_filtered_corrected_normalized.csv")
    )

    sys.stdout = sys.__stdout__
    log.close()


if __name__ == "__main__":
    main()
