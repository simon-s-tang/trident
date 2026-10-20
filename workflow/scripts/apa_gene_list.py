"""
Write the reliably expressed gene table for apatools from gene counts
(featureCounts matrix or per-sample htseq-count files).

Genes sharing a symbol are summed; genes with a mean count per sample above
--min-mean are kept. apatools uses only the gene names and sample columns.
"""
import argparse
import os
import re
import sys

import pandas as pd


# Gene symbols from the GTF
def gene_names(gtf):
    """gene_id -> gene_name, from the GTF's gene lines."""
    gene_id = re.compile(r'gene_id "([^"]+)"')
    gene_name = re.compile(r'gene_name "([^"]+)"')
    names = {}
    with open(gtf) as handle:
        for line in handle:
            fields = line.split("\t", 8)
            if len(fields) < 9 or fields[2] != "gene":
                continue
            gid = gene_id.search(fields[8]).group(1)
            name = gene_name.search(fields[8])
            names[gid] = name.group(1) if name else gid
    return pd.Series(names)


# Count readers
def read_featurecounts(path, samples):
    """featureCounts matrix -> gene_id x sample. Its columns are BAM paths named
    <sample>_Aligned..., so the sample is the file name before '_Aligned'."""
    counts = pd.read_csv(path, sep="\t", comment="#", index_col=0).iloc[:, 5:]
    counts.columns = [os.path.basename(c).split("_Aligned")[0] for c in counts.columns]
    missing = set(samples) - set(counts.columns)
    if missing:
        sys.exit(f"Samples not in {path}: {sorted(missing)}")
    return counts[samples]


def read_htseq(paths, samples):
    """htseq-count files -> gene_id x sample. One file per sample, in the same
    order as --samples; the special __no_feature etc. rows are dropped."""
    if len(paths) != len(samples):
        sys.exit(f"{len(paths)} htseq files for {len(samples)} samples")
    columns = {}
    for path, sample in zip(paths, samples):
        table = pd.read_csv(path, sep="\t", header=None, index_col=0)
        table = table[~table.index.str.startswith("__")]
        columns[sample] = table.iloc[:, -1]
    return pd.DataFrame(columns)


# Command line
def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--format", choices=["featurecounts", "htseq"], required=True)
    parser.add_argument("--counts", nargs="+", required=True,
                        help="the featureCounts matrix, or one htseq file per sample")
    parser.add_argument("--gtf", required=True, help="GTF the counts were made with")
    parser.add_argument("--samples", nargs="+", required=True,
                        help="sample names, in the order of the htseq files")
    parser.add_argument("--min-mean", type=float, required=True,
                        help="keep genes whose mean count per sample is above this")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if args.format == "featurecounts":
        counts = read_featurecounts(args.counts[0], args.samples)
    else:
        counts = read_htseq(args.counts, args.samples)

    names = gene_names(args.gtf)
    symbols = names.reindex(counts.index).fillna(pd.Series(counts.index, index=counts.index))
    by_symbol = counts.groupby(symbols.values).sum()

    kept = by_symbol[by_symbol.mean(axis=1) > args.min_mean]
    kept.to_csv(args.out)
    print(
        f"{len(counts)} genes ({len(by_symbol)} symbols); {len(kept)} with a mean "
        f"count per sample above {args.min_mean} written to {args.out}"
    )


if __name__ == "__main__":
    main()
