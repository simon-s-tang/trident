"""
Build the small test dataset in .test/data from real reference files.

The reference is real: a slice of the GRCh38 genome, the GENCODE annotation of
the genes that lie fully inside it, and the matching part of the 3' UTR
reference (isoforms and windows; the pipeline builds their kallisto index).
The reads are
simulated from those genes, so the data can be shared and the expected
signal is known:

- 4 samples in 2 conditions (control, treated); sampleA1 has two lanes;
- reverse-stranded (dUTP) paired-end reads, from the annotated transcripts
  (with condition-dependent transcript usage, which gives rMATS events) and
  from the 3' UTR isoforms (with the isoform usage reversed in the treated
  samples for a third of the multi-isoform genes);
- a few percent of short fragments carry adapter read-through, for fastp.

Run it with the trident environment (it needs pysam and numpy):

    python .test/make_test_data.py \
        --genome GRCh38.primary_assembly.genome.fa \
        --gtf gencode.v49.annotation.gtf \
        --utr-isoforms GRCh38_PolyASite_Ensembl_RefSeq_Filt_Reformat_full.gtf \
        --utr-windows GRCh38_NoOverlap_full.gtf
"""
import argparse
import gzip
import os
import re

import numpy as np
import pysam

CONTIG = "test_chr12"
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
SAMPLES = {
    "sampleA1": ("control", ["L001", "L002"]),
    "sampleA2": ("control", ["L001"]),
    "sampleB1": ("treated", ["L001"]),
    "sampleB2": ("treated", ["L001"]),
}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def attribute(attributes, key):
    match = re.search(rf'{key} "([^"]+)"', attributes)
    return match.group(1) if match else None


# Reference
def read_annotation(gtf, chrom, start, end):
    """Gene, transcript and exon lines of the genes fully inside chrom:start-end."""
    genes, lines = {}, []
    with open(gtf) as handle:
        for line in handle:
            if not line.startswith(chrom + "\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] == "gene" and int(fields[3]) >= start and int(fields[4]) <= end:
                genes[attribute(fields[8], "gene_id")] = attribute(fields[8], "gene_name")
            if fields[2] in ("gene", "transcript", "exon"):
                lines.append(fields)
    return genes, [f for f in lines if attribute(f[8], "gene_id") in genes]


def read_utr_gtf(gtf, chrom, start, end, gene_names):
    """3' UTR GTF lines inside the region whose gene symbol is a kept gene."""
    kept = []
    with open(gtf) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != chrom:
                continue
            if int(fields[3]) < start or int(fields[4]) > end:
                continue
            if attribute(fields[8], "gene_id") in gene_names:
                kept.append(fields)
    return kept


def shift(fields, offset, contig):
    shifted = list(fields)
    shifted[0] = contig
    shifted[3] = str(int(fields[3]) - offset)
    shifted[4] = str(int(fields[4]) - offset)
    return shifted


def write_gtf(path, lines):
    with open(path, "w") as out:
        for fields in lines:
            out.write("\t".join(fields) + "\n")


def interval_seq(genome, fields):
    seq = genome[int(fields[3]) - 1 : int(fields[4])]
    return revcomp(seq) if fields[6] == "-" else seq


# Molecules to simulate reads from
def transcripts(annotation, genome):
    """Spliced sequences of each gene's transcripts: gene -> [(id, seq)]."""
    exons = {}
    for fields in annotation:
        if fields[2] == "exon":
            exons.setdefault(attribute(fields[8], "transcript_id"), []).append(fields)
    by_gene = {}
    for transcript_id, parts in exons.items():
        parts.sort(key=lambda f: int(f[3]))
        seq = "".join(genome[int(f[3]) - 1 : int(f[4])] for f in parts)
        if parts[0][6] == "-":
            seq = revcomp(seq)
        by_gene.setdefault(attribute(parts[0][8], "gene_id"), []).append((transcript_id, seq))
    return by_gene


def simulate_sample(rng, molecules, weights, n_pairs, read_length, error_rate):
    """Read pairs from molecules, in proportion to weight x length."""
    lengths = np.array([len(seq) for seq in molecules])
    p = weights * np.where(lengths >= 150, lengths, 0)
    p = p / p.sum()
    picks = rng.choice(len(molecules), size=n_pairs, p=p)
    pairs = []
    for i in picks:
        seq = molecules[i]
        if rng.random() < 0.03:
            frag_len = int(rng.integers(40, read_length))
        else:
            frag_len = int(np.clip(rng.normal(220, 40), 100, len(seq)))
        start = int(rng.integers(0, len(seq) - frag_len + 1))
        fragment = seq[start : start + frag_len]
        # dUTP: read 1 is antisense to the transcript, read 2 sense
        r1 = (revcomp(fragment) + ADAPTER_R1 + "A" * read_length)[:read_length]
        r2 = (fragment + ADAPTER_R2 + "A" * read_length)[:read_length]
        pairs.append((add_errors(rng, r1, error_rate), add_errors(rng, r2, error_rate)))
    return pairs


def add_errors(rng, read, error_rate):
    read = list(read)
    for pos in np.flatnonzero(rng.random(len(read)) < error_rate):
        read[pos] = "ACGT"[(("ACGT".find(read[pos]) + int(rng.integers(1, 4))) % 4)]
    return "".join(read)


def write_fastq(path, name, reads):
    with gzip.open(path, "wt", compresslevel=9) as out:
        for i, read in enumerate(reads):
            out.write(f"@{name}:{i}\n{read}\n+\n{'F' * len(read)}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--gtf", required=True, help="GENCODE annotation")
    parser.add_argument("--utr-isoforms", required=True, help="3' UTR isoforms (gtf_no_focus)")
    parser.add_argument("--utr-windows", required=True, help="3' UTR windows (gtf_focus)")
    parser.add_argument("--region", default="chr12:56000001-58000000")
    parser.add_argument("--pairs", type=int, default=15000, help="read pairs per sample")
    parser.add_argument("--read-length", type=int, default=75)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "data"))
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    chrom, coords = args.region.split(":")
    start, end = (int(x) for x in coords.split("-"))
    offset = start - 1
    ref_dir = os.path.join(args.outdir, "ref")
    fastq_dir = os.path.join(args.outdir, "fastq")
    os.makedirs(ref_dir, exist_ok=True)
    os.makedirs(fastq_dir, exist_ok=True)

    # Genome slice
    genome = pysam.FastaFile(args.genome).fetch(chrom, start - 1, end).upper()
    with open(os.path.join(ref_dir, "genome.fa"), "w") as out:
        out.write(f">{CONTIG}\n")
        for i in range(0, len(genome), 60):
            out.write(genome[i : i + 60] + "\n")

    # Annotation of the genes fully inside the slice
    genes, annotation = read_annotation(args.gtf, chrom, start, end)
    annotation = [shift(f, offset, CONTIG) for f in annotation]
    write_gtf(os.path.join(ref_dir, "annotation.gtf"), annotation)

    # 3' UTR reference of those genes (Ensembl-style contig names, e.g. "12")
    utr_chrom = chrom.removeprefix("chr")
    gene_names = set(genes.values())
    isoforms = [shift(f, offset, CONTIG) for f in read_utr_gtf(args.utr_isoforms, utr_chrom, start, end, gene_names)]
    windows = [shift(f, offset, CONTIG) for f in read_utr_gtf(args.utr_windows, utr_chrom, start, end, gene_names)]
    window_ids = {attribute(f[8], "transcript_id") for f in windows}
    isoforms = [f for f in isoforms if attribute(f[8], "transcript_id") in window_ids]
    write_gtf(os.path.join(ref_dir, "utr_isoforms.gtf"), isoforms)
    write_gtf(os.path.join(ref_dir, "utr_windows.gtf"), windows)

    # Molecules: annotated transcripts, and 3' UTR isoforms (proximal first)
    by_gene = transcripts(annotation, genome)
    utr_by_gene = {}
    for fields in sorted(isoforms, key=lambda f: int(f[4]) - int(f[3])):
        utr_by_gene.setdefault(attribute(fields[8], "gene_id"), []).append(interval_seq(genome, fields))

    molecules, kind = [], []
    gene_level = {}
    for gene_id, gene_transcripts in by_gene.items():
        gene_level[gene_id] = rng.lognormal(0, 1.5)
        usage = {c: rng.dirichlet(np.ones(len(gene_transcripts))) for c in ("control", "treated")}
        for j, (_, seq) in enumerate(gene_transcripts):
            molecules.append(seq)
            kind.append(("tx", gene_id, j, usage))
    shifted_genes = set()
    for name, seqs in utr_by_gene.items():
        gene_id = next(g for g, n in genes.items() if n == name)
        base = rng.dirichlet(np.ones(len(seqs)) * 2)
        if len(seqs) > 1 and rng.random() < 1 / 3:
            shifted_genes.add(name)
        for j, seq in enumerate(seqs):
            molecules.append(seq)
            kind.append(("utr", gene_id, j, (name, base)))

    def weights(condition):
        w = np.empty(len(molecules))
        for i, (k, gene_id, j, info) in enumerate(kind):
            level = gene_level.get(gene_id, 1.0)
            if k == "tx":
                w[i] = level * info[condition][j]
            else:
                name, base = info
                share = base.copy()
                if condition == "treated" and name in shifted_genes:
                    share = share[::-1]
                w[i] = 3 * level * share[j]
        return w

    # Reads, sample sheets
    samples_tsv = ["sample_name\tcondition"]
    units_tsv = ["sample_name\tunit_name\tfq1\tfq2"]
    for sample, (condition, lanes) in SAMPLES.items():
        pairs = simulate_sample(rng, molecules, weights(condition) * rng.lognormal(0, 0.2, len(molecules)),
                                args.pairs, args.read_length, error_rate=0.002)
        samples_tsv.append(f"{sample}\t{condition}")
        for k, lane in enumerate(lanes):
            lane_pairs = pairs[k :: len(lanes)]
            fq = [os.path.join(fastq_dir, f"{sample}_{lane}_R{r}.fastq.gz") for r in (1, 2)]
            write_fastq(fq[0], f"{sample}_{lane}", [p[0] for p in lane_pairs])
            write_fastq(fq[1], f"{sample}_{lane}", [p[1] for p in lane_pairs])
            units_tsv.append(f"{sample}\t{lane}\tdata/fastq/{os.path.basename(fq[0])}\tdata/fastq/{os.path.basename(fq[1])}")

    config_dir = os.path.join(os.path.dirname(os.path.abspath(args.outdir)), "config")
    os.makedirs(config_dir, exist_ok=True)
    with open(os.path.join(config_dir, "samples.tsv"), "w") as out:
        out.write("\n".join(samples_tsv) + "\n")
    with open(os.path.join(config_dir, "units.tsv"), "w") as out:
        out.write("\n".join(units_tsv) + "\n")

    print(f"{len(genes)} genes, {len(by_gene)} with exons; {len(utr_by_gene)} with 3' UTR isoforms "
          f"({sum(len(s) > 1 for s in utr_by_gene.values())} with several); "
          f"isoform usage reversed in treated for: {', '.join(sorted(shifted_genes))}")


if __name__ == "__main__":
    main()
