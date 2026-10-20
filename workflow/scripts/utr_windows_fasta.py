"""
Write the sequences of the 3' UTR windows (apa_reference.gtf_focus) as FASTA,
for the kallisto index.

Each window is named by its transcript_id; windows on the minus strand are
reverse-complemented, and a window with several exons is joined in transcript
order. The windows' contig names are matched to the genome's, also when only
one of them uses the "chr" prefix (e.g. Ensembl "1" and GENCODE "chr1", "MT"
and "chrM"). The genome is read one contig at a time, without an index.
"""
import argparse
import re
import sys

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def genome_contig_to_gtf(name, gtf_contigs):
    """The GTF contig name a genome contig corresponds to, or None."""
    candidates = [name, name.removeprefix("chr")]
    if name.startswith("chr"):
        candidates.append("MT" if name == "chrM" else None)
    else:
        candidates += [f"chr{name}", "chrM" if name == "MT" else None]
    return next((c for c in candidates if c in gtf_contigs), None)


def read_windows(gtf):
    """transcript_id -> (contig, strand, [(start, end), ...]), in GTF order."""
    transcript_id = re.compile(r'transcript_id "([^"]+)"')
    windows = {}
    with open(gtf) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            tid = transcript_id.search(fields[8]).group(1)
            contig, strand, parts = windows.setdefault(tid, (fields[0], fields[6], []))
            parts.append((int(fields[3]), int(fields[4])))
    return windows


def read_genome(path):
    """(name, sequence) for each contig of a FASTA file."""
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.rstrip("\n"))
    if name is not None:
        yield name, "".join(chunks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True, help="genome FASTA")
    parser.add_argument("--gtf", required=True, help="3' UTR windows GTF (gtf_focus)")
    parser.add_argument("--out", required=True, help="output FASTA")
    args = parser.parse_args()

    windows = read_windows(args.gtf)
    gtf_contigs = {contig for contig, _, _ in windows.values()}
    print(f"{len(windows)} windows on {len(gtf_contigs)} contigs in {args.gtf}")

    sequences = {}
    found = set()
    for name, seq in read_genome(args.genome):
        contig = genome_contig_to_gtf(name, gtf_contigs)
        if contig is None:
            continue
        found.add(contig)
        for tid, (wcontig, strand, parts) in windows.items():
            if wcontig != contig:
                continue
            window = "".join(seq[start - 1 : end] for start, end in sorted(parts)).upper()
            sequences[tid] = window.translate(COMPLEMENT)[::-1] if strand == "-" else window

    missing = gtf_contigs - found
    if missing:
        sys.exit(f"Contigs of {args.gtf} not found in {args.genome}: {', '.join(sorted(missing))}")

    with open(args.out, "w") as out:
        for tid in windows:
            out.write(f">{tid}\n{sequences[tid]}\n")
    print(f"Wrote {len(sequences)} window sequences to {args.out}")


if __name__ == "__main__":
    main()
