"""
Write units.tsv and samples.tsv for a dataset from a folder of paired-end FASTQs.

--pattern is a regex matched against R1 file names, with named groups:
  sample  sample name (all units of a sample are aligned together)
  unit    sequencing unit, e.g. the lane
  read    the read-1 token (e.g. R1); the R2 file has 2 in place of its 1

--metadata adds per-sample columns (e.g. condition) from a CSV/TSV whose
--metadata-id column equals the sample name or starts with "<sample>_".

Example (files like SAMPLE1_L001_R1.fastq.gz):
  python workflow/scripts/make_sample_sheets.py --fastq-dir /path/to/fastq \
    --pattern '^(?P<sample>.+?)_(?P<unit>L\d+)_(?P<read>R1)\.fastq\.gz$' \
    --metadata metadata.csv --metadata-id sample_id --out-dir config/my_dataset
"""
import argparse
import os
import re
import sys

import pandas as pd


# One row per R1 file, with its R2 mate
def find_units(fastq_dir, pattern):
    """One row per R1 file matching `pattern`, with its R2 mate."""
    regex = re.compile(pattern)
    missing_groups = {"sample", "unit", "read"} - set(regex.groupindex)
    if missing_groups:
        sys.exit(f"--pattern needs the named groups: {', '.join(sorted(missing_groups))}")

    rows = []
    for name in sorted(os.listdir(fastq_dir)):
        match = regex.match(name)
        if not match:
            continue
        start, end = match.span("read")
        r2_name = name[:start] + match.group("read").replace("1", "2") + name[end:]
        r2_path = os.path.join(fastq_dir, r2_name)
        if not os.path.isfile(r2_path):
            sys.exit(f"No R2 file for {name} (looked for {r2_name})")
        rows.append(
            {
                "sample_name": match.group("sample"),
                "unit_name": match.group("unit"),
                "fq1": os.path.join(fastq_dir, name),
                "fq2": r2_path,
            }
        )

    if not rows:
        sys.exit(f"No file in {fastq_dir} matches --pattern")
    units = pd.DataFrame(rows).sort_values(["sample_name", "unit_name"])
    duplicated = units.duplicated(["sample_name", "unit_name"], keep=False)
    if duplicated.any():
        sys.exit(
            "These sample/unit pairs match more than one R1 file:\n"
            + units.loc[duplicated].to_string(index=False)
        )
    return units


# Optional per-sample metadata columns
def add_metadata(samples, metadata_file, id_column):
    """Join per-sample columns from `metadata_file` onto `samples`."""
    sep = "," if metadata_file.endswith(".csv") else "\t"
    metadata = pd.read_csv(metadata_file, sep=sep, dtype=str)

    def sample_for(metadata_id):
        hits = [
            s for s in samples["sample_name"]
            if metadata_id == s or metadata_id.startswith(s + "_")
        ]
        if len(hits) > 1:
            sys.exit(f"Metadata ID {metadata_id} matches several samples: {hits}")
        return hits[0] if hits else None

    metadata["sample_name"] = metadata[id_column].map(sample_for)
    unmatched = metadata.loc[metadata["sample_name"].isna(), id_column]
    if len(unmatched):
        print(f"Note: {len(unmatched)} metadata rows match no sample: {list(unmatched)}")
    metadata = metadata.dropna(subset=["sample_name"]).drop(columns=id_column)
    if metadata["sample_name"].duplicated().any():
        sys.exit("Some samples match more than one metadata row.")

    samples = samples.merge(metadata, on="sample_name", how="left")
    no_metadata = samples.loc[samples.drop(columns="sample_name").isna().all(axis=1)]
    if len(no_metadata):
        print(f"Note: no metadata for samples {list(no_metadata['sample_name'])}")
    return samples


# Command line
def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--fastq-dir", required=True, help="folder with the FASTQ files")
    parser.add_argument("--pattern", required=True, help="regex for R1 files, see above")
    parser.add_argument("--out-dir", required=True, help="where to write the two sheets")
    parser.add_argument("--metadata", help="optional CSV/TSV of per-sample columns")
    parser.add_argument("--metadata-id", default="sample_id", help="its sample ID column")
    args = parser.parse_args()

    units = find_units(os.path.abspath(args.fastq_dir), args.pattern)
    samples = units[["sample_name"]].drop_duplicates()
    if args.metadata:
        samples = add_metadata(samples, args.metadata, args.metadata_id)

    os.makedirs(args.out_dir, exist_ok=True)
    units.to_csv(os.path.join(args.out_dir, "units.tsv"), sep="\t", index=False)
    samples.to_csv(os.path.join(args.out_dir, "samples.tsv"), sep="\t", index=False)
    print(
        f"Wrote {len(samples)} samples and {len(units)} units to {args.out_dir}"
    )


if __name__ == "__main__":
    main()
