# Shared setup: config, sample sheets, output folders, helpers, targets


# Config
if "samples" not in config:
    sys.exit(
        "No dataset config given. Run with "
        "--configfile config/<dataset>/config.yaml (see README.md)."
    )


validate(config, schema="../schemas/config.schema.yaml")


# Sample and unit sheets
for sheet in ["samples", "units"]:
    if not os.path.isfile(config[sheet]):
        sys.exit(f"The {sheet} sheet {config[sheet]} does not exist.")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype=str)
    .set_index("sample_name", drop=False)
    .sort_index()
)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], sep="\t", dtype=str)
validate(units, schema="../schemas/units.schema.yaml")
units = units.set_index(["sample_name", "unit_name"], drop=False).sort_index()

unmatched = set(samples["sample_name"]) ^ set(units["sample_name"])
if unmatched:
    sys.exit(
        "These samples are in only one of the samples and units sheets: "
        + ", ".join(sorted(unmatched))
    )


# Output folders, one per module, under the working directory (--directory)
preprocess_dir = "results/preprocess"
trimmed_dir = f"{preprocess_dir}/trimmed_fastq"
fastp_dir = f"{preprocess_dir}/fastp"
fastqc_raw_dir = f"{preprocess_dir}/fastqc_raw"
fastqc_trimmed_dir = f"{preprocess_dir}/fastqc_trimmed"

align_dir = "results/align"
star_pass1_dir = f"{align_dir}/star_pass1"
star_sj_merged_dir = f"{align_dir}/star_sj_merged"
star_pass2_dir = f"{align_dir}/star_pass2"

expression_dir = "results/expression"
bam_namesorted_dir = f"{expression_dir}/bam_namesorted"
featurecounts_dir = f"{expression_dir}/featurecounts"
htseq_dir = f"{expression_dir}/htseq"

splicing_dir = "results/splicing"
rmats_prep_dir = f"{splicing_dir}/rmats_prep"
rmats_dir = f"{splicing_dir}/rmats"

apa_dir = "results/apa"
kallisto_dir = f"{apa_dir}/kallisto"
apatools_dir = f"{apa_dir}/apatools"

multiqc_report = "results/multiqc_report.html"


# Files shared between modules
star_bam = f"{star_pass2_dir}/{{sample}}_Aligned.sortedByCoord.out.bam"
namesorted_bam = f"{bam_namesorted_dir}/{{sample}}_Aligned.sortedByName.out.bam"
featurecounts_file = f"{featurecounts_dir}/feature_counts.txt"
htseq_file = f"{htseq_dir}/{{sample}}_counts.txt"
apa_gene_list_file = f"{apa_dir}/reliably_expressed_genes.csv"


wildcard_constraints:
    sample="|".join(map(re.escape, samples["sample_name"])),
    unit="|".join(map(re.escape, units["unit_name"])),
    direction="R1|R2",


# Strandedness, translated into each tool's flag
STRAND_FLAGS = {
    "reverse": {
        "featurecounts": 2,
        "htseq": "reverse",
        "rmats": "fr-firststrand",
        "kallisto": "--rf-stranded",
    },
    "forward": {
        "featurecounts": 1,
        "htseq": "yes",
        "rmats": "fr-secondstrand",
        "kallisto": "--fr-stranded",
    },
    "unstranded": {
        "featurecounts": 0,
        "htseq": "no",
        "rmats": "fr-unstranded",
        "kallisto": "",
    },
}
strand = STRAND_FLAGS[config["strandedness"]]


# Input functions
FASTQ_COLUMN = {"R1": "fq1", "R2": "fq2"}


def raw_fastq(wildcards):
    """The raw FASTQ of one unit and read direction."""
    unit = units.loc[(wildcards.sample, wildcards.unit)]
    return unit[FASTQ_COLUMN[wildcards.direction]]


def raw_fastqs(direction):
    """All raw FASTQs of a sample for one read direction, one per unit."""

    def raw_fastqs_of_sample(wildcards):
        return list(units.loc[wildcards.sample, FASTQ_COLUMN[direction]])

    return raw_fastqs_of_sample


def trimmed_fastqs(direction):
    """All trimmed FASTQs of a sample for one read direction, one per unit."""

    def trimmed_fastqs_of_sample(wildcards):
        return expand(
            f"{trimmed_dir}/{wildcards.sample}_{{unit}}_{direction}.fastq.gz",
            unit=units.loc[wildcards.sample, "unit_name"],
        )

    return trimmed_fastqs_of_sample


def mapped_fastqs(direction):
    """The FASTQs that STAR and kallisto read: trimmed or raw (run_trimming)."""
    return branch(
        config["run_trimming"],
        then=trimmed_fastqs(direction),
        otherwise=raw_fastqs(direction),
    )


# Resource functions
def star_mem_mb(wildcards, attempt):
    """STAR's memory: star_mapping.mem_mb, times the attempt number on retries."""
    return config["star_mapping"]["mem_mb"] * attempt


# Parameter functions
def output_dir(wildcards, output):
    """The folder of the rule's first output."""
    return os.path.dirname(output[0])


def star_prefix(wildcards, output):
    """STAR's --outFileNamePrefix: <output folder>/<sample>_."""
    return os.path.join(os.path.dirname(output[0]), f"{wildcards.sample}_")


def star_reads(wildcards, input):
    """STAR's --readFilesIn: comma-separated R1 files, then R2 files."""
    return " ".join(",".join(files) for files in [input.fastq1, input.fastq2])


def star_read_groups(wildcards):
    """STAR's --outSAMattrRGline: one read group per unit of the sample."""
    sample = wildcards.sample
    return " , ".join(
        f"ID:{sample}_{unit}\tSM:{sample}\tPU:{unit}\tLB:{sample}\tPL:ILLUMINA"
        for unit in units.loc[sample, "unit_name"]
    )


def kallisto_reads(wildcards, input):
    """kallisto's read files: R1 and R2 of each unit, in pairs."""
    return " ".join(f"{r1} {r2}" for r1, r2 in zip(input.fastq1, input.fastq2))


def rmats_prep_workdir(wildcards, output):
    """Scratch folder for one sample's rMATS prep, next to its output."""
    return os.path.join(os.path.dirname(output.rmats), f"{wildcards.sample}_work")


def rmats_post_tmpdir(wildcards, output):
    """Scratch folder for rMATS post, next to the rMATS output folder."""
    return os.path.join(os.path.dirname(os.path.dirname(output.b1)), "rmats_post_tmp")


def comma_joined_bams(wildcards, input):
    """rMATS's BAM list: the BAMs of all samples, comma-separated."""
    return ",".join(input.bams)


def lowercase_bool(value):
    """A config boolean as the true/false the APA scripts expect."""
    return str(value).lower()


# Per-unit lists for expand(..., zip)
sample_unit_pairs = units.index.tolist()

unit_samples = [s for s, u in sample_unit_pairs]
unit_names = [u for s, u in sample_unit_pairs]

read_samples = [s for s, u in sample_unit_pairs for _ in ["R1", "R2"]]
read_units = [u for s, u in sample_unit_pairs for _ in ["R1", "R2"]]
read_directions = ["R1", "R2"] * len(sample_unit_pairs)


# Final targets of each module
expression_targets = [
    featurecounts_file,
    *collect(htseq_file, sample=samples["sample_name"]),
    multiqc_report,
]

RMATS_EVENTS = ["SE", "MXE", "A3SS", "A5SS", "RI"]
RMATS_COUNTS = ["JC", "JCEC"]

splicing_targets = [
    *collect(
        f"{rmats_dir}/{{event}}.MATS.{{count}}.txt",
        event=RMATS_EVENTS,
        count=RMATS_COUNTS,
    ),
    f"{rmats_dir}/sample_order.txt",
]

apa_targets = [
    f"{apatools_dir}/table_of_pairs.csv",
    f"{apatools_dir}/df_counts_filtered_corrected_normalized.csv",
]

all_targets = []
if config["modules"]["expression"]:
    all_targets += expression_targets
if config["modules"]["splicing"]:
    all_targets += splicing_targets
if config["modules"]["apa"]:
    all_targets += apa_targets
    if config["apa_quantify"]["group_column"] not in samples.columns:
        sys.exit(
            f"apa_quantify.group_column is {config['apa_quantify']['group_column']}, "
            f"but {config['samples']} has no such column."
        )
