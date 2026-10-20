# Gene expression: counts and the dataset-wide MultiQC report


# Name-sorted BAMs for counting
rule bam_sort:
    input:
        bam=star_bam,
    output:
        bam=namesorted_bam,
    log:
        "logs/expression/bam_sort_{sample}.log",
    benchmark:
        "benchmarks/expression/bam_sort_{sample}.txt"
    conda:
        "../envs/trident.yaml"
    threads: config["bam_sort"]["threads"]
    resources:
        mem_mb=config["bam_sort"]["mem_mb"],
        runtime=config["bam_sort"]["runtime"],
    shell:
        "(samtools --version | head -1; "
        "samtools sort -n -@ {threads} {input.bam} -o {output.bam}) > {log} 2>&1"


# featureCounts: gene x sample matrix (-p --countReadPairs counts fragments)
rule featurecounts:
    input:
        bams=collect(namesorted_bam, sample=samples["sample_name"]),
        gtf=ancient(config["gtf"]),
    output:
        counts=featurecounts_file,
        summary=report(
            f"{featurecounts_file}.summary",
            caption="../report/featurecounts.rst",
            category="Expression",
        ),
    log:
        "logs/expression/featurecounts.log",
    benchmark:
        "benchmarks/expression/featurecounts.txt"
    conda:
        "../envs/trident.yaml"
    threads: config["featurecounts"]["threads"]
    resources:
        mem_mb=config["featurecounts"]["mem_mb"],
        runtime=config["featurecounts"]["runtime"],
    params:
        strand=strand["featurecounts"],
        feature_type=config["featurecounts"]["feature_type"],
        attribute=config["featurecounts"]["attribute"],
        extra=config["featurecounts"]["extra"],
    shell:
        "featureCounts -a {input.gtf} -s {params.strand} -p --countReadPairs "
        "-t {params.feature_type} -T {threads} -g {params.attribute} {params.extra} "
        "-o {output.counts} {input.bams} > {log} 2>&1"


# htseq-count: per-sample gene counts (input must be name-sorted)
rule htseq_count:
    input:
        bam=namesorted_bam,
        gtf=ancient(config["gtf"]),
    output:
        counts=htseq_file,
    log:
        "logs/expression/htseq_{sample}.log",
    benchmark:
        "benchmarks/expression/htseq_{sample}.txt"
    conda:
        "../envs/trident.yaml"
    resources:
        mem_mb=config["htseq_count"]["mem_mb"],
        runtime=config["htseq_count"]["runtime"],
    params:
        strand=strand["htseq"],
        feature_type=config["htseq_count"]["feature_type"],
        attribute=config["htseq_count"]["attribute"],
        mode=config["htseq_count"]["mode"],
        nonunique=config["htseq_count"]["nonunique"],
        extra=config["htseq_count"]["extra"],
    shell:
        "htseq-count -f bam -r name -s {params.strand} -t {params.feature_type} "
        "-i {params.attribute} -m {params.mode} --nonunique={params.nonunique} "
        "{params.extra} {input.bam} {input.gtf} > {output.counts} 2> {log}"


# MultiQC over every QC and count file of the run
rule multiqc:
    input:
        fastqc_raw=collect(
            f"{fastqc_raw_dir}/{{sample}}_{{unit}}_{{direction}}_fastqc.zip",
            zip,
            sample=read_samples,
            unit=read_units,
            direction=read_directions,
        ),
        fastp=collect(
            f"{fastp_dir}/{{sample}}_{{unit}}_fastp.json",
            zip,
            sample=unit_samples,
            unit=unit_names,
        ),
        fastqc_trimmed=collect(
            f"{fastqc_trimmed_dir}/{{sample}}_{{unit}}_{{direction}}_fastqc.zip",
            zip,
            sample=read_samples,
            unit=read_units,
            direction=read_directions,
        ),
        star=collect(
            f"{star_pass2_dir}/{{sample}}_Log.final.out", sample=samples["sample_name"]
        ),
        featurecounts=f"{featurecounts_file}.summary",
        htseq=collect(htseq_file, sample=samples["sample_name"]),
        kallisto=collect(
            "logs/apa/kallisto_{sample}.log",
            sample=samples["sample_name"] if config["modules"]["apa"] else [],
        ),
    output:
        report=report(
            multiqc_report,
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/expression/multiqc.log",
    conda:
        "../envs/trident.yaml"
    resources:
        mem_mb=config["multiqc"]["mem_mb"],
        runtime=config["multiqc"]["runtime"],
    params:
        outdir=output_dir,
        extra=config["multiqc"]["extra"],
    shell:
        "multiqc {input} --outdir {params.outdir} --filename multiqc_report.html "
        "--force {params.extra} > {log} 2>&1"
