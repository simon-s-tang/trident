# Alignment: STAR index and two-pass alignment


# STAR index
rule star_index:
    input:
        genome=ancient(config["genome"]),
        gtf=ancient(config["gtf"]),
    output:
        index=directory(config["star_index_dir"]),
    log:
        "logs/align/star_index.log",
    benchmark:
        "benchmarks/align/star_index.txt"
    conda:
        "../envs/trident.yaml"
    threads: config["star_index"]["threads"]
    resources:
        mem_mb=config["star_index"]["mem_mb"],
        runtime=config["star_index"]["runtime"],
    params:
        sjdb_overhang=config["read_length"] - 1,
        extra=config["star_index"]["extra"],
    shell:
        "(STAR --version; "
        "STAR --runMode genomeGenerate "
        "--runThreadN {threads} "
        "--genomeDir {output.index} "
        "--genomeFastaFiles {input.genome} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.sjdb_overhang} "
        "--outTmpDir {resources.tmpdir}/star_index_$RANDOM "
        "{params.extra}) > {log} 2>&1"


# Pass 1: junction discovery only (no BAM written)
rule star_pass1:
    input:
        index=config["star_index_dir"],
        fastq1=mapped_fastqs("R1"),
        fastq2=mapped_fastqs("R2"),
    output:
        sj=f"{star_pass1_dir}/{{sample}}_SJ.out.tab",
    log:
        "logs/align/star_pass1_{sample}.log",
    benchmark:
        "benchmarks/align/star_pass1_{sample}.txt"
    retries: config["star_mapping"]["retries"]
    conda:
        "../envs/trident.yaml"
    threads: config["star_mapping"]["threads"]
    resources:
        mem_mb=star_mem_mb,
        runtime=config["star_mapping"]["runtime"],
    params:
        prefix=star_prefix,
        reads=star_reads,
        extra=config["star_mapping"]["extra_pass1"],
    shell:
        "(STAR --version; "
        "STAR --genomeDir {input.index} "
        "--readFilesIn {params.reads} "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype None --readFilesCommand zcat "
        "--runThreadN {threads} "
        "{params.extra}) > {log} 2>&1"


# Pool pass-1 junctions across samples
rule sj_merge:
    input:
        sjs=collect(
            f"{star_pass1_dir}/{{sample}}_SJ.out.tab", sample=samples["sample_name"]
        ),
    output:
        merged_sjs=f"{star_sj_merged_dir}/SJ.out.pass1_merged.tab",
    log:
        "logs/align/sj_merge.log",
    conda:
        "../envs/trident.yaml"
    resources:
        mem_mb=config["sj_merge"]["mem_mb"],
        runtime=config["sj_merge"]["runtime"],
    params:
        min_unique_reads=config["sj_merge"]["min_unique_reads"],
    shell:
        "(cat {input.sjs} | awk '$7 >= {params.min_unique_reads}' | cut -f1-4 "
        "| sort -u > {output.merged_sjs}) 2> {log}"


# Pass 2: final coordinate-sorted BAMs, with pooled junctions
rule star_pass2:
    input:
        index=config["star_index_dir"],
        fastq1=mapped_fastqs("R1"),
        fastq2=mapped_fastqs("R2"),
        merged_sjs=f"{star_sj_merged_dir}/SJ.out.pass1_merged.tab",
    output:
        bam=star_bam,
        log_final=f"{star_pass2_dir}/{{sample}}_Log.final.out",
    log:
        "logs/align/star_pass2_{sample}.log",
    benchmark:
        "benchmarks/align/star_pass2_{sample}.txt"
    retries: config["star_mapping"]["retries"]
    conda:
        "../envs/trident.yaml"
    threads: config["star_mapping"]["threads"]
    resources:
        mem_mb=star_mem_mb,
        runtime=config["star_mapping"]["runtime"],
    params:
        prefix=star_prefix,
        reads=star_reads,
        read_groups=star_read_groups,
        out_sam_unmapped=config["star_mapping"]["out_sam_unmapped"],
        limit_sjdb_insert_nsj=config["star_mapping"]["limit_sjdb_insert_nsj"],
        extra=config["star_mapping"]["extra_pass2"],
    shell:
        "(STAR --version; "
        "STAR --genomeDir {input.index} "
        "--readFilesIn {params.reads} "
        "--outSAMattrRGline {params.read_groups} "
        "--outFileNamePrefix {params.prefix} "
        "--outSAMtype BAM SortedByCoordinate --readFilesCommand zcat "
        "--outSAMunmapped {params.out_sam_unmapped} "
        "--runThreadN {threads} "
        "--limitSjdbInsertNsj {params.limit_sjdb_insert_nsj} "
        "--sjdbFileChrStartEnd {input.merged_sjs} "
        "{params.extra}) > {log} 2>&1"
