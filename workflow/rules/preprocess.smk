# Preprocessing: adapter trimming (fastp) and read QC (FastQC)


# Adapter trimming, per sample and unit
rule fastp:
    input:
        r1=lookup(
            within=units,
            query="sample_name == '{sample}' & unit_name == '{unit}'",
            cols="fq1",
        ),
        r2=lookup(
            within=units,
            query="sample_name == '{sample}' & unit_name == '{unit}'",
            cols="fq2",
        ),
    output:
        r1=f"{trimmed_dir}/{{sample}}_{{unit}}_R1.fastq.gz",
        r2=f"{trimmed_dir}/{{sample}}_{{unit}}_R2.fastq.gz",
        json=f"{fastp_dir}/{{sample}}_{{unit}}_fastp.json",
        html=f"{fastp_dir}/{{sample}}_{{unit}}_fastp.html",
    log:
        "logs/preprocess/fastp_{sample}_{unit}.log",
    benchmark:
        "benchmarks/preprocess/fastp_{sample}_{unit}.txt"
    conda:
        "../envs/trident.yaml"
    threads: config["fastp"]["threads"]
    resources:
        mem_mb=config["fastp"]["mem_mb"],
        runtime=config["fastp"]["runtime"],
    params:
        adapter_r1=config["adapters"]["r1"],
        adapter_r2=config["adapters"]["r2"],
        extra=config["fastp"]["extra"],
    shell:
        "(fastp --version; "
        "fastp "
        "--in1 {input.r1} --in2 {input.r2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "--adapter_sequence {params.adapter_r1} "
        "--adapter_sequence_r2 {params.adapter_r2} "
        "{params.extra} "
        "--json {output.json} --html {output.html} "
        "--thread {threads}) > {log} 2>&1"


# FastQC on raw reads
rule fastqc_raw:
    input:
        fastq=raw_fastq,
    output:
        zip=f"{fastqc_raw_dir}/{{sample}}_{{unit}}_{{direction}}_fastqc.zip",
        html=f"{fastqc_raw_dir}/{{sample}}_{{unit}}_{{direction}}_fastqc.html",
    log:
        "logs/preprocess/fastqc_raw_{sample}_{unit}_{direction}.log",
    benchmark:
        "benchmarks/preprocess/fastqc_raw_{sample}_{unit}_{direction}.txt"
    conda:
        "../envs/trident.yaml"
    threads: config["fastqc"]["threads"]
    resources:
        mem_mb=config["fastqc"]["mem_mb"],
        runtime=config["fastqc"]["runtime"],
    params:
        outdir=output_dir,
    shell:
        "(fastqc --version; "
        "zcat {input.fastq} | fastqc -o {params.outdir} -t {threads} "
        "stdin:{wildcards.sample}_{wildcards.unit}_{wildcards.direction}) > {log} 2>&1"


# FastQC on trimmed reads
rule fastqc_trimmed:
    input:
        fastq=f"{trimmed_dir}/{{sample}}_{{unit}}_{{direction}}.fastq.gz",
    output:
        zip=f"{fastqc_trimmed_dir}/{{sample}}_{{unit}}_{{direction}}_fastqc.zip",
        html=f"{fastqc_trimmed_dir}/{{sample}}_{{unit}}_{{direction}}_fastqc.html",
    log:
        "logs/preprocess/fastqc_trimmed_{sample}_{unit}_{direction}.log",
    benchmark:
        "benchmarks/preprocess/fastqc_trimmed_{sample}_{unit}_{direction}.txt"
    conda:
        "../envs/trident.yaml"
    threads: config["fastqc"]["threads"]
    resources:
        mem_mb=config["fastqc"]["mem_mb"],
        runtime=config["fastqc"]["runtime"],
    params:
        outdir=output_dir,
    shell:
        "(fastqc --version; "
        "zcat {input.fastq} | fastqc -o {params.outdir} -t {threads} "
        "stdin:{wildcards.sample}_{wildcards.unit}_{wildcards.direction}) > {log} 2>&1"
