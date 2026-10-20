# 3' UTR isoform (APA) quantification: kallisto + apatools (no statistical test)


# kallisto index of the 3' UTR windows; built only if it does not exist
rule kallisto_index:
    input:
        genome=ancient(config["genome"]),
        gtf_focus=ancient(config["apa_reference"]["gtf_focus"]),
    output:
        index=config["apa_reference"]["kallisto_index"],
        fasta=temp(f"{apa_dir}/utr_windows.fa"),
    log:
        "logs/apa/kallisto_index.log",
    benchmark:
        "benchmarks/apa/kallisto_index.txt"
    conda:
        "../envs/trident.yaml"
    threads: config["kallisto_index"]["threads"]
    resources:
        mem_mb=config["kallisto_index"]["mem_mb"],
        runtime=config["kallisto_index"]["runtime"],
    params:
        script=os.path.join(workflow.basedir, "scripts", "utr_windows_fasta.py"),
        extra=config["kallisto_index"]["extra"],
    shell:
        "(python {params.script} --genome {input.genome} --gtf {input.gtf_focus} "
        "--out {output.fasta} && "
        "kallisto index -i {output.index} -t {threads} {params.extra} {output.fasta}) "
        "> {log} 2>&1"


# kallisto, per sample, on the same reads STAR maps
rule kallisto_quant:
    input:
        index=ancient(config["apa_reference"]["kallisto_index"]),
        fastq1=mapped_fastqs("R1"),
        fastq2=mapped_fastqs("R2"),
    output:
        abundance=f"{kallisto_dir}/{{sample}}/abundance.tsv",
        run_info=f"{kallisto_dir}/{{sample}}/run_info.json",
    log:
        "logs/apa/kallisto_{sample}.log",
    benchmark:
        "benchmarks/apa/kallisto_{sample}.txt"
    conda:
        "../envs/trident.yaml"
    threads: config["kallisto_quant"]["threads"]
    resources:
        mem_mb=config["kallisto_quant"]["mem_mb"],
        runtime=config["kallisto_quant"]["runtime"],
    params:
        outdir=output_dir,
        reads=kallisto_reads,
        strand=strand["kallisto"],
        extra=config["kallisto_quant"]["extra"],
    shell:
        "kallisto quant --index {input.index} -o {params.outdir} -t {threads} "
        "{params.strand} {params.extra} {params.reads} 2> {log}"


# Reliably expressed genes, from featureCounts or htseq counts
rule apa_gene_list:
    input:
        counts=branch(
            config["apa_gene_list"]["counts"] == "featurecounts",
            then=featurecounts_file,
            otherwise=collect(htseq_file, sample=samples["sample_name"]),
        ),
        gtf=ancient(config["gtf"]),
    output:
        genes=apa_gene_list_file,
    log:
        "logs/apa/apa_gene_list.log",
    benchmark:
        "benchmarks/apa/apa_gene_list.txt"
    conda:
        "../envs/trident.yaml"
    resources:
        mem_mb=config["apa_gene_list"]["mem_mb"],
        runtime=config["apa_gene_list"]["runtime"],
    params:
        script=os.path.join(workflow.basedir, "scripts", "apa_gene_list.py"),
        format=config["apa_gene_list"]["counts"],
        samples=list(samples["sample_name"]),
        min_mean=config["apa_gene_list"]["min_mean"],
    shell:
        "python {params.script} "
        "--format {params.format} --counts {input.counts} "
        "--gtf {input.gtf} --samples {params.samples} --min-mean {params.min_mean} "
        "--out {output.genes} > {log} 2>&1"


# apatools: filtered isoform counts and the table of proximal/distal pairs
rule apa_quantify:
    input:
        kallisto=collect(
            f"{kallisto_dir}/{{sample}}/abundance.tsv", sample=samples["sample_name"]
        ),
        genes=apa_gene_list_file,
        samples=config["samples"],
        gtf_focus=ancient(config["apa_reference"]["gtf_focus"]),
        gtf_no_focus=ancient(config["apa_reference"]["gtf_no_focus"]),
    output:
        table_of_pairs=report(
            f"{apatools_dir}/table_of_pairs.csv",
            caption="../report/table_of_pairs.rst",
            category="APA",
        ),
        counts=f"{apatools_dir}/df_counts_filtered_corrected_normalized.csv",
        figures=report(
            directory(f"{apatools_dir}/Figures"),
            patterns=[f"{{name}}.{config['apa_quantify']['figures_extension']}"],
            caption="../report/apa_figures.rst",
            category="APA",
        ),
    log:
        "logs/apa/apa_quantify.log",
    benchmark:
        "benchmarks/apa/apa_quantify.txt"
    conda:
        "../envs/trident.yaml"
    resources:
        mem_mb=config["apa_quantify"]["mem_mb"],
        runtime=config["apa_quantify"]["runtime"],
    params:
        script=os.path.join(workflow.basedir, "scripts", "apa_quantify.py"),
        outdir=output_dir,
        group_column=config["apa_quantify"]["group_column"],
        focus_nt=config["apa_quantify"]["focus_nt"],
        minimal_required_length=config["apa_quantify"]["minimal_required_length"],
        quantile_threshold=config["apa_quantify"]["quantile_threshold"],
        background_quantile=lowercase_bool(
            config["apa_quantify"]["background_quantile"]
        ),
        pct_in_treatment=config["apa_quantify"]["pct_in_treatment"],
        quantile_normalization=lowercase_bool(
            config["apa_quantify"]["quantile_normalization"]
        ),
        figures_extension=config["apa_quantify"]["figures_extension"],
    shell:
        "python {params.script} "
        "--kallisto {input.kallisto} --genes {input.genes} "
        "--samples {input.samples} --group-column {params.group_column} "
        "--gtf-focus {input.gtf_focus} --gtf-no-focus {input.gtf_no_focus} "
        "--focus-nt {params.focus_nt} "
        "--minimal-required-length {params.minimal_required_length} "
        "--quantile-threshold {params.quantile_threshold} "
        "--background-quantile {params.background_quantile} "
        "--pct-in-treatment {params.pct_in_treatment} "
        "--quantile-normalization {params.quantile_normalization} "
        "--figures-extension {params.figures_extension} "
        "--outdir {params.outdir} > {log} 2>&1"
