# Alternative splicing quantification with rMATS


# rMATS prep, per sample; output renamed from rMATS's timestamped name
rule rmats_prep:
    input:
        bam=star_bam,
        gtf=ancient(config["gtf"]),
    output:
        rmats=f"{rmats_prep_dir}/{{sample}}.rmats",
        read_outcomes=f"{rmats_prep_dir}/{{sample}}_read_outcomes_by_bam.txt",
    log:
        "logs/splicing/rmats_prep_{sample}.log",
    benchmark:
        "benchmarks/splicing/rmats_prep_{sample}.txt"
    conda:
        "../envs/trident.yaml"
    threads: config["rmats_prep"]["threads"]
    resources:
        mem_mb=config["rmats_prep"]["mem_mb"],
        runtime=config["rmats_prep"]["runtime"],
    params:
        workdir=rmats_prep_workdir,
        read_length=config["read_length"],
        lib_type=strand["rmats"],
        options=config["rmats"]["options"],
    shell:
        """
        rm -rf {params.workdir} && mkdir -p {params.workdir}/od
        printf '%s' {input.bam} >{params.workdir}/b1.txt
        rmats.py --b1 {params.workdir}/b1.txt --gtf {input.gtf} -t paired \
            --readLength {params.read_length} --libType {params.lib_type} {params.options} \
            --nthread {threads} --od {params.workdir}/od --tmp {params.workdir} \
            --task prep >{log} 2>&1
        mv {params.workdir}/*_0.rmats {output.rmats}
        mv {params.workdir}/*_read_outcomes_by_bam.txt {output.read_outcomes}
        rm -rf {params.workdir}
        """


# rMATS post over all samples as one group (--statoff); records the sample order
rule rmats_post:
    input:
        rmats=collect(
            f"{rmats_prep_dir}/{{sample}}.rmats", sample=samples["sample_name"]
        ),
        bams=collect(star_bam, sample=samples["sample_name"]),
        gtf=ancient(config["gtf"]),
    output:
        b1=f"{rmats_dir}/b1.txt",
        sample_order=f"{rmats_dir}/sample_order.txt",
        tables=collect(
            f"{rmats_dir}/{{event}}.MATS.{{count}}.txt",
            event=RMATS_EVENTS,
            count=RMATS_COUNTS,
        ),
        summary=report(
            f"{rmats_dir}/summary.txt",
            caption="../report/rmats.rst",
            category="Splicing",
        ),
    log:
        "logs/splicing/rmats_post.log",
    benchmark:
        "benchmarks/splicing/rmats_post.txt"
    conda:
        "../envs/trident.yaml"
    threads: config["rmats_post"]["threads"]
    resources:
        mem_mb=config["rmats_post"]["mem_mb"],
        runtime=config["rmats_post"]["runtime"],
    params:
        bams=comma_joined_bams,
        sample_names=list(samples["sample_name"]),
        tmp=rmats_post_tmpdir,
        outdir=output_dir,
        read_length=config["read_length"],
        lib_type=strand["rmats"],
        options=config["rmats"]["options"],
        extra=config["rmats_post"]["extra"],
    shell:
        """
        rm -rf {params.tmp} && mkdir -p {params.tmp}
        for f in {input.rmats}; do ln -s "$(readlink -f "$f")" {params.tmp}/; done
        printf '%s' {params.bams} >{output.b1}
        printf '%s\\n' {params.sample_names} >{output.sample_order}
        rmats.py --b1 {output.b1} --gtf {input.gtf} -t paired \
            --readLength {params.read_length} --libType {params.lib_type} {params.options} \
            --nthread {threads} --od {params.outdir} --tmp {params.tmp} \
            --task post --statoff {params.extra} >{log} 2>&1
        rm -rf {params.tmp}
        """
