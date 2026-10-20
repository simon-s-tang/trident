TRIDENT quantified gene expression, alternative splicing and alternative
polyadenylation for {{ snakemake.config["samples"] }}. Paired-end reads were
trimmed with fastp (``run_trimming``: {{ snakemake.config["run_trimming"] }}),
aligned with STAR in two passes, and counted with featureCounts and htseq-count;
splicing events were quantified with rMATS (``--statoff``) and 3' UTR isoforms
with kallisto and apatools. No differential tests were run.
