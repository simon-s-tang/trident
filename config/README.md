# Configuration

**The files you edit:** your dataset config and its two sample sheets, copied
from `config/example_dataset/` to `config/<dataset>/`. Leave the defaults file,
`config/config.yaml`, and everything in `workflow/` unchanged: any default can
be overridden from the dataset config, and your settings then survive updates
to TRIDENT.

Settings come in two files:

| File | What it sets | Defaults |
|---|---|---|
| `config/config.yaml`: the defaults (don't edit) | every tool option; each rule's threads, memory (`mem_mb`, MB) and run time (`runtime`, minutes); and which modules run (`modules`: `expression`, `splicing`, `apa`) | yes |
| `config/<dataset>/config.yaml` (yours) | the dataset: reference genome and GTF, STAR index folder, 3′ UTR reference (APA only), sample sheets, read length, strandedness, adapters, whether to trim. It can also override any key of `config/config.yaml`. | none: every dataset sets its own |

Copy `config/example_dataset/` to `config/<dataset>/` to start. Unknown keys are
rejected, so a typo in an override stops the run with an error. Use absolute
paths in a dataset config: relative paths are read from the run folder.

## Run folder

Each run writes into the folder given with `--directory`:
`results/` (outputs), `logs/`, `benchmarks/` and `.snakemake/`. Use one run
folder per dataset.

## Threads, memory and run time

Each rule's section of `config/config.yaml` sets its `threads`, `mem_mb` and
`runtime`. To change them for one dataset, override them in its dataset
config, e.g.

```yaml
star_mapping:
  mem_mb: 64000
```

To change them for one run only, use the command line:

```bash
snakemake ... --set-threads star_pass2=16 --set-resources star_pass2:mem_mb=64000
```

STAR's memory grows with each retry (`star_mapping.retries`, default 3: 48 GB,
then 96, 144, 192 GB), so an out-of-memory STAR job is resubmitted bigger.

## Sample sheets

- `samples.tsv`: `sample_name`, plus any metadata columns. The APA module
  groups samples by the column set in `apa_quantify.group_column` (default
  `condition`).
- `units.tsv`: `sample_name`, `unit_name`, `fq1`, `fq2`, with one row per
  sample and sequencing unit (lane). Units of the same sample are aligned
  together.

`workflow/scripts/make_sample_sheets.py` builds both from a FASTQ folder,
using a regex for the file names (`--help` shows an example).

## 3′ UTR reference (APA module)

`apa_reference` needs two GTFs: `gtf_no_focus` (the 3′ UTR isoforms) and
`gtf_focus` (the windows that are counted). `kallisto_index` is the path of the
windows' kallisto index; TRIDENT builds it there from `gtf_focus` and the
genome if it does not exist yet. See the README for what the reference is. If
you don't have one, set `modules: {apa: false}`.
