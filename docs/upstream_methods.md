# Data generation and verification

This repository verifies archived data and reproduces their numerical summaries.
The source data were generated through the following stages.

| Stage | Inputs | Verification provided here |
|---|---|---|
| Haploid joint calling | 169 BAMs, FSP237 reference sequence and calling intervals | Genotype representation checked in the exported catalogue and panel. |
| Genome-wide observation measurement | BAMs, indexes, reference and annotation | Archived observation summaries checked; depth scans are not repeated. |
| Host-reaction linkage | Phenotype workbook and sample identifiers | Exported phenotype cells and sample joins checked. |
| Common observation masks | Per-position depth records and coding intervals | Archived BED intervals recounted and intersected. |
| Variant catalogue and selected panel | Joint calls, masks and retained-call criteria | Catalogue partitions, selected genotypes and shared denominators recounted. |
| Stratified read audit | 142 BAMs, reference and selected targets | Saved pileups independently parsed. |
| Follow-up audit of 425 sites | 142 BAMs, reference and fixed targets | Saved pileups independently parsed; carrier-specific support summarized. |
| Resource assembly | Completed data tables and audit records | File manifests and cross-file consistency checked. |

## Calling and observation rules

Original joint calling used bcftools 1.20 / HTSlib 1.20, `mpileup --ignore-RG -q 20 -Q 20
-d 1000 --skip-any-set 3844 -a FORMAT/DP,FORMAT/AD -Ou`, followed by `call -m --ploidy 1
-A -v -V indels --annotate GQ -Ob`. The BAM list, reference and 24 interval files were
provided as `-b`, `-f` and `-R`. Chunk files were concatenated in reference order.
The archived joint-calling source agrees with these header commands and requires
metadata from the preceding pilot. The article uses the raw joint calls followed
by observation and call-retention rules, rather than the caller's provisional panels.

Genome-wide observation measurement used samtools depth 1.17 at MQ/BQ 0/0, 20/20 and 40/20, with flag mask 3844.
Unknown mapping qualities (255) and absent sequence/quality strings were excluded.
Overlapping mates were counted separately; BAQ was not applied in the depth measurements.
The genomic mask for the selected panel required DP >= 5 in >=128 of 142 isolates
at MQ40/BQ20. The allele calls retain the original calling filters.

A genotype was retained when haploid GT was 0 or 1, DP >= 5, GQ >= 20, site QUAL >= 30,
FILTER was PASS or dot, and the called-allele AD was >=90% of max(DP, AD_REF + AD_ALT).
Dot FILTER means an unfiltered record, not a recorded PASS. No maximum depth was imposed.
All 493,398 emitted biallelic SNPs remain in the catalogue; unretained calls remain dot.

The read audit used samtools 1.17 at MQ20/BQ20 without BAQ, MQ40/BQ20 without BAQ and
MQ40/BQ20 with recalculated BAQ. Support needs >=5 observations and >=90% REF or ALT
of all printed observations. The samtools default overlap handling differs from depth.
No maximum pileup depth was imposed. The chosen alleles are read-evidence labels, not
new GQ values or a replacement VCF. Sites without an emitted variant record are not
certified reference by any of these exported variant tables.

## Reproducing data generation

Rebuilding these inputs requires the original reads or alignments, exact reference
and annotation, depth records, masks, source code and configuration files. Input
identities and relative paths need to accompany those files. The verification
workflow here starts with the completed resource export and audit archives.

The catalogue also uses the position-observation records, selected genotype
matrix and sample/host joins. These dependencies are needed together for a
raw-data rebuild. Source identifiers are listed in the
[archive label glossary](internal_labels.md).

## Input integrity

Full verification uses three archives: the resource export, the completed
follow-up audit and the observation masks. The mask check is reported separately
from the historical follow-up status. Original files and identifiers remain
unchanged so their recorded hashes can be checked. Runtime caches (`.pyc`, `.pyo`
and `__pycache__`) are excluded from source manifests.

SAMtools documentation:
[depth](https://www.htslib.org/doc/1.17/samtools-depth.html) and
[mpileup](https://www.htslib.org/doc/1.17/samtools-mpileup.html).
