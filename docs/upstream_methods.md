# Upstream calculations and exact scope

The portable driver replays an already generated resource. The following upstream stages
made the inputs. They are documented here so the release is not mistaken for a raw-read pipeline.

| Stage | Executed source retained in author archive | External input | Re-executed by this code? |
|---|---|---|---|
| Haploid joint calling | Original CODE5 source, header and retained command record | 169 BAMs; exact FSP237 FASTA | No |
| Observation atlas | production_sources/observation_atlas | 169 BAMs and indexes; reference; supplied annotation | No |
| Gene/host linkage and sample summaries | Original SD04 and SD05 snapshots | Completed atlas; phenotype workbook; reference annotations | No original atlas regeneration |
| Common observation masks | Original SD06 snapshot and recovered mask files | Frozen SD03 support arrays and SD04 intervals | Archived BED intervals recounted separately |
| Site observation, selected matrix and host/pair joins | Original SD07, SD08 and SD09 snapshots | Joint calls; SD06 masks; sample/host joins | Exported final matrix recounted |
| Catalogue and partitions | production_sources/catalogue | Joint VCF; frozen common-position mask; sample order | No raw-VCF scan |
| Targeted read extraction | production_sources/read_support | 142 BAMs, reference and 660 fixed targets | No BAM extraction |
| Release assembly | production_sources/release_assembly | Completed catalogue and read audit | Previously run on Ceres |
| Targeted audit of 425 sites | Completed SD13 v0.13.1 code and REVIEW run | 142 BAMs, reference and 425 frozen targets | Archived pileups recounted; BAM extraction not repeated |
| Portable verification and manuscript tables | This repository | Frozen resource export, including pileups | Yes |

Original joint calling used bcftools 1.20 / HTSlib 1.20, `mpileup --ignore-RG -q 20 -Q 20
-d 1000 --skip-any-set 3844 -a FORMAT/DP,FORMAT/AD -Ou`, followed by `call -m --ploidy 1
-A -v -V indels --annotate GQ -Ob`. The BAM list, reference and 24 interval files were
provided as `-b`, `-f` and `-R`. Chunk files were concatenated in reference order.
The recovered original CODE5 v0.5.0 source agrees with these header commands.
Its startup also requires frozen CODE4B metadata; the original source snapshot is
provenance, not a self-contained raw-input rebuild. Its provisional base/strict
panels are not the final selected panel described here; the paper uses the raw
joint calls and subsequent observation/call-retention rules.

The atlas used samtools depth 1.17 at MQ/BQ 0/0, 20/20 and 40/20, with flag mask 3844.
Unknown mapping qualities (255) and absent sequence/quality strings were excluded.
Overlapping mates were counted separately; BAQ was not applied in the depth atlas.
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

## Items to complete before raw-input public reproduction

Deposit or link the raw data, exact reference and annotation, masks and atlas arrays;
replace account-specific paths in the original launch configurations with a portable
manifest; run that rebuild against the original outputs and record the comparison.
This is an explicit pending scope extension. The current entry point already works
for the exported resource and is not labeled an end-to-end raw-read rebuild.

Documentation: https://www.htslib.org/doc/1.17/samtools-depth.html and
https://www.htslib.org/doc/1.17/samtools-mpileup.html .

## Upstream dependencies and release separation

The catalogue stage depends on SD07 (position observation and mask membership),
SD08 (frozen selected matrix) and SD09 (sample/host joins and pair summaries).
These inputs must accompany their source snapshots for a future raw-to-resource
rebuild. The atlas also pins earlier pilot and annotation-audit metadata.
Original source is retained under private provenance with original names and
hashes; obsolete launchers are not invoked by this repository.

The full article replay requires the exported resource, the completed
SD13 v0.13.1 REVIEW directory and the recovered original common-mask subset.
The mask check is a separate optional `--maskdata` step: its result does not
replace the original SD13 `MASK_AUDIT_PENDING` status or claim a new depth scan. It never launches the failed v0.13.0 package and
never validates Python runtime caches. Keep `.pyc`, `.pyo` and `__pycache__` out
of release manifests. Exact source and missing-asset inventories accompany the
private author review package.
