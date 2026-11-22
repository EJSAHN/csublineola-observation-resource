# Archive label glossary

Original data files contain stage identifiers that are needed to match records
and verify their checksums. The documentation and command examples use the
analysis names below. The stored keys remain unchanged.

| Recorded label | Meaning |
|---|---|
| CODE4B | Pilot metadata required by the original joint-calling workflow. |
| CODE5 | Haploid joint calling across all 169 alignments. |
| SD03 | Genome-wide observation measurements. |
| SD04 / SD05 | Gene and host-reaction linkage and sample summaries. |
| SD06 | Common observation masks. |
| SD07 | Position observation and mask membership. |
| SD08 / SD09 | Selected genotype matrix and sample/host joins. |
| SD10 | Variant catalogue construction. |
| SD11 | Stratified read-support audit. |
| SD12 | Resource assembly. |
| SD13 | Follow-up read-support audit. |
| SD10_71 | Zero-based row 71 of the emitted catalogue, not a genomic coordinate. |
| NO_ATLAS_SAMPLE / phenotype_only_no_CODE5_sample | Phenotype record without a supplied alignment in this resource. |
| legacy_sparse_27 | The 27 supplied sparse alignments. |
| legacy_other | The 20 primary isolates outside the FSP411–532 identifier range. |
| new_411_532 | The 122 primary isolates labeled FSP411–532. |
| A90Q30 | Called-allele support at least 90% of max(DP, sumAD) and site QUAL at least 30; other depth and genotype-quality criteria also apply. |
| FROZEN_COMMON_CORE | Observation-mask identifier; read its metadata for depth, quality, sample count and genomic scope. |
| genotype_accuracy_assessed_in_SD06 | Historical bookkeeping field; it contains no measured genotype-accuracy estimate. |
| MASK_AUDIT_PENDING | Mask check was pending in the original follow-up archive; the later archived-mask check has a separate report. |

Join biological records by FSP identifier and genomic coordinates. The stage
labels do not establish sequencing strategy, population membership or collection
date.

## Compatibility names

| Earlier name | Current interface |
|---|---|
| `--sd13` / PowerShell `-SD13` | `--followup` / `-Followup` |
| `--maskdata` / PowerShell `-MaskData` | `--masks` / `-Masks` |
| `run_ceres.sh` | `run_slurm.sh`, with account and partition supplied at submission. |
| `scripts/replay_sd13.py --run ...` | `scripts/verify_followup.py --data ...` |

The original `verify_sd11.py` and `verify_sd13.py` filenames, source descriptions
and identity JSON files are retained to preserve comparisons with the executed
parsers. They are implementations called by the public commands, not additional
analysis steps users need to run. Existing JSON keys such as `sd13` and
`manuscript_summary`, output filenames and the `RESOURCE_REPLAY=PASS` status line
also remain stable for downstream scripts.
