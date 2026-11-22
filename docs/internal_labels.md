# Identifiers in the resource and frozen validators

Identifiers used by existing hashes, VCFs and joins are retained. Reader-facing
labels explain their meaning without asserting unconfirmed sequencing batches.

| Recorded label | Meaning |
|---|---|
| SD10 | Catalogue construction stage |
| SD11 | Fixed-site read-evidence stage |
| SD12 | Resource assembly stage |
| CODE5 | Haploid joint calling across all 169 alignments |
| SD10_71 | Zero-based row 71 of the emitted biallelic catalogue, not a coordinate |
| NO_ATLAS_SAMPLE / phenotype_only_no_CODE5_sample | Phenotype record without a supplied alignment in this resource |
| legacy_sparse_27 | Supplied sparse-alignment subset, 27 isolates |
| legacy_other | Remaining 20 primary isolates outside the FSP411–532 label range |
| new_411_532 | 122 primary isolates labeled FSP411–532 |
| A90Q30 | Called-allele support at least 90% of max(DP,sumAD), site QUAL at least30; other DP/GQ criteria also apply |
| FROZEN_COMMON_CORE | An observation-mask label; the exact depth, quality, sample count and genomic scope must be read with its metadata |
| genotype_accuracy_assessed_in_SD06 | Historical bookkeeping field, not a measured accuracy or a validation certificate |

A reader should join by FSP identifier and genomic coordinates. Do not infer library
strategy, biological populations or time of collection from these source labels.
The original validator files and historical identifiers are immutable provenance.
Any future public alias must have an explicit mapping and freshly generated manifest;
it must not silently replace keys inside an already hashed export.
