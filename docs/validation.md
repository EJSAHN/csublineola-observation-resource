# Numerical validation and its scope

Code version 0.3.2 preserves the numerical scripts and tests from version 0.3.1.
The changes concern repository documentation, licensing and the restored CI/ignore files.
The measurements below describe the recorded validation of the underlying replay code;
they are not evidence that every upstream analysis was rerun for this maintenance version.

The frozen-export replay ran on Linux / Python 3.12.14 and passed in 43.344 seconds on the authoring host. This is a measured local time, not a Ceres queue estimate. Eighteen semantic, interval and manifest-integrity unit tests passed. This includes all inputs, the final SD13 result-hash recheck and the recovered-mask interval replay.

The replay verified 471 source files, 1,225 phenotype cells, 493,398 catalogue records, 43,814 selected TSV/VCF records, 6,221,588 genotype cells, 345,250 missing genotypes and 281,160 support labels. An independent parser checks the archived pileups; original BAM extraction was not rerun here.

A separate standard-library implementation, summarize_panel_structure.py, recounts per-isolate ALT and missingness, carrier concentration, minor allele counts, transition/transversion counts, pairwise differences with shared denominators, and audit support by carrier group. Two runs produced identical seven output files. It does not import the production calling or audit engines.

The updated Figure 1 and Figure 2 were generated from the recount tables with Matplotlib 3.10.8. PNG byte identity was checked in that environment. Other plotting versions may produce byte differences; numeric source-table checks remain authoritative. Vector PDF timestamps are not equality targets.

No full-BAM processing, new joint variant calling, independent sequencing validation or wet-lab replication was performed by this code. Ceres and Windows wrappers are supplied; this updated portable driver has not been run by the user on those environments. SD13 completed the targeted read audit and recorded its mask lookup as pending. The exact original SD06 genome/CDS masks were subsequently recovered, verified against the original SD07 contract, and independently recounted through a separate --maskdata step. The historical SD13 status is retained.

The added SD13 wrapper verifies 60 REVIEW metadata entries, 43 source entries and
4,415 result entries; source entries overlap the REVIEW manifest and are not
additional unique files. Its byte-identical independent parser recounts all
181,050 archived support cells (4,888,350 compared fields) with zero mismatches.
The five targeted-summary outputs match the separately generated audit tables
byte for byte. The parser runs on temporary copies because the original parser
writes its verification report. Original BAM extraction was not repeated.

Tests reject missing inputs, changed hashes, unlisted source files, unsafe
manifest paths and cache entries in source manifests. Automatically generated
runtime caches are not release inputs. Tests use temporary synthetic files.

The archived-mask verifier checks five original resource files against the frozen
SD07 contract identity, recounts BED spans and independently intersects the genome
mask with the union of 42,025 CDS intervals. It reports 55,005,039 / 71,016,121
reference ACGT bp and 18,812,552 / 19,682,398 CDS ACGT bp. No BAMs, FASTA or GFF
were reread in this interval replay. Four added tests check interval unions,
half-open intersections, missing inputs and a changed contract hash.
