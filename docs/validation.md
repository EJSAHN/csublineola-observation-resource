# Validation

## Current checks

Automated tests cover phenotype and genotype semantics, pileup parsing, interval
operations, manifest integrity and command compatibility. Test data are synthetic
and contain no research records. The workflow also checks that the original
verifier files retain their recorded SHA-256 hashes.
Repository attributes preserve LF line endings when cloning on Windows, keeping
the verifier files identical across platforms.

The public command names use `--followup` and `--masks`; the earlier options remain
valid. All nine numerical scripts from v0.3.2 are unchanged. The follow-up command
calls the same implementation and retains its output schema.

The updated command and the archived v0.3.2 command were run on the same three
input archives on Linux with Python 3.12.14. Both passed. Of 40 output files,
39 were byte-identical; `run_status.json` differed only in elapsed time
(42.555 versus 42.084 seconds). This comparison verifies the command update
against the archived calculations; it is not a new BAM-based analysis.

## Recorded resource verification

The original full-data verification completed on Linux with Python 3.12.14 in
43.344 seconds. Runtime depends on the machine and storage. The 18 tests supplied
with v0.3.2 passed.

The resource checks covered 471 source files, 1,225 phenotype cells, 493,398
catalogue records, 43,814 selected TSV/VCF records, 6,221,588 genotype cells,
345,250 missing genotypes and 281,160 read-support labels. A separate parser
recounted the saved pileups.

`summarize_panel_structure.py` independently recounts per-isolate alternate and
missing calls, carrier concentration, minor-allele counts, transition/transversion
counts and pairwise differences using shared observation denominators. Two runs
produced the same seven output files. Figure 1 and Figure 2 were generated with
Matplotlib 3.10.8; their numerical tables provide the comparison when rendering
versions produce different image bytes.

## Follow-up read audit

The follow-up verifier checks 60 metadata entries, 43 source entries and 4,415
result entries. The source entries overlap the metadata manifest and are not
additional unique files. The original parser recounts 181,050 support cells
(4,888,350 compared fields), with zero mismatches in the recorded verification.
The five targeted-summary outputs match the separately generated audit tables
byte for byte.

The parser writes reports on temporary copies. Input manifests are checked before
and after execution. Tests reject missing inputs, changed checksums, unlisted
files, paths outside the archive and runtime caches listed as source files.

## Observation masks

The mask verifier checks five resource files against the recorded contract,
recounts BED spans and intersects the genome mask with the union of 42,025
coding-sequence intervals. Recorded counts are 55,005,039 / 71,016,121 reference
A/C/G/T bases and 18,812,552 / 19,682,398 coding-sequence A/C/G/T bases.

This interval check uses the archived masks and their recorded reference
identities. It does not rescan the BAMs, FASTA or annotation. The original
follow-up archive recorded its mask check as pending; the later mask check is
reported separately and leaves that historical status intact.

## Execution scope

The workflow verifies saved records, pileups and intervals. It does not repeat
read alignment, variant calling or laboratory assays. Command-wrapper tests
check argument forwarding and failure propagation; they do not establish
successful execution on a particular Slurm installation.
