# C. sublineola observation resource

Version 0.3.1. Resource-replay code for the *Colletotrichum sublineola* host-reaction and
genomic-observation data descriptor.

The code reproduces the sample inventories, host-reaction summaries, SNP-panel accounting
and read-support summaries reported in the article. The resource concerns fungal isolates
and their reactions on sorghum host differentials.

## What the code runs

1. Verifies every file against the resource manifest.
2. Checks all phenotype cells, the complete biallelic catalogue and the selected TSV/VCF
   genotypes, preserving missing values.
3. Parses the archived pileups with a separate parser and reproduces read-support classifications.
4. Generates the numerical source tables used in the article.
5. Recounts alternate-allele concentration, minor-allele counts, per-isolate missingness,
   transition/transversion counts, pairwise comparison denominators and carrier-specific read support.
6. With `--sd13`, verifies the completed follow-up read audit (internal run label SD13) and
   parses its 181,050 support records before generating the follow-up audit table.
7. With `--maskdata`, verifies the archived observation masks against their recorded checksums
   and recounts genome-wide and coding-sequence intervals.

Requirements: Python 3.9 or newer, standard library only. No network download or package
installation is needed. Do not use Python's `-O` option. Use a new output directory for each run.

```bash
python -m unittest discover -s tests -v
python run.py --data /path/to/release_candidate --sd13 /path/to/completed_SD13_run --maskdata /path/to/recovered_mask_data --out /path/to/new_replay
```

Windows PowerShell:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run.ps1 -Data 'D:\projects\csublineola_resource\release_candidate' -Out 'D:\projects\csublineola_resource\replay_001' -Python 'python' -SD13 'D:\projects\csublineola_resource\completed_SD13_run' -MaskData 'D:\projects\csublineola_resource\recovered_mask_data'
```

The Python executable shown is an example; an activated Python 3.9+ environment also works.
The PowerShell and Slurm wrappers call the same Python entry point; they have been syntax
checked but not yet run on Windows or under Slurm.

On a Slurm cluster, submit the wrapper rather than running on a login node:

```bash
sbatch run_ceres.sh /absolute/path/to/this/repository /absolute/path/to/release_candidate /absolute/path/to/new_replay /absolute/path/to/completed_SD13_run /absolute/path/to/recovered_mask_data
```

The wrapper requests two CPUs, 8 GB and 30 minutes. Its first argument is the repository path,
because Slurm runs a spooled copy of the script. Change the account directive for another allocation.

## Inputs and outputs

| Input | Option | Contents |
|---|---|---|
| Resource export | `--data` (required) | The intact `release_candidate` directory: sample records, catalogue, selected panel, stratified read audit, validation and metadata. |
| Follow-up read audit | `--sd13` (optional) | The completed run directory with `code/`, `results/`, `runner_status.json` and three SHA-256 manifests. Required to reproduce the follow-up audit table. Keep it intact. |
| Observation masks | `--maskdata` (optional) | Genome-wide and coding-sequence mask BED files, coding-sequence intervals, mask summaries and the recorded contract that pins their checksums. Required for the mask-size check. |

None of these inputs is included in this repository. Their public repository records,
DOIs and licenses are assigned by the data owners.

Outputs are `verification/`, `tables/`, `panel_structure/`, optional `targeted_audit/`,
optional `recovered_mask/` and `run_status.json`. The follow-up audit table is written to
`targeted_audit/paper_summary/`. A failed check stops the run. A `PASS` confirms the archived
resource and the stated calculations; it does not certify every genomic call.

## Scope

The code reproduces the article's descriptive tables from archived data and recounts the
original pileups with a separate parser. It does not regenerate the observation atlas or the
joint variant calls from BAM files, realign reads or repeat laboratory assays. The upstream
workflows that produced the archived data require the original alignments, reference,
annotation and run metadata; they are described in `docs/upstream_methods.md`.

The export verifiers are byte-identical to the checkers used when the resource was assembled;
their hashes are recorded in `docs/original_verifier_hashes.json`. The follow-up audit parser
is byte-identical to the parser in the completed run (`docs/sd13_verifier_identity.json`). It
runs on a temporary copy, and input manifests are checked before and after the replay.
Summary code is covered by unit tests; see `docs/validation.md`.

No manuscript, correspondence, unpublished reference, phenotype table or genotype matrix is
included here. The continuous-integration workflow uses synthetic test data only; a full-data
check requires the separate resource export.

## Figures

After a successful full-data replay, figures can be regenerated with Python 3.10+ and the
plotting library version recorded in `requirements-figures.txt`:

```bash
python scripts/plot_tables.py --tables /path/to/new_replay/tables --out /path/to/new_figures
```

This optional step writes both article figures as PNG and vector PDF. Figure labels and
plotting ranges are specific to this resource.

## Data semantics

R, resistant host reaction; S, susceptible host reaction; V, variable recorded outcomes.
V is neither a missing value nor a probability of disease. Haploid genotype 0 is the reference
allele, 1 is the alternate allele and a dot is an unretained or missing call. The selected
panel combines a position-level observation mask with a retained-call threshold; it is not a
new genotype callset at mapping quality 40. See the article, `docs/upstream_methods.md` and
the data dictionary in the data release.

## Variation in the selected panel

At 43,389 of 43,814 selected sites, alternate calls occur only in eight isolates (FSP348,
FSP350, FSP351, FSP353, FSP357, FSP360, FSP363 and FSP369). This grouping is defined in the
summary script for reporting; it is not a taxonomic or population assignment. The other
134 isolates carry alternate calls at 425 selected sites: 421 are polymorphic among their
retained calls, and at four sites all retained calls are alternate. In the follow-up audit,
under mapping quality 40 and base quality 20 with recalculated base alignment quality,
1,934 of 1,979 baseline alternate calls in these 134 isolates retain alternate support,
including 1,453 of 1,489 at the 421 polymorphic sites. These figures describe read support
in the same alignments, not independent sequencing validation.

Internal identifiers are kept as join keys because existing checksums and VCF records depend
on them; `docs/internal_labels.md` explains each one. Join records by FSP identifier and
genomic coordinates.

## Release status

Repository URL, archived code DOI, data DOI and code license are pending; see
`LICENSE_STATUS.md` before the first public release.
