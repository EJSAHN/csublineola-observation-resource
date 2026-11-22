# C. sublineola observation resource

Verification code for the *Colletotrichum sublineola* host-reaction and
genomic-observation resource. It checks archived data, recounts read support and
observation masks, and reproduces the article's tables and figures.

[![Code archive v0.3.2](https://zenodo.org/badge/DOI/10.5281/zenodo.24432601.svg)](https://doi.org/10.5281/zenodo.24432601)

The archived code release is [v0.3.2](https://github.com/EJSAHN/csublineola-observation-resource/tree/v0.3.2).
This branch includes subsequent documentation and command-interface updates;
see [CHANGELOG.md](CHANGELOG.md). The archived release remains available for
reproducing that version.

## Quick start

Python 3.9 or newer is required. Verification uses the Python standard library.
Run the synthetic tests without downloading any research data:

```bash
python -m unittest discover -s tests -v
```

For the full set of numerical tables, provide the resource export, follow-up
read-audit archive and observation masks:

```bash
python run.py --data /path/to/resource --followup /path/to/followup_audit --masks /path/to/observation_masks --out /path/to/new_results
```

Use an unused output directory and run Python without `-O`, which disables some
checks. Inputs are checked before use and remain unchanged.

### Windows PowerShell

Activate a Python 3.9+ environment, then run:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run.ps1 `
  -Data 'D:\projects\csublineola_resource\resource' `
  -Followup 'D:\projects\csublineola_resource\followup_audit' `
  -Masks 'D:\projects\csublineola_resource\observation_masks' `
  -Out 'D:\projects\csublineola_resource\results_001' -Python 'python'
```

### Slurm

Use your allocation and partition names:

```bash
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION run_slurm.sh \
  /absolute/path/to/repository /absolute/path/to/resource \
  /absolute/path/to/new_results /absolute/path/to/followup_audit \
  /absolute/path/to/observation_masks
```

The script requests two CPUs, 8 GB and 30 minutes. The repository path is required
because Slurm executes a spooled copy of the script. Activate the required Python
environment before submission.

## Inputs

| Input | Option | Contents |
|---|---|---|
| Resource export | `--data` | Sample records, variant catalogue, selected genotypes, stratified read-audit pileups and checksum manifests. |
| Follow-up read audit | `--followup` | Completed audit archive with `code/`, `results/`, `runner_status.json` and three SHA-256 manifests. Required for the follow-up audit table. |
| Observation masks | `--masks` | Genome-wide and coding-sequence BED masks, coding-sequence intervals, summaries and the checksum contract. Required for the mask-size check. |

Only `--data` is required for a basic run. Include both other inputs to reproduce
all numerical tables. Directory names can be chosen freely; preserve the files
and relative paths inside each archive.

Research data are distributed separately. Download links for these three inputs
are not yet available in this repository. Access and reuse terms will accompany
the data records supplied by the data owners.

## Checks and outputs

The workflow verifies checksums, phenotype cells, the biallelic catalogue and
agreement between the selected TSV and VCF genotypes. A separate parser recounts
saved pileups and their read-support classifications. Further summaries describe
allele counts, missingness, pairwise differences and support within carrier groups.

| Output | Contents |
|---|---|
| `verification/` | Resource and pileup verification reports. |
| `tables/` | Numerical tables for the article and figure generation. |
| `panel_structure/` | Allele counts, missingness and pairwise summaries. |
| `targeted_audit/` | Follow-up audit verification; its summary tables are in `paper_summary/`. |
| `recovered_mask/` | Observation-mask checksums, interval counts and intersections. |
| `run_status.json` | Overall status, software environment and component results. |

The follow-up audit and mask directories are created when their respective inputs are supplied.
A failed check stops the run. These checks establish agreement with the archived
records; genotype accuracy requires additional evidence.

## Figures

Install the plotting dependencies from `requirements-figures.txt` in a Python
3.10+ environment, then use the tables produced above:

```bash
python -m pip install -r requirements-figures.txt
python scripts/plot_tables.py --tables /path/to/new_results/tables --out /path/to/new_figures
```

This writes both article figures as PNG and vector PDF. Labels and plotting
ranges are specific to this resource.

## Interpretation

R denotes a resistant host reaction, S a susceptible host reaction and V variable
recorded outcomes. V is neither a missing value nor a probability of disease.
For haploid genotypes, 0 is the reference allele, 1 the alternate allele and a dot
an unretained or missing call. The selected panel combines an observation mask
with a retained-call threshold; the allele calls retain their original calling
filters.

At 43,389 of 43,814 selected sites, alternate calls occur only in eight isolates
(FSP348, FSP350, FSP351, FSP353, FSP357, FSP360, FSP363 and FSP369). The remaining
134 isolates carry alternate calls at 425 sites: 421 are polymorphic among retained
calls, and four have only alternate retained calls. Under mapping quality 40,
base quality 20 and recalculated base alignment quality, the follow-up audit
supports 1,934 of 1,979 baseline alternate calls in these 134 isolates, including
1,453 of 1,489 at the 421 polymorphic sites. These are read-support results from
the same supplied alignments. The eight-isolate grouping is descriptive and is
not a species or population assignment.

## Methods and validation

The workflow starts from archived records and pileups. Joint calling and
genome-wide observation measurement are described in
[upstream methods](docs/upstream_methods.md). Their BAM-based calculations are
outside this repository's verification workflow.

The original parsers and numerical scripts are preserved, with verifier
identities recorded in [original_verifier_hashes.json](docs/original_verifier_hashes.json)
and [sd13_verifier_identity.json](docs/sd13_verifier_identity.json). The follow-up
parser writes its report on temporary copies. The workflow checks input manifests
again afterwards. See [validation](docs/validation.md) for results and the scope of
the tests.

Existing commands still work: `--sd13` is an alias for `--followup`, `--maskdata`
for `--masks`, and the PowerShell parameters `-SD13` and `-MaskData` are accepted.
`run_ceres.sh` forwards to `run_slurm.sh`; pass your Slurm account and partition at
submission. Stored identifiers and output keys are explained in the
[archive label glossary](docs/internal_labels.md).

## Citation and license

The archived v0.3.2 code is available at
[doi:10.5281/zenodo.24432601](https://doi.org/10.5281/zenodo.24432601).
When using a later commit, record its commit identifier as well.

Original code and documentation are released under
[CC0 1.0 Universal](LICENSE). Separate data, reference assemblies, annotations and
third-party software retain their own terms; see [NOTICE.md](NOTICE.md).
