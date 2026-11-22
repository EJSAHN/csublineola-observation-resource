# Changelog

## 0.3.4 — 2026-11-21

- Replace links to the withdrawn v0.3.2 Zenodo record with the preserved source
  commit, keeping the earlier verification baseline accessible.
- Link the README badge to the Zenodo version collection and explain how to cite
  the DOI for a specific release.
- Include verified Zenodo title and creator metadata in the release.
- Leave numerical scripts, tests, command wrappers and verifier hashes unchanged
  from v0.3.3.

## 0.3.3 — 2026-11-21

- Use `--followup` and `--masks` in commands and examples; keep the previous options
  and PowerShell parameters as compatibility aliases.
- Add `scripts/verify_followup.py` as a descriptive entry point for the existing
  follow-up verifier.
- Add `run_slurm.sh` with allocation and partition supplied at submission;
  retain `run_ceres.sh` as a compatibility entry point.
- Replace development-stage wording in the documentation with descriptions of
  the data and analyses, and document the v0.3.2 source baseline.
- Preserve all nine existing numerical scripts and their recorded verifier hashes.
- Add command-compatibility and wrapper tests, including a Windows PowerShell
  job in continuous integration. Compare a full verification run with v0.3.2:
  all numerical outputs match, with only elapsed time changing in the status file.
- Preserve LF line endings during checkout so Windows Git settings do not alter
  the byte identities of the archived verifiers.

## 0.3.2 — 2026-11-21

- Restore the automated code checks and repository ignore rules.
- Run Python 3.9 and 3.12 checks on Ubuntu 22.04.
- Add the CC0 1.0 dedication and clarify code and data licensing.
- Document input archives, commands and verification scope.
- Preserve the numerical scripts, tests, wrappers and verifier identities from 0.3.1.

Source baseline: [commit bd07093bd81c](https://github.com/EJSAHN/csublineola-observation-resource/tree/bd07093bd81c684dcb1078b973e179d400ee398f).
The original v0.3.2 Zenodo record was withdrawn; the source remains in the Git
history and is the baseline for the comparisons described above.
