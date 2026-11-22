"""Check public command compatibility, wrapper failures and verifier identity."""
import contextlib
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('resource_command', ROOT / 'run.py')
COMMAND = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMMAND)


class CommandChecks(unittest.TestCase):
    def test_old_and_new_options_select_the_same_inputs(self):
        common = ['--data', 'resource with spaces', '--out', 'new results']
        parser = COMMAND.argument_parser()
        old = parser.parse_args(common + ['--sd13', 'followup archive', '--maskdata', 'masks'])
        new = parser.parse_args(common + ['--followup', 'followup archive', '--masks', 'masks'])
        self.assertEqual(vars(old), vars(new))
        self.assertEqual(new.followup, Path('followup archive'))
        self.assertEqual(new.masks, Path('masks'))

    def test_optional_archives_can_be_omitted(self):
        args = COMMAND.argument_parser().parse_args(['--data', 'resource', '--out', 'new'])
        self.assertIsNone(args.followup)
        self.assertIsNone(args.masks)

    def test_conflicting_aliases_are_rejected(self):
        for pair in (['--followup', 'a', '--sd13', 'b'], ['--masks', 'a', '--maskdata', 'b']):
            with self.subTest(pair=pair), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as result:
                    COMMAND.argument_parser().parse_args(['--data', 'x', '--out', 'y'] + pair)
                self.assertEqual(result.exception.code, 2)

    def test_followup_entrypoint_uses_the_original_implementation(self):
        import replay_sd13
        import verify_followup
        self.assertIs(verify_followup.verify, replay_sd13.replay)

    def test_original_verifier_hashes_match(self):
        entries = json.loads((ROOT / 'docs/original_verifier_hashes.json').read_text())
        for entry in entries:
            with self.subTest(file=entry['public_file']):
                self.assertEqual(hashlib.sha256((ROOT / entry['public_file']).read_bytes()).hexdigest(), entry['sha256'])
        followup = json.loads((ROOT / 'docs/sd13_verifier_identity.json').read_text())
        self.assertEqual(hashlib.sha256((ROOT / followup['portable_file']).read_bytes()).hexdigest(), followup['sha256'])


@unittest.skipIf(os.name == 'nt' or shutil.which('bash') is None, 'Requires a POSIX shell')
class SlurmWrapperChecks(unittest.TestCase):
    def invoke(self, script, code=0, followup='followup archive', masks='mask archive'):
        with tempfile.TemporaryDirectory() as tmp:
            temp = Path(tmp)
            fake_bin = temp / 'bin'
            fake_bin.mkdir()
            executable = fake_bin / 'python3'
            executable.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\nexit "$TEST_EXIT"\n')
            executable.chmod(0o755)
            env = dict(os.environ, PATH=str(fake_bin) + os.pathsep + os.environ.get('PATH', ''), TEST_EXIT=str(code))
            return subprocess.run(['bash', str(ROOT / script), str(ROOT), 'resource archive', 'new results', followup, masks], env=env, text=True, capture_output=True)

    def test_both_wrappers_forward_paths_with_spaces(self):
        expected = [str(ROOT / 'run.py'), '--data', 'resource archive', '--out', 'new results', '--followup', 'followup archive', '--masks', 'mask archive']
        for script in ('run_slurm.sh', 'run_ceres.sh'):
            with self.subTest(script=script):
                result = self.invoke(script)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.splitlines(), expected)

    def test_wrappers_propagate_python_failure(self):
        for script in ('run_slurm.sh', 'run_ceres.sh'):
            with self.subTest(script=script):
                self.assertEqual(self.invoke(script, code=17).returncode, 17)

    def test_masks_can_be_used_without_followup(self):
        result = self.invoke('run_slurm.sh', followup='')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), [str(ROOT / 'run.py'), '--data', 'resource archive', '--out', 'new results', '--masks', 'mask archive'])


POWERSHELL = shutil.which('powershell') or shutil.which('pwsh')


@unittest.skipUnless(POWERSHELL, 'Requires PowerShell')
class PowerShellWrapperChecks(unittest.TestCase):
    def invoke(self, legacy=False, code=0):
        with tempfile.TemporaryDirectory(prefix='command check ') as tmp:
            temp = Path(tmp)
            shutil.copyfile(ROOT / 'run.ps1', temp / 'run.ps1')
            (temp / 'run.py').write_text('import json, os, sys\nprint(json.dumps(sys.argv[1:]))\nsys.exit(int(os.environ["TEST_EXIT"]))\n')
            args = [POWERSHELL, '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', str(temp / 'run.ps1'), '-Data', 'resource archive', '-Out', 'new results', '-Python', sys.executable, '-SD13' if legacy else '-Followup', 'followup archive', '-MaskData' if legacy else '-Masks', 'mask archive']
            return subprocess.run(args, env=dict(os.environ, TEST_EXIT=str(code)), text=True, capture_output=True)

    def test_new_and_legacy_parameters_forward_the_same_paths(self):
        expected = ['--data', 'resource archive', '--out', 'new results', '--followup', 'followup archive', '--masks', 'mask archive']
        for legacy in (False, True):
            with self.subTest(legacy=legacy):
                result = self.invoke(legacy=legacy)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(json.loads(result.stdout), expected)

    def test_python_failure_fails_the_wrapper(self):
        self.assertNotEqual(self.invoke(code=17).returncode, 0)


if __name__ == '__main__':
    unittest.main()
