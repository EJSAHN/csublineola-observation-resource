#!/usr/bin/env python3
"""Verify a follow-up read-audit archive and recount its saved pileups.

The implementation and output schema are shared with the original command so
that archived checksums and existing downstream readers continue to work.
"""
import argparse
import json
from pathlib import Path

from replay_sd13 import replay as verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True, type=Path,
                        help='Completed follow-up read-audit archive')
    parser.add_argument('--out', required=True, type=Path,
                        help='New results directory')
    args = parser.parse_args()
    print(json.dumps(verify(args.data, args.out), indent=2))


if __name__ == '__main__':
    main()
