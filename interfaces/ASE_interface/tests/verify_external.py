"""Prove the candidate lane imports the installed wheel, not a checkout."""

import argparse
from importlib import metadata
from pathlib import Path


def verify_locations(imported, installed, bundled):
    imported = Path(imported).resolve()
    installed = Path(installed).resolve()
    bundled = Path(bundled).resolve()
    if imported.is_relative_to(bundled):
        raise AssertionError(f'source shadow: {imported} is under {bundled}')
    if imported != installed:
        raise AssertionError(f'import is not from installed distribution: {imported} != {installed}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bundled-root', type=Path, required=True)
    parser.add_argument('--candidate-root', type=Path, required=True)
    args = parser.parse_args()

    import abacuslite

    distribution = metadata.distribution('abacuslite')
    installed = distribution.locate_file('abacuslite/__init__.py')
    imported = Path(abacuslite.__file__)
    verify_locations(imported, installed, args.bundled_root)
    if imported.resolve().is_relative_to(args.candidate_root.resolve()):
        raise AssertionError(f'candidate source shadow: {imported}')
    if not installed.is_file():
        raise AssertionError(f'installed distribution lacks abacuslite package: {installed}')
    print(f'abacuslite distribution: {distribution.metadata["Name"]} {distribution.version}')
    print(f'abacuslite import: {imported.resolve()}')
    print(f'abacuslite distribution path: {installed.resolve()}')


if __name__ == '__main__':
    main()
