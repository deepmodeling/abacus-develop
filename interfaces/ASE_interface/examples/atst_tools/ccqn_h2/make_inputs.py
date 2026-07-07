from pathlib import Path

from ase import Atoms
from ase.io import write


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"


def main():
    INPUTS.mkdir(exist_ok=True)
    atoms = Atoms(
        "H2",
        positions=[(4.0, 4.0, 3.50), (4.0, 4.0, 4.50)],
        cell=[8.0, 8.0, 8.0],
        pbc=True,
    )
    write(INPUTS / "ts_guess.extxyz", atoms)
    print(f"Wrote {INPUTS / 'ts_guess.extxyz'}")


if __name__ == "__main__":
    main()
