from pathlib import Path

from ase.build import bulk
from ase.io import write


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"


def main():
    INPUTS.mkdir(exist_ok=True)
    atoms = bulk("Si", "diamond", a=5.43, cubic=True)
    write(INPUTS / "init.traj", atoms)
    print(f"Wrote {INPUTS / 'init.traj'}")


if __name__ == "__main__":
    main()
