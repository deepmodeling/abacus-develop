from pathlib import Path

from ase.build import bulk
from ase.io import write
from ase.mep import NEB


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"


def main():
    INPUTS.mkdir(exist_ok=True)
    initial = bulk("Si", "diamond", a=5.43)
    final = initial.copy()
    final.positions[1, 0] += 0.35

    images = [initial.copy()]
    images.extend(initial.copy() for _ in range(3))
    images.append(final.copy())
    neb = NEB(images, method="improvedtangent")
    neb.interpolate(method="linear")

    write(INPUTS / "init.extxyz", initial)
    write(INPUTS / "final.extxyz", final)
    write(INPUTS / "init_neb_chain.traj", images)
    print(f"Wrote {INPUTS / 'init_neb_chain.traj'}")


if __name__ == "__main__":
    main()
