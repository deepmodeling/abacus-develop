"""Compare matching unprofiled DeltaSpin benchmark runs."""

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load(label: str) -> dict:
    result = json.loads((ROOT / "results" / label / "result.json").read_text())
    assert result["return_code"] == 0 and result["runtime_complete"]
    assert result["threads"] > 0 and result["ranks"] == 1
    assert not result["profiled"]
    return result


def compare(reference: dict, candidate: dict) -> dict:
    assert reference["case"] == candidate["case"]
    assert reference["threads"] == candidate["threads"]
    assert reference["ranks"] == candidate["ranks"]
    assert reference["step_limit_override"] == candidate["step_limit_override"]
    assert reference["input_sha256"] == candidate["input_sha256"]
    assert reference["asset_sha256"] == candidate["asset_sha256"]
    changes = {"energy_eV": abs(reference["energy_eV"] - candidate["energy_eV"])}
    for field in ("forces_eV_per_A", "stress_kbar", "moments_uB", "magnetic_forces_eV_per_uB"):
        assert len(reference[field]) == len(candidate[field]) > 0
        assert all(len(x) == len(y) == 3 for x, y in zip(reference[field], candidate[field]))
        changes[field] = max(abs(a - b) for x, y in zip(reference[field], candidate[field])
                             for a, b in zip(x, y))
    return {
        "reference_seconds": reference["elapsed_seconds"],
        "candidate_seconds": candidate["elapsed_seconds"],
        "speedup": reference["elapsed_seconds"] / candidate["elapsed_seconds"],
        "reference_scf_steps": reference["scf_steps"],
        "candidate_scf_steps": candidate["scf_steps"],
        "reference_converged": reference["native_scf_converged"],
        "candidate_converged": candidate["native_scf_converged"],
        "maximum_absolute_changes": changes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, help="reference label under results/")
    parser.add_argument("--candidate", required=True, help="candidate label under results/")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    reference = load(args.reference)
    candidate = load(args.candidate)
    result = compare(reference, candidate)
    result.update(
        case=reference["case"],
        atom_count=len(reference["forces_eV_per_A"]),
        threads=reference["threads"],
        ranks=reference["ranks"],
        reference_label=args.reference,
        candidate_label=args.candidate,
        reference_peak_gpu_memory_MiB=reference["sampled_peak_gpu_memory_MiB"],
        candidate_peak_gpu_memory_MiB=candidate["sampled_peak_gpu_memory_MiB"],
        step_limit_override=reference["step_limit_override"],
        reference_spin_rms_uB=reference["spin_rms_uB"],
        candidate_spin_rms_uB=candidate["spin_rms_uB"],
    )
    output = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as stream:
            stream.write(output)
    print(output, end="")


if __name__ == "__main__":
    main()
