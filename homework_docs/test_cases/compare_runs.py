#!/usr/bin/env python3
"""Compare original vs rerun results, focusing on cold-start (r1) behavior."""
import os, re, glob

OLD_DIR = "/home/aunixt/abacus-develop/homework_docs/test_cases/results"
NEW_DIR = "/home/aunixt/abacus-develop/homework_docs/test_cases/results_rerun"
CASES = ["gaas_tiny", "gaas_small", "gaas_medium"]
CONFIGS = [
    ("1", "1", "serial"),
    ("1", "4", "omp4"),
    ("1", "8", "omp8"),
    ("2", "2", "mix_np2_omp2"),
    ("2", "4", "mix_np2_omp4"),
    ("4", "2", "mix_np4_omp2"),
]

def get_wall(path):
    """Extract TOTAL Time from stdout."""
    try:
        with open(path) as f:
            for line in f:
                if "TOTAL  Time" in line:
                    parts = line.strip().split()
                    return int(parts[-1])
    except: pass
    return None

def get_timers(path):
    """Extract specific function timers from TIME STATISTICS."""
    timers = {}
    try:
        with open(path) as f:
            lines = f.readlines()
        start = None
        for i, line in enumerate(lines):
            if "TIME STATISTICS" in line:
                start = i
                break
        if start is None:
            return timers
        for line in lines[start:]:
            parts = line.strip().split()
            if len(parts) >= 5:
                name = parts[1]
                time_val = parts[2]
                if name in ["gatherp_scatters", "gathers_scatterp", "real2recip",
                           "recip2real", "hPsi", "diag_once", "hamilt2rho_single",
                           "before_scf", "initialize_psi", "diag_subspace",
                           "veff_pw", "nonlocal_pw"]:
                    try:
                        timers[name] = float(time_val)
                    except ValueError:
                        pass
    except: pass
    return timers

print("=" * 100)
print(f"{'Case':15s} {'Config':20s} {'Wall OLD':>10s} {'Wall NEW':>10s} {'Δwall':>8s} {'gatherp OLD':>12s} {'gatherp NEW':>12s} {'Δgp':>8s} {'gathers OLD':>12s} {'gathers NEW':>12s}")
print("=" * 100)

for case in CASES:
    for np, omp, label in CONFIGS:
        for rep in ["r1", "r2", "r3"]:
            fname = f"np{np}_omp{omp}_{rep}.stdout"
            old_path = os.path.join(OLD_DIR, case, fname)
            new_path = os.path.join(NEW_DIR, case, fname)

            old_wall = get_wall(old_path)
            new_wall = get_wall(new_path)
            old_t = get_timers(old_path)
            new_t = get_timers(new_path)

            old_gp = old_t.get("gatherp_scatters", 0)
            new_gp = new_t.get("gatherp_scatters", 0)
            old_gs = old_t.get("gathers_scatterp", 0)
            new_gs = new_t.get("gathers_scatterp", 0)

            if old_wall is not None and new_wall is not None:
                dwall = new_wall - old_wall
                gp_str = f"{old_gp:.2f}s" if old_gp else "—"
                ngp_str = f"{new_gp:.2f}s" if new_gp else "—"
                gs_str = f"{old_gs:.2f}s" if old_gs else "—"
                ngs_str = f"{new_gs:.2f}s" if new_gs else "—"
                gp_diff = new_gp - old_gp if old_gp and new_gp else 0
                print(f"{case:15s} {label+'_'+rep:20s} {old_wall:>4d}s {new_wall:>4d}s {dwall:>+4d}s "
                      f"{gp_str:>12s} {ngp_str:>12s} {gp_diff:>+7.2f}s "
                      f"{gs_str:>12s} {ngs_str:>12s}")

    # Blank line between cases
    print()

# Summary: cold-start ratio comparison
print("\n" + "=" * 80)
print("COLD-START RATIO (r1/r2 wall) COMPARISON")
print("=" * 80)
print(f"{'Case':15s} {'Config':20s} {'OLD r1/r2':>10s} {'NEW r1/r2':>10s} {'Δ':>8s}")
print("-" * 70)

for case in CASES:
    for np, omp, label in CONFIGS:
        old_r1 = get_wall(os.path.join(OLD_DIR, case, f"np{np}_omp{omp}_r1.stdout"))
        old_r2 = get_wall(os.path.join(OLD_DIR, case, f"np{np}_omp{omp}_r2.stdout"))
        new_r1 = get_wall(os.path.join(NEW_DIR, case, f"np{np}_omp{omp}_r1.stdout"))
        new_r2 = get_wall(os.path.join(NEW_DIR, case, f"np{np}_omp{omp}_r2.stdout"))

        if old_r1 and old_r2 and new_r1 and new_r2:
            old_ratio = old_r1 / old_r2
            new_ratio = new_r1 / new_r2
            delta = new_ratio - old_ratio
            print(f"{case:15s} {label:20s} {old_ratio:>6.2f}×  {new_ratio:>6.2f}×  {delta:>+6.2f}")
