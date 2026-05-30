#!/usr/bin/env python3
"""Extract timer statistics from ABACUS benchmark output files."""
import os
import re
import sys
from collections import defaultdict

def parse_timer_table(stdout_path):
    """Parse TIME STATISTICS table from stdout file."""
    with open(stdout_path, 'r') as f:
        content = f.read()

    if 'TIME STATISTICS' not in content:
        return None

    # Extract the TIME STATISTICS table section
    match = re.search(r'TIME STATISTICS\n-+\n(.*?)-+\n', content, re.DOTALL)
    if not match:
        return None

    table_text = match.group(1)

    # Parse header
    lines = table_text.strip().split('\n')
    if len(lines) < 2:
        return None

    # data lines are from line 2 onwards (skip the CLASS_NAME header and separator)
    timers = {}
    for line in lines[1:]:  # skip header
        line = line.strip()
        if not line:
            continue
        # Format: CLASS_NAME NAME TIME/s CALLS AVG/s PER/%
        parts = line.split()
        if len(parts) >= 6:
            class_name = parts[0]
            name = parts[1]
            time_s = float(parts[2])
            calls = int(parts[3])
            avg_s = float(parts[4])
            per_pct = float(parts[5])
            key = f"{class_name}::{name}"
            timers[key] = {
                'time_s': time_s,
                'calls': calls,
                'avg_s': avg_s,
                'per_pct': per_pct
            }

    # Also extract total time
    total_match = re.search(r'TOTAL  Time  : (\d+)', content)
    wall_time = int(total_match.group(1)) if total_match else None

    # Extract number of SCF iterations
    scf_iters = len(re.findall(r'#ELEC ITER#', content))

    return {
        'timers': timers,
        'wall_time': wall_time,
        'scf_iters': scf_iters
    }

def compute_stats(values):
    """Compute average and stddev."""
    if not values:
        return None, None
    avg = sum(values) / len(values)
    if len(values) > 1:
        var = sum((v - avg)**2 for v in values) / (len(values) - 1)
        std = var ** 0.5
    else:
        std = 0.0
    return avg, std

def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else "results"

    # Target timer keys for Workflow C
    target_keys = [
        "PW_Basis_K::gatherp_scatters",
        "PW_Basis_K::gathers_scatterp",
        "PW_Basis_K::real2recip",
        "PW_Basis_K::recip2real",
        "PW_Basis::collect_local_pw",
        "PW_Basis::collect_uniqgg",
        "PW_Basis_K::setupIndGk",
        "PW_Basis_K::collect_local_pw",
        "ESolver_KS_PW::before_scf",
        "ESolver_KS_PW::hamilt2rho_single",
        "HSolverPW::solve",
        "Operator::hPsi",
        "total",
    ]

    # Parse run name: np{N}_omp{M}_r{R}.stdout
    run_pattern = re.compile(r'np(\d+)_omp(\d+)_r(\d+)\.stdout$')

    for case_name in sorted(os.listdir(results_dir)):
        case_dir = os.path.join(results_dir, case_name)
        if not os.path.isdir(case_dir):
            continue

        # Group runs by np x omp
        configs = defaultdict(list)
        for fname in sorted(os.listdir(case_dir)):
            m = run_pattern.match(fname)
            if not m:
                continue
            np, omp, rep = int(m.group(1)), int(m.group(2)), int(m.group(3))
            fpath = os.path.join(case_dir, fname)
            result = parse_timer_table(fpath)
            if result:
                configs[(np, omp)].append(result)

        print(f"\n{'='*80}")
        print(f"Case: {case_name}")
        print(f"{'='*80}")

        for (np, omp) in sorted(configs.keys()):
            runs = configs[(np, omp)]
            n_runs = len(runs)

            print(f"\n  Config: np={np}, omp={omp} ({n_runs} runs)")
            print(f"  {'Timer':<40} {'Avg/s':>8} {'Std':>8} {'Calls':>8} {'%':>6}")
            print(f"  {'-'*40} {'-'*8} {'-'*8} {'-'*8} {'-'*6}")

            wall_times = [r['wall_time'] for r in runs if r['wall_time'] is not None]
            if wall_times:
                avg_wt, std_wt = compute_stats(wall_times)
                print(f"  {'[Wall Time]':<40} {avg_wt:>8.1f} {std_wt:>8.1f}")

            for key in target_keys:
                values = []
                calls_list = []
                pct_list = []
                for r in runs:
                    if key in r['timers']:
                        values.append(r['timers'][key]['time_s'])
                        calls_list.append(r['timers'][key]['calls'])
                        pct_list.append(r['timers'][key]['per_pct'])

                if values:
                    avg_t, std_t = compute_stats(values)
                    avg_c = sum(calls_list) / len(calls_list)
                    avg_p = sum(pct_list) / len(pct_list)
                    # Shorten key for display
                    short_key = key.split('::')[-1] if '::' in key else key
                    cls = key.split('::')[0] if '::' in key else ''
                    display = f"{cls}::{short_key}" if cls else short_key
                    if len(display) > 40:
                        display = display[:39]
                    print(f"  {display:<40} {avg_t:>8.3f} {std_t:>8.3f} {avg_c:>8.0f} {avg_p:>5.1f}%")

        # Summary table
        print(f"\n  --- Summary: gather/scatter as % of FFT ---")
        for (np, omp) in sorted(configs.keys()):
            runs = configs[(np, omp)]
            for r in runs:
                timers = r['timers']
                gs_key = "PW_Basis_K::gathers_scatterp"
                gp_key = "PW_Basis_K::gatherp_scatters"
                r2r_key = "PW_Basis_K::recip2real"
                r2c_key = "PW_Basis_K::real2recip"

                gp = timers.get(gp_key, {}).get('time_s', 0)
                gs = timers.get(gs_key, {}).get('time_s', 0)
                r2c = timers.get(r2c_key, {}).get('time_s', 0)
                r2r = timers.get(r2r_key, {}).get('time_s', 0)

                gp_pct = (gp / r2c * 100) if r2c > 0 else 0
                gs_pct = (gs / r2r * 100) if r2r > 0 else 0
                break  # just first rep for summary
            print(f"  np={np} omp={omp}: gatherp={gp_pct:.0f}% of real2recip, gathers={gs_pct:.0f}% of recip2real")

if __name__ == '__main__':
    main()
