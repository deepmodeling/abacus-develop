#!/usr/bin/env python3
"""
Example: LCAO workflow with breakpoint support

This example demonstrates how to use the LCAOWorkflow class to run
LCAO calculations with Python-controlled SCF and breakpoint support.

Usage:
    python lcao_workflow_example.py

Requirements:
    - pyabacus with ESolver support
    - Input files (INPUT, STRU, KPT, etc.) in current directory
"""

import numpy as np
from pathlib import Path


def example_basic_scf():
    """
    Basic SCF calculation example.

    Shows how to run a simple SCF calculation and get results.
    """
    from pyabacus.esolver import LCAOWorkflow

    print("=" * 60)
    print("Example 1: Basic SCF Calculation")
    print("=" * 60)

    # Initialize workflow
    workflow = LCAOWorkflow("./", gamma_only=True)
    workflow.initialize()

    # Run SCF
    result = workflow.run_scf(max_iter=100)

    # Print results
    print(result.summary())
    print(f"\nEnergy breakdown:")
    for key, value in result.energy.to_dict().items():
        print(f"  {key}: {value:.8f} Ry")


def example_with_callbacks():
    """
    SCF calculation with callbacks example.

    Shows how to register callbacks to monitor SCF progress
    and inspect state at breakpoints.
    """
    from pyabacus.esolver import LCAOWorkflow

    print("\n" + "=" * 60)
    print("Example 2: SCF with Callbacks")
    print("=" * 60)

    # Initialize workflow
    workflow = LCAOWorkflow("./", gamma_only=True)
    workflow.initialize()

    # Define callback for each iteration
    def print_iteration_info(wf, iter_num):
        energy = wf.energy
        drho = wf.drho
        print(f"  Iter {iter_num:3d}: E = {energy.etot:16.8f} Ry, drho = {drho:.2e}")

    # Define callback for breakpoint before after_scf
    def save_final_state(wf):
        print("\n[Breakpoint] Before after_scf - saving state...")

        # Get charge density
        charge = wf.charge
        if charge.rho.size > 0:
            print(f"  Charge density shape: {charge.rho.shape}")
            print(f"  Total charge: {charge.total_charge():.6f}")
            # Save to file
            np.save("charge_density.npy", charge.rho)
            print("  Saved charge density to charge_density.npy")

        # Get energy
        energy = wf.energy
        print(f"  Total energy: {energy.etot:.8f} Ry")

        # Get Hamiltonian (if available)
        hamiltonian = wf.hamiltonian
        if hamiltonian.nbasis > 0:
            print(f"  Number of basis functions: {hamiltonian.nbasis}")
            print(f"  Number of k-points: {hamiltonian.nks}")

        print("[Breakpoint] State inspection complete\n")

    # Register callbacks
    workflow.register_callback('after_iter', print_iteration_info)
    workflow.register_callback('before_after_scf', save_final_state)

    # Run SCF
    print("\nStarting SCF iterations:")
    result = workflow.run_scf(max_iter=100)

    print(f"\nFinal result: {'Converged' if result.converged else 'Not converged'}")


def example_manual_control():
    """
    Manual SCF control example.

    Shows how to manually control the SCF loop for maximum flexibility.
    """
    from pyabacus.esolver import LCAOWorkflow

    print("\n" + "=" * 60)
    print("Example 3: Manual SCF Control")
    print("=" * 60)

    # Initialize workflow
    workflow = LCAOWorkflow("./", gamma_only=True)
    workflow.initialize()

    # Manual SCF control
    workflow.before_scf(istep=0)

    print("\nManual SCF loop:")
    max_iter = 100
    for iter_num in range(1, max_iter + 1):
        # Run single iteration
        workflow.run_scf_step(iter_num)

        # Get current state
        energy = workflow.energy
        drho = workflow.drho

        print(f"  Iter {iter_num}: E = {energy.etot:.8f} Ry")

        # Custom convergence check or early termination
        if workflow.is_converged:
            print(f"\n  Converged at iteration {iter_num}")
            break

        # Example: Custom breakpoint at iteration 5
        if iter_num == 5:
            print("\n  [Custom breakpoint at iter 5]")
            print(f"    Current energy: {energy.etot:.8f} Ry")
            print(f"    Current drho: {drho:.2e}")
            # Could save intermediate state here

    # Inspect state before finalization
    print("\n[Before after_scf]")
    charge = workflow.charge
    hamiltonian = workflow.hamiltonian
    print(f"  Charge nspin: {charge.nspin}")
    print(f"  Hamiltonian nbasis: {hamiltonian.nbasis}")

    # Finalize
    workflow.after_scf(istep=0)
    print("\nSCF completed.")


def example_multi_k():
    """
    Multi-k calculation example.

    Shows how to run calculations with multiple k-points.
    """
    from pyabacus.esolver import LCAOWorkflow

    print("\n" + "=" * 60)
    print("Example 4: Multi-k Calculation")
    print("=" * 60)

    # Initialize workflow with multi-k
    workflow = LCAOWorkflow("./", gamma_only=False)
    workflow.initialize()

    # Run SCF
    result = workflow.run_scf(max_iter=100)

    print(f"\nNumber of k-points: {workflow.nks}")
    print(f"Number of bands: {workflow.nbands}")

    # Access k-point specific data
    for ik in range(min(workflow.nks, 3)):  # Show first 3 k-points
        kvec = workflow.get_kvec(ik)
        eigenvalues = workflow.get_eigenvalues(ik)
        print(f"\nK-point {ik}: ({kvec[0]:.4f}, {kvec[1]:.4f}, {kvec[2]:.4f})")
        if eigenvalues.size > 0:
            print(f"  Eigenvalues (first 5): {eigenvalues[:5]}")


def example_data_extraction():
    """
    Data extraction example.

    Shows how to extract various data for post-processing.
    """
    from pyabacus.esolver import LCAOWorkflow

    print("\n" + "=" * 60)
    print("Example 5: Data Extraction")
    print("=" * 60)

    workflow = LCAOWorkflow("./", gamma_only=True)
    workflow.initialize()

    # Run SCF
    result = workflow.run_scf(max_iter=100)

    # Extract data
    print("\n1. Energy Data:")
    energy = result.energy
    print(f"   Total energy: {energy.etot:.8f} Ry ({energy.etot * 13.6057:.8f} eV)")
    energy_ev = energy.to_eV()
    print(f"   Band energy: {energy_ev.eband:.8f} eV")

    print("\n2. Charge Density:")
    if result.charge is not None and result.charge.rho.size > 0:
        charge = result.charge
        print(f"   Shape: {charge.rho.shape}")
        print(f"   Min/Max: {charge.rho.min():.6f} / {charge.rho.max():.6f}")

    print("\n3. Hamiltonian Matrices:")
    hamiltonian = workflow.hamiltonian
    if hamiltonian.nbasis > 0:
        print(f"   Number of basis: {hamiltonian.nbasis}")
        print(f"   Number of k-points: {hamiltonian.nks}")
        if len(hamiltonian.Hk) > 0:
            print(f"   H(k=0) shape: {hamiltonian.Hk[0].shape}")

    print("\n4. Density Matrix:")
    dm = workflow.density_matrix
    if dm.nks > 0:
        print(f"   DM dimensions: {dm.nrow} x {dm.ncol}")
        print(f"   Number of k-points: {dm.nks}")


def main():
    """Run all examples."""
    print("PyABACUS LCAO Workflow Examples")
    print("================================\n")

    # Check if input files exist
    input_file = Path("INPUT")
    if not input_file.exists():
        print("Note: INPUT file not found in current directory.")
        print("These examples require ABACUS input files (INPUT, STRU, etc.)")
        print("Please run from a directory with valid input files.\n")
        print("Showing example code structure only...\n")

        # Show code structure without running
        import inspect
        for func in [example_basic_scf, example_with_callbacks,
                     example_manual_control, example_data_extraction]:
            print(f"\n{'=' * 60}")
            print(f"Function: {func.__name__}")
            print("=" * 60)
            print(func.__doc__)
        return

    # Run examples
    try:
        example_basic_scf()
    except Exception as e:
        print(f"Example 1 failed: {e}")

    try:
        example_with_callbacks()
    except Exception as e:
        print(f"Example 2 failed: {e}")

    try:
        example_manual_control()
    except Exception as e:
        print(f"Example 3 failed: {e}")

    try:
        example_data_extraction()
    except Exception as e:
        print(f"Example 5 failed: {e}")


if __name__ == "__main__":
    main()
