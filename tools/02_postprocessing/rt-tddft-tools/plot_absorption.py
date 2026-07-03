"""
Plot absorption spectrum from ABACUS RT-TDDFT output.
Author: Taoni Bao
Contact: baotaoni@pku.edu.cn

Usage Example:
python plot_absorption.py \
    --direc 0 1 \
    --efield_path ./OUT.ABACUS/efield_0.txt ./OUT.ABACUS/efield_1.txt \
    --response_source dipole \
    --spectrum_type epsilon \
    --signal_path ./OUT.ABACUS/dipole_s1.txt \
    --td_dt 0.002 \
    --step_end 5000 \
    --material_name 'C$_6$H$_6$' \
    --energy_range 0.0 20.0 \
    --wavelength_range 80.0 200.0 \
    --volume 11697.413291555

Parameter Descriptions:
Required Parameters:
  --direc               List of field directions: 0=x, 1=y, 2=z. Matches 1-to-1 with --efield_path.
  --efield_path         List of paths to E-field files. Matches 1-to-1 with --direc.

Optional Parameters:
  --response_source     Source type for external response. Choices: ['dipole', 'current']. Default: 'dipole'.
  --spectrum_type       Target spectrum property to calculate. Choices: ['epsilon', 'sigma']. Default: 'epsilon'.
  --signal_path         Path to the dipole or current file. If omitted, default path is auto-selected based on --response_source.
  --td_dt               Time step size in femtoseconds (fs). Default: 0.002.
  --material_name       Name of the material for figure titles (supports LaTeX syntax). Default: 'C$_6$H$_6$'.
  --step_start          The starting step index for data analysis. Default: 0.
  --step_end            The ending step index for data analysis. Default: 5000.
  --pad_factor          Zero-padding factor multiplied with original signal length for FFT. Default: 10.
  --decay_beta          Exponential damping factor windowing parameter. Default: 5.0.
  --remove_dc           Whether to subtract the mean value to remove the zero-frequency DC component. Default: False.
  --energy_range        X-axis rendering boundaries for energy domain plots [Emin, Emax] in eV. Default: [0.0, 20.0].
  --wavelength_range    X-axis rendering boundaries for wavelength domain plots [WLmin, WLmax] in nm. Default: [80.0, 200.0].
  --ylim_energy         Y-axis viewport clipping limits for energy domain plots [Ymin, Ymax]. Default: None (auto).
  --ylim_wavelength     Y-axis viewport clipping limits for wavelength domain plots [Ymin, Ymax]. Default: None (auto).
  --volume              Cell volume in atomic units (Bohr^3). Only applied to dipole response source. Default: 1.0.
"""

import argparse
from typing import List, Tuple, Optional
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
from scipy.signal import hilbert as hilbert_scipy

# Matplotlib global settings for robust rendering across different environments
plt.rcParams.update({
    "text.usetex": False,
    "legend.fontsize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
})

# ==========================================
# Physical Constants & Unit Conversion Factors
# ==========================================
FREQ2EV = sc.h / sc.eV * 1e15  # Conversion: 1/fs -> eV
# Conversion: V/Å -> a.u.
V_ANGSTROM_TO_AU_EFIELD = (sc.e / (4 * sc.pi * sc.epsilon_0 * sc.physical_constants['Bohr radius'][0]**2) * 1e-10)
AU_TIME = sc.physical_constants['atomic unit of time'][0]
FS_TO_AU_TIME = 1e-15 / AU_TIME  # Conversion: fs -> a.u. of time


# ======================
# Data Loading Utilities
# ======================

def load_signal(filename: str, step_start: int, step_end: int, remove_dc: bool) -> np.ndarray:
    """Load polarization response data (dipole/current). Strict checking is enforced."""
    try:
        data = np.loadtxt(filename)
    except OSError as e:
        raise RuntimeError(f"Error loading signal file {filename}: {e}")

    available_rows = data.shape[0]

    # Enforce strict validation: signal file length must cover the requested step range
    if available_rows < step_end:
        raise ValueError(
            f"Execution halted: Signal file '{filename}' contains only {available_rows} rows, "
            f"which is less than the requested --step_end ({step_end}). "
            f"Please reduce --step_end or ensure the simulation completed fully."
        )

    # Slice data within the targeted window and extract direction components
    sliced = data[step_start:step_end, 1:].T  # shape: (3, n_steps)

    if remove_dc:
        for i in range(3):
            sliced[i] -= np.mean(sliced[i])
    else:
        for i in range(3):
            sliced[i] -= sliced[i, 0]
    return sliced


def load_efield(file_groups: List[List[str]], step_start: int, step_end: int, dt: float, remove_dc: bool) -> np.ndarray:
    """Load and sum E-field profiles, mapping them exactly to absolute simulation time."""
    n_steps = step_end - step_start
    efield = np.zeros((3, n_steps))

    for i, files in enumerate(file_groups):
        for f in files:
            try:
                arr = np.loadtxt(f)
            except OSError as e:
                raise RuntimeError(f"Error loading efield file {f}: {e}")

            if arr.size == 0:
                continue

            # Ensure 2D shape for single-row files to prevent index out of bounds
            if arr.ndim == 1:
                arr = arr.reshape(1, -1)

            times_fs = arr[:, 0]
            fields = arr[:, 1]

            # Map absolute physical time to a 0-based global step index
            # E.g., t = 0.006 fs with dt = 0.002 fs -> global_step = 2
            abs_step_indices = np.round(times_fs / dt).astype(int) - 1

            # Calculate relative indices within the current analysis window [step_start, step_end)
            rel_indices = abs_step_indices - step_start

            # Create boolean mask to filter out data outside the requested window
            valid_mask = (rel_indices >= 0) & (rel_indices < n_steps)

            valid_rel_indices = rel_indices[valid_mask]
            valid_fields = fields[valid_mask]

            # Map valid fields into the pre-allocated zero array.
            # Time slots before e-field starts or after it ends remain naturally zero.
            efield[i, valid_rel_indices] += valid_fields / V_ANGSTROM_TO_AU_EFIELD

    if remove_dc and n_steps > 0:
        for i in range(3):
            # Only subtract DC component if the channel contains a non-zero signal
            if np.any(efield[i] != 0):
                efield[i] -= np.mean(efield[i])

    return efield


def validate_inputs(direc: List[int], efield_path: List[str]) -> None:
    if len(direc) != len(efield_path):
        raise ValueError(f"Input mismatch: Number of --direc ({len(direc)}) != --efield_path ({len(efield_path)})")
    if not all(d in (0, 1, 2) for d in direc):
        raise ValueError("Invalid direction: --direc values must be chosen from 0 (x), 1 (y), or 2 (z)")


def group_files_by_direction(direc: List[int], paths: List[str]) -> List[List[str]]:
    groups = [[], [], []]
    for d, p in zip(direc, paths):
        groups[d].append(p)
    return groups


def zero_pad_and_mask(signal: np.ndarray, dt: float, pad_factor: int = 10, decay_beta: float = 5.0) -> np.ndarray:
    n = len(signal)
    t = np.arange(n) * dt
    mask = np.exp(-decay_beta * t / t[-1]) if t[-1] != 0 else np.ones_like(t)
    padded = np.zeros(n * pad_factor)
    padded[:n] = signal * mask
    return padded


# ======================
# Core Calculator Class
# ======================

class AbsorptionSpectrum:
    def __init__(self, dipole: np.ndarray, efield: np.ndarray, dt: float,
                 response_source: str, spectrum_type: str, pad_factor: int = 10,
                 decay_beta: float = 5.0, volume: float = 1.0):
        self.dipole = dipole
        self.efield = efield
        self.dt = dt
        self.dt_au = self.dt * FS_TO_AU_TIME
        self.response_source = response_source
        self.spectrum_type = spectrum_type
        self.pad_factor = pad_factor
        self.decay_beta = decay_beta
        self.volume = volume
        if self.response_source == "dipole" and self.volume <= 0:
            raise ValueError("volume must be > 0 when --response_source is 'dipole'")
        self.N_orig = dipole.shape[1]
        self.N_pad = self.N_orig * pad_factor

    def compute_complex_alpha(self, direction: int) -> Tuple[np.ndarray, np.ndarray]:
        """Compute the complex response function."""
        d_pad = zero_pad_and_mask(self.dipole[direction], self.dt_au, self.pad_factor, self.decay_beta)
        e_pad = zero_pad_and_mask(self.efield[direction], self.dt_au, self.pad_factor, self.decay_beta)

        # Convert FFT output to physics convention (e^{iwt}) by applying complex conjugate
        d_fft = np.conj(np.fft.fft(d_pad))
        e_fft = np.conj(np.fft.fft(e_pad))

        # Safe complex division using boolean mask to completely bypass complex ufunc mask bugs
        ratio = np.zeros_like(d_fft, dtype=complex)
        non_zero_mask = (e_fft != 0)
        ratio[non_zero_mask] = d_fft[non_zero_mask] / e_fft[non_zero_mask]

        freqs_full = np.fft.fftfreq(self.N_pad, d=self.dt_au)
        omega = 2 * np.pi * freqs_full

        if self.response_source == "dipole":
            if self.spectrum_type == "epsilon":
                alpha_real = 1 + (4 * np.pi * ratio.real) / self.volume
                alpha_imag = (4 * np.pi * ratio.imag) / self.volume
                return alpha_real, alpha_imag
            elif self.spectrum_type == "sigma":
                alpha_real = (omega * ratio.imag) / self.volume
                alpha_imag = (-omega * ratio.real) / self.volume
                return alpha_real, alpha_imag

        elif self.response_source == "current":
            if self.spectrum_type == "sigma":
                return ratio.real, ratio.imag
            elif self.spectrum_type == "epsilon":
                alpha_real = np.zeros_like(ratio.real)
                alpha_imag = np.zeros_like(ratio.real)
                np.divide(4 * np.pi * ratio.imag, omega, out=alpha_real, where=omega != 0)
                alpha_real = 1 - alpha_real
                np.divide(4 * np.pi * ratio.real, omega, out=alpha_imag, where=omega != 0)
                return alpha_real, alpha_imag

        raise ValueError("Invalid setup: Check your response_source and spectrum_type parameters.")

    def get_positive_alpha(self, direction: int) -> Tuple[np.ndarray, np.ndarray]:
        alpha_real, alpha_imag = self.compute_complex_alpha(direction)
        half = self.N_pad // 2
        return alpha_real[:half], alpha_imag[:half]

    def get_energy_axis(self) -> np.ndarray:
        freqs = np.fft.fftfreq(self.N_pad, d=self.dt)[:self.N_pad//2]
        return freqs * FREQ2EV


# ==========================
# Specialized Plotting Tools
# ==========================

def calculate_y_limits(x_vals: np.ndarray, y_arrays: List[np.ndarray], x_range: Tuple[float, float],
                       percentile: float = 99.5, ignore_zero_energy: float = 0.2) -> Optional[Tuple[float, float]]:
    """
    Calculate automatic Y-axis limits using percentile filtering.
    A low-energy threshold is utilized to bypass divergences or numerical artifacts near zero frequency.
    """
    safe_x_min = max(x_range[0], ignore_zero_energy)
    mask = (x_vals >= safe_x_min) & (x_vals <= x_range[1])

    # Fallback to the original requested range if the exclusion window covers the entire domain
    if not np.any(mask):
        mask = (x_vals >= x_range[0]) & (x_vals <= x_range[1])
    if not np.any(mask):
        return (0, 1)

    valid_y = np.concatenate([arr[mask] for arr in y_arrays])
    valid_y = valid_y[np.isfinite(valid_y)]

    if len(valid_y) == 0:
        return (0, 1)

    y_low = np.percentile(valid_y, 100 - percentile)
    y_high = np.percentile(valid_y, percentile)
    margin = 0.1 * (y_high - y_low) if y_high != y_low else 0.1 * max(abs(y_high), 1.0)

    return y_low - margin, y_high + margin


def plot_time_series(ax, time: np.ndarray, data: np.ndarray, directions: List[int], labels: List[str], ylabel: str):
    for d in directions:
        ax.plot(time, data[d], label=labels[d])
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel(ylabel)
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.6)


def plot_fft(ax, energies: np.ndarray, signals: List[np.ndarray], directions: List[int], x_range: Tuple[float, float]):
    labels = {0: "$x$", 1: "$y$", 2: "$z$"}
    for d, sig in zip(directions, signals):
        sig_pos = sig[:len(energies)]
        ax.plot(energies, sig_pos.real, '-', label=f'Re[{labels[d]}]')
        ax.plot(energies, sig_pos.imag, '--', label=f'Im[{labels[d]}]')
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("FFT Amplitude")
    ax.set_xlim(x_range)
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.6)


def generate_triple_plots(x_vals: np.ndarray, alphas_r: List[np.ndarray], alphas_i: List[np.ndarray],
                          directions: List[int], xlabel: str, x_range: Tuple[float, float],
                          ylim: Optional[Tuple[float, float]], title_base: str, filename_base: str):
    """Generate and save three standalone layout variants: Mixed components, Real-only, and Imag-only."""
    axis_labels = {0: "$x$", 1: "$y$", 2: "$z$"}
    is_energy = "Energy" in xlabel
    zero_cutoff = 0.2 if is_energy else 0.0

    plot_configurations = [
        {"suffix": "", "title_append": "", "arrays": alphas_r + alphas_i, "plot_r": True, "plot_i": True},
        {"suffix": "_real", "title_append": " (Real Part)", "arrays": alphas_r, "plot_r": True, "plot_i": False},
        {"suffix": "_imag", "title_append": " (Imaginary Part)", "arrays": alphas_i, "plot_r": False, "plot_i": True},
    ]

    for config in plot_configurations:
        fig, ax = plt.subplots(figsize=(10, 6))
        for idx, d in enumerate(directions):
            base_label = axis_labels[d]
            if config["plot_r"]:
                label_str = f"Re[{base_label}]" if config["plot_i"] else base_label
                ax.plot(x_vals, alphas_r[idx], '-', label=label_str)
            if config["plot_i"]:
                label_str = f"Im[{base_label}]" if config["plot_r"] else base_label
                style = '--' if config["plot_r"] else '-'
                ax.plot(x_vals, alphas_i[idx], style, label=label_str)

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Intensity")
        ax.set_xlim(x_range)

        current_ylim = ylim or calculate_y_limits(x_vals, config["arrays"], x_range, ignore_zero_energy=zero_cutoff)
        if current_ylim:
            ax.set_ylim(current_ylim)

        ax.legend()
        ax.set_title(f"{title_base}{config['title_append']}", fontsize=18, y=1.02)
        ax.grid(True, linestyle=':', alpha=0.6)
        fig.savefig(f"{filename_base}{config['suffix']}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def export_spectrum_data(filename: str, energies: np.ndarray, wavelengths: np.ndarray,
                         data_arrays: List[np.ndarray], component_labels: List[str]):
    """Export calculated physical properties to a structured text file."""
    data_out = np.column_stack([energies, wavelengths] + data_arrays)
    header_str = " ".join(["Energy(eV)", "Wavelength(nm)"] + component_labels)
    np.savetxt(filename, data_out, header=header_str)


def apply_fermi_dirac_taper(energies: np.ndarray, signal: np.ndarray, cutoff_ev: float = 40.0, smear_ev: float = 2.0, baseline: float = 0.0) -> np.ndarray:
    """
    Apply a smooth Fermi-Dirac-like tapering window in the frequency domain.
    Forces the optical response to smoothly decay to a specified vacuum baseline at high energies.
    """
    energy_abs = np.abs(energies)
    taper_window = 1.0 / (1.0 + np.exp((energy_abs - cutoff_ev) / smear_ev))
    return (signal - baseline) * taper_window + baseline


def kk_hilbert_transform(signal_pos: np.ndarray, omega: np.ndarray, is_odd: bool = True) -> np.ndarray:
    N = len(omega)
    if N < 2:
        return np.zeros_like(omega)

    # Construct the full symmetric frequency axis array assuming uniform spacing
    omega_full = np.concatenate((-omega[-1:0:-1], omega))
    signal_full = np.zeros_like(omega_full)

    # Reconstruct the full signal utilizing symmetry properties
    signal_full[N - 1:] = signal_pos
    if is_odd:
        signal_full[:N - 1] = -signal_pos[1:][::-1]
    else:
        signal_full[:N - 1] = signal_pos[1:][::-1]

    # Extend signal with zero-padding to completely eliminate periodic boundary wrap-around
    pad_length = len(signal_full) * 10
    analytic = hilbert_scipy(signal_full, N=pad_length)

    # Extract imaginary part corresponding to the Hilbert transform and truncate padding
    hilbert_trans = analytic.imag[:len(signal_full)]
    return hilbert_trans[N - 1:]


# ======================
# Main Execution Entry
# ======================

def main():
    parser = argparse.ArgumentParser(description="Script for processing ABACUS RT-TDDFT absorption spectra.")

    # Operational configuration switches
    parser.add_argument("--response_source", choices=["dipole", "current"], default="dipole",
                        help="Select polarization physical driving mechanism source. Choices: ['dipole', 'current'].")
    parser.add_argument("--spectrum_type", choices=["epsilon", "sigma"], default="epsilon",
                        help="Select optical calculation properties parameter. Choices: ['epsilon', 'sigma'].")

    # Input file specifications
    parser.add_argument("--signal_path", type=str, default=None,
                        help="File path to polarization input data. Configured dynamically matching response source if omitted.")
    parser.add_argument("--efield_path", type=str, nargs="+", required=True,
                        help="File paths to electric field source files.")

    # Signal windowing parameters
    parser.add_argument("--td_dt", type=float, default=0.002, help="Simulation numerical time step in fs.")
    parser.add_argument("--direc", type=int, nargs="+", required=True,
                        help="Field simulation spatial directions list: 0=x, 1=y, 2=z.")
    parser.add_argument("--material_name", type=str, default=r"C$_6$H$_6$",
                        help="Target description context used inside labels.")
    parser.add_argument("--step_start", type=int, default=0)
    parser.add_argument("--step_end", type=int, default=5000)
    parser.add_argument("--pad_factor", type=int, default=10,
                        help="Zero padding coefficient indicator for FFT scaling.")
    parser.add_argument("--decay_beta", type=float, default=5.0,
                        help="Damping window filter boundary argument parameter.")
    parser.add_argument("--remove_dc", type=bool, default=False,
                        help="Enable filtering out low-frequency background signals.")
    parser.add_argument("--volume", type=float, default=1.0,
                        help="Cell volume in atomic units (Bohr^3). Only applied to dipole response source.")

    # Visual canvas bounding limits properties
    parser.add_argument("--energy_range", type=float, nargs=2,
                        default=[0.0, 20.0], help="X-axis display window bounds inside Energy plots.")
    parser.add_argument("--wavelength_range", type=float, nargs=2,
                        default=[80.0, 200.0], help="X-axis display window bounds inside Wavelength plots.")
    parser.add_argument("--ylim_energy", type=float, nargs=2, default=None,
                        help="Y-axis bounding box limits overrides inside Energy plots.")
    parser.add_argument("--ylim_wavelength", type=float, nargs=2, default=None,
                        help="Y-axis bounding box limits overrides inside Wavelength plots.")

    args = parser.parse_args()

    if args.step_start < 0:
        raise ValueError(f"--step_start must be >= 0 (got {args.step_start})")
    if args.step_end <= args.step_start:
        raise ValueError(f"--step_end ({args.step_end}) must be > --step_start ({args.step_start})")

    # Intelligent fallback tracking path assignments matching runtime selections
    if args.signal_path is None:
        if args.response_source == "dipole":
            args.signal_path = "./OUT.ABACUS/dipole_s1.txt"
        else:
            args.signal_path = "./OUT.ABACUS/current_tot.txt"

    print(f"--> Selected physical response: [{args.response_source}]")
    print(f"--> Selected spectrum variant:  [{args.spectrum_type}]")
    print(f"--> Target parsing data path:    {args.signal_path}")

    validate_inputs(args.direc, args.efield_path)
    Efile_groups = group_files_by_direction(args.direc, args.efield_path)

    # Load arrays via distinct specific mechanisms
    signal_data = load_signal(args.signal_path, args.step_start, args.step_end, args.remove_dc)
    # Load electric fields utilizing absolute time step mapping
    efield_data = load_efield(Efile_groups, args.step_start, args.step_end, args.td_dt, args.remove_dc)

    # Initialize calculation core workspace
    absorber = AbsorptionSpectrum(
        dipole=signal_data, efield=efield_data, dt=args.td_dt,
        response_source=args.response_source, spectrum_type=args.spectrum_type,
        pad_factor=args.pad_factor, decay_beta=args.decay_beta, volume=args.volume
    )

    energies = absorber.get_energy_axis()
    with np.errstate(divide='ignore'):
        wavelengths = np.where(energies > 0, (sc.h * sc.c / sc.eV * 1e9) / energies, np.inf)

    directions = sorted(set(args.direc))

    # Compile array structures across active components
    alphas_real = []
    alphas_imag = []
    for d in range(3):
        r, i = absorber.get_positive_alpha(d)
        alphas_real.append(r)
        alphas_imag.append(i)

    # Export structured raw metrics matrices to external text files
    export_spectrum_data(
        f"{args.spectrum_type}_dat_real.txt", energies, wavelengths, alphas_real,
        ["alpha_real_X", "alpha_real_Y", "alpha_real_Z"]
    )
    export_spectrum_data(
        f"{args.spectrum_type}_dat_imag.txt", energies, wavelengths, alphas_imag,
        ["alpha_imag_X", "alpha_imag_Y", "alpha_imag_Z"]
    )

    time = np.arange(args.step_end - args.step_start) * args.td_dt

    # Adjust terminology and identifiers map dynamically based on driving response origin
    if args.response_source == "dipole":
        time_ylabel, time_legends, file_prefix = "Dipole (a.u.)", ["$P_x$", "$P_y$", "$P_z$"], "dipole"
    else:
        time_ylabel, time_legends, file_prefix = "Current (a.u.)", ["$J_x$", "$J_y$", "$J_z$"], "current"

    # Figure 1: Time domain response analysis
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_time_series(ax, time, signal_data, directions, time_legends, time_ylabel)
    ax.set_xlim(0, args.step_end * args.td_dt)
    fig.savefig(f"{file_prefix}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Figure 2 & 3 series: Generate combinations of graphs automatically across domains
    title_label = r"$\epsilon$" if args.spectrum_type == "epsilon" else r"$\sigma$"
    title_string = rf"{args.material_name} {title_label} Spectrum"

    # Render all energy component variants
    generate_triple_plots(energies, [alphas_real[d] for d in directions], [alphas_imag[d] for d in directions],
                          directions, "Energy (eV)", tuple(args.energy_range), getattr(args, "ylim_energy", None),
                          title_string, args.spectrum_type)

    # Render all wavelength component variants
    generate_triple_plots(wavelengths, [alphas_real[d] for d in directions], [alphas_imag[d] for d in directions],
                          directions, "Wavelength (nm)", tuple(
                              args.wavelength_range), getattr(args, "ylim_wavelength", None),
                          title_string, f"{args.spectrum_type}_wavelength")

    # Figure 4: Incident field fourier response mapping
    fig, ax = plt.subplots(figsize=(10, 6))
    efield_ffts = [np.fft.fft(zero_pad_and_mask(efield_data[d], args.td_dt, args.pad_factor)) for d in directions]
    plot_fft(ax, energies, efield_ffts, directions, tuple(args.energy_range))
    fig.savefig("efield_fourier.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Figure 5: Polarization domain tracking fourier mapping
    fig, ax = plt.subplots(figsize=(10, 6))
    signal_ffts = [np.fft.fft(zero_pad_and_mask(signal_data[d], args.td_dt, args.pad_factor)) for d in directions]
    plot_fft(ax, energies, signal_ffts, directions, tuple(args.energy_range))
    fig.savefig(f"{file_prefix}_fourier.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Figure 6: Incident electric field input time tracking
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_time_series(ax, time, efield_data, directions, ["$E_x$", "$E_y$", "$E_z$"], "Electric Field (a.u.)")
    ax.set_xlim(0, 1.0)
    fig.savefig("efield_time.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Figure 7: Kramers-Kronig mathematical consistency checks validation workspace
    n_dirs = len(directions)
    fig_kk, axes_kk = plt.subplots(n_dirs, 2, figsize=(14, 5 * n_dirs), squeeze=False)

    # Initialize a dedicated calculator to strictly enforce dielectric function (epsilon) evaluation
    kk_absorber = AbsorptionSpectrum(
        dipole=signal_data, efield=efield_data, dt=args.td_dt,
        response_source=args.response_source, spectrum_type="epsilon",
        pad_factor=args.pad_factor, decay_beta=args.decay_beta, volume=args.volume
    )

    for idx, d in enumerate(directions):
        eps_real, eps_imag = kk_absorber.get_positive_alpha(d)

        # Isolate the high-frequency numerical smearing strictly to the KK verification step
        eps_real_smooth = apply_fermi_dirac_taper(energies, eps_real, baseline=1.0)
        eps_imag_smooth = apply_fermi_dirac_taper(energies, eps_imag, baseline=0.0)

        # Perform Kramers-Kronig relations analysis
        eps_real_KK = 1.0 - kk_hilbert_transform(eps_imag_smooth, energies, is_odd=True)
        eps_imag_KK = kk_hilbert_transform(eps_real_smooth - 1.0, energies, is_odd=False)

        # Real domain validation subplot plotting
        ax1 = axes_kk[idx, 0]
        ax1.plot(energies, eps_real, 'b-', label=r"$\mathrm{Re}[\epsilon(\omega)]$ (TDDFT)")
        ax1.plot(energies, eps_real_KK, 'r--', label=r"$1 - \mathcal{H}[\mathrm{Im}(\epsilon)](\omega)$")
        ax1.set_xlabel("Energy (eV)")
        ax1.set_ylabel(r"Re[$\epsilon(\omega)$]")
        ax1.set_xlim(tuple(args.energy_range))

        ylim_real = calculate_y_limits(energies, [eps_real, eps_real_KK], tuple(args.energy_range))
        if ylim_real:
            ax1.set_ylim(ylim_real)

        ax1.legend()
        ax1.grid(True, linestyle=':', alpha=0.6)

        # Imaginary domain validation subplot plotting
        ax2 = axes_kk[idx, 1]
        ax2.plot(energies, eps_imag, 'g-', label=r"$\mathrm{Im}[\epsilon(\omega)]$ (TDDFT)")
        ax2.plot(energies, eps_imag_KK, 'm--', label=r"$\mathcal{H}[\mathrm{Re}(\epsilon) - 1](\omega)$")
        ax2.set_xlabel("Energy (eV)")
        ax2.set_ylabel(r"Im[$\epsilon(\omega)$]")
        ax2.set_xlim(tuple(args.energy_range))

        ylim_imag = calculate_y_limits(energies, [eps_imag, eps_imag_KK], tuple(args.energy_range))
        if ylim_imag:
            ax2.set_ylim(ylim_imag)

        ax2.legend()
        ax2.grid(True, linestyle=':', alpha=0.6)

    fig_kk.tight_layout()
    fig_kk.savefig("kk_check.png", dpi=300, bbox_inches="tight")
    plt.close(fig_kk)

    print(">>> Success: All data files and multi-component plots saved successfully.")


if __name__ == "__main__":
    main()
