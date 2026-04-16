"""
Plot absorption spectrum from ABACUS RT-TDDFT output.
Author: Taoni Bao
Contact: baotaoni@pku.edu.cn

Supports:
- Multi-directional dipole & E-field files
- Automatic padding & FFT
- Solid vs molecule mode (via config)
- Multiple output plots and data files

Usage:
  python plot_absorption.py --material_name "Si" --step_end 5000 --efield_path Ex.txt Ey.txt --direc 0 1 --system_type "dipole_sigma"
"""

import argparse
from typing import List, Tuple, Optional
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

# LaTeX rendering
plt.rcParams.update({
    "text.usetex": False,
    "legend.fontsize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
})

# Conversion: 1/fs → eV
FREQ2EV = sc.h / sc.eV * 1e15
# Conversion: V/Å → a.u. of E-field, 1 a.u. = 51.422 V/Å
V_ANGSTROM_TO_AU_EFIELD = (sc.e / (4 * sc.pi * sc.epsilon_0 * sc.physical_constants['Bohr radius'][0]**2) * 1e-10)
# Conversion: fs → a.u. of time, 1 fs = 41.34137 a.u.
AU_TIME = sc.physical_constants['atomic unit of time'][0]
FS_TO_AU_TIME = 1e-15 / AU_TIME

# ======================
# Utility Functions
# ======================

def validate_inputs(direc: List[int], efield_path: List[str]) -> None:
    if len(direc) != len(efield_path):
        raise ValueError(f"Number of --direc ({len(direc)}) != --efield_path ({len(efield_path)})")
    if not all(d in (0, 1, 2) for d in direc):
        raise ValueError("--direc values must be 0 (x), 1 (y), or 2 (z)")

def group_files_by_direction(direc: List[int], paths: List[str]) -> List[List[str]]:
    """Return [[x_files], [y_files], [z_files]]"""
    groups = [[], [], []]
    for d, p in zip(direc, paths):
        groups[d].append(p)
    return groups

def zero_pad_and_mask(signal: np.ndarray, dt: float, pad_factor: int = 10, decay_beta: float = 5.0) -> np.ndarray:
    """Apply exponential decay mask and zero-padding for FFT."""
    n = len(signal)
    t = np.arange(n) * dt
    mask = np.exp(-decay_beta * t / t[-1]) if t[-1] != 0 else np.ones_like(t)
    padded = np.zeros(n * pad_factor)
    padded[:n] = signal * mask
    return padded

def energy_from_fft(N: int, dt: float) -> np.ndarray:
    """Return positive-frequency energy array in eV."""
    freqs = np.fft.fftfreq(N, d=dt)[:N//2]  # in 1/fs
    return freqs * FREQ2EV

# ======================
# Data Loaders
# ======================

def load_dipole(filename: str, step_start: int, step_end: int, remove_dc: bool) -> np.ndarray:
    """Load current/dipole. Optionally remove DC component."""
    data = np.loadtxt(filename)
    end = min(step_end or len(data), data.shape[0])
    sliced = data[step_start:end, 1:].T  # shape: (3, N)

    if remove_dc:
        # Remove DC (zero-frequency) component by subtracting mean
        for i in range(3):
            sliced[i] -= np.mean(sliced[i])
    else:
        # Subtract initial value
        for i in range(3):
            sliced[i] -= sliced[i, 0]
    return sliced

def load_efield(file_groups: List[List[str]], step_start: int, step_end: int, remove_dc: bool) -> np.ndarray:
    """Load and sum E-field files per direction. Convert V/Å → a.u. Optionally remove DC component."""
    n_steps = step_end - step_start
    efield = np.zeros((3, n_steps))
    for i, files in enumerate(file_groups):
        for f in files:
            arr = np.loadtxt(f)
            efield[i] += arr[step_start:step_end, 1] / V_ANGSTROM_TO_AU_EFIELD  # V/Å → a.u.
        # Remove DC (zero-frequency) component by subtracting mean
        if remove_dc and n_steps > 0:
            efield[i] -= np.mean(efield[i])
    return efield

# ======================
# Core Calculator
# ======================

class AbsorptionSpectrum:
    def __init__(
        self,
        dipole: np.ndarray,
        efield: np.ndarray,
        dt: float,
        system_type: str,
        pad_factor: int = 10,
        decay_beta: float = 5.0,
    ):
        self.dipole = dipole
        self.efield = efield
        self.dt = dt
        self.dt_au = self.dt * FS_TO_AU_TIME
        self.system_type = system_type
        self.pad_factor = pad_factor
        self.decay_beta = decay_beta
        self.N_orig = dipole.shape[1]
        self.N_pad = self.N_orig * pad_factor

    def compute_alpha(self, direction: int) -> np.ndarray:
        """Compute alpha(ω) for a given direction (0,1,2). Returns full FFT array."""
        d_pad = zero_pad_and_mask(self.dipole[direction], self.dt_au, self.pad_factor, self.decay_beta)
        e_pad = zero_pad_and_mask(self.efield[direction], self.dt_au, self.pad_factor, self.decay_beta)

        d_fft = np.fft.fft(d_pad)
        e_fft = np.fft.fft(e_pad)

        ratio = np.divide(d_fft, e_fft, out=np.zeros_like(d_fft), where=e_fft != 0)

        if self.system_type == "dipole_epsilon":
            return np.abs(ratio.imag)
        elif self.system_type == "dipole_sigma":
            # Full frequency array (length = N_pad)
            freqs_full = np.fft.fftfreq(self.N_pad, d=self.dt_au)  # f in a.u.^{-1}
            omega = 2 * np.pi * freqs_full                         # ω in a.u.^{-1}
            alpha = np.multiply(ratio.imag, omega)
            return np.abs(alpha)
        elif self.system_type == "current_epsilon":
            # Full frequency array (length = N_pad)
            freqs_full = np.fft.fftfreq(self.N_pad, d=self.dt_au)  # f in a.u.^{-1}
            omega = 2 * np.pi * freqs_full                         # ω in a.u.^{-1}
            alpha = np.divide(ratio.real, omega, out=np.zeros_like(ratio.real), where=omega != 0)
            return np.abs(alpha)
        elif self.system_type == "current_sigma":
            return np.abs(ratio.real)
        else:
            raise ValueError("system_type must be 'dipole_epsilon', 'dipole_sigma', 'current_epsilon', or 'current_sigma'")

    def get_positive_alpha(self, direction: int) -> np.ndarray:
        """Return alpha for positive frequencies only."""
        return self.compute_alpha(direction)[:self.N_pad // 2]

    def get_energy_axis(self) -> np.ndarray:
        return energy_from_fft(self.N_pad, self.dt)

# ======================
# Plotting Utilities
# ======================

def plot_time_series(ax, time: np.ndarray, data: np.ndarray, directions: List[int], labels: List[str], ylabel: str):
    for d in directions:
        ax.plot(time, data[d], label=labels[d])
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel(ylabel)
    ax.legend()

def plot_absorption(ax, x_vals: np.ndarray, alphas: List[np.ndarray], directions: List[int], xlabel: str, x_range: Optional[Tuple] = None):
    labels = {0: "$x$", 1: "$y$", 2: "$z$"}
    for d, alpha in zip(directions, alphas):
        ax.plot(x_vals, alpha, label=labels[d])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Intensity")
    ax.legend()
    if x_range:
        ax.set_xlim(x_range)
        mask = (x_vals >= min(x_range)) & (x_vals <= max(x_range))
        if np.any(mask):
            ax.set_ylim(0, 1.2 * max(max(a[mask]) for a in alphas if len(a) == len(x_vals)))

def plot_fft(ax, energies: np.ndarray, signals: List[np.ndarray], directions: List[int], data_type: str):
    """Plot both real and imaginary parts of FFT for given signals."""
    labels = {0: "$x$", 1: "$y$", 2: "$z$"}
    for d, sig in zip(directions, signals):
        sig_pos = sig[:len(energies)]  # positive frequency part
        ax.plot(energies, sig_pos.real, '-', label=f'Re[{labels[d]}]')
        ax.plot(energies, sig_pos.imag, '--', label=f'Im[{labels[d]}]')
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("FFT Amplitude")
    ax.set_xlim(0, 20)
    ax.legend()

# ======================
# Kramers-Kronig Relations
# ======================

def compute_response_function(dipole_fft: np.ndarray, efield_fft: np.ndarray) -> np.ndarray:
    """Compute complex response χ(ω) = D(ω) / E(ω)"""
    return np.divide(dipole_fft, efield_fft, out=np.zeros_like(dipole_fft), where=efield_fft != 0)

def auto_ylim(data: np.ndarray, x: np.ndarray, x_min: float = 0.5, x_max: float = 20.0, percentile: float = 99.5) -> Tuple[float, float]:
    """
    Automatically determine y-limits by:
      1. Restricting to energy range [x_min, x_max]
      2. Using percentile to ignore outliers
    """
    mask = (x >= x_min) & (x <= x_max)
    if not np.any(mask):
        return (np.nanmin(data), np.nanmax(data))
    y_valid = data[mask]
    y_valid = y_valid[np.isfinite(y_valid)]
    if len(y_valid) == 0:
        return (0, 1)

    # Remove extreme outliers via percentile
    y_low = np.percentile(y_valid, 100 - percentile)
    y_high = np.percentile(y_valid, percentile)
    margin = 0.1 * (y_high - y_low) if y_high != y_low else 0.1 * max(abs(y_high), 1)
    return y_low - margin, y_high + margin

def kk_hilbert_transform(signal_pos: np.ndarray, omega: np.ndarray, is_odd: bool = True) -> np.ndarray:
    """
    Compute the Hilbert transform of a function defined on positive frequencies.

    Parameters:
        signal_pos: values of the function on omega > 0 (array of length N)
        omega: positive frequency grid (must be uniformly spaced, starting near 0)
        is_odd: if True, extend as odd function (for χ''); if False, extend as even (for χ')

    Returns:
        ht_pos: Hilbert transform on the same positive omega grid
    """
    from scipy.interpolate import interp1d
    from scipy.signal import hilbert as hilbert_scipy

    N = len(omega)
    if N < 2:
        return np.zeros_like(omega)

    domega = omega[1] - omega[0]
    omega_max = omega[-1]

    # Full frequency grid (symmetric around 0)
    omega_full = np.linspace(-omega_max, omega_max, 2 * N - 1)

    # Interpolator for positive side
    interp = interp1d(omega, signal_pos, kind='linear', bounds_error=False, fill_value=0.0)

    # Extend to negative frequencies
    signal_full = np.zeros_like(omega_full)
    pos_mask = omega_full >= 0
    signal_full[pos_mask] = interp(omega_full[pos_mask])
    if is_odd:
        signal_full[~pos_mask] = -interp(-omega_full[~pos_mask])  # odd extension
    else:
        signal_full[~pos_mask] = interp(-omega_full[~pos_mask])   # even extension

    # Compute analytic signal → Hilbert transform is imaginary part
    analytic = hilbert_scipy(signal_full)
    ht_full = analytic.imag  # this is H[signal_full]

    # Extract positive frequencies (center at index N-1)
    ht_pos = ht_full[N - 1:]

    return ht_pos

# ======================
# Main Entry Point
# ======================

def main():
    parser = argparse.ArgumentParser(description="Plot RT-TDDFT absorption spectrum from ABACUS output.")
    parser.add_argument("--td_dt", type=float, default=0.00484, help="Time step in fs")
    parser.add_argument("--direc", type=int, nargs="+", required=True, help="Field directions: 0=x,1=y,2=z")
    parser.add_argument("--efield_path", type=str, nargs="+", required=True, help="E-field file paths")
    parser.add_argument("--dipolefile", type=str, default="./OUT.ABACUS/dipole_s1.txt")
    # parser.add_argument("--dipolefile", type=str, default="./OUT.ABACUS/dipole_s1.txt")
    parser.add_argument("--material_name", type=str,
                        default=r"CH$_3$COOH", help="Name for plot title")
    parser.add_argument("--step_start", type=int, default=0)
    parser.add_argument("--step_end", type=int, default=8000)
    parser.add_argument("--system_type", choices=["dipole_epsilon", "dipole_sigma", "current_epsilon", "current_sigma"],
                        default="dipole_sigma", help="Type of system and response for absorption calculation")
    parser.add_argument("--pad_factor", type=int, default=10, help="Zero-padding factor for FFT")
    parser.add_argument("--decay_beta", type=float, default=5, help="Exponential decay parameter for windowing")
    parser.add_argument("--remove_dc", type=bool, default=False,
                        help="Whether to remove DC component from time signals")

    args = parser.parse_args()

    validate_inputs(args.direc, args.efield_path)
    Efile_groups = group_files_by_direction(args.direc, args.efield_path)

    # Load data
    dipole_data = load_dipole(args.dipolefile, args.step_start, args.step_end, args.remove_dc)
    efield_data = load_efield(Efile_groups, args.step_start, args.step_end, args.remove_dc)

    # Compute absorption
    absorber = AbsorptionSpectrum(
        dipole=dipole_data,
        efield=efield_data,
        dt=args.td_dt,
        system_type=args.system_type,
        pad_factor=args.pad_factor,
        decay_beta=args.decay_beta,
    )

    energies = absorber.get_energy_axis()

    with np.errstate(divide='ignore'):
        wavelengths = np.where(energies > 0, (sc.h * sc.c / sc.eV * 1e9) / energies, np.inf)

    directions = sorted(set(args.direc))  # unique directions to plot

    # --- Save data ---
    alphas = [absorber.get_positive_alpha(d) for d in range(3)]
    data_out = np.column_stack([energies, wavelengths] + alphas)
    np.savetxt("abs_dat.txt", data_out, header="Energy(eV) Wavelength(nm) alpha_X alpha_Y alpha_Z")

    # --- Plotting ---
    time = np.arange(args.step_end - args.step_start) * args.td_dt

    # 1. Dipole vs time
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_time_series(ax, time, dipole_data, directions, ["$J_x$", "$J_y$", "$J_z$"], "Current (a.u.)")
    ax.set_xlim(0, args.step_end * args.td_dt)
    fig.savefig("dipole.png", dpi=300, bbox_inches="tight")

    # 2. Absorption spectrum
    fig, ax = plt.subplots(figsize=(10, 6))
    alphas_to_plot = [absorber.get_positive_alpha(d) for d in directions]
    xmax = 20
    # plot_absorption(ax, energies, alphas_to_plot, directions, x_range=(5, xmax))
    plot_absorption(ax, energies, alphas_to_plot, directions, xlabel="Energy (eV)", x_range=(4, 20))
    # Automatically set y-limits
    all_alpha = np.concatenate([a for a in alphas_to_plot if len(a) == len(energies)])
    if len(all_alpha) > 0:
        ymin, ymax = auto_ylim(all_alpha, energies, x_min=0.5, x_max=xmax)
        ax.set_ylim(0, ymax)
    ax.set_ylim(0, 100)
    title = rf"{args.material_name} Absorption Spectrum"
    ax.set_title(title, fontsize=20, y=1.02)
    fig.savefig("abs.png", dpi=300, bbox_inches="tight")

    fig, ax = plt.subplots(figsize=(10, 6))
    alphas_to_plot = [absorber.get_positive_alpha(d) for d in directions]

    # Set wavelength range for plotting
    wl_min, wl_max = 80, 200
    plot_absorption(ax, wavelengths, alphas_to_plot, directions, xlabel="Wavelength (nm)", x_range=(wl_min, wl_max))

    all_alpha = np.concatenate([a for a in alphas_to_plot if len(a) == len(wavelengths)])
    if len(all_alpha) > 0:
        ymin, ymax = auto_ylim(all_alpha, wavelengths, x_min=wl_min, x_max=wl_max)
        ax.set_ylim(0, ymax)

    title = rf"{args.material_name} Absorption Spectrum"
    ax.set_title(title, fontsize=20, y=1.02)
    fig.savefig("abs_wavelength.png", dpi=300, bbox_inches="tight")

    # 3. E-field Fourier
    fig, ax = plt.subplots(figsize=(10, 6))
    efield_ffts = [np.fft.fft(zero_pad_and_mask(efield_data[d], args.td_dt, args.pad_factor)) for d in directions]
    plot_fft(ax, energies, efield_ffts, directions, "efield")
    fig.savefig("efield_fourier.png", dpi=300, bbox_inches="tight")

    # 4. Dipole Fourier
    fig, ax = plt.subplots(figsize=(10, 6))
    dipole_ffts = [np.fft.fft(zero_pad_and_mask(dipole_data[d], args.td_dt, args.pad_factor)) for d in directions]
    plot_fft(ax, energies, dipole_ffts, directions, "dipole")
    fig.savefig("dipole_fourier.png", dpi=300, bbox_inches="tight")

    # 5. E-field vs time
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_time_series(ax, time, efield_data, directions, ["$E_x$", "$E_y$", "$E_z$"], "Electric Field (a.u.)")
    ax.set_xlim(0, 1)
    fig.savefig("efield_time.png", dpi=300, bbox_inches="tight")

    # 6. Kramers-Kronig Validation
    n_dirs = len(directions)
    fig_kk, axes_kk = plt.subplots(n_dirs, 2, figsize=(14, 5 * n_dirs), squeeze=False)

    for idx, d in enumerate(directions):
        # --- Get complex response ---
        d_pad = zero_pad_and_mask(dipole_data[d], args.td_dt, args.pad_factor, args.decay_beta)
        e_pad = zero_pad_and_mask(efield_data[d], args.td_dt, args.pad_factor, args.decay_beta)
        d_fft = np.fft.fft(d_pad)
        e_fft = np.fft.fft(e_pad)
        chi = compute_response_function(d_fft, e_fft)
        chi_pos = chi[:len(energies)]

        chi_real = chi_pos.real
        chi_imag = chi_pos.imag
        omega = energies  # in eV

        # --- Reconstruct real part from imag ---
        chi_real_KK = kk_hilbert_transform(chi_imag, omega, is_odd=True)  # χ'' is odd → use odd extension

        # --- Reconstruct imag part from real ---
        chi_imag_KK = -kk_hilbert_transform(chi_real, omega, is_odd=False)  # χ' is even → use even extension

        # --- Plot: Real part comparison ---
        ax1 = axes_kk[idx, 0]
        ax1.plot(omega, chi_real, 'b-', lw=1.5, label=r"$\chi'(\omega)$ (FFT)")
        ax1.plot(omega, chi_real_KK, 'r--', lw=1.5, label=r"$\mathcal{H}[\chi''](\omega)$")
        ax1.set_xlabel("Energy (eV)")
        ax1.set_ylabel(r"$\chi'(\omega)$")
        ax1.set_title(f"Direction {['$x$','$y$','$z$'][d]}: Real part KK check", fontsize=18, y=1.02)
        ax1.legend()
        ax1.grid(True, linestyle=':', alpha=0.6)
        ax1.set_xlim(0, min(20, omega[-1]))
        # Automatically set y-limits
        ymin, ymax = auto_ylim(chi_real, omega, x_min=0.5, x_max=20.0)
        ymin_kk, ymax_kk = auto_ylim(chi_real_KK, omega, x_min=0.5, x_max=20.0)
        ax1.set_ylim(-max(abs(ymin), abs(ymin_kk)), max(ymax, ymax_kk))

        # --- Plot: Imaginary part comparison ---
        ax2 = axes_kk[idx, 1]
        ax2.plot(omega, chi_imag, 'g-', lw=1.5, label=r"$\chi''(\omega)$ (FFT)")
        ax2.plot(omega, chi_imag_KK, 'm--', lw=1.5, label=r"$-\mathcal{H}[\chi'](\omega)$")
        ax2.set_xlabel("Energy (eV)")
        ax2.set_ylabel(r"$\chi''(\omega)$")
        ax2.set_title(f"Direction {['$x$','$y$','$z$'][d]}: Imaginary part KK check", fontsize=18, y=1.02)
        ax2.legend()
        ax2.grid(True, linestyle=':', alpha=0.6)
        ax2.set_xlim(0, min(20, omega[-1]))
        # Automatically set y-limits
        ymin2, ymax2 = auto_ylim(chi_imag, omega, x_min=0.5, x_max=20.0)
        ymin2_kk, ymax2_kk = auto_ylim(chi_imag_KK, omega, x_min=0.5, x_max=20.0)
        ax2.set_ylim(-max(abs(ymin2), abs(ymin2_kk)), max(ymax2, ymax2_kk))

    fig_kk.tight_layout()
    fig_kk.savefig("kk_check.png", dpi=300, bbox_inches="tight")

    print("Plots and data saved successfully.")

if __name__ == "__main__":
    main()
