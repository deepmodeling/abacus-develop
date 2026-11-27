#!/usr/bin/env python3
"""Plot an ABACUS exciton-density slice.

The data text format is deliberately compact and self-describing:
  For example:
    # ABACUS_EXCITON_SLICE 1
    # state 0 energy_Ry 0.2265 density_kind conditional_elec
    # fixed_particle hole position_bohr 2.5 2.5 2.5
    # plane ca slice_pos_bohr 0
    # grid 199 199 u_range -5 6 v_range -5 6
    # bvk 7 7
    # lattice_bohr
    # a 0 5.15 5.15
    # b 5.15 0 5.15
    # c 5.15 5.15 0
    # axes_bohr
    # u 5.15 5.15 0
    # v 0 5.15 5.15
    # atoms 2
    # Si 0 0 0
    # Si 2.575 2.575 2.575
    # data
    <nu rows, each containing nv density values>

Usage:
    python plot_exciton_slice.py FILE.dat [-o FILE.png]
"""

import argparse
from pathlib import Path

import numpy as np


_MAGIC = "ABACUS_EXCITON_SLICE"


def _vector(values):
    return np.asarray([float(value) for value in values], dtype=float)


def _parse_v1(comments):
    if comments[0].split() != [_MAGIC, "1"]:
        raise ValueError("Unsupported exciton slice format: " + comments[0])

    meta = {"format_version": 1, "atoms": []}
    atom_lines_remaining = 0
    for line in comments[1:]:
        fields = line.split()
        if not fields:
            continue
        if atom_lines_remaining:
            if len(fields) != 4:
                raise ValueError("Invalid atom record: " + line)
            meta["atoms"].append({"label": fields[0], "pos": _vector(fields[1:])})
            atom_lines_remaining -= 1
            continue

        key = fields[0]
        if key == "state":
            if len(fields) != 6 or fields[2] != "energy_Ry" or fields[4] != "density_kind":
                raise ValueError("Invalid state record: " + line)
            meta["state"] = int(fields[1])
            meta["energy_Ry"] = float(fields[3])
            meta["density_kind"] = fields[5]
        elif key == "fixed_particle":
            meta["fixed_particle"] = fields[1]
            if len(fields) == 6 and fields[2] == "position_bohr":
                meta["hole"] = _vector(fields[3:6])
            elif len(fields) != 2:
                raise ValueError("Invalid fixed-particle record: " + line)
        elif key == "plane":
            if len(fields) != 4 or fields[2] != "slice_pos_bohr":
                raise ValueError("Invalid plane record: " + line)
            meta["plane"] = fields[1]
            meta["slice_pos"] = float(fields[3])
        elif key == "grid":
            if (len(fields) != 9 or fields[3] != "u_range"
                    or fields[6] != "v_range"):
                raise ValueError("Invalid grid record: " + line)
            meta["nu"], meta["nv"] = int(fields[1]), int(fields[2])
            meta["u_min"], meta["u_max"] = float(fields[4]), float(fields[5])
            meta["v_min"], meta["v_max"] = float(fields[7]), float(fields[8])
        elif key == "bvk":
            if len(fields) != 3:
                raise ValueError("Invalid BvK record: " + line)
            meta["nk_u"], meta["nk_v"] = int(fields[1]), int(fields[2])
        elif key in ("a", "b", "c", "u", "v"):
            if len(fields) != 4:
                raise ValueError("Invalid vector record: " + line)
            meta[key if key in ("a", "b", "c") else key + "_vec"] = _vector(fields[1:])
        elif key == "atoms":
            if len(fields) != 2:
                raise ValueError("Invalid atoms record: " + line)
            meta["n_atoms"] = int(fields[1])
            atom_lines_remaining = meta["n_atoms"]
        elif key not in ("lattice_bohr", "axes_bohr", "data"):
            raise ValueError("Unknown format-v1 header record: " + line)

    if atom_lines_remaining:
        raise ValueError("Slice header ended before all atom records were read")
    return meta

def _validate(meta, data, filename):
    required = (
        "state", "energy_Ry", "density_kind", "fixed_particle", "plane",
        "slice_pos", "nu", "nv", "nk_u", "nk_v", "u_min", "u_max",
        "v_min", "v_max", "u_vec", "v_vec", "a", "b", "c", "atoms",
    )
    missing = [key for key in required if key not in meta]
    if missing:
        raise ValueError("{} is missing metadata: {}".format(
            filename, ", ".join(missing)))
    if meta["plane"] not in ("ab", "bc", "ca"):
        raise ValueError("Unsupported slice plane: " + meta["plane"])
    if data.shape != (meta["nu"], meta["nv"]):
        raise ValueError(
            "{} declares a {}x{} grid but contains {}".format(
                filename, meta["nu"], meta["nv"], data.shape))
    if not np.all(np.isfinite(data)):
        raise ValueError("{} contains non-finite density values".format(filename))
    if np.any(data < 0.0):
        raise ValueError("{} contains negative density values".format(filename))
    if np.linalg.norm(meta["u_vec"]) == 0 or np.linalg.norm(meta["v_vec"]) == 0:
        raise ValueError("{} contains a zero slice vector".format(filename))


def read_slice(filename):
    """Read metadata and the density matrix from an old or format-v1 file."""
    comments = []
    with open(filename, encoding="utf-8") as stream:
        for raw_line in stream:
            stripped = raw_line.strip()
            if stripped.startswith("#"):
                comments.append(stripped[1:].strip())
    if not comments:
        raise ValueError("{} has no exciton slice header".format(filename))

    if comments[0].startswith(_MAGIC):
        meta = _parse_v1(comments)
    else:
        raise ValueError("{} has no format-v1 header".format(filename))
    data = np.loadtxt(filename, comments="#", dtype=float)
    data = np.atleast_2d(data)
    _validate(meta, data, filename)
    return meta, data

def _density_labels(kind):
    labels = {
        "conditional_elec": (r"$|\psi_{\rm cond}|^2$", "Conditional Electron Density"),
        "conditional_hole": (r"$|\psi_{\rm cond}|^2$", "Conditional Hole Density"),
        "average_hole": (r"$\rho_{\rm h}^{\rm avg}$", "Average Hole Density"),
        "average_elec": (r"$\rho_{\rm e}^{\rm avg}$", "Average Electron Density"),
    }
    return labels.get(kind, (kind, kind.replace("_", " ").title()))


def plot(data, meta, output):
    """Render the density in Cartesian coordinates for a possibly skew plane."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.lines import Line2D
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    u_vec, v_vec = meta["u_vec"], meta["v_vec"]
    u_length = np.linalg.norm(u_vec)
    v_length = np.linalg.norm(v_vec)
    u_hat = u_vec / u_length
    v_perp = v_vec - np.dot(v_vec, u_hat) * u_hat
    v_perp_length = np.linalg.norm(v_perp)
    if v_perp_length == 0:
        raise ValueError("Slice vectors are collinear")
    v_perp_hat = v_perp / v_perp_length
    v_along_u = np.dot(v_vec, u_hat)
    angle = np.degrees(np.arccos(np.clip(
        np.dot(u_vec, v_vec) / (u_length * v_length), -1.0, 1.0)))

    u_edges = np.linspace(meta["u_min"], meta["u_max"], meta["nu"])
    v_edges = np.linspace(meta["v_min"], meta["v_max"], meta["nv"])
    du = u_edges[1] - u_edges[0]
    dv = v_edges[1] - v_edges[0]
    u_edges = np.r_[u_edges - 0.5 * du, u_edges[-1] + 0.5 * du]
    v_edges = np.r_[v_edges - 0.5 * dv, v_edges[-1] + 0.5 * dv]
    x_edges = u_edges[:, None] * u_length + v_edges[None, :] * v_along_u
    y_edges = np.broadcast_to(
        v_edges[None, :] * v_perp_length, x_edges.shape)

    maximum = data.max()
    if maximum <= 0:
        raise ValueError("Density data contains no positive values")
    positive = data[data > 0]
    minimum = max(positive.min(), maximum * 1.0e-6)

    fig, axis = plt.subplots(figsize=(12, 10))
    image = axis.pcolormesh(
        x_edges, y_edges, data, cmap="hot",
        norm=LogNorm(vmin=minimum, vmax=maximum),
        shading="flat", rasterized=True)
    color_label, title_kind = _density_labels(meta["density_kind"])
    divider = make_axes_locatable(axis)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    fig.colorbar(image, cax=cax, label=color_label)

    def project(position):
        return np.dot(position, u_hat), np.dot(position, v_perp_hat)

    for index in range(int(np.floor(meta["u_min"])),
                       int(np.ceil(meta["u_max"])) + 1):
        start = index * u_vec + meta["v_min"] * v_vec
        end = index * u_vec + meta["v_max"] * v_vec
        axis.plot(*zip(project(start), project(end)),
                  color="black", ls="--", lw=0.6, alpha=0.55)
    for index in range(int(np.floor(meta["v_min"])),
                       int(np.ceil(meta["v_max"])) + 1):
        start = meta["u_min"] * u_vec + index * v_vec
        end = meta["u_max"] * u_vec + index * v_vec
        axis.plot(*zip(project(start), project(end)),
                  color="black", ls="--", lw=0.6, alpha=0.55)

    bvk_u0 = -(meta["nk_u"] // 2)
    bvk_v0 = -(meta["nk_v"] // 2)
    for index in (bvk_u0, bvk_u0 + meta["nk_u"]):
        axis.plot(*zip(
            project(index * u_vec + meta["v_min"] * v_vec),
            project(index * u_vec + meta["v_max"] * v_vec)),
            color="black", lw=2.8, alpha=0.94)
    for index in (bvk_v0, bvk_v0 + meta["nk_v"]):
        axis.plot(*zip(
            project(meta["u_min"] * u_vec + index * v_vec),
            project(meta["u_max"] * u_vec + index * v_vec)),
            color="black", lw=2.8, alpha=0.94)

    perpendicular = (
        meta["c"] if meta["plane"] == "ab"
        else meta["a"] if meta["plane"] == "bc" else meta["b"])
    perpendicular_length = np.linalg.norm(perpendicular)
    perpendicular_hat = perpendicular / perpendicular_length
    labels = sorted({atom["label"] for atom in meta["atoms"]})
    colors = plt.cm.tab10(np.linspace(0, 1, max(1, len(labels))))
    label_colors = dict(zip(labels, colors))
    for atom in meta["atoms"]:
        home_x, home_y = project(atom["pos"])
        depth = np.dot(atom["pos"], perpendicular_hat) - meta["slice_pos"]
        depth = ((depth + perpendicular_length / 2) % perpendicular_length
                 - perpendicular_length / 2)
        alpha = max(0.2, 1.0 - abs(depth) / (perpendicular_length * 0.5))
        for cell_u in range(int(np.floor(meta["u_min"])),
                            int(np.ceil(meta["u_max"]))):
            for cell_v in range(int(np.floor(meta["v_min"])),
                                int(np.ceil(meta["v_max"]))):
                axis.plot(
                    home_x + cell_u * u_length + cell_v * v_along_u,
                    home_y + cell_v * v_perp_length,
                    "o", color=label_colors[atom["label"]], ms=6,
                    alpha=alpha, mec="white", mew=0.5)

    handles = []
    if "hole" in meta:
        fixed_x, fixed_y = project(meta["hole"])
        axis.plot(
            fixed_x, fixed_y, "o", color="lime", ms=10,
            mec="black", mew=1.5, zorder=10)
        handles.append(Line2D(
            [0], [0], marker="o", color="white", markerfacecolor="lime",
            markersize=8, markeredgecolor="black", markeredgewidth=1.5,
            label="Fixed {} (home cell)".format(meta["fixed_particle"])))
    for label in labels:
        handles.append(Line2D(
            [0], [0], marker="o", color="white",
            markerfacecolor=label_colors[label], markersize=7,
            markeredgecolor="white", markeredgewidth=0.5, label=label))
    if handles:
        axis.legend(loc="upper right", fontsize=9, handles=handles)

    axis.set_xlabel("Distance along u (|u| = {:.2f} Bohr)".format(u_length))
    axis.set_ylabel(
        "Distance perpendicular to u (|v| = {:.2f} Bohr, angle = {:.0f} deg)"
        .format(v_length, angle))
    axis.set_aspect("equal")
    fixed_text = (
        "{} at ({:.2f}, {:.2f}, {:.2f}) Bohr | ".format(
            meta["fixed_particle"].title(), *meta["hole"])
        if "hole" in meta else "")
    axis.set_title(
        "{} | State {} | {}-plane\n"
        "{}BvK {}x{} | grid {}x{}".format(
            title_kind, meta["state"], meta["plane"], fixed_text,
            meta["nk_u"], meta["nk_v"], meta["nu"], meta["nv"]))
    fig.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_file", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    args = parser.parse_args()

    meta, data = read_slice(args.data_file)
    output = args.output or args.data_file.with_suffix(".png")
    plot(data, meta, output)
    print("{} | format v{} | grid {} | BvK {}x{} | saved {}".format(
        meta["density_kind"], meta["format_version"], data.shape,
        meta["nk_u"], meta["nk_v"], output))


if __name__ == "__main__":
    main()
