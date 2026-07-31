#!/usr/bin/env python3
"""Plot an ABACUS exciton-density slice.

The current text format is deliberately compact and self-describing::

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

Files written before format version 1 are still accepted.

Usage:
    python plot_cond_silce.py FILE.dat [-o FILE.png]
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


def _parse_legacy(comments):
    raw = {}
    for line in comments:
        fields = line.split(None, 1)
        if len(fields) == 2:
            raw[fields[0]] = fields[1]

    integer_keys = (
        "state", "grid_nu", "grid_nv", "BvK_u", "BvK_v",
        "u_min", "u_max", "v_min", "v_max", "n_atoms", "has_hole",
    )
    for key in integer_keys:
        if key in raw:
            raw[key] = int(raw[key])

    float_keys = (
        "energy_Ry", "slice_pos", "hole_fix_x", "hole_fix_y", "hole_fix_z",
        "u_vec_x", "u_vec_y", "u_vec_z", "v_vec_x", "v_vec_y", "v_vec_z",
        "cell_a_x", "cell_a_y", "cell_a_z", "cell_b_x", "cell_b_y", "cell_b_z",
        "cell_c_x", "cell_c_y", "cell_c_z",
    )
    for key in float_keys:
        if key in raw:
            raw[key] = float(raw[key])

    meta = {
        "format_version": 0,
        "state": raw["state"],
        "energy_Ry": raw["energy_Ry"],
        "density_kind": raw.get("density_kind", "conditional_elec"),
        "fixed_particle": raw.get("fixed_particle", "hole"),
        "plane": raw.get("plane", "ab"),
        "slice_pos": raw.get("slice_pos", 0.0),
        "nu": raw["grid_nu"],
        "nv": raw["grid_nv"],
        "nk_u": raw.get("BvK_u", 1),
        "nk_v": raw.get("BvK_v", 1),
        "u_min": raw["u_min"],
        "u_max": raw["u_max"],
        "v_min": raw["v_min"],
        "v_max": raw["v_max"],
        "u_vec": _vector([raw["u_vec_x"], raw["u_vec_y"], raw["u_vec_z"]]),
        "v_vec": _vector([raw["v_vec_x"], raw["v_vec_y"], raw["v_vec_z"]]),
        "a": _vector([raw["cell_a_x"], raw["cell_a_y"], raw["cell_a_z"]]),
        "b": _vector([raw["cell_b_x"], raw["cell_b_y"], raw["cell_b_z"]]),
        "c": _vector([raw["cell_c_x"], raw["cell_c_y"], raw["cell_c_z"]]),
        "atoms": [],
    }
    has_fixed = raw.get(
        "has_hole", 1 if meta["density_kind"].startswith("conditional") else 0)
    if has_fixed:
        meta["hole"] = _vector([
            raw["hole_fix_x"], raw["hole_fix_y"], raw["hole_fix_z"]])
    for index in range(raw["n_atoms"]):
        meta["atoms"].append({
            "label": raw["atom_label_{}".format(index)],
            "pos": _vector([
                raw["atom_x_{}".format(index)],
                raw["atom_y_{}".format(index)],
                raw["atom_z_{}".format(index)],
            ]),
        })
    meta["n_atoms"] = len(meta["atoms"])
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

    meta = (_parse_v1(comments) if comments[0].startswith(_MAGIC)
            else _parse_legacy(comments))
    data = np.loadtxt(filename, comments="#", dtype=float)
    data = np.atleast_2d(data)
    _validate(meta, data, filename)
    return meta, data


def parse(filename):
    """Compatibility wrapper returning only metadata."""
    return read_slice(filename)[0]


def load_data(filename):
    """Compatibility wrapper returning only the density matrix."""
    return read_slice(filename)[1]


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
        x_edges, y_edges, data.T, cmap="hot",
        norm=LogNorm(vmin=minimum, vmax=maximum),
        shading="flat", rasterized=True)
    color_label, title_kind = _density_labels(meta["density_kind"])
    fig.colorbar(image, ax=axis, label=color_label, shrink=0.78)

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
        "{} | State {} | {}-plane\n{}E = {:.6f} Ry = {:.3f} eV\n"
        "BvK {}x{} | grid {}x{}".format(
            title_kind, meta["state"], meta["plane"], fixed_text,
            meta["energy_Ry"], meta["energy_Ry"] * 13.605693,
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
