"""scf 多帧累积下 read_results 必须返回当前结构帧的力（spec R1/R2，legacyio 主路径）。

Task 3 GREEN：`AbacusTemplate.read_results` 对 calculation=='scf' 按坐标选择"当前结构帧"
（core.py::_select_scf_frame_for_structure），非 scf（relax/md）保持 [-1] 末帧语义。
金样 multiframe_scf_trial_last/running_scf.log 共 4 帧：帧 0（ION=1）= INIT（== 当前 STRU），
帧 1-2 为位移/trial 试探结构，末帧（ION=4）= P0 重算（≠ STRU）。

Task 4 扩充（回归矩阵）：容差两侧、md/relax/sella/ccqn 五路径、双后端（legacyio/latestio）、
真实非 scf 金样（multiframe_md_legacy / multiframe_md_latest）、latestio eig_occ.txt 一次性消费（P5）。
"""
import re
import shutil
from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.units import Bohr

from abacuslite.core import (
    AbacusTemplate,
    _select_scf_frame_for_structure,
)
from abacuslite.io.generalio import read_stru
from abacuslite.io.legacyio import read_abacus_out
from abacuslite import core as abacus_core


TESTFILES = (
    Path(__file__).resolve().parents[1] / "fixtures"
)
GOLD = TESTFILES / "multiframe_scf_trial_last"
LTS_ROUNDED_CELL = TESTFILES / "lts_scf_cell_rounding"
LOG = GOLD / "running_scf.log"

# 帧 0（INIT，== 当前 STRU）的力；末帧（P0 重算）力不同，用于断言"返回当前结构帧而非末帧"
STRU_FRAME_FORCES = np.array(
    [[-0.5, -0.01, 0.0], [0.5, 0.01, 0.0], [0.0, 0.0, 0.0]]
)
LAST_FRAME_FORCES = np.array(
    [[-0.45, -0.014, 0.0], [0.45, 0.014, 0.0], [0.0, 0.0, 0.0]]
)

# 真实 md 金样末帧力（multiframe_md_legacy=首个 MD 步/唯一帧；multiframe_md_latest=MD_dump 第 2 帧）
MD_LEGACY_FORCES = np.array(
    [[1.0712942948, 2.1781413019, 1.2152201133],
     [-1.0712942948, -2.1781413019, -1.2152201133]]
)
MD_LATEST_FORCES = np.array(
    [[0.1191299488, 0.2678988387, 1.5033416041],
     [-0.1191299488, -0.2678988387, -1.5033416041]]
)
MD_LEGACY_ENERGY = -1940.5509086572
MD_LATEST_ENERGY = -1940.6311063727


@pytest.fixture(autouse=True)
def _pin_legacyio_backend():
    """本模块以 legacyio 为帧选择主路径（spec R1/R2），测试间钉住 __LEGACYIO__=True，
    避免其他测试模块（如 test_abacuslite_profile 的 beta 版本用例）改写该模块级全局
    导致的 CI 顺序依赖失败。latestio 专项用例自行 monkeypatch 为 False 覆盖。"""
    original = abacus_core.__LEGACYIO__
    abacus_core.__LEGACYIO__ = True
    yield
    abacus_core.__LEGACYIO__ = original

# tauc_/taud_ 原子坐标行（x y z mag vx vy vz）：前 3 个浮点为坐标
_COORD_LINE = re.compile(r"^(taud_\S+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)(.*)$")


def _stru_positions_cartesian(stru, cell):
    """read_stru 的 DIRECT 分数坐标 -> Cartesian Å（与 running log 同源晶胞换算）。

    STRU 与 running log 描述同一晶胞，直接用帧的 cell（Å）换算（ase 3.28 无 abacus
    STRU 读取器，故用 read_stru）。与 core.py 帧选择采用同一坐标系/容差口径（spec P3 ③）。
    """
    frac = np.array(
        [a["coord"] for sp in stru["species"] for a in sp["atom"]], dtype=float
    )
    return frac @ np.asarray(cell)


def _cartesian_stru_text():
    """把金样 DIRECT STRU 变换为同结构的 Cartesian STRU（write_stru 默认即 Cartesian）。

    体系与金样一致：H(0.4,0.5,0.5)/H(0.6,0.5,0.5)/Au(0.5,0.5,0.1)（Direct）⇔
    H(4,5,5)/H(6,5,5)/Au(5,5,1)（Cartesian Å，晶胞 10 Å 立方）。
    """
    return (
        (GOLD / "STRU")
        .read_text()
        .replace("Direct", "Cartesian")
        .replace("0.4000000000 0.5000000000 0.5000000000",
                 "4.0000000000 5.0000000000 5.0000000000")
        .replace("0.6000000000 0.5000000000 0.5000000000",
                 "6.0000000000 5.0000000000 5.0000000000")
        .replace("0.5000000000 0.5000000000 0.1000000000",
                 "5.0000000000 5.0000000000 1.0000000000")
    )


def _cartesian_log_text(alat=1.0):
    """Express the same 10-Angstrom cell and Cartesian positions in units of alat."""
    lines = []
    for line in (GOLD / "running_scf.log").read_text().splitlines(keepends=True):
        if re.match(r"^DIRECT COORDINATES\s*$", line.rstrip("\n")):
            lines.append("CARTESIAN COORDINATES\n")
            continue
        m = _COORD_LINE.match(line.rstrip("\n"))
        if m:
            x, y, z = (float(m.group(i)) * 10.0 / alat for i in (2, 3, 4))
            lines.append(f"{m.group(1)} {x:20.10f} {y:20.10f} {z:20.10f}{m.group(5)}\n")
            continue
        lines.append(line)
    return ("".join(lines)
            .replace("(Angstrom) = 1.0", f"(Angstrom) = {alat}")
            .replace("(Bohr) = 1.889726", f"(Bohr) = {alat / Bohr}")
            .replace("+10", f"+{10.0 / alat}"))


def _log_with_extra_atom_per_frame():
    """Corrupt the raw parser stream without changing its force tables."""
    lines = []
    for line in (GOLD / "running_scf.log").read_text().splitlines():
        if "TOTAL ATOM NUMBER = 3" in line:
            line = line.replace("TOTAL ATOM NUMBER = 3", "TOTAL ATOM NUMBER = 4")
        lines.append(line)
        if line.lstrip().startswith("taud_Au1"):
            lines.append(
                "taud_X1            0.2000000000        0.3000000000        "
                "0.4000000000  0.0000        0.0000000000        "
                "0.0000000000        0.0000000000"
            )
    return "\n".join(lines) + "\n"


def _log_with_translated_first_frame():
    """Translate each raw first-frame atom by a different integer direct vector."""
    shifts = {
        "taud_H1": np.array([1, -2, 0]),
        "taud_H2": np.array([-1, 0, 3]),
        "taud_Au1": np.array([0, 1, -1]),
    }
    lines = []
    in_first_frame = False
    translated = 0
    for line in (GOLD / "running_scf.log").read_text().splitlines():
        if line.strip() == "DIRECT COORDINATES" and not in_first_frame:
            in_first_frame = True
        match = _COORD_LINE.match(line)
        if in_first_frame and match and match.group(1) in shifts:
            coordinates = np.array(
                [float(match.group(index)) for index in (2, 3, 4)]
            ) + shifts[match.group(1)]
            line = (
                f"{match.group(1)} {coordinates[0]:20.10f} "
                f"{coordinates[1]:20.10f} {coordinates[2]:20.10f}{match.group(5)}"
            )
            translated += 1
            if translated == 3:
                in_first_frame = False
        lines.append(line)
    return "\n".join(lines) + "\n"


def _template(calc="scf", atomorder=None):
    """构造最小 AbacusTemplate：仅 read_results 依赖的字段（suffix/calculation/atomorder）。"""
    template = AbacusTemplate()
    template.suffix = "ABACUS"
    template.calculation = calc
    template.atomorder = [0, 1, 2] if atomorder is None else atomorder
    return template


def _rounded_cell_run(tmp_path):
    """Copy the attributed one-frame LTS output without modifying its geometry."""
    shutil.copytree(LTS_ROUNDED_CELL, tmp_path, dirs_exist_ok=True)
    template = _template(atomorder=[0])
    return template, tmp_path


@pytest.mark.parametrize("json_contents", [None, "{invalid JSON"])
def test_lts_rounded_log_cell_matches_without_abacus_json(tmp_path, json_contents):
    """Six-significant-digit lattice rounding stays inside the default tolerance.

    abacus.json is not consulted at all: the recomputed frame matches the STRU
    once the tolerance covers the log's own printing resolution.
    """
    template, directory = _rounded_cell_run(tmp_path)
    if json_contents is not None:
        (directory / "abacus.json").write_text(json_contents)

    results = template.read_results(directory)
    assert results["energy"] == pytest.approx(-13.5410411603)
    np.testing.assert_allclose(results["forces"], [[0, 0, 0]], atol=1e-12)


def test_recorded_lts_rounding_exceeds_old_tolerance(tmp_path):
    """The old 1e-5 bound rejects a valid saved run; retain the 1e-4 default."""
    _, directory = _rounded_cell_run(tmp_path)
    log = directory / "OUT.ABACUS" / "running_scf.log"
    frames = read_abacus_out(log)
    stru = read_stru(directory / "STRU")
    expected_cell = np.asarray(stru["lat"]["vec"]) * stru["lat"]["const"] * Bohr
    difference = np.max(np.abs(np.asarray(frames[0].cell) - expected_cell))
    assert difference == pytest.approx(1.65221e-5, abs=1e-12)
    with pytest.raises(RuntimeError):
        _select_scf_frame_for_structure(
            frames, log, directory, [0], atol=1e-5, legacyio=True
        )
    assert _select_scf_frame_for_structure(
        frames, log, directory, [0], legacyio=True
    ) is frames[0]


@pytest.mark.parametrize("scale", [1.234564321, 1.2345652])
def test_lts_rounded_cell_nonunit_alat_and_nonorthogonal_vectors(tmp_path, scale):
    """The precise Cartesian header, rather than its rounded Angstrom line, sets alat."""
    template, directory = _rounded_cell_run(tmp_path)
    # LTS prints the rounded Angstrom line with its own 0.529177 constant;
    # at the second value this differs from rounding with ASE's Bohr value.
    lts_printed_scale = (scale / Bohr) * 0.529177
    vectors = np.array([[10.0000165221, 0.25, 0],
                        [0, 11.000002299, 0.5],
                        [0.125, 0, 11.9999880759]])
    stru = directory / "STRU"
    text = stru.read_text().replace(
        "1.8897261258369282", str(scale / Bohr)
    ).replace(
        "10.000016522100001 0.0 0.0\n"
        "0.0 11.000002299 0.0\n"
        "0.0 0.0 11.999988075900001",
        "\n".join(" ".join(map(str, row)) for row in vectors),
    )
    stru.write_text(text)

    log_path = directory / "OUT.ABACUS" / "running_scf.log"
    lines = log_path.read_text().splitlines()
    header = next(i for i, line in enumerate(lines)
                  if line.strip() ==
                  "Lattice vectors: (Cartesian coordinate: in unit of a_0)")
    lines[header + 1:header + 4] = [
        " ".join(f"{value:+.6g}" for value in row) for row in vectors
    ]
    lines = [line.replace("lattice constant (Angstrom) = 1",
                          f"lattice constant (Angstrom) = {lts_printed_scale:.6g}")
             .replace("1.889726125837 Bohr ).", f"{scale / Bohr:.12f} Bohr ).")
             for line in lines]
    log_path.write_text("\n".join(lines) + "\n")

    assert abs(10 * float(f"{lts_printed_scale:.6g}") - vectors[0, 0] * scale) > 1e-5
    results = template.read_results(directory)
    assert results["energy"] == pytest.approx(-13.5410411603)


def test_lts_rounded_cell_displaced_log_frame_fails_closed(tmp_path):
    """Print rounding is absorbed; a real geometry difference is still rejected."""
    template, directory = _rounded_cell_run(tmp_path)
    log = directory / "OUT.ABACUS" / "running_scf.log"
    frames = read_abacus_out(log)
    assert _select_scf_frame_for_structure(
        frames, log, directory, [0], legacyio=True
    ) is frames[0]

    displaced = frames[0].copy()
    displaced.positions = displaced.positions + np.array([0.0, 0.0, 0.05])
    with pytest.raises(RuntimeError, match="fail-closed"):
        _select_scf_frame_for_structure(
            [displaced], log, directory, [0], legacyio=True
        )


def test_write_input_backs_up_previous_abacus_json(tmp_path):
    path = tmp_path / "abacus.json"
    path.write_text('{"old": true}')
    atoms = Atoms("H", positions=[[1, 2, 3]], cell=[10, 11, 12], pbc=True)
    AbacusTemplate().write_input(
        SimpleNamespace(pseudo_dir=None, orbital_dir=None), tmp_path, atoms,
        parameters={"calculation": "scf"}, properties=["energy"],
    )
    assert not path.exists()
    assert (tmp_path / "abacus.json.bak.0").read_text() == '{"old": true}'


def _read_results_fixture(tmp_path, calc="scf", stru_text=None, log_text=None):
    """把金样布置成 read_results 期望的目录布局（STRU + OUT.ABACUS/running_<calc>.log）。

    非 scf 时把同一 running_scf.log 复制为 running_<calc>.log，仅用于验证非 scf 分支
    保持 [-1] 末帧语义（不做坐标匹配）；stru_text/log_text 可覆盖 STRU 与日志内容。
    """
    (tmp_path / "STRU").write_text(
        stru_text if stru_text is not None else (GOLD / "STRU").read_text(),
        encoding="utf-8",
    )
    out = tmp_path / "OUT.ABACUS"
    out.mkdir()
    (out / f"running_{calc}.log").write_text(
        log_text if log_text is not None else (GOLD / "running_scf.log").read_text(),
        encoding="utf-8",
    )
    return tmp_path


def _md_fixture_legacy(tmp_path):
    """布置真实 legacy md 金样：OUT.ABACUS/running_md.log（multiframe_md_legacy，真实 ABACUS md 日志）。

    刻意不写 STRU：原生 md 路径不做坐标匹配，缺 STRU 也必须正常返回（不 fail-closed）。
    """
    out = tmp_path / "OUT.ABACUS"
    out.mkdir()
    shutil.copy(TESTFILES / "multiframe_md_legacy" / "running_md.log", out / "running_md.log")
    return tmp_path


def _md_fixture_latest(tmp_path):
    """布置 latestio MD_dump 金样：OUT.ABACUS/running_md.log + MD_dump + eig_occ.txt。"""
    out = tmp_path / "OUT.ABACUS"
    out.mkdir()
    for name in ("running_md.log", "MD_dump", "eig_occ.txt"):
        shutil.copy(TESTFILES / "multiframe_md_latest" / name, out / name)
    return tmp_path


class _SyntheticFrame:
    """Small parser-shaped frame used to exercise fail-closed validation paths."""

    def __init__(self, symbols, positions, cell):
        self._symbols = list(symbols)
        self.positions = np.asarray(positions)
        self.cell = np.asarray(cell)

    def get_chemical_symbols(self):
        return list(self._symbols)


def _synthetic_stru_text(symbols, coordinates, cell, coordinate_type="Direct"):
    """Build a minimal valid STRU fixture with a hand-selected cell and positions."""
    symbols = list(symbols)
    coordinates = np.asarray(coordinates, dtype=float)
    cell = np.asarray(cell, dtype=float)
    species = list(dict.fromkeys(symbols))
    lines = [
        "ATOMIC_SPECIES",
        *[f"{symbol} 1.0 {symbol}.upf" for symbol in species],
        "",
        "LATTICE_CONSTANT",
        f"{1.0 / Bohr:.16g}",
        "",
        "LATTICE_VECTORS",
        *[" ".join(f"{value:.16g}" for value in vector) for vector in cell],
        "",
        "ATOMIC_POSITIONS",
        coordinate_type,
        "",
    ]
    for symbol in species:
        lines.extend([symbol, "0.0"])
        indices = [index for index, value in enumerate(symbols) if value == symbol]
        lines.append(str(len(indices)))
        for index in indices:
            lines.append(
                " ".join(f"{value:.16g}" for value in coordinates[index])
                + " m 1 1 1"
            )
    return "\n".join(lines) + "\n"


def _synthetic_direct_case(tmp_path, symbols=("H", "He")):
    """Return a small non-orthogonal direct-coordinate matcher fixture."""
    cell = np.array(
        [[4.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.2, 0.4, 5.0]],
        dtype=float,
    )
    frac = np.array([[0.1, 0.2, 0.3], [0.65, 0.7, 0.8]], dtype=float)
    stru = tmp_path / "STRU"
    stru.write_text(_synthetic_stru_text(symbols, frac, cell), encoding="utf-8")
    return stru, cell, frac, list(symbols)


def test_golden_last_frame_differs_from_current_structure():
    """金样语义：末帧为试探/P0 重算结构，与当前 STRU 不一致（RED 场景前提）。"""
    frames = read_abacus_out(LOG, sort_atoms_with=None)
    stru = read_stru(GOLD / "STRU")
    last_cart = frames[-1].positions @ np.asarray(frames[-1].cell)
    current_cart = _stru_positions_cartesian(stru, frames[-1].cell)
    assert not np.allclose(last_cart, current_cart, atol=1e-4)


def test_golden_first_frame_is_current_structure():
    """金样语义：帧 0（INIT）== 当前 STRU，是 read_results 应返回的帧。"""
    frames = read_abacus_out(LOG, sort_atoms_with=None)
    stru = read_stru(GOLD / "STRU")
    first_cart = frames[0].positions @ np.asarray(frames[0].cell)
    current_cart = _stru_positions_cartesian(stru, frames[0].cell)
    assert np.allclose(first_cart, current_cart, atol=1e-4)


def test_multiframe_scf_read_results_returns_current_structure_force(tmp_path):
    """scf 多帧累积：read_results 返回当前结构帧（帧 0）的力，而非末帧试探结构的力。"""
    directory = _read_results_fixture(tmp_path)
    results = _template().read_results(directory)

    assert not np.allclose(STRU_FRAME_FORCES, LAST_FRAME_FORCES)  # 语义：两帧力不同
    assert np.allclose(results["forces"], STRU_FRAME_FORCES, atol=1e-12)


def test_scf_read_results_preserves_backend_neutral_efermi(tmp_path):
    """The Fermi level survives the file-I/O boundary under the CP result key."""
    directory = _read_results_fixture(tmp_path)
    results = _template().read_results(directory)

    # The selected parser value in this nspin=1 log is -2.04 eV.
    assert results["efermi"] == pytest.approx(-2.04)
    assert results["fermi_level"] == pytest.approx(results["efermi"])


def test_scf_entropy_metadata_follows_selected_frame(tmp_path):
    """Select energy metadata from the matching geometry, not the final trial."""
    log = (GOLD / "running_scf.log").read_text()
    for energy, entropy in [("-2898.0000000000", "-0.1"), ("-2898.4000000000", "-0.2"),
                            ("-2898.8000000000", "-0.3"), ("-2898.6000000000", "-0.4")]:
        log = re.sub(r"( E_KohnSham[^\n]*" + re.escape(energy) + r"[^\n]*\n)",
                     r"\1 E_entropy(-TS) -0.01 " + entropy + "\n", log)
    results = _template().read_results(_read_results_fixture(tmp_path, log_text=log))
    assert results["energy"] == pytest.approx(-2898.0)
    assert results["free_energy"] == pytest.approx(-2898.0)
    assert results["energy_without_entropy"] == pytest.approx(-2897.9)


def test_missing_entropy_is_not_fabricated(tmp_path):
    """Absent entropy remains unknown while standard energies are preserved."""
    results = _template().read_results(_read_results_fixture(tmp_path))
    assert results["free_energy"] == results["energy"]
    assert "energy_without_entropy" not in results


def test_scf_read_results_fails_closed_when_no_frame_matches(tmp_path):
    """无匹配帧：STRU 与 running log 全部帧不一致时 fail-closed（含 log 路径/帧数/差异摘要）。"""
    stru_text = (
        (GOLD / "STRU")
        .read_text()
        .replace("0.4000000000 0.5000000000 0.5000000000",
                 "0.4300000000 0.5000000000 0.5000000000")
        .replace("0.6000000000 0.5000000000 0.5000000000",
                 "0.5700000000 0.5000000000 0.5000000000")
    )
    directory = _read_results_fixture(tmp_path, stru_text=stru_text)
    with pytest.raises(RuntimeError) as excinfo:
        _template().read_results(directory)
    msg = str(excinfo.value)
    assert "running_scf.log" in msg
    assert "4 frames" in msg
    assert "fail-closed" in msg


def test_native_relax_multiframe_keeps_last_frame_semantics(tmp_path):
    """relax 五路径之一（原生单次调用）：ABACUS 内部 relax 一次调用在 running_relax.log 累积
    多帧，read_results 保持 [-1] 末帧语义、不做坐标匹配（spec R2/P2）。"""
    directory = _read_results_fixture(tmp_path, calc="relax")
    results = _template(calc="relax").read_results(directory)

    assert not np.allclose(STRU_FRAME_FORCES, LAST_FRAME_FORCES)  # 语义：两帧力不同
    assert np.allclose(results["forces"], LAST_FRAME_FORCES, atol=1e-12)


def test_frame_selection_tolerance_both_sides(tmp_path):
    """绝对 Å 容差两侧（spec P3 ③）：容差内近帧命中；容差外无匹配帧 fail-closed。

    同源数据对照：金样帧 0（== 当前 STRU，Cartesian Å 最大差异≈4e-7）在默认 atol=1e-4 内命中；
    帧 1（位移帧，差异 0.2 Å）在 atol=1e-4 外 fail-closed；放宽 atol=0.5 后同一帧命中。
    """
    directory = _read_results_fixture(tmp_path)
    log = directory / "OUT.ABACUS" / "running_scf.log"
    frames = read_abacus_out(log, sort_atoms_with=[0, 1, 2])
    expected = _stru_positions_cartesian(read_stru(directory / "STRU"), frames[0].cell)

    # 容差内：返回"当前结构帧"（帧 0），而非末帧（试探结构）
    selected = _select_scf_frame_for_structure(frames, log, directory, [0, 1, 2], atol=1e-4)
    assert selected is frames[0]
    assert np.allclose(selected.positions @ np.asarray(selected.cell), expected, atol=1e-4)

    # 容差外：仅位移帧（0.2 Å > 1e-4）→ 无匹配帧 fail-closed（异常含日志路径/帧数/差异摘要）
    with pytest.raises(RuntimeError) as excinfo:
        _select_scf_frame_for_structure([frames[1]], log, directory, [0, 1, 2], atol=1e-4)
    msg = str(excinfo.value)
    assert "fail-closed" in msg
    assert "running_scf.log" in msg
    assert "1 frame" in msg

    # 对照：放宽容差（0.5 Å）后同一位移帧命中 → 容差两侧语义
    loose = _select_scf_frame_for_structure([frames[1]], log, directory, [0, 1, 2], atol=0.5)
    assert loose is frames[1]


def test_native_md_legacy_keeps_last_frame_semantics_real_gold(tmp_path):
    """真实 legacy md 金样（multiframe_md_legacy，真实 ABACUS md 日志的首个完整 MD 步）：
    原生 md 走 [-1] 末帧语义、不做坐标匹配——目录无 STRU 也正常返回力/能量（不 fail-closed）。"""
    directory = _md_fixture_legacy(tmp_path)
    results = _template(calc="md", atomorder=[0, 1]).read_results(directory)
    assert np.allclose(results["forces"], MD_LEGACY_FORCES, atol=1e-12)
    assert results["energy"] == pytest.approx(MD_LEGACY_ENERGY)
    assert results["free_energy"] == pytest.approx(MD_LEGACY_ENERGY)
    assert results["energy_without_entropy"] == pytest.approx(-1940.3851061129, abs=1e-9)


def test_native_md_latestio_md_dump_keeps_last_frame_semantics(tmp_path, monkeypatch):
    """latestio MD_dump 路径（ABACUS ≥3.11，md 轨迹在 MD_dump 而非 running log）：
    原生 md 走 [-1] 末帧语义（MD_dump 2 帧取最后一帧的力）、不做坐标匹配——无 STRU 正常返回。"""
    from abacuslite import core as abacus_core

    monkeypatch.setattr(abacus_core, "__LEGACYIO__", False)
    directory = _md_fixture_latest(tmp_path)
    results = _template(calc="md", atomorder=[0, 1]).read_results(directory)
    assert np.allclose(results["forces"], MD_LATEST_FORCES, atol=1e-12)
    assert results["energy"] == pytest.approx(MD_LATEST_ENERGY)
    assert results["free_energy"] == pytest.approx(MD_LATEST_ENERGY)
    assert results["energy_without_entropy"] == pytest.approx(-1940.4234898941, abs=1e-9)


def test_latestio_md_eig_occ_repeatable_read(tmp_path, monkeypatch):
    """Reading MD results twice must preserve eig_occ.txt and its contents."""
    from abacuslite import core as abacus_core

    monkeypatch.setattr(abacus_core, "__LEGACYIO__", False)
    directory = _md_fixture_latest(tmp_path)
    eig_occ = directory / "OUT.ABACUS" / "eig_occ.txt"
    gold_eig_occ = (TESTFILES / "multiframe_md_latest" / "eig_occ.txt").read_text(encoding="utf-8")

    results = _template(calc="md", atomorder=[0, 1]).read_results(directory)
    assert eig_occ.read_text() == gold_eig_occ

    results2 = _template(calc="md", atomorder=[0, 1]).read_results(directory)
    assert np.allclose(results2["forces"], MD_LATEST_FORCES, atol=1e-12)


def test_scf_frame_selection_reused_by_relax_ase_sella_ccqn(tmp_path):
    """relax-ASE / sella / ccqn 以 calculation='scf' 驱动，无专门分支，复用 scf 帧选择。

    ABACUS 无 sella/ccqn 计算类型；relax-ASE 逐步、sella 试探步、ccqn 自循环均在每个试探点
    以 calculation='scf' 调用 ABACUS（多帧累积 running_scf.log），与既有
    test_multiframe_scf_read_results_returns_current_structure_force 同夹具同路径。本冒烟锁定
    该复用语义（sella 2 步停根因回归钉）：返回当前结构帧力，而非末帧试探结构力。
    """
    directory = _read_results_fixture(tmp_path)
    results = _template(calc="scf").read_results(directory)
    assert not np.allclose(STRU_FRAME_FORCES, LAST_FRAME_FORCES)  # 语义：两帧力不同
    assert np.allclose(results["forces"], STRU_FRAME_FORCES, atol=1e-12)


def test_scf_frame_selection_respects_atomorder_revmap(tmp_path):
    """统一排列：非平凡 atomorder（revmap）下 STRU 与帧仍按同一映射对齐（spec P3 ①）。"""
    directory = _read_results_fixture(tmp_path)
    # ASE 序 [H, Au, H] 的 revmap：species 序 [H, H, Au] -> ASE 序
    results = _template(atomorder=[0, 2, 1]).read_results(directory)

    assert np.allclose(results["forces"][1], 0.0, atol=1e-12)  # 中间原子为 Au，力为 0
    assert np.allclose(results["forces"][0], [-0.5, -0.01, 0.0], atol=1e-12)
    assert np.allclose(results["forces"][2], [0.5, 0.01, 0.0], atol=1e-12)


def test_scf_frame_selection_uses_log_coordinate_system_not_stru(tmp_path):
    """F1：帧侧坐标系以 running log 实际坐标头为准，而非 STRU coord_type。

    STRU 为 Cartesian（write_stru 默认），金样日志为 DIRECT（分数帧）——旧实现按 STRU
    判 is_direct=False 会把分数坐标当 Cartesian Å，合法 scf 假 fail-closed；修复后正常
    返回当前结构帧的力。
    """
    directory = _read_results_fixture(tmp_path, stru_text=_cartesian_stru_text())
    results = _template().read_results(directory)
    assert np.allclose(results["forces"], STRU_FRAME_FORCES, atol=1e-12)


def test_scf_frame_selection_cartesian_log_branch(tmp_path):
    """Cartesian log coordinates coincide with Angstrom when alat is 1 Angstrom."""
    directory = _read_results_fixture(tmp_path, log_text=_cartesian_log_text())
    results = _template().read_results(directory)
    assert np.allclose(results["forces"], STRU_FRAME_FORCES, atol=1e-12)


@pytest.mark.parametrize("legacyio", [True, False])
def test_cartesian_frame_matching_scales_lattice_constant(tmp_path, legacyio):
    """Both log readers expose Cartesian coordinates in units of alat."""
    from abacuslite.io import legacyio as legacy, latestio as latest

    directory = _read_results_fixture(tmp_path, log_text=_cartesian_log_text(alat=2.0))
    log = directory / "OUT.ABACUS" / "running_scf.log"
    reader = legacy if legacyio else latest
    trajectory = reader.read_traj_from_running_log(log)
    frames = [Atoms("H2Au", positions=item["coords"], cell=item["cell"], pbc=True)
              for item in trajectory]
    assert trajectory[0]["alat_in_angstrom"] == 2.0
    assert np.allclose(frames[0].positions[0], [2.0, 2.5, 2.5])
    selected = _select_scf_frame_for_structure(
        frames, log, directory, [0, 1, 2], legacyio=legacyio
    )
    assert selected is frames[0]
    if legacyio:
        assert np.allclose(_template().read_results(directory)["forces"], STRU_FRAME_FORCES)


def test_scf_read_results_without_atomorder_uses_identity(tmp_path):
    """F2：atomorder=None（未走 write_input 直接 read_results）时按 identity 守卫，
    scf 帧选择正常工作而非裸 TypeError（与 read_abacus_out sort_atoms_with=None 语义一致）。"""
    directory = _read_results_fixture(tmp_path)
    template = _template()
    template.atomorder = None
    results = template.read_results(directory)
    assert np.allclose(results["forces"], STRU_FRAME_FORCES, atol=1e-12)


def test_scf_selector_accepts_independent_integer_lattice_translations_and_latest_match(tmp_path):
    """PBC identity uses per-atom integer lattice translations and reverse latest-match order."""
    stru, cell, frac, symbols = _synthetic_direct_case(tmp_path)
    shifts = np.array([[1, -2, 0], [-1, 0, 3]], dtype=float)
    latest_shifts = np.array([[-2, 1, 1], [3, -1, -2]], dtype=float)
    frames = [
        Atoms(symbols=symbols, positions=frac, cell=cell, pbc=True),
        Atoms(symbols=symbols, positions=frac + shifts, cell=cell, pbc=True),
        Atoms(symbols=symbols, positions=frac + latest_shifts, cell=cell, pbc=True),
    ]

    selected = _select_scf_frame_for_structure(
        frames, LOG, tmp_path, [0, 1], atol=1e-10
    )

    assert selected is frames[-1]


def test_scf_selector_accepts_cartesian_integer_lattice_translations(tmp_path):
    """Cartesian running logs use the same fractional residual test after unit normalization."""
    _, cell, frac, symbols = _synthetic_direct_case(tmp_path)
    cart = frac @ cell
    stru = tmp_path / "STRU"
    stru.write_text(
        _synthetic_stru_text(symbols, cart, cell, coordinate_type="Cartesian"),
        encoding="utf-8",
    )
    log = tmp_path / "running_cartesian.log"
    log.write_text(_cartesian_log_text(), encoding="utf-8")
    shifts = np.array([[1, -2, 0], [-1, 0, 3]], dtype=float)
    frame = Atoms(
        symbols=symbols,
        positions=cart + shifts @ cell,
        cell=cell,
        pbc=True,
    )

    assert _select_scf_frame_for_structure(
        [frame], log, tmp_path, [0, 1], atol=1e-10
    ) is frame


@pytest.mark.parametrize("delta", [9e-6, 11e-6, 5e-5])
@pytest.mark.parametrize("field", ["cell", "positions"])
def test_default_scf_identity_tolerance_covers_log_print_resolution(
    tmp_path, delta, field
):
    """Default 1e-4 A absorbs six-significant-digit lattice printing.

    ``5e-5`` is the half-ulp a component in ``[10, 100)`` Angstrom can pick up,
    so the former ``1e-5`` default could never match a rounded but correct cell.
    """
    _, cell, frac, symbols = _synthetic_direct_case(tmp_path)
    frame_cell, frame_frac = cell.copy(), frac.copy()
    if field == "cell":
        frame_cell[0, 0] += delta
    else:
        frame_frac[0] += np.array([delta, 0, 0]) @ np.linalg.inv(cell)
    frame = Atoms(symbols=symbols, positions=frame_frac, cell=frame_cell, pbc=True)
    assert _select_scf_frame_for_structure([frame], LOG, tmp_path, [0, 1]) is frame


@pytest.mark.parametrize("field", ["cell", "positions"])
def test_default_scf_identity_tolerance_rejects_beyond_print_resolution(
    tmp_path, field
):
    """A difference above the print resolution still fails closed."""
    delta = 2e-4
    _, cell, frac, symbols = _synthetic_direct_case(tmp_path)
    frame_cell, frame_frac = cell.copy(), frac.copy()
    if field == "cell":
        frame_cell[0, 0] += delta
    else:
        frame_frac[0] += np.array([delta, 0, 0]) @ np.linalg.inv(cell)
    frame = Atoms(symbols=symbols, positions=frame_frac, cell=frame_cell, pbc=True)
    with pytest.raises(RuntimeError):
        _select_scf_frame_for_structure([frame], LOG, tmp_path, [0, 1])
    assert _select_scf_frame_for_structure(
        [frame], LOG, tmp_path, [0, 1], atol=1e-3
    ) is frame


def test_latest_coordinate_parser_uses_periodic_matcher(tmp_path, monkeypatch):
    """The latest running-log parser's Direct header feeds the same matcher contract."""
    from abacuslite.io.latestio import (
        read_traj_from_running_log,
    )

    monkeypatch.setattr(abacus_core, "__LEGACYIO__", False)
    log = Path(__file__).resolve().parents[2] / "abacuslite/io/testfiles/lcao-symm1-nspin1-multik-scf_"
    trajectory = read_traj_from_running_log(log)
    data = trajectory[0]
    cell = np.asarray(data["cell"], dtype=float)
    frac = np.asarray(data["coords"], dtype=float)
    symbols = ["Zn", "O"]
    (tmp_path / "STRU").write_text(
        _synthetic_stru_text(symbols, frac, cell), encoding="utf-8"
    )
    frame = Atoms(
        symbols=symbols,
        positions=frac + np.array([[1, 0, -2], [-1, 2, 1]], dtype=float),
        cell=cell,
        pbc=True,
    )

    assert _select_scf_frame_for_structure(
        [frame], log, tmp_path, [0, 1], atol=1e-8
    ) is frame


def test_scf_selector_uses_fractional_residual_for_triclinic_cell(tmp_path):
    """A triclinic cell must compare lattice coordinates, not Cartesian-axis modulo."""
    stru, cell, frac, symbols = _synthetic_direct_case(tmp_path)
    shift = np.array([[1, -1, 2], [-2, 3, -1]], dtype=float)
    frame = Atoms(symbols=symbols, positions=frac + shift, cell=cell, pbc=True)

    assert _select_scf_frame_for_structure(
        [frame], LOG, stru.parent, [0, 1], atol=1e-10
    ) is frame


def test_scf_selector_skips_valid_position_and_cell_mismatches(tmp_path):
    """Finite valid trial frames may mismatch while an earlier valid frame remains selectable."""
    _, cell, frac, symbols = _synthetic_direct_case(tmp_path)
    matching = Atoms(symbols=symbols, positions=frac, cell=cell, pbc=True)
    different_cell = cell.copy()
    different_cell[0, 0] += 0.25
    valid_mismatch = Atoms(
        symbols=symbols,
        positions=frac + np.array([[0.17, 0.0, 0.0], [0.0, -0.23, 0.0]]),
        cell=different_cell,
        pbc=True,
    )

    assert _select_scf_frame_for_structure(
        [matching, valid_mismatch], LOG, tmp_path, [0, 1], atol=1e-8
    ) is matching


@pytest.mark.parametrize(
    ("kind", "expected_keyword"),
    [
        ("positions-shape", "shape"),
        ("positions-finite", "finite"),
        ("cell-shape", "cell"),
        ("cell-finite", "finite"),
        ("cell-singular", "singular"),
        ("atom-count", "atom"),
        ("symbols", "symbol"),
    ],
)
def test_scf_selector_rejects_corrupt_frame_even_after_earlier_match(
    tmp_path, kind, expected_keyword
):
    """A corrupt later frame must fail the whole SCF instead of hiding behind an earlier match."""
    _, cell, frac, symbols = _synthetic_direct_case(tmp_path)
    matching = Atoms(symbols=symbols, positions=frac, cell=cell, pbc=True)
    corrupt_positions = frac.copy()
    corrupt_cell = cell.copy()
    corrupt_symbols = symbols
    if kind == "positions-shape":
        corrupt_positions = np.zeros((2, 2))
    elif kind == "positions-finite":
        corrupt_positions[0, 0] = np.nan
    elif kind == "cell-shape":
        corrupt_cell = np.zeros((3, 2))
    elif kind == "cell-finite":
        corrupt_cell[1, 1] = np.inf
    elif kind == "cell-singular":
        corrupt_cell[1] = corrupt_cell[0]
    elif kind == "atom-count":
        corrupt_positions = np.vstack([frac, [[0.2, 0.3, 0.4]]])
        corrupt_symbols = [*symbols, "H"]
    elif kind == "symbols":
        corrupt_symbols = ["H", "H"]
    corrupt = _SyntheticFrame(corrupt_symbols, corrupt_positions, corrupt_cell)

    with pytest.raises(RuntimeError) as excinfo:
        _select_scf_frame_for_structure(
            [matching, corrupt], LOG, tmp_path, [0, 1], atol=1e-8
        )
    message = str(excinfo.value).lower()
    assert "frame 1" in message
    assert expected_keyword in message
    assert "fail-closed" in message


def test_scf_selector_does_not_mutate_candidate_atoms_or_stru(tmp_path):
    """Normalization must leave positions, cell, forces, and STRU bytes untouched."""
    _, cell, frac, symbols = _synthetic_direct_case(tmp_path)
    frame = Atoms(symbols=symbols, positions=frac, cell=cell, pbc=True)
    forces = np.array([[1.2, -0.3, 0.7], [-1.2, 0.3, -0.7]])
    frame.calc = SinglePointCalculator(frame, energy=-2.0, forces=forces)
    positions_before = frame.positions.copy()
    cell_before = np.asarray(frame.cell).copy()
    forces_before = frame.calc.results["forces"].copy()
    stru_before = (tmp_path / "STRU").read_bytes()

    assert _select_scf_frame_for_structure(
        [frame], LOG, tmp_path, [0, 1], atol=1e-10
    ) is frame
    np.testing.assert_array_equal(frame.positions, positions_before)
    np.testing.assert_array_equal(np.asarray(frame.cell), cell_before)
    np.testing.assert_array_equal(frame.calc.results["forces"], forces_before)
    assert (tmp_path / "STRU").read_bytes() == stru_before


def test_scf_read_results_rejects_invalid_atomorder_before_parser_reorder(tmp_path):
    """Out-of-range atomorder is a bounded RuntimeError, never a parser IndexError."""
    directory = _read_results_fixture(tmp_path)
    template = _template(atomorder=[0, 1, 99])

    with pytest.raises(RuntimeError) as excinfo:
        template.read_results(directory)
    message = str(excinfo.value)
    assert "atomorder" in message
    assert "STRU" in message
    assert "99" in message


def test_scf_read_results_rejects_raw_extra_atom_before_parser_reorder(tmp_path):
    """Parser sorting must not discard an extra raw atom before SCF validation."""
    directory = _read_results_fixture(tmp_path, log_text=_log_with_extra_atom_per_frame())

    with pytest.raises(RuntimeError) as excinfo:
        _template().read_results(directory)
    message = str(excinfo.value)
    assert "running_scf.log" in message
    assert "frame 0" in message
    assert "atom count" in message
    assert "fail-closed" in message


def test_scf_read_results_rejects_stru_declared_atom_count_mismatch(tmp_path):
    """A STRU species declaration must agree with its actual atom records."""
    stru_text = (GOLD / "STRU").read_text(encoding="utf-8").replace(
        "Au\n0.0\n1\n0.5000000000 0.5000000000 0.1000000000",
        "Au\n0.0\n2\n0.5000000000 0.5000000000 0.1000000000",
        1,
    )
    directory = _read_results_fixture(tmp_path, stru_text=stru_text)

    with pytest.raises(RuntimeError) as excinfo:
        _template().read_results(directory)
    message = str(excinfo.value).lower()
    assert "stru" in message
    assert "au" in message
    assert "declared" in message
    assert "fail-closed" in message


def test_scf_read_results_accepts_translated_real_parser_frame(tmp_path):
    """The real legacy parser and atom reorder preserve a periodic translated SCF frame."""
    directory = _read_results_fixture(tmp_path, log_text=_log_with_translated_first_frame())
    results = _template(atomorder=[0, 2, 1]).read_results(directory)

    np.testing.assert_allclose(
        results["forces"], STRU_FRAME_FORCES[[0, 2, 1]], atol=1e-12
    )


def test_write_input_tracks_custom_stru_file_for_scf_selection(tmp_path):
    """write_input's non-default STRU name is the file consumed by SCF identity checks."""
    atoms = Atoms(
        symbols=["H", "Au", "H"],
        positions=[[4.0, 5.0, 5.0], [5.0, 5.0, 1.0], [6.0, 5.0, 5.0]],
        cell=[10.0, 10.0, 10.0],
        pbc=True,
    )
    template = AbacusTemplate()
    template.write_input(
        SimpleNamespace(pseudo_dir=None, orbital_dir=None),
        tmp_path,
        atoms,
        parameters={"calculation": "scf", "stru_file": "CURRENT.STRU"},
        properties=["energy", "forces"],
    )
    out = tmp_path / "OUT.ABACUS"
    out.mkdir()
    shutil.copy(GOLD / "running_scf.log", out / "running_scf.log")
    # A misleading default STRU proves read_results did not hard-code its filename.
    (tmp_path / "STRU").write_text(
        (GOLD / "STRU").read_text().replace(
            "0.4000000000 0.5000000000 0.5000000000",
            "0.4300000000 0.5000000000 0.5000000000",
        ),
        encoding="utf-8",
    )

    assert template.stru_file == "CURRENT.STRU"
    results = template.read_results(tmp_path)
    np.testing.assert_allclose(
        results["forces"], STRU_FRAME_FORCES[[0, 2, 1]], atol=1e-12
    )
