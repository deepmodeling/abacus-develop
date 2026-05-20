# NSCF 测试集分析与设计文档

## 1. 背景

tests/17_DS_DFTU 测试集中包含 6 个 NSCF (non-self-consistent field) 相关算例，覆盖了 DeltaSpin 和 DFT+U 在非自洽计算模式下的功能验证。这些算例需要预收敛的电荷密度文件（`autotest-CHARGE-DENSITY.restart`）和 DFT+U 相关的 onsite 矩阵文件（`onsite.dm`）。

### 1.1 现有 NSCF 算例清单

| # | 算例 | 基组 | 功能 | 电荷密度大小 |
|---|------|------|------|-------------|
| 55 | PW_DS_NSCF_S4_XY | PW | DeltaSpin + nscf | 247K |
| 60 | PW_DFTU_DS_NSCF_Band_XY | PW | DFT+U + DeltaSpin + band | 85K |
| 61 | LCAO_DS_NSCF_S4_XY | LCAO | DeltaSpin + nscf | 30K |
| 62 | LCAO_DFTU_NSCF_Band_XY | LCAO | DFT+U + band (实际是 scf, scf_nmax=1) | 30K |
| 63 | LCAO_DFTU_DS_NSCF_Band_XY | LCAO | DFT+U + DeltaSpin + band | 30K |
| 64 | PW_DFTU_NSCF_Band_XY | PW | DFT+U + band | 85K |

### 1.2 当前问题

1. **文件体积**：PW 算例的电荷密度文件较大（85K-247K），随测试系统增大可能进一步增长
2. **维护成本**：预收敛文件需要与代码版本同步更新，增加维护负担
3. **CI 稳定性**：部分算例因收敛性问题被禁用（如 62、63）
4. **依赖管理**：虽然当前算例自包含文件，但未来扩展可能引入跨测试依赖

## 2. NSCF 核心功能分析

### 2.1 NSCF 计算的关键逻辑

NSCF 模式 (`calculation = "nscf"`) 与 SCF 模式的主要区别：

| 功能点 | SCF | NSCF |
|--------|-----|------|
| 电荷密度更新 | 迭代更新 | 跳过 (`skip_charge = true`) |
| 哈密顿量构建 | 基于自洽电荷 | 基于读取的电荷 |
| 本征值求解 | 迭代直至收敛 | 单次对角化 |
| DeltaSpin lambda 优化 | 在 SCF 循环内 | 在固定电荷下独立迭代 |
| 能带输出 | 可选 | 主要用途之一 |

### 2.2 需要测试的核心功能模块

```
NSCF 核心功能
├── 电荷密度读取与初始化
│   ├── PW: 平面波系数读取 (charge density in PW basis)
│   └── LCAO: 密度矩阵读取 (density matrix in LCAO basis)
├── 固定电荷下的哈密顿量构建
│   ├── DFT+U 有效势 (eff_pot_pw_index / cal_v_of_u)
│   └── DeltaSpin 约束势 (delta_hcc / pauli_to_moment)
├── DeltaSpin lambda 迭代优化
│   ├── 固定电荷下的 lambda 收敛
│   ├── sc_thr 收敛判断
│   └── direction_only 模式
├── 能带结构计算
│   ├── 高对称 k 点路径处理
│   └── band.txt 输出格式
└── onsite 矩阵处理
    ├── onsite.dm 读取与验证
    └── PW 投影子索引映射
```

## 3. 测试方案设计

### 3.1 方案对比

| 方案 | 优点 | 缺点 | 适用场景 |
|------|------|------|---------|
| **A. 纯单元测试** | 快速、独立、无需大文件 | 覆盖范围有限，难以测试端到端流程 | 核心算法逻辑验证 |
| **B. SCF+NSCF 工作流** | 完整覆盖、自生成数据 | 耗时较长、依赖 SCF 稳定性 | 集成测试、回归测试 |
| **C. 混合方案 (推荐)** | 兼顾速度与覆盖率 | 设计复杂度较高 | CI + 开发者测试 |

### 3.2 推荐方案：分层测试策略

```
┌─────────────────────────────────────────────────────────────┐
│ Layer 1: 单元测试 (Unit Tests)                              │
│ - 目标: 核心算法逻辑验证                                     │
│ - 位置: source/*/test/ 目录下                               │
│ - 数据: 手工构造的小规模测试数据                              │
│ - 耗时: < 1s per test                                     │
├─────────────────────────────────────────────────────────────┤
│ Layer 2: 组件测试 (Component Tests)                          │
│ - 目标: 模块间交互验证                                       │
│ - 位置: tests/17_DS_DFTU/ 目录下                            │
│ - 数据: 小系统预收敛文件 (< 50K)                            │
│ - 耗时: < 5s per test                                     │
├─────────────────────────────────────────────────────────────┤
│ Layer 3: 集成测试 (Integration Tests)                        │
│ - 目标: 完整工作流验证                                       │
│ - 位置: tests/17_DS_DFTU/ 目录下                            │
│ - 数据: SCF 生成电荷密度 → NSCF 读取                        │
│ - 耗时: 10-60s per test                                   │
└─────────────────────────────────────────────────────────────┘
```

### 3.3 各层详细设计

#### Layer 1: 单元测试 (新增)

**建议新增测试文件：**

1. **`nscf_charge_read_test.cpp`** — 电荷密度读取逻辑
   ```
   测试用例:
   - PWChargeReadTest: 验证平面波电荷系数读取
   - LCAODensityMatrixReadTest: 验证密度矩阵读取
   - ChargeSkipTest: 验证 skip_charge 标志正确跳过电荷更新
   ```

2. **`nscf_hamilt_build_test.cpp`** — 固定电荷下哈密顿量构建
   ```
   测试用例:
   - FixedChargeHamiltTest: 验证电荷固定时哈密顿量构建正确
   - DFTU_NSCF_PotentialTest: 验证 NSCF 模式下 DFT+U 有效势
   - DeltaSpin_NSCF_ConstraintTest: 验证 NSCF 模式下约束势
   ```

3. **`nscf_band_output_test.cpp`** — 能带输出验证
   ```
   测试用例:
   - BandDataLayoutTest: 验证能带数据布局
   - BandKPointMappingTest: 验证 k 点路径映射
   - BandOutputFormatTest: 验证 band.txt 格式
   ```

#### Layer 2: 组件测试 (改造现有)

**改造策略：**

1. **减小现有算例规模**：
   - 使用更小的基组（降低 ecutwfc）
   - 减少 k 点数量
   - 使用更简单的晶体结构

2. **文件压缩**：
   - 对电荷密度文件进行压缩存储
   - CI 运行时解压

3. **共享预收敛数据**：
   - 建立 `common_charge_data/` 目录
   - 多个算例共享同一 SCF 收敛结果

#### Layer 3: 集成测试 (工作流)

**工作流设计：**

```
SCF 阶段                          NSCF 阶段
┌─────────────┐                  ┌─────────────┐
│ 1. SCF 计算  │ ── charge ──→   │ 4. 读取电荷  │
│    (小系统)  │                  │    密度文件   │
└─────────────┘                  └─────────────┘
         │                                │
         ↓                                ↓
┌─────────────┐                  ┌─────────────┐
│ 2. 验证收敛  │                  │ 5. NSCF 计算 │
│    (能量/密度)│                  │    (固定电荷) │
└─────────────┘                  └─────────────┘
         │                                │
         ↓                                ↓
┌─────────────┐                  ┌─────────────┐
│ 3. 保存电荷  │                  │ 6. 验证输出  │
│    密度文件  │                  │    (能带/磁矩) │
└─────────────┘                  └─────────────┘
```

**工作流实现选项：**

| 选项 | 实现方式 | 复杂度 |
|------|---------|--------|
| **A. CTest fixture** | 使用 CTest fixture 功能定义依赖 | 中 |
| **B. 脚本编排** | Python/Bash 脚本串联 SCF→NSCF | 低 |
| **C. 单算例两阶段** | INPUT 中设置 `calculation = scf` 后手动切换 | 低 |

**推荐选项 B**：使用 Python 脚本编排工作流，灵活性最高。

### 3.4 具体实施建议

#### 阶段一：单元测试覆盖核心逻辑 (1-2 周)

1. 创建 `source/source_esolver/test/nscf_charge_read_test.cpp`
2. 创建 `source/source_hsolver/test/nscf_hamilt_build_test.cpp`
3. 创建 `source/source_io/test/nscf_band_output_test.cpp`
4. 集成到现有 CMake 测试系统

#### 阶段二：组件测试优化 (1 周)

1. 审查现有 55/60/61/62/63/64 算例
2. 缩小算例规模（ecutwfc, k 点数）
3. 建立共享电荷密度目录
4. 更新 CASES_CPU.txt 启用被禁用的算例

#### 阶段三：集成工作流测试 (2-3 周)

1. 设计 Python 工作流编排脚本
2. 创建 2-3 个典型 SCF→NSCF 工作流测试
3. 集成到 CI 系统（可选，仅在 nightly 运行）

## 4. 代码修改建议

### 4.1 NSCF 相关代码可测试性改进

```cpp
// 当前代码 (esolver_ks_pw.cpp:208)
bool skip_charge = PARAM.inp.calculation == "nscf" ? true : false;

// 建议：提取为可测试函数
namespace nscf_utils {
    bool should_skip_charge_update(const std::string& calculation);
    bool should_read_charge_from_file(const std::string& init_chg);
    void validate_nscf_input(const Input_Param& inp);
}
```

### 4.2 电荷密度读取接口抽象

```cpp
// 建议：抽象电荷密度读取接口
class ChargeDensityReader {
public:
    virtual ~ChargeDensityReader() = default;
    virtual void read_charge(const std::string& dir, 
                             const std::string& suffix) = 0;
    virtual bool is_valid() const = 0;
};

// PW 实现
class PWChargeDensityReader : public ChargeDensityReader { ... };

// LCAO 实现
class LCAOChargeDensityReader : public ChargeDensityReader { ... };
```

### 4.3 工作流测试辅助函数

```python
# tests/17_DS_DFTU/workflow_helpers.py
def run_scf_then_nscf(scf_input, nscf_input, system_dir):
    """Run SCF calculation, then use charge density for NSCF."""
    # 1. Run SCF
    scf_result = run_abacus(scf_input, system_dir)
    assert scf_result.converged
    
    # 2. Copy charge density
    copy_charge_density(scf_result.output_dir, nscf_input.read_dir)
    
    # 3. Run NSCF
    nscf_result = run_abacus(nscf_input, system_dir)
    
    # 4. Validate
    validate_nscf_output(nscf_result, expected_values)
    return nscf_result
```

## 5. 总结与建议

### 5.1 核心结论

| 问题 | 结论 |
|------|------|
| 电荷密度文件是否过大？ | 当前 30K-247K 尚可接受，但随系统增大会成问题 |
| 能否用单元测试覆盖？ | **部分可以**：核心算法逻辑可单元测试，端到端流程需集成测试 |
| SCF+NSCF 工作流是否必要？ | **是**：用于验证完整数据流和数值一致性 |
| 推荐方案？ | **分层测试策略**：单元测试 (快速) + 组件测试 (中速) + 工作流测试 (完整) |

### 5.2 优先级排序

1. **高优先级**：单元测试覆盖 NSCF 核心逻辑（charge read, hamilt build, band output）
2. **中优先级**：优化现有组件测试算例规模，建立共享数据目录
3. **低优先级**：SCF+NSCF 工作流集成测试（可作为 nightly CI）

### 5.3 预期收益

| 指标 | 当前 | 实施后 |
|------|------|--------|
| 单元测试覆盖率 | ~60% | ~85% |
| CI 测试耗时 | 2-5 min | 1-2 min (单元) + 5-10 min (集成) |
| 维护成本 | 高 (大文件同步) | 中 (小文件 + 自动生成) |
| 问题定位速度 | 慢 (集成测试失败难定位) | 快 (单元测试快速定位) |

## 6. 参考资料

- ABACUS 用户手册: https://abacus.deepmodeling.com/en/latest/advanced/opt-latt.html
- NSCF 计算文档: https://abacus.deepmodeling.com/en/latest/quick_start/calculation.html
- DeltaSpin 模块文档: source/source_lcao/module_deltaspin/README.md
- DFT+U 模块文档: source/source_lcao/module_dftu/README.md
