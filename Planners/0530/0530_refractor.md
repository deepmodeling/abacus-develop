# ABACUS 重构 PR —— 三人分工方案

## 零冲突保证

四个重构点修改的文件集合经逐文件核对，**完全不相交**。三人可在不同分支并行工作，不会产生 git merge 冲突。

---

## 分工总览

| 人   | 重构编号 | 内容                           | 涉及文件数   | 难度  | 分支名                           |
| ---- | -------- | ------------------------------ | ------------ | ----- | -------------------------------- |
| A    | #5 + #9  | 积分器工厂 + DP/NEP 缓冲区复用 | 6            | 低-中 | `refactor/md-factory-and-buffer` |
| B    | #7       | 统计纯函数                     | 6（+1 新建） | 中    | `refactor/md-pure-statistics`    |
| C    | #12      | 测试夹具                       | 8（+1 新建） | 低    | `refactor/md-test-fixtures`      |

---

## 人 A：`#5 积分器工厂` + `#9 缓冲区复用`

| 文件                                    | 操作                                  | 参考章节（`代码修改与重构.md`） |
| --------------------------------------- | ------------------------------------- | ------------------------------- |
| `source/source_md/run_md.cpp`           | 修改 — 抽取工厂函数 + unique_ptr      | §3.3                            |
| `source/source_md/run_md.h`             | 修改 — 前向声明 MD_base               | §3.4                            |
| `source/source_esolver/esolver_dp.h`    | 修改 — 新增 cell/coord/f/v 缓冲区成员 | §2.3                            |
| `source/source_esolver/esolver_dp.cpp`  | 修改 — 复用成员缓冲区                 | §2.4                            |
| `source/source_esolver/esolver_nep.h`   | 修改 — 新增 cell/coord 缓冲区成员     | §2.5                            |
| `source/source_esolver/esolver_nep.cpp` | 修改 — 复用成员缓冲区                 | §2.6                            |

**提交信息**：

```
refactor: extract MD integrator factory and reuse DP/NEP buffers

- Replace manual new/delete in md_line() with factory function + unique_ptr
- Move DP/NEP per-step cell/coord/f/v allocations to member buffers
  allocated once in before_all_runners()

Co-Authored-By: ... <...>
```

**建议两个 commit**：#5 和 #9 分别提交，共用一个分支。

---

## 人 B：`#7 抽取温度、动能、应力统计为纯计算组件`

| 文件                               | 操作                                    | 参考章节（`代码修改与重构.md`） |
| ---------------------------------- | --------------------------------------- | ------------------------------- |
| `source/source_md/md_statistics.h` | **新建**                                | §4.3                            |
| `source/source_md/md_func.h`       | 修改 — 新增纯函数声明                   | §4.4                            |
| `source/source_md/md_func.cpp`     | 修改 — 新增纯函数实现，旧接口改为包装器 | §4.5                            |
| `source/source_md/nhchain.cpp`     | 修改 — 调用点改用纯函数                 | §4.6                            |
| `source/source_md/msst.cpp`        | 修改 — 调用点改用纯函数                 | §4.7                            |
| `source/source_md/fire.cpp`        | 可选 — `calc_fire_projection` 纯函数    | §4.8                            |

**提交信息**：

```
refactor: extract MD statistics as pure functions

- Add md_statistics.h with MDKineticState, MDStressState, FIREProjection
- Add calc_kinetic_state() and calc_stress_state() pure functions
- Old interfaces (current_temp, compute_stress) kept as wrappers
- Update nhchain.cpp and msst.cpp call sites

Co-Authored-By: ... <...>
```

**注意事项**：

1. `calc_stress_state` 中 `state.t_vector.create(3,3)` 之后直接做 `+=`，需确认 `create` 是否清零。如果不清零，需在 create 后手动把 t_vector 元素置零。
2. `fire.cpp` 的修改是可选的，稳妥起见可以先不改，留在后续 PR。

---

## 人 C：`#12 建立测试夹具/benchmark 夹具`

| 文件                                      | 操作                                           | 参考章节（`代码修改与重构.md`） |
| ----------------------------------------- | ---------------------------------------------- | ------------------------------- |
| `source/source_md/test/md_test_fixture.h` | **新建**                                       | §1.3                            |
| `source/source_md/test/verlet_test.cpp`   | 修改 — 继承 MdIntegratorFixture\<Verlet\>      | §1.4                            |
| `source/source_md/test/fire_test.cpp`     | 修改 — 继承 MdIntegratorFixture\<FIRE\>        | §1.5                            |
| `source/source_md/test/langevin_test.cpp` | 修改 — 继承 MdIntegratorFixture\<Langevin\>    | §1.6                            |
| `source/source_md/test/msst_test.cpp`     | 修改 — 继承 MdIntegratorFixture\<MSST\>        | §1.7                            |
| `source/source_md/test/nhchain_test.cpp`  | 修改 — 继承 MdTestBase，手动创建 Nose_Hoover   | §1.8                            |
| `source/source_md/test/md_func_test.cpp`  | 修改 — 继承 MdFuncTestFixture，数组改用 vector | §1.9                            |
| `source/source_md/test/lj_pot_test.cpp`   | 修改 — 继承 LjPotTestFixture                   | §1.10                           |
| `source/source_md/test/CMakeLists.txt`    | 修改 — 移除重复的 `global_variable.cpp` 行     | §1.11                           |

**提交信息**：

```
refactor: introduce MD test fixtures to reduce duplication

- Add md_test_fixture.h with MdTestBase, MdIntegratorFixture<T>,
  MdFuncTestFixture, LjPotTestFixture
- Replace manual new/delete in all MD test SetUp/TearDown with RAII
- Remove duplicate dependency in test CMakeLists.txt

Co-Authored-By: ... <...>
```

**注意事项**：

1. `nhchain_test.cpp` 的 NHC 测试没有继承 `MdIntegratorFixture<Nose_Hoover>`，而是直接继承 `MdTestBase` 再手动创建积分器，因为 NHC 需要在 SetUp 中额外设置 NPT 参数。这在 `代码修改与重构.md` §1.8 的代码里已经体现了。
2. 各测试的 TEST_F 主体代码**完全不变**，只有 `dynamic_cast<XXX*>(mdrun)` 改成 `dynamic_cast<XXX*>(mdrun.get())`。
3. 做完后运行 `ctest --test-dir build -R 'MODULE_MD' --output-on-failure` 确认所有测试通过。

---

## Git 操作流程

### 前置准备（一人做）

```bash
# 在 GitHub 网页上 Fork https://github.com/deepmodeling/abacus-develop
# 然后：
git clone https://github.com/<fork-user>/abacus-develop.git abacus-pr
cd abacus-pr
git remote add upstream https://github.com/deepmodeling/abacus-develop.git
```

三人可以 clone 同一个 fork 到各自目录，或共用同一个 clone 开不同分支。

### 各自操作（三人完全并行）

```bash
cd abacus-pr

# 1. 从最新 main 出发
git checkout main
git pull upstream main

# 2. 创建独立分支
git checkout -b <你的分支名>

# 3. 修改文件（参考上方分工表和"代码修改与重构.md"对应章节）

# 4. 构建验证
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build -j$(nproc)

# 5. 运行测试
ctest --test-dir build -R 'MODULE_MD' --output-on-failure

# 6. 提交并推送
git add <修改的文件>
git commit -m "<提交信息>"
git push -u origin <你的分支名>
```

### 在 GitHub 创建 PR

推送后在 GitHub 网页上：

- 进入你的 fork 页面
- 点 "Compare & pull request"
- base 选 `deepmodeling/abacus-develop` 的 `main` 分支
- head 选你刚推的分支
- PR 标题用提交信息的第一行
- 正文写清楚改动目的和范围

---

## 建议 PR 顺序

1. **C 先提** — 测试夹具最安全，reviewer 最容易通过
2. **B 接着提** — 统计纯函数依赖旧接口包装器兼容
3. **A 最后提** — 工厂和缓冲区分属两个独立改动，合在一个 PR 也可以

三人也可以同时提（文件不冲突），但分先后能避免 reviewer 困惑。

---

## 额外可做（进度快时）

- **#11 FIRE 参数覆盖 bugfix**：`fire.cpp` 构造函数里 `force_thr = 1e-3` 覆盖了用户输入。这是一行修复，可以分给人 B 顺便做，单独一个 commit。只改 `fire.cpp`，不与其他人的改动冲突。
- **#8 ESolver 结果与日志分离**：改 `esolver_dp.cpp`、`esolver_nep.cpp`，和人 A 的 #9 有同一文件冲突，不建议今晚做。