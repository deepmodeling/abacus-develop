# ABACUS CI/CD 提速方案（可落地版）

> 配套文档：《CI速度差异分析.md》（同目录，含全部实测数据）
> 原则：每个方案给出**具体文件改动 + 预期收益 + 风险 + 验证方法**，可直接拆成 PR。
> 已验证的关键事实：
> - ccache 已通过 `CMakeLists.txt:299-303` 接入全部语言（含 CUDA），装了就用，无需改 CMake
> - CMake 只用 `git log -1` 取 commit 信息（`CMakeLists.txt:163-174`），**不需要 fetch-depth: 0**
> - `Autotest.sh -n` 是每个用例的 MPI 进程数，不是用例并行数
> - 上游 cuda.yml 仍是 `-j4` 单 job；你的 f3e844ff7（nproc+matrix）只在本地/已关闭的 #7713

---

## 方案总览（按 ROI 排序）

| 方案 | 改动量 | 预期收益 | 风险 | 建议批次 |
|---|---|---|---|---|
| B. CUDA 编译提速 | cuda.yml 2 行 | Build 33min → 8~15min | 极低 | **PR-1（立即）** |
| D3. checkout 浅克隆 | test.yml 1 行 | 冷节点省 ~12min | 极低 | **PR-1（立即）** |
| C. ccache 可靠化 | 2 个 yml 数行 + 1 个新 yml | 消除 2.7 vs 37.5min 波动 | 低 | PR-2 |
| A. 镜像稳定化 | 2 个 Dockerfile 删段落 | init 8~78min → ~1min | 中（需维护者共识） | PR-3 |
| E. toolchain 归零 | 随 A 一起做 | 冷节点省 3~10min | 低 | 随 PR-3 |
| D. 测试并行化 | test.yml 重构 / cuda.yml | CPU 18min→~6min | 中 | PR-4 |

预期总效果：**test.yml 24~100min → 12~20min；cuda.yml 53~131min → 15~25min。**

---

## 方案 B：CUDA 编译提速（PR-1，最先做）

### B1. 裁剪 CUDA 架构（7 种 → 1~2 种）

**问题**：`CMakeLists.txt:445-476` 未指定时，CUDA 12.2 下默认编 60/70/75/80/86/89/90 共 7 种架构，61 个 `.cu` 每个编 7 遍。CI 只需要能在 CI 的那块 GPU 上跑。

**改动**（`.github/workflows/cuda.yml` Configure & Build 步骤）：

```yaml
      - name: Configure & Build
        run: |
          nvidia-smi
          source toolchain/install/setup
          rm -rf build
          cmake -B build -G Ninja \
            -DUSE_CUDA=ON \
            -DBUILD_TESTING=ON \
            -DENABLE_FLOAT_FFTW=ON \
            -DCMAKE_CUDA_ARCHITECTURES=80 \   # ← 新增：按 CI 实际 GPU 填
            -DCMAKE_INSTALL_PREFIX="${GITHUB_WORKSPACE}/install"
          cmake --build build -j $(nproc)      # ← B2：-j4 改掉
          cmake --install build
```

**架构怎么填**：打开最近一次 CUDA run 的 Configure & Build 日志，第一行 `nvidia-smi` 输出就是卡型：

| GPU | Compute Capability | 填 |
|---|---|---|
| V100 | 7.0 | `70` |
| T4 | 7.5 | `75` |
| A100 | 8.0 | `80` |
| A10 / RTX 30 系 | 8.6 | `86` |
| L40 / RTX 40 系 | 8.9 | `89` |
| H100 | 9.0 | `90` |

池子里若有多种卡，填 2 个（如 `80;90`），仍比 7 种省 ~70% nvcc 工作量。**不要**用 `native`——你本地 matrix 版把 build 和 test 拆成了不同 job，可能落在不同型号的机器上。

**预期收益**：33min → 8~15min（-70%~-85% nvcc 时间）。测试有效性不受影响（二进制只需在 CI GPU 上可运行）。

### B2. 放开并行度 `-j4` → `-j $(nproc)`

**问题**：`-j4` 自 #4032 时代写死。GPU 节点通常 16+ 核。

**注意**：必须和 B1 一起做。7 架构 × 高并行 → nvcc 单进程 2~4GB 内存，容易 OOM（你之前关心的 OOM 问题就在这里）；架构裁到 1~2 种后放开并行是安全的。你本地 commit f3e844ff7 已包含此项，可直接摘出来。

---

## 方案 D3：checkout 浅克隆（并入 PR-1）

**问题**：test.yml `fetch-depth: 0` 全历史克隆，冷节点实测 11.9min。已验证 CMake 只执行 `git log -1`，不需要历史。

**改动**（`.github/workflows/test.yml`）：

```yaml
      - name: Checkout repository
        uses: actions/checkout@v7
        with:
          fetch-depth: 1        # ← 0 改 1
          submodules: 'false'
```

**风险排查**：`git submodule update --init --recursive` 不受 depth 影响；`version_check.yml` 是独立 workflow，如它依赖 tag 历史则保持它自己的 depth 不变（与本 PR 无关）。

---

## 方案 C：ccache 可靠化（PR-2，消除最大波动源）

### C1. 显式容量 + 可观测（改动极小，先做）

**问题**：两个 workflow 只挂载 `/tmp/ccache`，从未设置容量（ccache 4.x 默认约 10GB，7 架构 CUDA 产物 + 日常 churn 有驱逐压力），也从不打印命中率——慢的时候无法自证。

**改动**（test.yml 和 cuda.yml 的 Build 步骤前后各加一步）：

```yaml
      - name: Setup ccache
        run: |
          ccache --max-size=30G
          ccache --zero-stats
          ccache -s

      # ... Configure / Build ...

      - name: ccache statistics
        if: always()
        run: ccache -s
```

**收益**：①容量兜底；②每次运行的日志里直接看到 `cache hit rate`，下次再出现"2min vs 37min"一眼定位是缓存问题还是机器问题。这也是给维护者讲故事的证据。

### C2. 确认动态池的缓存持久化（需要 runner 管理员配合）

**问题**：6/17 #7476 切到动态 runner 池（K8s Pod），`/tmp/ccache` 是宿主机目录还是 emptyDir 决定了缓存能不能跨 run 存活。实测同一 PR 越推越快，说明**部分节点**是能保住的，但新扩的节点是冷的。

**行动**（在 issue/PR 里 @ 维护者确认两件事）：
1. `/tmp/ccache` 是否 hostPath/PV？多节点池是否每节点独立一份？（独立 → 命中率随节点数稀释）
2. 若是 emptyDir，改为 hostPath 或 PVC；
3. 终极方案：ccache 4.x 支持二级远程缓存，池内共享一份：
   ```yaml
   env:
     CCACHE_REMOTE_STORAGE: "http://内网缓存服务|layout=bazel"   # 或 redis:/s3
     CCACHE_REMOTE_ONLY: "false"   # 本地+远程双层
   ```
   这样新扩的冷节点也能命中远程缓存，**直接消除"节点运气"因子**。

### C3. cache-warmer：给池子预热（消除"合并风暴后第一批 PR 变慢"）

**问题**：大重构合入 develop 后，全池缓存失效，之后第一批 PR 的 CI 全部变慢（7/29-31 的 347 文件/125 头文件风暴就是案例）。

**方案**：新增 `.github/workflows/cache_warmer.yml`，在 develop 有合并时自动"只编不测"，让池子缓存始终贴着最新 develop：

```yaml
name: Cache Warmer
on:
  push:
    branches: [develop]
  workflow_dispatch:

concurrency:
  group: cache-warmer
  cancel-in-progress: true

jobs:
  warm:
    runs-on: X64
    if: github.repository_owner == 'deepmodeling'
    container:
      image: ghcr.io/deepmodeling/abacus-gnu   # 方案 A 落地后换稳定 tag
      volumes:
        - /tmp/ccache:/github/home/.ccache
    steps:
      - uses: actions/checkout@v7
        with:
          fetch-depth: 1
          submodules: recursive
      - name: Build only (prime ccache)
        run: |
          sudo apt-get update && sudo apt-get install -y ccache ninja-build gfortran
          cmake -B build -G Ninja -DBUILD_TESTING=ON
          cmake --build build -j $(nproc)
          ccache -s
```

GPU 池同理可加一条 `runs-on: gpu` 的 warm job（USE_CUDA=ON）。成本是 runner 时间，收益是所有 PR 的 Build 稳定在"只差 PR 自己改动"的水平。

---

## 方案 A：镜像稳定化（PR-3，收益最大，需维护者共识）

### 问题回顾

- `devcontainer.yml`：**每次 push 到 develop** 都重建并推送 `abacus-{gnu,intel,cuda}:latest`；
- `Dockerfile.gnu:31-41` / `Dockerfile.cuda:27-36`：cache-bust（`ADD .../develop /dev/null`）+ **在镜像里完整编译一份 ABACUS**；
- 结果：镜像几乎天天变 → 动态 runner 天天重拉（init 实测 8~78min）；
- 而 CI 从不用镜像里那份二进制（每个 PR 都 checkout 自己重编）——**镜像内编译对 CI 是纯浪费**。

### 改法 A1（推荐：删，不加新东西，最容易过审）

从 `Dockerfile.gnu` 删 31-41 行、`Dockerfile.cuda` 删 27-36 行（cache-bust + 镜像内编译段落）。效果：

- 镜像内容只在**依赖本身变化**（改 Dockerfile）时才变 → devcontainer 每次重建产出的镜像层摘要不变 → runner 端 `docker pull` 变 no-op → **Initialize containers 稳定 ~0-1min**；
- 镜像体积：gnu 减 ~1GB、cuda 减数 GB（abacus 构建产物 + 源码）→ 真需要重拉时也更快；
- 符合 AGENTS.md「复用现有 Docker 资产、不新增容器」——这是做减法。

**需要处理的一个副作用**：这两个镜像的定位本来是"给用户体验 ABACUS"（Dockerfile 注释写明 aimed for evaluating ABACUS），删掉内置二进制后，想直接 `docker run` 用 abacus 的用户会受影响。PR 里提供两个选项让维护者选：
- 选项 1：接受——evaluating 用户改用 release 二进制/自行编译，README/文档更新一句；
- 选项 2：把"含 abacus 的 eval 镜像"拆成单独的 `Dockerfile.eval`（CI 用纯净 deps 镜像，eval 镜像照旧天天新）——代价是多一个 Dockerfile，需要在 PR 里按 AGENTS.md 要求写明理由。

### 改法 A2（可叠加：拉镜像换国内源）

`Dockerfile.gnu` 注释里已写明镜像同时发布在 `registry.dp.tech/deepmodeling/abacus-*`，devcontainer.yml 也推 `dp-harbor-registry.us-east-1.cr.aliyuncs.com`。runner 在国内，却让它们从 ghcr.io 拉。

```yaml
    container:
      image: registry.dp.tech/deepmodeling/abacus-cuda   # 原: ghcr.io/deepmodeling/abacus-cuda
```

一行改动。先验证匿名可拉（找台机器 `docker pull registry.dp.tech/deepmodeling/abacus-gnu:latest` 试试），拉不通就用 dp-harbor 地址（需配 docker login，workflow 里用 secrets）。

---

## 方案 E：toolchain 时间归零（随 PR-3 的镜像一起改）

**问题**：7/1 #7449 起每次运行跑完整 `install_stage4.sh`（dftd4/cereal/rapidjson/**libtorch ~200MB 下载**/libnpy/libri/libcomm/nep）。热节点 20 秒，冷节点数分钟；且 **cuda.yml 的构建根本没开 MLALGO，libtorch 白下载**。

**改法**：
1. 把 stage4 产物预装进 deps 镜像（Dockerfile 里加一段 stage4 安装，libtorch gnu 镜像已有）；
2. CI 的 toolchain 步骤保留（幂等，检测到已装会秒过）——热路径不变，冷路径归零；
3. cuda.yml 如确认不需要 libtorch/LibRI，toolchain 步骤可改为只装必需项（`install_dftd4.sh` 等单独调用），或 configure 时显式关掉对应 feature 避免误装。

---

## 方案 D：测试提速（PR-4，改动最大放最后）

### D1. test.yml：26 个串行步骤 → 分组并行

**现状**：单 job 里 16 个 Module_* 单测 + 10 个集成套件串行跑 ~18min（实测：01_PW 4.4、08_EXX 3.0、03_NAO_multik 2.6、09_DeePKS 2.1…）。

**轻量版（diff 小，推荐先上）**：16 个 Module_* 步骤合并成一个并行调用：

```yaml
      - name: All Module Unittests
        env:
          GTEST_COLOR: 'yes'
          OMP_NUM_THREADS: '2'
        run: |
          ctest --test-dir build -j $(nproc) --timeout 1700 \
            -R 'MODULE_' -E PERF_MODULE_HSOLVER_KERNELS
```

集成套件（每个内部 Autotest.sh 串行跑用例）保持分步即可。**收益：单测 ~4.5min → ~1min。**

**完整版（参照你 cuda.yml 的 matrix 改法）**：build job 产出 `install/` artifact → 10 个集成套件 + 1 个单测 job 并行下载运行。18min → ~5min，但 yml 改动大，建议等 PR-1~3 合入后再提。

### D2. 01_PW GPU 套件瘦身（cuda.yml）

**问题**：#7690 加的 73 个用例 `-n 1` 串行，+4.8min/次，且用例是为 CPU 写的。

**改法**：
1. `-n 1` → `-n 2`（与其他 GPU 套件一致）；
2. `CASES_GPU.txt` 拆两份：`CASES_GPU_SMOKE.txt`（~20 个代表用例，PR 跑）+ 全量留给定时任务：

```yaml
on:
  pull_request:
  schedule:
    - cron: '0 18 * * *'   # 每晚全量（UTC）
  workflow_dispatch:

jobs:
  gpu-test:
    steps:
      - name: Select case list
        run: |
          if [ "${{ github.event_name }}" = "pull_request" ]; then
            echo "CASES=CASES_GPU_SMOKE.txt" >> "$GITHUB_ENV"
          else
            echo "CASES=CASES_GPU.txt" >> "$GITHUB_ENV"
          fi
      - name: Test 01_PW on GPU
        run: |
          cd tests/01_PW && bash ../integrate/Autotest.sh -n 2 -a abacus -f "$CASES"
```

---

## 落地路线图（给导师汇报版）

| 批次 | 内容 | 文件 | 预期收益 | 过审难度 |
|---|---|---|---|---|
| **PR-1** | B1 架构裁剪 + B2 nproc + D3 浅克隆 | cuda.yml、test.yml | CUDA Build 33→8~15min；冷节点 checkout 省 12min | ★ 极易（几行配置） |
| **PR-2** | C1 ccache 容量+统计；C3 cache-warmer | 2 个 yml + 1 个新 yml | 消除 2.7 vs 37.5min 波动；风暴后自动回暖 | ★★ 易（纯新增步骤） |
| **PR-3** | A1 Dockerfile 删 cache-bust+镜像内编译（+E stage4 预装 +A2 国内源） | 2 个 Dockerfile、devcontainer.yml、2 个 yml | init 8~78→~1min；toolchain 归零 | ★★★ 中（改镜像定位，需共识） |
| **PR-4** | D1 单测合并并行 + D2 GPU 套件瘦身（之后视情况上完整 matrix） | test.yml、cuda.yml、CASES_GPU | CPU 测试 18→~6min；GPU 测试省 ~3min | ★★ 易~中 |

**验证方法（每个 PR 必做，写进 PR 描述，符合 AGENTS.md「报告确切验证」要求）**：
1. 从 fork 提交后观察本 PR 的 CI 运行时长对比（注意：fork PR 改 workflow 文件首次需维护者 approve，先发 PR 后 @ 维护者）；
2. 贴 `ccache -s` 的 hit rate 前后对比；
3. 贴 Initialize containers / Build 步骤耗时的前后截图；
4. 确认 GPU 测试全部通过（架构裁剪后二进制必须能在 CI 卡上跑——这也是 nvidia-smi 要先看的原因）。

**风险提示（写进 PR）**：
- B1 若填错架构 → GPU 测试全挂（立即可见，回退一行即可）；
- B2 不做 B1 直接加并行 → nvcc 内存 ×7 架构可能 OOM；
- A1 影响"用镜像体验 abacus"的用户路径，PR 里给维护者两个选项；
- C3 cache-warmer 占用 runner 时间，建议 `concurrency` 防堆积。

---

## 附：为什么这些方案能成立（事实锚点）

| 方案依赖的事实 | 出处 |
|---|---|
| ccache 自动接入（含 CUDA） | CMakeLists.txt:299-303 |
| 构建只需 git log -1 | CMakeLists.txt:163-174 |
| 7 架构默认列表 | CMakeLists.txt:445-476 |
| -j4 写死 | 上游 cuda.yml `cmake --build build -j4` |
| 镜像每次合并重建 + cache-bust + 镜像内编译 | devcontainer.yml；Dockerfile.gnu:31-41；Dockerfile.cuda:27-36 |
| 国内镜像源已存在 | Dockerfile.gnu:4 注释；devcontainer.yml:39 |
| 01_PW GPU -n1 串行 73 例 | #7690（4b6ce1c63） |
| init 8~78min、Build 2.7 vs 37.5min 实测 | GitHub API runs/jobs 数据（见分析文档） |
