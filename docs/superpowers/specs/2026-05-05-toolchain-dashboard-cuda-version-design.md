# Toolchain 启动看板输出 CUDA 版本

## 背景

toolchain 的启动看板（`ui_show_summary()`）目前会探测并输出服务器 GPU 型号，但不输出 CUDA 相关信息。

## 目标

- 在启动看板的 “System Information” 区域中，增加 CUDA 相关信息输出
- 不引入对 `nvcc` 的依赖
- 在无 NVIDIA 环境下保持可运行，并明确显示不可用状态

## 方案

- 复用 `nvidia-smi`：
  - `NVIDIA Driver`：使用 `nvidia-smi --query-gpu=driver_version --format=csv,noheader`，取第一行
  - `CUDA Version`：解析 `nvidia-smi` banner 行中的 `CUDA Version: X.Y`
- 当 `nvidia-smi` 不存在或解析失败时：
  - `NVIDIA Driver` 与 `CUDA Version` 均输出 `unavailable`

## 输出格式

在现有 `GPU:` 行附近增加两行：

- `NVIDIA Driver: ...`
- `CUDA Version: ...`

并保持树形前缀（`├─` / `└─`）的连贯性。

