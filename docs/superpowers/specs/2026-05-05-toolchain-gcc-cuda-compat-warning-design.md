# Toolchain 启动检查：GCC 与 CUDA 版本兼容性风险提示

## 背景

当系统 GCC 版本较新（例如 14+）而 CUDA 版本较旧（例如 < 13）时，可能存在编译器与 CUDA 工具链/驱动之间的兼容性风险，需要在启动检查阶段提示用户。

## 目标

- 在 toolchain 的启动检查（System Validation）中加入风险提示
- 逻辑与现有 GCC 版本检查保持一致
- 仅在 `with_gcc != __INSTALL__`（使用系统 GCC/自定义 GCC 路径）时触发该检查

## 触发条件

- `gcc_major >= 14`
- 且可从 `nvidia-smi` 解析到 `CUDA Version`，并满足 `cuda_major < 13`

## 实现策略

- 放置在 `config_validator.sh` 的 `validate_system_requirements()` 中，复用已获取的 `gcc_version/gcc_major`
- CUDA 版本解析来源为 `nvidia-smi` banner 行的 `CUDA Version: X.Y`
- 当无法获取/解析 CUDA major 时不提示，避免误报

## 输出

- 以 warning 形式输出，提示 “GCC>=14 且 CUDA<13 可能存在不兼容风险”

