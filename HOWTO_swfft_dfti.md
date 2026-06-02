# CPE-DFTI plane-wave sticks FFT (Sunway swFFT) — factory backend

Branch `feat/swfft-dfti` (off `develop`), single commit. Offloads `FFT_CPU`'s
local 1D sticks FFTs to the Sunway CPEs via the swFFT (xMath-SACA) DFTI API,
packaged as a **separate FFT backend selected through the `FFT_Bundle` factory**
— `FFT_CPU` itself carries no DFTI `#ifdef`. No box/2DECOMP, no multi-process.

Compiled only on Sunway (`USE_SWDFTI`, default ON under `USE_SW`). On any other
build the backend file isn't compiled and the factory uses plain `FFT_CPU`, so
the result is **byte-identical to develop** (verified: `fft_cpu.cpp`/`fft_bundle.cpp`
compile clean at `USE_SW=OFF`; `fft_swdfti` excluded from baseline; CMake clean).

## What changed
- **`module_fft/fft_swdfti.{h,cpp}`** (new): `FFT_SWDFTI<double> : public FFT_CPU<double>`
  overrides `fftzfor/fftzbac` (batched 1D-z on CPE) and `fftxyfor/fftxybac`
  (strided 1D-x on CPE; y on FFTW — CPE loses on the small per-slice batch), and
  `setupFFT` (builds the DFTI descriptors after the base FFTW plans). Non-xprime /
  disabled cases delegate to `FFT_CPU`. `DftiInitAthread` spawns the CPEs.
- **`fft_bundle.cpp`**: factory — device `"cpu"` (double) builds `FFT_SWDFTI` when
  `__SWDFTI`, else `FFT_CPU`. The only backend-selection point.
- **`fft_cpu.h`**: members `private` → `protected` (subclass reuses plans/dims).
- **`CMakeLists.txt`**: `USE_SWDFTI` option → `__SWDFTI`; `-mieee` (IEEE FP, via
  `CheckCXXCompilerFlag`) under `USE_SW`; link the objcopy-isolated
  `libswfft_xmath_iso.a` in place of raw `${SW_MATH}/libswfft.a`.
- **`module_pw/CMakeLists.txt`**: compiles `fft_swdfti.cpp` under `USE_SWDFTI`.

## Performance (4GaAs, ecut60, 54³)
sticks + CPE-DFTI vs baseline FFTW sticks: `veff_pw` 1.7–1.8×, scales with np
(np 2/4/6 total ≈ 55/32/24 s). Energy bit-identical.

## Build on the Sunway machine
1. **Isolate the xMath swfft symbols** (mandatory — bundled `fftw_*` would hijack
   ABACUS's FFTW and break the density FFT → E_Hartree=0 / −4079 eV):
   ```bash
   SW_MATH=/usr/sw/yyzlib/xMath-SACA
   swnm $SW_MATH/libswfft.a | grep -E ' [TDBW] (d?fftw_|fftwf_)' \
     | awk '{print $3, "swfftpriv_"$3}' > fftw_iso.map
   swobjcopy --redefine-syms=fftw_iso.map \
     $SW_MATH/libswfft.a source/source_base/module_fft/libswfft_xmath_iso.a
   ```
   Use `swobjcopy`/`swnm` (sw_64), NOT the host tools.
2. `cmake -S . -B build -DUSE_SW=ON && cmake --build build -j`
   (`-DUSE_SWDFTI=OFF` → plain FFTW sticks.)
3. Run via direct bsub:
   `bsub -b -q q_swhnu -n <np> -cgsp 64 -share_size 4096 -host_stack 128 -o run.out ./abacus`

## Runtime toggle
- `ABACUS_NO_DFTI=1` — disable DFTI at runtime (FFT_SWDFTI falls back to FFTW), for A/B.

## Parked WIP
Earlier `cpu3d-experiment` uncommitted work is in `git stash@{0}`:
`git checkout cpu3d-experiment && git stash pop`.
