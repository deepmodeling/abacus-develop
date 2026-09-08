# DFT-D3 data generator

`generate_d3_data.py` converts immutable numerical data from the pinned
s-dftd3 v1.5.0 source revision into compact C++ include files used by ABACUS.
It is a maintainer tool only: s-dftd3 and mctc-lib are not ABACUS build-time or
run-time dependencies.

The generator verifies every input checksum before parsing it. To regenerate
or check the committed files, provide the s-dftd3 source tree and the mctc-lib
source tree used by that release:

```sh
python3 tools/dftd3/generate_d3_data.py \
  --s-dftd3-root /path/to/simple-dftd3-1.5.0 \
  --mctc-lib-root /path/to/mctc-lib \
  --output-root .

python3 tools/dftd3/generate_d3_data.py \
  --s-dftd3-root /path/to/simple-dftd3-1.5.0 \
  --mctc-lib-root /path/to/mctc-lib \
  --output-root . --check
```

Updating the numerical specification is an explicit review event: update the
pinned revision, input checksums, generated files, and s-dftd3 oracle tests in
the same change.
