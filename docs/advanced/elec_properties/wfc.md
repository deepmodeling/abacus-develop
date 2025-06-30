# Extracting Wave Functions

ABACUS is able to output electron wave functions in both PW and LCAO basis calculations. One can find the examples in [examples/11_wfc](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/11_wfc).

## Wave Function in G-Space

To output wave functions in G-space, add one of the following keywords to the `INPUT` file while performing SCF calculation:
- **PW basis**: Set [`out_wfc_pw`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-pw) to `1`. Output file format: `wfs[spin]k[kpoint]_pw.txt`, where `[spin]` is the spin channel index, and `[kpoint]` the k-point index.

- **LCAO basis**: Set [`out_wfc_lcao`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-lcao) to `1`.  
  - **Multi-k calculations**: Generates multiple files `wfs[spin]k[kpoint]_nao.txt`.
  - **Gamma-only calculations**: `wfs[spin]_nao.txt` instead.

## Wave Function in Real Space

One can also choose to output real-space wave functions with the keyword [`out_wfc_norm`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-norm) or [`out_wfc_re_im`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-re-im).

Notice: When the [`basis_type`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#basis-type) is `lcao`, only `get_wf` [`calculation`](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#calculation) is effective. An example is [examples/11_wfc/lcao_ienvelope_Si2](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/11_wfc/lcao_ienvelope_Si2).