# Extracting Wave Functions

ABACUS is able to output electron wave functions in both PW and LCAO basis calculations. One can find the examples in [examples/wfc](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/wfc).

## wave function in G space
For the wave function in G space, one only needs to do a ground-state energy calculation with one additional keyword in the INPUT file: '***[out_wfc_pw](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-pw)***' for PW basis calculation, and '***[out_wfc_lcao](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out-wfc-lcao)***' for LCAO basis calculation.
In the PW basis case, the wave function is output in a file called `WAVEFUNC${k}.txt`, where `${k}` is the index of K point. \
In the LCAO basis case, several `LOWF_K_${k}.dat` files will be output in multi-k calculation and `LOWF_GAMMA_S1.dat` in gamma-only calculation. 

## wave function in real space

One can also choose to output real-space wave function in PW basis calculation with the key word ***[out_wfc_r](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#out_wfc_r)***.

After calculation, an additional directory named `wfc_realspace` will appear in the `OUT.${system}` directory.

Notice: when the ***[basis_type](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#basis_type)*** is `lcao`, only `get_wf` ***[calculation](https://abacus-rtd.readthedocs.io/en/latest/advanced/input_files/input-main.html#calculation)*** is effective. An example is [examples/wfc/lcao_ienvelope_Si2](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/wfc/lcao_ienvelope_Si2). 