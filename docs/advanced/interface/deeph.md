# DeepH

[DeepH](https://doi.org/10.1038/s43588-022-00265-6) applies meaching learning to predict the Hamiltonian in atomic basis representation. For such purpose, DeepH uses the Hamiltonian and overlap matrices from DFT calculations. Here we introduce how to extract relevant information from ABACUS for the purpose of DeepH training and prediction.

Detailed instructions on installing and running DeepH can be found on its official [website](https://deeph-pack.deepmodeling.com/en/latest/#deeph). An [example](https://deeph-pack.deepmodeling.com/en/latest/demo/demo3.html) for using DeepH with ABACUS is also provided.

Here I intend not to repeat information from the above sources, but to add some minor details related to the setting of ABACUS `INPUT` files.

> Note: Use the LCAO basis for DeepH-related calculations

As mentioned in the README.md file in the above-mentioned example, there are two stages where users need to run ABACUS calculations.

The first stage is during the data preparation phase, where we need to run a series of SCF calculations and output the Hamiltonian and overlap matrices. For such purpose, one needs to add the following line in the `INPUT` file:

```
out_mat_hs2 1
```

Files named data-HR-sparse_SPIN`${x}`.csr and data-SR-sparse_SPIN`${x}`.csr will be generated, which contain the Hamiltonian and overlap matrices respectively in csr format. `${x}` takes value of 0 or 1, based on the spin component. More details on this keyword can be found in the [list of input keywords](../input_files/input-main.md#out_mat_hs2).

The second stage is during the inference phase. After DeepH training completes, we can apply the model to predict the Hamiltonian on other systems. For that purpose, we also need the overlap matrices from the new systems, but no SCF calculation is required.

For that purpose, in `INPUT` file we need to make the following specification of the keyword `calculation`:

```
calculation get_S
```

A file named `SR.csr` will be generated in the working directory, which contains the overlap matrix.
