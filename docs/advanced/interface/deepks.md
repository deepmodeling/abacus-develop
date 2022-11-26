# DeePKS

[DeePKS](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00872) is a machine-learning aided density funcitonal model that fits the energy difference between highly accurate but computationally demanding method and effcient but less accurate method via neural-network. As such, the trained DeePKS model can provide highly accurate energetics (and forces) with relatively low computational cost, and can therefore act as a bridge to connect expensive quantum mechanic data and machine-learning-based potentials. While the original framework of DeePKS is for molecular systems, please refer to this [reference](https://arxiv.org/abs/2206.10093) for the application of DeePKS in periodic systems. 

Detailed instructions on installing and running DeePKS can be found on this [website](https://deepks-kit.readthedocs.io/en/latest/index.html). An [example](https://github.com/deepmodeling/deepks-kit/tree/abacus/examples/water_single_lda2pbe_abacus) for training DeePKS model with ABACUS is also provided. The DeePKS-related keywords in `INPUT` file can be found [here](http://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#deepks).

> Note: Use the LCAO basis for DeePKS-related calculations


