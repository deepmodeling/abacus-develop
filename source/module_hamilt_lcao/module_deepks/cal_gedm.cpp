/// cal_gedm : calculates d(E_delta)/d(pdm)
///       this is the term V(D) that enters the expression H_V_delta = |alpha>V(D)<alpha|
///       caculated using torch::autograd::grad
/// check_gedm : prints gedm for checking

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

inline void generate_py_files(const int lmaxd, const int nmaxd, const std::string& out_dir)
{

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    std::ofstream ofs("cal_gedm.py");
    ofs << "import torch" << std::endl;
    ofs << "import numpy as np" << std::endl << std::endl;
    ofs << "import sys" << std::endl;

    ofs << "from deepks.scf.enn.scf import BasisInfo" << std::endl;
    ofs << "from deepks.iterate.template_abacus import t_make_pdm" << std::endl;
    ofs << "from deepks.utils import load_yaml" << std::endl << std::endl;

    ofs << "basis = load_yaml('basis.yaml')['proj_basis']" << std::endl;
    ofs << "model = torch.jit.load(sys.argv[1])" << std::endl;
    ofs << "dm_eig = np.expand_dims(np.load('" << out_dir << "dm_eig.npy'),0)" << std::endl;
    ofs << "dm_eig = torch.tensor(dm_eig, "
           "dtype=torch.float64,requires_grad=True)"
        << std::endl
        << std::endl;

    ofs << "dm_flat,basis_info = t_make_pdm(dm_eig,basis)" << std::endl;
    ofs << "ec = model(dm_flat.double())" << std::endl;
    ofs << "gedm = "
           "torch.autograd.grad(ec,dm_eig,grad_outputs=torch.ones_like(ec))[0]"
        << std::endl
        << std::endl;

    ofs << "np.save('ec.npy',ec.double().detach().numpy())" << std::endl;
    ofs << "np.save('gedm.npy',gedm.double().numpy())" << std::endl;
    ofs.close();

    ofs.open("basis.yaml");
    ofs << "proj_basis:" << std::endl;
    for (int l = 0; l < lmaxd + 1; l++)
    {
        ofs << "  - - " << l << std::endl;
        ofs << "    - [";
        for (int i = 0; i < nmaxd + 1; i++)
        {
            ofs << "0";
            if (i != nmaxd)
            {
                ofs << ", ";
            }
        }
        ofs << "]" << std::endl;
    }
}

void LCAO_Deepks::cal_gedm_equiv(const int nat, const std::vector<torch::Tensor>& descriptor)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_gedm_equiv");

    LCAO_deepks_io::save_npy_d(nat,
                               this->des_per_atom,
                               this->inlmax,
                               this->inl_l,
                               PARAM.inp.deepks_equiv,
                               descriptor,
                               PARAM.globalv.global_out_dir,
                               GlobalV::MY_RANK); // libnpy needed

    generate_py_files(this->lmaxd, this->nmaxd, PARAM.globalv.global_out_dir);

    if (GlobalV::MY_RANK == 0)
    {
        std::string cmd = "python cal_gedm.py " + PARAM.inp.deepks_model;
        int stat = std::system(cmd.c_str());
        assert(stat == 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    LCAO_deepks_io::load_npy_gedm(nat, this->des_per_atom, this->gedm, this->E_delta, GlobalV::MY_RANK);

    std::string cmd = "rm -f cal_gedm.py basis.yaml ec.npy gedm.npy";
    std::system(cmd.c_str());
}

// obtain from the machine learning model dE_delta/dDescriptor
void LCAO_Deepks::cal_gedm(const int nat, const std::vector<torch::Tensor>& descriptor)
{

    if (PARAM.inp.deepks_equiv)
    {
        this->cal_gedm_equiv(nat, descriptor);
        return;
    }

    ModuleBase::TITLE("LCAO_Deepks", "cal_gedm");

    // forward
    std::vector<torch::jit::IValue> inputs;

    // input_dim:(natom, des_per_atom)
    inputs.push_back(torch::cat(descriptor, 0).reshape({1, nat, this->des_per_atom}));
    std::vector<torch::Tensor> ec;
    ec.push_back(module.forward(inputs).toTensor()); // Hartree
    this->E_delta = ec[0].item<double>() * 2;        // Ry; *2 is for Hartree to Ry

    // cal gedm
    std::vector<torch::Tensor> gedm_shell;
    gedm_shell.push_back(torch::ones_like(ec[0]));
    this->gedm_tensor = torch::autograd::grad(ec,
                                              this->pdm,
                                              gedm_shell,
                                              /*retain_grad=*/true,
                                              /*create_graph=*/false,
                                              /*allow_unused=*/true);

    // gedm_tensor(Hartree) to gedm(Ry)
    for (int inl = 0; inl < inlmax; ++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        auto accessor = this->gedm_tensor[inl].accessor<double, 2>();
        for (int m1 = 0; m1 < nm; ++m1)
        {
            for (int m2 = 0; m2 < nm; ++m2)
            {
                int index = m1 * nm + m2;
                //*2 is for Hartree to Ry
                this->gedm[inl][index] = accessor[m1][m2] * 2;
            }
        }
    }
    return;
}

void LCAO_Deepks::check_gedm()
{
    std::ofstream ofs("gedm.dat");

    for (int inl = 0; inl < inlmax; inl++)
    {
        int nm = 2 * inl_l[inl] + 1;
        for (int m1 = 0; m1 < nm; ++m1)
        {
            for (int m2 = 0; m2 < nm; ++m2)
            {
                int index = m1 * nm + m2;
                //*2 is for Hartree to Ry
                ofs << this->gedm[inl][index] << " ";
            }
        }
        ofs << std::endl;
    }
}

#endif
