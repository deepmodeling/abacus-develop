// This file contains interfaces with libtorch,
// including loading of model and calculating gradients
// as well as subroutines that prints the results for checking

#ifdef __DEEPKS
#include "deepks_basic.h"

#include "module_base/timer.h"
#include "module_parameter/parameter.h"

// d(Descriptor) / d(projected density matrix)
// Dimension is different for each inl, so there's a vector of tensors
void DeePKS_domain::cal_gevdm(const int nat,
                              const int inlmax,
                              const int* inl_l,
                              const std::vector<torch::Tensor>& pdm,
                              std::vector<torch::Tensor>& gevdm)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_gevdm");
    ModuleBase::timer::tick("DeePKS_domain", "cal_gevdm");
    // cal gevdm(d(EigenValue(D))/dD)
    int nlmax = inlmax / nat;
    for (int nl = 0; nl < nlmax; ++nl)
    {
        std::vector<torch::Tensor> avmmv;
        for (int iat = 0; iat < nat; ++iat)
        {
            int inl = iat * nlmax + nl;
            int nm = 2 * inl_l[inl] + 1;
            // repeat each block for nm times in an additional dimension
            torch::Tensor tmp_x = pdm[inl].reshape({nm, nm}).unsqueeze(0).repeat({nm, 1, 1});
            // torch::Tensor tmp_y = std::get<0>(torch::symeig(tmp_x, true));
            torch::Tensor tmp_y = std::get<0>(torch::linalg::eigh(tmp_x, "U"));
            torch::Tensor tmp_yshell = torch::eye(nm, torch::TensorOptions().dtype(torch::kFloat64));
            std::vector<torch::Tensor> tmp_rpt; // repeated-pdm-tensor (x)
            std::vector<torch::Tensor> tmp_rdt; // repeated-d-tensor (y)
            std::vector<torch::Tensor> tmp_gst; // gvx-shell
            tmp_rpt.push_back(tmp_x);
            tmp_rdt.push_back(tmp_y);
            tmp_gst.push_back(tmp_yshell);
            std::vector<torch::Tensor> tmp_res;
            tmp_res = torch::autograd::grad(tmp_rdt,
                                            tmp_rpt,
                                            tmp_gst,
                                            false,
                                            false,
                                            /*allow_unused*/ true); // nm(v)**nm*nm
            avmmv.push_back(tmp_res[0]);
        }
        torch::Tensor avmm = torch::stack(avmmv, 0); // nat*nv**nm*nm
        gevdm.push_back(avmm);
    }
    assert(gevdm.size() == nlmax);
    ModuleBase::timer::tick("DeePKS_domain", "cal_gevdm");
    return;
}

void DeePKS_domain::load_model(const std::string& model_file, torch::jit::script::Module& model)
{
    ModuleBase::TITLE("DeePKS_domain", "load_model");
    ModuleBase::timer::tick("DeePKS_domain", "load_model");

    try
    {
        model = torch::jit::load(model_file);
    }
    catch (const c10::Error& e)
    {
        std::cerr << "error loading the model" << std::endl;
        ModuleBase::timer::tick("DeePKS_domain", "load_model");
        return;
    }
    ModuleBase::timer::tick("DeePKS_domain", "load_model");
    return;
}

inline void generate_py_files(const int lmaxd, const int nmaxd, const std::string& out_dir)
{

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    std::ofstream ofs("cal_edelta_gedm.py");
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

void DeePKS_domain::cal_edelta_gedm_equiv(const int nat,
                                          const int lmaxd,
                                          const int nmaxd,
                                          const int inlmax,
                                          const int des_per_atom,
                                          const int* inl_l,
                                          const std::vector<torch::Tensor>& descriptor,
                                          double** gedm,
                                          double& E_delta)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_edelta_gedm_equiv");
    ModuleBase::timer::tick("DeePKS_domain", "cal_edelta_gedm_equiv");

    LCAO_deepks_io::save_npy_d(nat,
                               des_per_atom,
                               inlmax,
                               inl_l,
                               PARAM.inp.deepks_equiv,
                               descriptor,
                               PARAM.globalv.global_out_dir,
                               GlobalV::MY_RANK); // libnpy needed

    generate_py_files(lmaxd, nmaxd, PARAM.globalv.global_out_dir);

    if (GlobalV::MY_RANK == 0)
    {
        std::string cmd = "python cal_edelta_gedm.py " + PARAM.inp.deepks_model;
        int stat = std::system(cmd.c_str());
        assert(stat == 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    LCAO_deepks_io::load_npy_gedm(nat, des_per_atom, gedm, E_delta, GlobalV::MY_RANK);

    std::string cmd = "rm -f cal_edelta_gedm.py basis.yaml ec.npy gedm.npy";
    std::system(cmd.c_str());

    ModuleBase::timer::tick("DeePKS_domain", "cal_edelta_gedm_equiv");
    return;
}

// obtain from the machine learning model dE_delta/dDescriptor
// E_delta is also calculated here
void DeePKS_domain::cal_edelta_gedm(const int nat,
                                    const int lmaxd,
                                    const int nmaxd,
                                    const int inlmax,
                                    const int des_per_atom,
                                    const int* inl_l,
                                    const std::vector<torch::Tensor>& descriptor,
                                    const std::vector<torch::Tensor>& pdm,
                                    torch::jit::script::Module& model_deepks,
                                    double** gedm,
                                    double& E_delta)
{
    if (PARAM.inp.deepks_equiv)
    {
        DeePKS_domain::cal_edelta_gedm_equiv(nat, lmaxd, nmaxd, inlmax, des_per_atom, inl_l, descriptor, gedm, E_delta);
        return;
    }
    ModuleBase::TITLE("DeePKS_domain", "cal_edelta_gedm");
    ModuleBase::timer::tick("DeePKS_domain", "cal_edelta_gedm");

    // forward
    std::vector<torch::jit::IValue> inputs;

    // input_dim:(natom, des_per_atom)
    inputs.push_back(torch::cat(descriptor, 0).reshape({1, nat, des_per_atom}));
    std::vector<torch::Tensor> ec;
    ec.push_back(model_deepks.forward(inputs).toTensor()); // Hartree
    E_delta = ec[0].item<double>() * 2;                    // Ry; *2 is for Hartree to Ry

    // cal gedm
    std::vector<torch::Tensor> gedm_shell;
    gedm_shell.push_back(torch::ones_like(ec[0]));
    std::vector<torch::Tensor> gedm_tensor = torch::autograd::grad(ec,
                                                                   pdm,
                                                                   gedm_shell,
                                                                   /*retain_grad=*/true,
                                                                   /*create_graph=*/false,
                                                                   /*allow_unused=*/true);

    // gedm_tensor(Hartree) to gedm(Ry)
    for (int inl = 0; inl < inlmax; ++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        auto accessor = gedm_tensor[inl].accessor<double, 2>();
        for (int m1 = 0; m1 < nm; ++m1)
        {
            for (int m2 = 0; m2 < nm; ++m2)
            {
                int index = m1 * nm + m2;
                gedm[inl][index] = accessor[m1][m2] * 2; //*2 is for Hartree to Ry
            }
        }
    }
    ModuleBase::timer::tick("DeePKS_domain", "cal_edelta_gedm");
    return;
}

void DeePKS_domain::check_gedm(const int inlmax, const int* inl_l, double** gedm)
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
                ofs << gedm[inl][index] << " ";
            }
        }
        ofs << std::endl;
    }
}

#endif
