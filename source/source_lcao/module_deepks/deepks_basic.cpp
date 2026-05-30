// This file contains interfaces with libtorch,
// including loading of model and calculating gradients
// as well as subroutines that prints the results for checking

#ifdef __MLALGO
#include "deepks_basic.h"

#include "source_base/atom_in.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#include <cmath>

#ifdef __MPI
#include <mpi.h>
#endif

// d(Descriptor) / d(projected density matrix)
// Dimension is different for each inl, so there's a vector of tensors
void DeePKS_domain::cal_gevdm(const int nat,
                              const DeePKS_Param& deepks_param,
                              const std::vector<torch::Tensor>& pdm,
                              std::vector<torch::Tensor>& gevdm)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_gevdm");
    ModuleBase::timer::start("DeePKS_domain", "cal_gevdm");
    // cal gevdm(d(EigenValue(D))/dD)
    int nlmax = deepks_param.inlmax / nat;
    for (int nl = 0; nl < nlmax; ++nl)
    {
        std::vector<torch::Tensor> avmmv;
        for (int iat = 0; iat < nat; ++iat)
        {
            int inl = iat * nlmax + nl;
            int nm = 2 * deepks_param.inl2l[inl] + 1;
            // repeat each block for nm times in an additional dimension
            torch::Tensor tmp_x = pdm[inl].reshape({nm, nm}).unsqueeze(0).repeat({nm, 1, 1});
            // torch::Tensor tmp_y = std::get<0>(torch::symeig(tmp_x, true));
            torch::Tensor tmp_y = std::get<0>(torch::linalg_eigh(tmp_x, "U"));
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
    ModuleBase::timer::end("DeePKS_domain", "cal_gevdm");
    return;
}

void DeePKS_domain::load_model(const std::string& model_file, torch::jit::script::Module& model)
{
    ModuleBase::TITLE("DeePKS_domain", "load_model");
    ModuleBase::timer::start("DeePKS_domain", "load_model");

    // check whether file exists
    std::ifstream ifs(model_file.c_str());
    if (!ifs)
    {
        ModuleBase::timer::end("DeePKS_domain", "load_model");
        ModuleBase::WARNING_QUIT("DeePKS_domain::load_model", "No model file named " + model_file + ", please check!");
        return;
    }
    ifs.close();
    try
    {
        model = torch::jit::load(model_file);
    }
    catch (const c10::Error& e)
    {
        std::cerr << "error loading the model" << std::endl;
        ModuleBase::timer::end("DeePKS_domain", "load_model");
        return;
    }
    ModuleBase::timer::end("DeePKS_domain", "load_model");
    return;
}

void DeePKS_domain::cal_edelta_gedm_equiv(const int nat,
                                          const DeePKS_Param& deepks_param,
                                          const std::vector<torch::Tensor>& descriptor,
                                          torch::jit::script::Module& model_deepks,
                                          double** gedm,
                                          double& E_delta,
                                          const int rank)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_edelta_gedm_equiv");
    ModuleBase::timer::start("DeePKS_domain", "cal_edelta_gedm_equiv");

    if (rank == 0)
    {
        const int basis_size
            = static_cast<int>(std::llround(std::sqrt(static_cast<double>(deepks_param.des_per_atom))));
        if (basis_size * basis_size != deepks_param.des_per_atom)
        {
            ModuleBase::WARNING_QUIT("DeePKS_domain::cal_edelta_gedm_equiv",
                                     "Invalid des_per_atom for equivariant DeePKS: it must be a perfect square.");
        }

        torch::Tensor dm_eig = torch::cat(descriptor, 0).reshape({1, nat, deepks_param.des_per_atom});
        dm_eig = dm_eig.to(torch::kFloat64).requires_grad_(true);
        torch::Tensor dm = dm_eig.reshape({1, nat, basis_size, basis_size});

        if (static_cast<int>(deepks_param.nchi_d_l.size()) != deepks_param.lmaxd + 1)
        {
            ModuleBase::WARNING_QUIT(
                "DeePKS_domain::cal_edelta_gedm_equiv",
                "Invalid nchi_d_l in DeePKS parameters: expected size lmaxd + 1 for equivariant shell construction.");
        }

        std::vector<torch::Tensor> ovlp_shells;
        int total_shells = 0;
        for (int l = 0; l <= deepks_param.lmaxd; ++l)
        {
            total_shells += deepks_param.nchi_d_l[l];
        }
        ovlp_shells.reserve(total_shells);
        int offset = 0;
        for (int l = 0; l <= deepks_param.lmaxd; ++l)
        {
            const int nm = 2 * l + 1;
            for (int n = 0; n < deepks_param.nchi_d_l[l]; ++n)
            {
                torch::Tensor po = torch::zeros({basis_size, 1, nm}, torch::TensorOptions().dtype(torch::kFloat64));
                auto accessor = po.accessor<double, 3>();
                for (int m = 0; m < nm; ++m)
                {
                    accessor[offset + m][0][m] = 1.0;
                }
                ovlp_shells.push_back(po);
                offset += nm;
            }
        }
        if (offset != basis_size)
        {
            ModuleBase::WARNING_QUIT("DeePKS_domain::cal_edelta_gedm_equiv",
                                     "Invalid shell layout: accumulated shell offset does not match basis size.");
        }

        std::vector<torch::Tensor> dm_flat;
        dm_flat.reserve(ovlp_shells.size());
        for (const auto& po : ovlp_shells)
        {
            // Equivalent to python:
            // torch.einsum('rap,...rs,saq->...apq', po, dm, po)
            torch::Tensor pdm_shell = torch::einsum("rap,...rs,saq->...apq", {po, dm, po});
            dm_flat.push_back(pdm_shell.squeeze(-3));
        }

        c10::List<torch::Tensor> model_input;
        for (const auto& pdm_shell : dm_flat)
        {
            model_input.push_back(pdm_shell);
        }

        std::vector<torch::jit::IValue> inputs;
        inputs.emplace_back(model_input);

        torch::Tensor ec;
        try
        {
            ec = model_deepks.forward(inputs).toTensor(); // Hartree
        }
        catch (const c10::Error& e)
        {
            ModuleBase::WARNING_QUIT("DeePKS_domain::cal_edelta_gedm_equiv",
                                     "Failed to evaluate equivariant DeePKS model in C++.");
            throw;
        }

        E_delta = ec.item<double>() * 2.0; // Hartree to Ry

        std::vector<torch::Tensor> grad_outputs{torch::ones_like(ec)};
        std::vector<torch::Tensor> grad_inputs{dm_eig};
        torch::Tensor gedm_tensor = torch::autograd::grad({ec}, grad_inputs, grad_outputs,
                                                           /*retain_graph=*/false,
                                                           /*create_graph=*/false,
                                                           /*allow_unused=*/false)[0];

        torch::Tensor gedm_nat = gedm_tensor.reshape({nat, deepks_param.des_per_atom});
        auto accessor = gedm_nat.accessor<double, 2>();
        for (int iat = 0; iat < nat; ++iat)
        {
            for (int ides = 0; ides < deepks_param.des_per_atom; ++ides)
            {
                gedm[iat][ides] = accessor[iat][ides] * 2.0; // Hartree to Ry
            }
        }
    }

#ifdef __MPI
    for (int iat = 0; iat < nat; ++iat)
    {
        MPI_Bcast(gedm[iat], deepks_param.des_per_atom, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&E_delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    ModuleBase::timer::end("DeePKS_domain", "cal_edelta_gedm_equiv");
    return;
}

// obtain from the machine learning model dE_delta/dDescriptor
// E_delta is also calculated here
void DeePKS_domain::cal_edelta_gedm(const int nat,
                                    const DeePKS_Param& deepks_param,
                                    const std::vector<torch::Tensor>& descriptor,
                                    const std::vector<torch::Tensor>& pdm,
                                    torch::jit::script::Module& model_deepks,
                                    double** gedm,
                                    double& E_delta)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_edelta_gedm");
    ModuleBase::timer::start("DeePKS_domain", "cal_edelta_gedm");

    // forward
    std::vector<torch::jit::IValue> inputs;

    // input_dim:(natom, des_per_atom)
    inputs.push_back(torch::cat(descriptor, 0).reshape({1, nat, deepks_param.des_per_atom}));
    std::vector<torch::Tensor> ec;
    try
    {
        ec.push_back(model_deepks.forward(inputs).toTensor()); // Hartree
    }
    catch (const c10::Error& e)
    {
        ModuleBase::WARNING_QUIT("DeePKS_domain::cal_edelta_gedm",
                                 "Please check whether the input shape required by model file matches the descriptor!");
        throw;
    }
    E_delta = ec[0].item<double>() * 2; // Ry; *2 is for Hartree to Ry

    // ec: [1, 1]
    ec[0].reshape({1, 1}).requires_grad_(true);

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
    for (int inl = 0; inl < deepks_param.inlmax; ++inl)
    {
        int nm = 2 * deepks_param.inl2l[inl] + 1;
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
    ModuleBase::timer::end("DeePKS_domain", "cal_edelta_gedm");
    return;
}

void DeePKS_domain::check_gedm(const DeePKS_Param& deepks_param, double** gedm)
{
    std::ofstream ofs("gedm.dat");

    for (int inl = 0; inl < deepks_param.inlmax; inl++)
    {
        int nm = 2 * deepks_param.inl2l[inl] + 1;
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

void DeePKS_domain::prepare_atom(const UnitCell& ucell, torch::Tensor& atom_out)
{
    int nat = ucell.nat;
    atom_out = torch::zeros({nat, 4}, torch::TensorOptions().dtype(torch::kFloat64));

    // get atom information
    atom_in AtomInfo;
    int index = 0;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            atom_out[index][0] = AtomInfo.atom_Z[ucell.atom_label[it]];

            // use bohr as unit
            atom_out[index][1] = ucell.atoms[it].tau[ia].x * ucell.lat0;
            atom_out[index][2] = ucell.atoms[it].tau[ia].y * ucell.lat0;
            atom_out[index][3] = ucell.atoms[it].tau[ia].z * ucell.lat0;
            index++;
        }
    }
}
void DeePKS_domain::prepare_box(const UnitCell& ucell, torch::Tensor& box_out)
{
    box_out = torch::zeros({9}, torch::TensorOptions().dtype(torch::kFloat64));

    // use bohr as unit
    box_out[0] = ucell.latvec.e11 * ucell.lat0;
    box_out[1] = ucell.latvec.e12 * ucell.lat0;
    box_out[2] = ucell.latvec.e13 * ucell.lat0;
    box_out[3] = ucell.latvec.e21 * ucell.lat0;
    box_out[4] = ucell.latvec.e22 * ucell.lat0;
    box_out[5] = ucell.latvec.e23 * ucell.lat0;
    box_out[6] = ucell.latvec.e31 * ucell.lat0;
    box_out[7] = ucell.latvec.e32 * ucell.lat0;
    box_out[8] = ucell.latvec.e33 * ucell.lat0;
}

#endif
