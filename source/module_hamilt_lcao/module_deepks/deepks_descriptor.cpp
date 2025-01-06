/// 1. cal_descriptor : obtains descriptors which are eigenvalues of pdm
///       by calling torch::linalg::eigh
/// 2. check_descriptor : prints descriptor for checking
/// 3. cal_descriptor_equiv : calculates descriptor in equivalent version

#ifdef __DEEPKS

#include "deepks_descriptor.h"

#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

void DeePKS_domain::cal_descriptor_equiv(const int nat,
                                         const int des_per_atom,
                                         const std::vector<torch::Tensor>& pdm,
                                         std::vector<torch::Tensor>& descriptor)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_descriptor_equiv");
    ModuleBase::timer::tick("DeePKS_domain", "cal_descriptor_equiv");

    assert(des_per_atom > 0);
    for (int iat = 0; iat < nat; iat++)
    {
        auto tmp = torch::zeros(des_per_atom, torch::kFloat64);
        std::memcpy(tmp.data_ptr(), pdm[iat].data_ptr<double>(), sizeof(double) * tmp.numel());
        descriptor.push_back(tmp);
    }

    ModuleBase::timer::tick("DeePKS_domain", "cal_descriptor_equiv");
}

// calculates descriptors from projected density matrices
void DeePKS_domain::cal_descriptor(const int nat,
                                   const int inlmax,
                                   const int* inl_l,
                                   const std::vector<torch::Tensor>& pdm,
                                   std::vector<torch::Tensor>& descriptor,
                                   const int des_per_atom = -1)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_descriptor");
    ModuleBase::timer::tick("DeePKS_domain", "cal_descriptor");

    if (PARAM.inp.deepks_equiv)
    {
        DeePKS_domain::cal_descriptor_equiv(nat, des_per_atom, pdm, descriptor);
        return;
    }

    for (int inl = 0; inl < inlmax; ++inl)
    {
        const int nm = 2 * inl_l[inl] + 1;
        pdm[inl].requires_grad_(true);
        descriptor.push_back(torch::ones({nm}, torch::requires_grad(true)));
    }

    // cal descriptor
    for (int inl = 0; inl < inlmax; ++inl)
    {
        torch::Tensor vd;
        std::tuple<torch::Tensor, torch::Tensor> d_v(descriptor[inl], vd);
        // d_v = torch::symeig(pdm[inl], /*eigenvalues=*/true,
        // /*upper=*/true);
        d_v = torch::linalg::eigh(pdm[inl], /*uplo*/ "U");
        descriptor[inl] = std::get<0>(d_v);
    }
    ModuleBase::timer::tick("DeePKS_domain", "cal_descriptor");
    return;
}

void DeePKS_domain::check_descriptor(const int inlmax,
                                     const int des_per_atom,
                                     const int* inl_l,
                                     const UnitCell& ucell,
                                     const std::string& out_dir,
                                     const std::vector<torch::Tensor>& descriptor)
{
    ModuleBase::TITLE("DeePKS_domain", "check_descriptor");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    // mohan updated 2024-07-25
    std::string file = out_dir + "deepks_desc.dat";

    std::ofstream ofs(file.c_str());
    ofs << std::setprecision(10);
    if (!PARAM.inp.deepks_equiv)
    {
        for (int it = 0; it < ucell.ntype; it++)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                int iat = ucell.itia2iat(it, ia);
                ofs << ucell.atoms[it].label << " atom_index " << ia + 1 << " n_descriptor " << des_per_atom
                    << std::endl;
                int id = 0;
                for (int inl = 0; inl < inlmax / ucell.nat; inl++)
                {
                    int nm = 2 * inl_l[inl] + 1;
                    const int ind = iat * inlmax / ucell.nat + inl;
                    auto accessor = descriptor[ind].accessor<double, 1>();
                    for (int im = 0; im < nm; im++)
                    {
                        ofs << accessor[im] << " ";
                        if (id % 8 == 7)
                        {
                            ofs << std::endl;
                        }
                        id++;
                    }
                }
                ofs << std::endl << std::endl;
            }
        }
    }
    else
    {
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            const int it = ucell.iat2it[iat];
            ofs << ucell.atoms[it].label << " atom_index " << iat + 1 << " n_descriptor " << des_per_atom << std::endl;
            auto accessor = descriptor[iat].accessor<double, 1>();
            for (int i = 0; i < des_per_atom; i++)
            {
                ofs << accessor[i] << " ";
                if (i % 8 == 7)
                {
                    ofs << std::endl;
                }
            }
            ofs << std::endl << std::endl;
        }
    }
    return;
}

#endif
