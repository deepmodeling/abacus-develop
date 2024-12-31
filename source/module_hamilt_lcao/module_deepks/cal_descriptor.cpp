/// 1. cal_descriptor : obtains descriptors which are eigenvalues of pdm
///       by calling torch::linalg::eigh
/// 2. check_descriptor : prints descriptor for checking

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

void LCAO_Deepks::cal_descriptor_equiv(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_descriptor_equiv");
    ModuleBase::timer::tick("LCAO_Deepks", "cal_descriptor_equiv");

    // a rather unnecessary way of writing this, but I'll do it for now
    if (!this->d_tensor.empty())
    {
        this->d_tensor.erase(this->d_tensor.begin(), this->d_tensor.end());
    }

    for (int iat = 0; iat < nat; iat++)
    {
        auto tmp = torch::zeros(des_per_atom, torch::kFloat64);
        std::memcpy(tmp.data_ptr(), pdm[iat].data_ptr<double>(), sizeof(double) * tmp.numel());
        this->d_tensor.push_back(tmp);
    }

    ModuleBase::timer::tick("LCAO_Deepks", "cal_descriptor_equiv");
}

// calculates descriptors from projected density matrices
void LCAO_Deepks::cal_descriptor(const int nat)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_descriptor");
    ModuleBase::timer::tick("LCAO_Deepks", "cal_descriptor");

    if (PARAM.inp.deepks_equiv)
    {
        this->cal_descriptor_equiv(nat);
        return;
    }

    // init d_tensor
    // if d_tensor is not empty, clear it !!
    if (!this->d_tensor.empty())
    {
        this->d_tensor.erase(this->d_tensor.begin(), this->d_tensor.end());
    }

    for (int inl = 0; inl < this->inlmax; ++inl)
    {
        const int nm = 2 * inl_l[inl] + 1;
        this->pdm[inl].requires_grad_(true);
        this->d_tensor.push_back(torch::ones({nm}, torch::requires_grad(true)));
    }

    // cal d_tensor
    for (int inl = 0; inl < inlmax; ++inl)
    {
        torch::Tensor vd;
        std::tuple<torch::Tensor, torch::Tensor> d_v(this->d_tensor[inl], vd);
        // d_v = torch::symeig(pdm[inl], /*eigenvalues=*/true,
        // /*upper=*/true);
        d_v = torch::linalg::eigh(pdm[inl], /*uplo*/ "U");
        d_tensor[inl] = std::get<0>(d_v);
    }
    ModuleBase::timer::tick("LCAO_Deepks", "cal_descriptor");
    return;
}

void LCAO_Deepks::check_descriptor(const UnitCell& ucell, const std::string& out_dir)
{
    ModuleBase::TITLE("LCAO_Deepks", "check_descriptor");

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
                ofs << ucell.atoms[it].label << " atom_index " << ia + 1 << " n_descriptor " << this->des_per_atom
                    << std::endl;
                int id = 0;
                for (int inl = 0; inl < inlmax / ucell.nat; inl++)
                {
                    int nm = 2 * inl_l[inl] + 1;
                    for (int im = 0; im < nm; im++)
                    {
                        const int ind = iat * inlmax / ucell.nat + inl;
                        ofs << d_tensor[ind].index({im}).item().toDouble() << " ";

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
            ofs << ucell.atoms[it].label << " atom_index " << iat + 1 << " n_descriptor " << this->des_per_atom
                << std::endl;
            for (int i = 0; i < this->des_per_atom; i++)
            {
                ofs << this->pdm[iat][i].item<double>() << " ";
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
