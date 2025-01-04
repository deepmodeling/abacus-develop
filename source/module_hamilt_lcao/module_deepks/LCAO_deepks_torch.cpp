// This file contains interfaces with libtorch,
// including loading of model and calculating gradients
// as well as subroutines that prints the results for checking

// The file contains 3 subroutines:
//  cal_gvepsl : gvepsl is used for training with stress label, which is derivative of
//        descriptors wrt strain tensor, calculated by
//        d(des)/d\epsilon_{ab} = d(pdm)/d\epsilon_{ab} * d(des)/d(pdm) = gdmepsl * gvdm
//        using einsum
//  cal_gevdm : d(des)/d(pdm)
//        calculated using torch::autograd::grad
//  load_model : loads model for applying V_delta

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"

// calculates stress of descriptors from gradient of projected density matrices
// gv_epsl:d(d)/d\epsilon_{\alpha\beta}, [natom][6][des_per_atom]
void LCAO_Deepks::cal_gvepsl(const int nat,
                             const std::vector<torch::Tensor>& gevdm,
                             const torch::Tensor& gdmepsl,
                             torch::Tensor& gvepsl)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_gvepsl");
    // dD/d\epsilon_{\alpha\beta}, tensor vector form of gdmepsl
    std::vector<torch::Tensor> gdmepsl_vector;
    auto accessor = gdmepsl.accessor<double, 4>();
    if (GlobalV::MY_RANK == 0)
    {
        // make gdmx as tensor
        int nlmax = this->inlmax / nat;
        for (int nl = 0; nl < nlmax; ++nl)
        {
            std::vector<torch::Tensor> bmmv;
            for (int i = 0; i < 6; ++i)
            {
                std::vector<torch::Tensor> ammv;
                for (int iat = 0; iat < nat; ++iat)
                {
                    int inl = iat * nlmax + nl;
                    int nm = 2 * this->inl_l[inl] + 1;
                    std::vector<double> mmv;
                    for (int m1 = 0; m1 < nm; ++m1)
                    {
                        for (int m2 = 0; m2 < nm; ++m2)
                        {
                            mmv.push_back(accessor[i][inl][m1][m2]);
                        }
                    } // nm^2
                    torch::Tensor mm
                        = torch::tensor(mmv, torch::TensorOptions().dtype(torch::kFloat64)).reshape({nm, nm}); // nm*nm
                    ammv.push_back(mm);
                }
                torch::Tensor bmm = torch::stack(ammv, 0); // nat*nm*nm
                bmmv.push_back(bmm);
            }
            gdmepsl_vector.push_back(torch::stack(bmmv, 0)); // nbt*3*nat*nm*nm
        }
        assert(gdmepsl_vector.size() == nlmax);

        // einsum for each inl:
        // gdmepsl_vector : b:npol * a:inl(projector) * m:nm * n:nm
        // gevdm : a:inl * v:nm (descriptor) * m:nm (pdm, dim1) * n:nm
        // (pdm, dim2) gvepsl_vector : b:npol * a:inl(projector) *
        // m:nm(descriptor)
        std::vector<torch::Tensor> gvepsl_vector;
        for (int nl = 0; nl < nlmax; ++nl)
        {
            gvepsl_vector.push_back(at::einsum("bamn, avmn->bav", {gdmepsl_vector[nl], gevdm[nl]}));
        }

        // cat nv-> \sum_nl(nv) = \sum_nl(nm_nl)=des_per_atom
        // concatenate index a(inl) and m(nm)
        gvepsl = torch::cat(gvepsl_vector, -1);
        assert(gvepsl.size(0) == 6);
        assert(gvepsl.size(1) == nat);
        assert(gvepsl.size(2) == this->des_per_atom);
    }

    return;
}

// d(Descriptor) / d(projected density matrix)
// Dimension is different for each inl, so there's a vector of tensors
void LCAO_Deepks::cal_gevdm(const int nat, std::vector<torch::Tensor>& gevdm)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_gevdm");
    if (!gevdm.empty())
    {
        gevdm.erase(gevdm.begin(), gevdm.end());
    }
    // cal gevdm(d(EigenValue(D))/dD)
    int nlmax = inlmax / nat;
    for (int nl = 0; nl < nlmax; ++nl)
    {
        std::vector<torch::Tensor> avmmv;
        for (int iat = 0; iat < nat; ++iat)
        {
            int inl = iat * nlmax + nl;
            int nm = 2 * this->inl_l[inl] + 1;
            // repeat each block for nm times in an additional dimension
            torch::Tensor tmp_x = this->pdm[inl].reshape({nm, nm}).unsqueeze(0).repeat({nm, 1, 1});
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
    return;
}

void LCAO_Deepks::load_model(const std::string& deepks_model)
{
    ModuleBase::TITLE("LCAO_Deepks", "load_model");

    try
    {
        this->module = torch::jit::load(deepks_model);
    }
    catch (const c10::Error& e)

    {
        std::cerr << "error loading the model" << std::endl;
        return;
    }
    return;
}

#endif
