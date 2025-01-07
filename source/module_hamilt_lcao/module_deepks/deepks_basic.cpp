// This file contains interfaces with libtorch,
// including loading of model and calculating gradients
// as well as subroutines that prints the results for checking

#ifdef __DEEPKS
#include "deepks_basic.h"

// d(Descriptor) / d(projected density matrix)
// Dimension is different for each inl, so there's a vector of tensors
void DeePKS_domain::cal_gevdm(const int nat,
                              const int inlmax,
                              const int* inl_l,
                              const std::vector<torch::Tensor>& pdm,
                              std::vector<torch::Tensor>& gevdm)
{
    ModuleBase::TITLE("DeePKS_domain", "cal_gevdm");
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
    return;
}

void DeePKS_domain::load_model(const std::string& model_file, torch::jit::script::Module& model)
{
    ModuleBase::TITLE("DeePKS_domain", "load_model");

    try
    {
        model = torch::jit::load(model_file);
    }
    catch (const c10::Error& e)
    {
        std::cerr << "error loading the model" << std::endl;
        return;
    }
    return;
}

#endif
