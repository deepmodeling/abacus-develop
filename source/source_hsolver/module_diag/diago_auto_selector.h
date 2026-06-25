#ifndef DIAGO_AUTO_SELECTOR_H_
#define DIAGO_AUTO_SELECTOR_H_

#include <cstdlib>
#include <sstream>
#include <string>

namespace hsolver
{

struct DiagoAutoSelectInput
{
    std::string current_method;
    std::string calculation;
    int nbands = 0;
    int nbasis = 0;
    int npw_total = 0;
    int nproc_in_pool = 1;
    int scf_iter = 1;
    bool gpu_device = false;
};

struct DiagoAutoSelectResult
{
    std::string method;
    std::string reason;
};

class DiagoAutoSelector
{
  public:
    static bool report_enabled()
    {
        return env_enabled("ABACUS_DIAGO_AUTO_REPORT") || auto_select_enabled();
    }

    static bool auto_select_enabled()
    {
        return env_enabled("ABACUS_DIAGO_AUTO_SELECT");
    }

    static DiagoAutoSelectResult recommend_pw(const DiagoAutoSelectInput& input)
    {
        DiagoAutoSelectResult result;
        result.method = input.current_method;

        const int nbands = input.nbands > 0 ? input.nbands : 1;
        const int basis = input.npw_total > 0 ? input.npw_total : input.nbasis;
        const double basis_per_band = static_cast<double>(basis > 0 ? basis : 1) / static_cast<double>(nbands);

        std::ostringstream reason;
        reason << "basis_per_band=" << basis_per_band
               << ", nbands=" << input.nbands
               << ", nproc_pool=" << input.nproc_in_pool
               << ", scf_iter=" << input.scf_iter
               << ", calculation=" << input.calculation
               << ", device=" << (input.gpu_device ? "GPU" : "CPU");

        if (input.gpu_device)
        {
            result.method = "bpcg";
            reason << "; recommend bpcg because block CG is the GPU-oriented iterative path";
        }
        else if (input.calculation == "nscf")
        {
            result.method = "dav";
            reason << "; recommend dav because nscf usually benefits from robust full convergence";
        }
        else if (input.nproc_in_pool > 1 && input.nbands >= 32)
        {
            result.method = "bpcg";
            reason << "; recommend bpcg because many bands with MPI can benefit from block operations";
        }
        else if (input.nbands >= 64)
        {
            result.method = "ppcg";
            reason << "; recommend ppcg because many bands make projected/block updates attractive";
        }
        else if (basis_per_band > 80.0 && input.scf_iter > 1)
        {
            result.method = "dav_subspace";
            reason << "; recommend dav_subspace for large PW subspaces after the initial SCF step";
        }
        else
        {
            result.method = "cg";
            reason << "; recommend cg as the conservative default for small or early PW solves";
        }

        result.reason = reason.str();
        return result;
    }

  private:
    static bool env_enabled(const char* name)
    {
        const char* value = std::getenv(name);
        return value != nullptr && value[0] != '\0' && value[0] != '0';
    }
};

} // namespace hsolver

#endif // DIAGO_AUTO_SELECTOR_H_
