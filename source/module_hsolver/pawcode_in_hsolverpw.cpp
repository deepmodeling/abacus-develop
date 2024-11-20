#include "hsolver_pw.h"

// #include "pawcode_in_hsolverpw.h"

#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_parameter/parameter.h"
#include "module_psi/psi.h"

#include <algorithm>
#include <vector>

#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#endif

#ifdef USE_PAW

namespace hsolver
{

template <typename T, typename Device>
void HSolverPW<T, Device>::paw_func_in_kloop(const int ik, ModulePW::PW_Basis_K* wfc_basis)
{
    if (PARAM.inp.use_paw)
    {
        const int npw = this->wfc_basis->npwk[ik];
        ModuleBase::Vector3<double>* _gk = new ModuleBase::Vector3<double>[npw];
        for (int ig = 0; ig < npw; ig++)
        {
            _gk[ig] = this->wfc_basis->getgpluskcar(ik, ig);
        }

        std::vector<double> kpt(3, 0);
        kpt[0] = this->wfc_basis->kvec_c[ik].x;
        kpt[1] = this->wfc_basis->kvec_c[ik].y;
        kpt[2] = this->wfc_basis->kvec_c[ik].z;

        double** kpg;
        double** gcar;
        kpg = new double*[npw];
        gcar = new double*[npw];
        for (int ipw = 0; ipw < npw; ipw++)
        {
            kpg[ipw] = new double[3];
            kpg[ipw][0] = _gk[ipw].x;
            kpg[ipw][1] = _gk[ipw].y;
            kpg[ipw][2] = _gk[ipw].z;

            gcar[ipw] = new double[3];
            gcar[ipw][0] = this->wfc_basis->getgcar(ik, ipw).x;
            gcar[ipw][1] = this->wfc_basis->getgcar(ik, ipw).y;
            gcar[ipw][2] = this->wfc_basis->getgcar(ik, ipw).z;
        }

        GlobalC::paw_cell.set_paw_k(npw,
                                    wfc_basis->npwk_max,
                                    kpt.data(),
                                    this->wfc_basis->get_ig2ix(ik).data(),
                                    this->wfc_basis->get_ig2iy(ik).data(),
                                    this->wfc_basis->get_ig2iz(ik).data(),
                                    (const double**)kpg,
                                    GlobalC::ucell.tpiba,
                                    (const double**)gcar);

        std::vector<double>().swap(kpt);
        for (int ipw = 0; ipw < npw; ipw++)
        {
            delete[] kpg[ipw];
            delete[] gcar[ipw];
        }
        delete[] kpg;
        delete[] gcar;

        GlobalC::paw_cell.get_vkb();

        GlobalC::paw_cell.set_currentk(ik);
    }
}


} // namespace hsolver

#endif