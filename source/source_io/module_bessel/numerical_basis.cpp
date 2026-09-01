#include "numerical_basis.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include "source_base/intarray.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/vector3.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "numerical_basis_jyjy.h"

#include <algorithm>
#include <cstring>
#include <functional>
#include <vector>
Numerical_Basis::Numerical_Basis()
{
}
Numerical_Basis::~Numerical_Basis()
{
}

//============================================================
// MEMBER FUNCTION :
// NAME : init
// DESCRIPTION : Two main functions:
// (1) start_from_file = true;
// Firstly, use check(1) to call bessel_basis.init
// to generate TableOne.
// Secondly readin C4 from file.
// Thirdly generate 3D atomic wfc in G space, put the
// results in psi.
//
// (2) If output overlap Q, start_from_file = false;
// Firstly, use check(0) to call bessel_basis,init
// to generate TableOne
// Secondly output overlap, use psi(evc) and jlq3d.
//============================================================
// void Numerical_Basis::start_from_file_k(const int& ik, ModuleBase::ComplexMatrix& psi, const Structure_Factor& sf,
//                                         const ModulePW::PW_Basis_K* wfcpw, const UnitCell& ucell)
// {
//     ModuleBase::TITLE("Numerical_Basis", "start_from_file_k");

//     if (!this->init_label)
//     {
//         // true stands for : start_from_file
//         this->bessel_basis.init(true, std::stod(PARAM.inp.bessel_nao_ecut), ucell.ntype, ucell.lmax,
//                                 PARAM.inp.bessel_nao_smooth, PARAM.inp.bessel_nao_sigma, PARAM.globalv.bessel_nao_rcut,
//                                 PARAM.inp.bessel_nao_tolerence, ucell);
//         this->mu_index = this->init_mu_index(ucell);
//         this->init_label = true;
//     }
//     this->numerical_atomic_wfc(ik, wfcpw, psi, sf, ucell);
// }

std::vector<ModuleBase::IntArray> Numerical_Basis::init_mu_index(const UnitCell& ucell)
{
    GlobalV::ofs_running << " Initialize the mu index" << std::endl;
    std::vector<ModuleBase::IntArray> mu_index_(ucell.ntype);

    int mu = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        mu_index_[it].create(ucell.atoms[it].na, ucell.atoms[it].nwl + 1, ucell.nmax,
                             2 * (ucell.atoms[it].nwl + 1) + 1); // m ==> 2*l+1

        mu_index_[it].zero_out();

        // mohan added 2021-01-03
        GlobalV::ofs_running << "Type " << it + 1 << " number_of_atoms " << ucell.atoms[it].na << " number_of_L "
                             << ucell.atoms[it].nwl + 1 << " number_of_n " << ucell.nmax << " number_of_m "
                             << 2 * (ucell.atoms[it].nwl + 1) + 1 << std::endl;

        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            for (int l = 0; l < ucell.atoms[it].nwl + 1; l++)
            {
                for (int n = 0; n < ucell.atoms[it].l_nchi[l]; n++)
                {
                    for (int m = 0; m < 2 * l + 1; m++)
                    {
                        mu_index_[it](ia, l, n, m) = mu;
                        mu++;
                    }
                }
            }
        }
    }
    return mu_index_;
}

void Numerical_Basis::numerical_atomic_wfc(const int& ik, const ModulePW::PW_Basis_K* wfcpw,
                                           ModuleBase::ComplexMatrix& psi, const Structure_Factor& sf,
                                           const UnitCell& ucell)
{
    ModuleBase::TITLE("Numerical_Basis", "numerical_atomic_wfc");
    const int np = wfcpw->npwk[ik];
    std::vector<ModuleBase::Vector3<double>> gk(np);
    for (int ig = 0; ig < np; ig++) {
        gk[ig] = wfcpw->getgpluskcar(ik, ig);
}

    const int total_lm = (ucell.lmax + 1) * (ucell.lmax + 1);
    ModuleBase::matrix ylm(total_lm, np);
    ModuleBase::YlmReal::Ylm_Real(total_lm, np, gk.data(), ylm);

    std::vector<double> flq(np);
    for (int it = 0; it < ucell.ntype; it++)
    {
        // OUT("it",it);
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            // OUT("ia",ia);
            std::complex<double>* sk = sf.get_sk(ik, it, ia, wfcpw);
            for (int l = 0; l < ucell.atoms[it].nwl + 1; l++)
            {
                // OUT("l",l);
                std::complex<double> lphase = pow(ModuleBase::IMAG_UNIT, l);
                for (int ic = 0; ic < ucell.atoms[it].l_nchi[l]; ic++)
                {
                    // OUT("ic",ic);
                    for (int ig = 0; ig < np; ig++)
                    {
                        flq[ig] = this->bessel_basis.Polynomial_Interpolation(it, l, ic, gk[ig].norm() * ucell.tpiba);
                    }

                    for (int m = 0; m < 2 * l + 1; m++)
                    {
                        // OUT("m",m);
                        const int lm = l * l + m;
                        for (int ig = 0; ig < np; ig++)
                        {
                            psi(this->mu_index[it](ia, l, ic, m), ig) = lphase * sk[ig] * ylm(lm, ig) * flq[ig];
                        }
                    }
                }
            }
            delete[] sk;
            sk = nullptr;
        }
    }
}
