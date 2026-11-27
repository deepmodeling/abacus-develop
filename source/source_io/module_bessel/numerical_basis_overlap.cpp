#include "numerical_basis.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/constants.h"
#include "source_base/global_variable.h"
#include "source_base/intarray.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/vector3.h"
#include "numerical_basis_jyjy.h"

#include <algorithm>
#include <cstring>
#include <functional>
#include <vector>

// The function is called in run_fp.cpp.
void Numerical_Basis::output_overlap(const psi::Psi<std::complex<double>>& psi,
                                     const Structure_Factor& sf,
                                     const K_Vectors& kv,
                                     const ModulePW::PW_Basis_K* wfcpw,
                                     const UnitCell& ucell,
                                     const int& index)
{
    ModuleBase::TITLE("Numerical_Basis", "output_overlap");
    ModuleBase::GlobalFunc::NEW_PART("Overlap Data For Spillage Minimization");
    const double bessel_nao_rcut = PARAM.inp.bessel_nao_rcuts[index];

    //---------------------------------------------------------
    // if the numerical_basis hasn't been initialized yet,
    // then we initial here.
    //---------------------------------------------------------
    if (!this->init_label)
    {
        // false stands for : 'Faln' is not used.
        this->bessel_basis.init(false, std::stod(PARAM.inp.bessel_nao_ecut), ucell.ntype, ucell.lmax,
                                PARAM.inp.bessel_nao_smooth, PARAM.inp.bessel_nao_sigma, bessel_nao_rcut,
                                PARAM.inp.bessel_nao_tolerence, ucell);
        this->mu_index = this->init_mu_index(ucell);
        this->init_label = true;
    }
    ModuleBase::GlobalFunc::MAKE_DIR(PARAM.inp.spillage_outdir);
    for (int derivative_order = 0; derivative_order <= 1; ++derivative_order) // Peize Lin add 2020.04.23
    {
        std::ofstream ofs;
        std::stringstream ss;
        ss << PARAM.inp.spillage_outdir << "/";

        if (PARAM.inp.bessel_nao_rcuts.size() > 1)
        {
            ss << "orb_matrix_rcut" << bessel_nao_rcut << "deriv";
        }
        else
        {
            ss << "orb_matrix.";
        } // to make it compatible with old version of orbital generation
        ss << derivative_order << ".dat";

        if (GlobalV::MY_RANK == 0)
        {
            ofs.open(ss.str().c_str());
        }

        // ALLOCATE MEMORY FOR THE OVERLAP MATRIX
        // OVERLAP : < J_mu | Psi >
        std::vector<ModuleBase::ComplexArray> overlap_Q(kv.get_nks());
        // OVERLAP : < J_mu | J_nu >
        std::vector<ModuleBase::ComplexArray> overlap_Sq(kv.get_nks());

        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of k points", kv.get_nks());
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of bands", PARAM.inp.nbands);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of local orbitals", PARAM.globalv.nlocal);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "number of eigenvalues of Jl(x)",
                                    this->bessel_basis.get_ecut_number());

        // CALCULATE THE OVERLAP MATRIX
        // nks now is the reduced k-points.
        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            const int npw = kv.ngk[ik];
            GlobalV::ofs_running << " --------------------------------------------------------" << std::endl;
            GlobalV::ofs_running << " Print the overlap matrixs Q and S for this kpoint" << std::endl;
            GlobalV::ofs_running << std::setw(8) << "ik" << std::setw(8) << "npw" << std::endl;
            GlobalV::ofs_running << std::setw(8) << ik + 1 << std::setw(8) << npw << std::endl;
            GlobalV::ofs_running << " --------------------------------------------------------" << std::endl;

            // search for all k-points.
            psi.fix_k(ik);
            overlap_Q[ik] = this->cal_overlap_Q(ik, npw, wfcpw, psi, static_cast<double>(derivative_order), sf, ucell);
            ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_overlap_Q");

            // (2) generate Sq matrix if necessary.
            if (PARAM.inp.out_spillage == 2)
            {
#ifndef __LCAO
                // compute <jY|jY> in plane-wave basis
                overlap_Sq[ik] = this->cal_overlap_Sq(ik, npw, static_cast<double>(derivative_order), sf, wfcpw, ucell);
#else
                // compute <jY|jY> with two-center integration
                assert(derivative_order == 0 || derivative_order == 1);
                char type = (derivative_order == 0) ? 'S' : 'T';
                std::vector<int> natom;
                std::vector<int> lmax;
                std::vector<std::vector<ModuleBase::Vector3<double>>> tau_cart;
                for (int it = 0; it < ucell.ntype; ++it)
                {
                    natom.push_back(ucell.atoms[it].na);
                    lmax.push_back(ucell.atoms[it].nwl);
                    tau_cart.emplace_back();

                    for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
                    {
                        tau_cart[it].push_back(ucell.atoms[it].tau[ia] * ucell.lat0);
                    }
                }

                overlap_Sq[ik] = NumericalBasis::cal_overlap_Sq(
                    type, ucell.lmaxmax, this->bessel_basis.get_ecut_number(), bessel_nao_rcut, tau_cart,
                    ucell.lat0 * ucell.latvec, NumericalBasis::indexgen(natom, lmax));
#endif
                ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "cal_overlap_Sq");
            }
        }

        const ModuleBase::matrix overlap_V
            = this->cal_overlap_V(wfcpw, psi, static_cast<double>(derivative_order), kv, ucell.tpiba);

        // ALTHOUGH THIS FUNCTION NAMES output_overlap, IT ACTUALLY OUTPUTS THE OVERLAP MATRIX HERE
#ifdef __MPI
        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            Parallel_Reduce::reduce_pool(overlap_Q[ik].ptr, overlap_Q[ik].getSize());
            // Parallel_Reduce::reduce_pool(overlap_Sq[ik].ptr, overlap_Sq[ik].getSize());
        }
        Parallel_Reduce::reduce_pool(overlap_V.c, overlap_V.nr * overlap_V.nc); // Peize Lin add 2020.04.23
#endif
        // exception handling following, for FileNotOpenFailure
        if (ofs.good()) {
            this->output_info(ofs, bessel_basis, kv, ucell); // header of orb_matrix* file
        } else {
            ModuleBase::WARNING_QUIT("Numerical_Basis", "Failed to open file for writing the overlap matrix.");
        }
        // because one stage of file io complete, re-check the file status.
        if (ofs.good()) {
            this->output_k(ofs, kv); // <WEIGHTS_OF_KPOINTS>...</WEIGHTS_OF_KPOINTS>
        } else {
            ModuleBase::WARNING_QUIT("Numerical_Basis", "Failed to write k-points to file.");
        }
        // because one stage of file io complete, re-check the file status.
        if (ofs.good()) {
            this->output_overlap_Q(ofs, overlap_Q, kv); // <OVERLAP_Q>...</OVERLAP_Q>
        } else {
            ModuleBase::WARNING_QUIT("Numerical_Basis", "Failed to write overlap Q to file.");
        }
        // because one stage of file io complete, re-check the file status.
        if (PARAM.inp.out_spillage == 2)
        {
            // caution: this is the largest matrix to be output, always flush
            if (ofs.good()) {
                this->output_overlap_Sq(ss.str(), ofs, overlap_Sq, kv); // <OVERLAP_Sq>...</OVERLAP_Sq>
            } else {
                ModuleBase::WARNING_QUIT("Numerical_Basis", "Failed to write overlap S to file.");
            }
        }
        // because one stage of file io complete, re-check the file status.
        if (ofs.good()) {
            this->output_overlap_V(ofs, overlap_V); // <OVERLAP_V>...</OVERLAP_V>
                                                    // Peize Lin add 2020.04.23
        } else {
            ModuleBase::WARNING_QUIT("Numerical_Basis", "Failed to write overlap V to file.");
        }
        if (GlobalV::MY_RANK == 0) {
            ofs.close();
        }
    }
    return;
}

ModuleBase::ComplexArray Numerical_Basis::cal_overlap_Q(const int& ik, const int& np, const ModulePW::PW_Basis_K* wfcpw,
                                                        const psi::Psi<std::complex<double>>& psi,
                                                        const double derivative_order, const Structure_Factor& sf,
                                                        const UnitCell& ucell) const
{
    ModuleBase::TITLE("Numerical_Basis", "cal_overlap_Q");
    ModuleBase::timer::start("Numerical_Basis", "cal_overlap_Q");

    GlobalV::ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS AND BLOCH WAVE FUNCTIONS"
                         << std::endl;
    GlobalV::ofs_running << " Q = < J_mu, q | Psi_n, k > " << std::endl;

    ModuleBase::ComplexArray overlap_Q(PARAM.inp.nbands, PARAM.globalv.nlocal, this->bessel_basis.get_ecut_number());
    overlap_Q.zero_out();

    const double normalization = (4 * ModuleBase::PI) / sqrt(ucell.omega); // Peize Lin add normalization 2015-12-29

    std::vector<ModuleBase::Vector3<double>> gk(np);
    for (int ig = 0; ig < np; ig++)
    {
        gk[ig] = wfcpw->getgpluskcar(ik, ig) * ucell.tpiba;
    }

    const std::vector<double> gpow = Numerical_Basis::cal_gpow(gk, derivative_order);

    const ModuleBase::realArray flq = this->cal_flq(gk, ucell.lmax);

    const ModuleBase::matrix ylm = Numerical_Basis::cal_ylm(gk, ucell.lmax);

    GlobalV::ofs_running << "\n " << std::setw(5) << "ik" << std::setw(8) << "Type1" << std::setw(8) << "Atom1"
                         << std::setw(8) << "L" << std::endl;

    for (int T = 0; T < ucell.ntype; T++)
    {
        // OUT("T",T);
        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            // OUT("I",I);
            std::complex<double>* sk = sf.get_sk(ik, T, I, wfcpw);
            for (int L = 0; L < ucell.atoms[T].nwl + 1; L++)
            {
                GlobalV::ofs_running << " " << std::setw(5) << ik + 1 << std::setw(8) << ucell.atoms[T].label
                                     << std::setw(8) << I + 1 << std::setw(8) << L << std::endl;
                // OUT("l",l);
                std::complex<double> lphase
                    = normalization * pow(ModuleBase::IMAG_UNIT, -L); // Peize Lin add normalization 2015-12-29
                for (int ie = 0; ie < this->bessel_basis.get_ecut_number(); ie++)
                {
                    const int N = 0;
                    assert(ucell.nmax == 1);
                    for (int m = 0; m < 2 * L + 1; m++)
                    {
                        const int lm = L * L + m;
                        for (int ib = 0; ib < PARAM.inp.nbands; ib++)
                        {
                            std::complex<double> overlap_tmp = ModuleBase::ZERO;
                            for (int ig = 0; ig < np; ig++)
                            {
                                const std::complex<double> local_tmp = lphase * sk[ig] * ylm(lm, ig) * flq(L, ie, ig)
                                                                       * gpow[ig]; // Peize Lin add for dpsi 2020.04.23
                                overlap_tmp += conj(local_tmp) * psi(ib, ig);      // psi is bloch orbitals
                            }
                            overlap_Q(ib, this->mu_index[T](I, L, N, m), ie) = overlap_tmp;
                        }
                    }
                } // end ie
            }     // end l
            delete[] sk;
            sk = nullptr;
        }
    }

    ModuleBase::timer::end("Numerical_Basis", "cal_overlap_Q");
    return overlap_Q;
}

ModuleBase::ComplexArray Numerical_Basis::cal_overlap_Sq(const int& ik, const int& np, const double derivative_order,
                                                         const Structure_Factor& sf, const ModulePW::PW_Basis_K* wfcpw,
                                                         const UnitCell& ucell) const
{
    ModuleBase::TITLE("Numerical_Basis", "cal_overlap_Sq");
    ModuleBase::timer::start("Numerical_Basis", "cal_overlap_Sq");

    GlobalV::ofs_running << " OUTPUT THE OVERLAP BETWEEN SPHERICAL BESSEL FUNCTIONS" << std::endl;
    GlobalV::ofs_running << " S = < J_mu,q1 | J_nu,q2 >" << std::endl;

    const int enumber = this->bessel_basis.get_ecut_number();
    ModuleBase::ComplexArray overlap_Sq(PARAM.globalv.nlocal, PARAM.globalv.nlocal, enumber, enumber);
    overlap_Sq.zero_out();

    const double normalization
        = (4 * ModuleBase::PI) * (4 * ModuleBase::PI) / ucell.omega; // Peize Lin add normalization 2015-12-29

    std::vector<ModuleBase::Vector3<double>> gk(np);
    for (int ig = 0; ig < np; ig++) {
        gk[ig] = wfcpw->getgpluskcar(ik, ig) * ucell.tpiba;
    }

    const std::vector<double> gpow = Numerical_Basis::cal_gpow(gk, derivative_order);

    const ModuleBase::realArray flq = this->cal_flq(gk, ucell.lmax);

    const ModuleBase::matrix ylm = Numerical_Basis::cal_ylm(gk, ucell.lmax);

    GlobalV::ofs_running << "\n " << std::setw(5) << "ik" << std::setw(8) << "Type1" << std::setw(8) << "Atom1"
                         << std::setw(8) << "L1" << std::setw(8) << "Type2" << std::setw(8) << "Atom2" << std::setw(8)
                         << "L2" << std::endl;

    for (int T1 = 0; T1 < ucell.ntype; T1++) // 1.1
    {
        for (int I1 = 0; I1 < ucell.atoms[T1].na; I1++) // 1.2
        {
            std::complex<double>* sk1 = sf.get_sk(ik, T1, I1, wfcpw);
            for (int T2 = 0; T2 < ucell.ntype; T2++) // 2.1
            {
                for (int I2 = 0; I2 < ucell.atoms[T2].na; I2++) // 2.2
                {
                    std::complex<double>* sk2 = sf.get_sk(ik, T2, I2, wfcpw);
                    for (int l1 = 0; l1 < ucell.atoms[T1].nwl + 1; l1++) // 1.3
                    {
                        const std::complex<double> lphase1
                            = normalization * pow(ModuleBase::IMAG_UNIT, l1); // Peize Lin add normalization 2015-12-29
                        for (int l2 = 0; l2 < ucell.atoms[T2].nwl + 1; l2++)  // 2.3
                        {
                            GlobalV::ofs_running << " " << std::setw(5) << ik + 1 << std::setw(8)
                                                 << ucell.atoms[T1].label << std::setw(8) << I1 + 1 << std::setw(8)
                                                 << l1 << std::setw(8) << ucell.atoms[T2].label << std::setw(8)
                                                 << I2 + 1 << std::setw(8) << l2 << std::setw(8) << std::endl;

                            const std::complex<double> lphase2 = pow(ModuleBase::IMAG_UNIT, l2);
                            for (int ic1 = 0; ic1 < ucell.nmax; ic1++) // 1.5
                            {
                                for (int ic2 = 0; ic2 < ucell.nmax; ic2++) // 2.5
                                {
                                    for (int m1 = 0; m1 < 2 * l1 + 1; m1++) // 1.6
                                    {
                                        const int lm1 = l1 * l1 + m1;
                                        const int iwt1 = this->mu_index[T1](I1, l1, ic1, m1);

                                        std::vector<std::complex<double>> about_ig1(np, std::complex<double>(0.0, 0.0));
                                        for (int ig = 0; ig < np; ig++) {
                                            about_ig1[ig] = conj(lphase1 * sk1[ig] * ylm(lm1, ig))
                                                            * gpow[ig]; // Peize Lin add for dpsi 2020.04.23
                                        }

                                        for (int m2 = 0; m2 < 2 * l2 + 1; m2++) // 2.6
                                        {
                                            const int lm2 = l2 * l2 + m2;
                                            const int iwt2 = this->mu_index[T2](I2, l2, ic2, m2);

                                            std::vector<std::complex<double>> about_ig2(np,
                                                                                        std::complex<double>(0.0, 0.0));
                                            for (int ig = 0; ig < np; ++ig) {
                                                about_ig2[ig] = lphase2 * sk2[ig] * ylm(lm2, ig) * about_ig1[ig];
                                            }

                                            ModuleBase::ComplexMatrix about_ig3_1(enumber, np);
                                            std::copy(&flq(l1, 0, 0), &flq(l1, 0, 0) + enumber * np, about_ig3_1.c);

                                            ModuleBase::ComplexMatrix about_ig3_2(enumber, np);
                                            for (int ie2 = 0; ie2 < enumber; ++ie2) {
                                                std::transform(&flq(l2, ie2, 0), &flq(l2, ie2, 0) + np,
                                                               about_ig2.data(), about_ig3_2.c + ie2 * np,
                                                               std::multiplies<std::complex<double>>());
                                            }

                                            BlasConnector::gemm('N', 'T', enumber, enumber, np, 1.0, about_ig3_1.c, np,
                                                                about_ig3_2.c, np, 1.0, &overlap_Sq(iwt1, iwt2, 0, 0),
                                                                enumber);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    delete[] sk2;
                    sk2 = nullptr;
                }
            }
            delete[] sk1;
            sk1 = nullptr;
        }
    }

    ModuleBase::timer::end("Numerical_Basis", "cal_overlap_Sq");
    return overlap_Sq;
}

// Peize Lin add for dpsi 2020.04.23
ModuleBase::matrix Numerical_Basis::cal_overlap_V(const ModulePW::PW_Basis_K* wfcpw,
                                                  const psi::Psi<std::complex<double>>& psi,
                                                  const double derivative_order, const K_Vectors& kv,
                                                  const double tpiba)
{
    ModuleBase::matrix overlap_V(kv.get_nks(), PARAM.inp.nbands);
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        std::vector<ModuleBase::Vector3<double>> gk(kv.ngk[ik]);
        for (int ig = 0; ig < gk.size(); ig++) {
            gk[ig] = wfcpw->getgpluskcar(ik, ig) * tpiba;
        }

        const std::vector<double> gpow = Numerical_Basis::cal_gpow(gk, derivative_order);

        for (int ib = 0; ib < PARAM.inp.nbands; ++ib) {
            for (int ig = 0; ig < kv.ngk[ik]; ++ig) {
                overlap_V(ik, ib) += norm(psi(ik, ib, ig)) * gpow[ig];
            }
        }
    }
    return overlap_V;
}

ModuleBase::realArray Numerical_Basis::cal_flq(const std::vector<ModuleBase::Vector3<double>>& gk,
                                               const int ucell_lmax) const
{
    const int np = gk.size();
    const int enumber = this->bessel_basis.get_ecut_number();

    // get flq(G) = \int f(r)jl(G*r) from interpolation table.
    ModuleBase::realArray flq(ucell_lmax + 1, enumber, np);
    for (int il = 0; il < ucell_lmax + 1; il++)
    {
        for (int ie = 0; ie < enumber; ie++)
        {
            for (int ig = 0; ig < np; ig++)
            {
                flq(il, ie, ig) = this->bessel_basis.Polynomial_Interpolation2(il, ie, gk[ig].norm());
            }
        }
    }
    return flq;
}

ModuleBase::matrix Numerical_Basis::cal_ylm(const std::vector<ModuleBase::Vector3<double>>& gk, const int ucell_lmax)
{
    const int total_lm = (ucell_lmax + 1) * (ucell_lmax + 1);
    ModuleBase::matrix ylm(total_lm, gk.size());
    ModuleBase::YlmReal::Ylm_Real(total_lm, gk.size(), gk.data(), ylm);
    return ylm;
}

std::vector<double> Numerical_Basis::cal_gpow(const std::vector<ModuleBase::Vector3<double>>& gk,
                                              const double derivative_order)
{
    constexpr double thr = 1E-12;
    std::vector<double> gpow(gk.size(), 0.0);
    for (int ig = 0; ig < gpow.size(); ++ig)
    {
        if (derivative_order >= 0)
        {
            gpow[ig] = std::pow(gk[ig].norm2(), derivative_order);
        }
        else
        {
            if (gk[ig].norm2() >= thr) {
                gpow[ig] = std::pow(gk[ig].norm2(), derivative_order);
            }
        }
    }
    return gpow;
}
