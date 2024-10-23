#include "elecond.h"

#include "module_parameter/parameter.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_elecstate/occupy.h"
#include "module_io/binstream.h"

#include <vector>

//------------------------------------------------------------------
// hbar = 6.62607015e-34/2pi
// e    = 1.6021766208e-19
// a    = 5.2917721067e-11
// m    = 9.1093897e-31
// 1 ha   = hbar^2/m/a^2/e  = 27.21136857564 eV
// 1 ry   = hbar^2/2m/a^2/e = 13.60568428782 eV = 2.17987092759e-18 J
// 1 t(ry^-1) = hbar/ry/e    = 4.837771834548454e-17 s
// factor = hbar*e^2/a^5/m^2*t^2  = 1.839939223835727e+07  (1 a.u. = 1.84e7 Sm^-1)
// 1 a.u. = factor*(2.17987092759e-18)^2/e^2 = 3.40599696130e+09 Wm^-1
// k = 1.380649e-23
// e/k = 11604.518026 , 1 eV = 11604.5 K
//------------------------------------------------------------------
#define TWOSQRT2LN2 2.354820045030949 // FWHM = 2sqrt(2ln2) * \sigma
#define FACTOR 1.839939223835727e7

EleCond::EleCond(UnitCell* p_ucell_in, K_Vectors* p_kv_in, elecstate::ElecState* p_elec_in,
                 ModulePW::PW_Basis_K* p_wfcpw_in, psi::Psi<std::complex<double>>* p_psi_in,
                 pseudopot_cell_vnl* p_ppcell_in)
{
    this->p_ppcell = p_ppcell_in;
    this->p_ucell = p_ucell_in;
    this->p_wfcpw = p_wfcpw_in;
    this->p_kv = p_kv_in;
    this->p_elec = p_elec_in;
    this->p_psi = p_psi_in;
}

void EleCond::KG(const int& smear_type, const double& fwhmin, const double& wcut, const double& dw_in,
                 const double& dt_in, const bool& nonlocal, ModuleBase::matrix& wg)
{
    //-----------------------------------------------------------
    //               KS conductivity
    //-----------------------------------------------------------
    std::cout << "Calculating conductivity..." << std::endl;
    int nw = ceil(wcut / dw_in);
    double dw = dw_in / ModuleBase::Ry_to_eV; // converge unit in eV to Ry
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    const double gamma = fwhmin / 2.0 / ModuleBase::Ry_to_eV;
    double dt = dt_in;           // unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    const double expfactor = 23; // exp(-23) = 1e-10
    int nt;                      // set nt empirically
    if (smear_type == 1)
    {
        nt = ceil(sqrt(2 * expfactor) / sigma / dt);
    }
    else if (smear_type == 2)
    {
        nt = ceil(expfactor / gamma / dt);
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW::calcondw", "smear_type should be 0 or 1");
    }
    std::cout << "nw: " << nw << " ; dw: " << dw * ModuleBase::Ry_to_eV << " eV" << std::endl;
    std::cout << "nt: " << nt << " ; dt: " << dt << " a.u.(ry^-1)" << std::endl;
    assert(nw >= 1);
    assert(nt >= 1);
    const int nk = this->p_kv->get_nks();

    std::vector<double> ct11(nt, 0);
    std::vector<double> ct12(nt, 0);
    std::vector<double> ct22(nt, 0);

    hamilt::Velocity velop(this->p_wfcpw, this->p_kv->isk.data(), this->p_ppcell, this->p_ucell, nonlocal);
    double decut = (wcut + fwhmin) / ModuleBase::Ry_to_eV;
    std::cout << "Recommended dt: " << 0.25 * M_PI / decut << " a.u." << std::endl;
    for (int ik = 0; ik < nk; ++ik)
    {
        velop.init(ik);
        jjresponse_ks(ik, nt, dt, decut, wg, velop, ct11.data(), ct12.data(), ct22.data());
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ct11.data(), nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct12.data(), nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct22.data(), nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    //------------------------------------------------------------------
    //                    Output
    //------------------------------------------------------------------
    if (GlobalV::MY_RANK == 0)
    {
        calcondw(nt, dt, smear_type, fwhmin, wcut, dw_in, ct11.data(), ct12.data(), ct22.data());
    }
}

void EleCond::jjresponse_ks(const int ik, const int nt, const double dt, const double decut, ModuleBase::matrix& wg,
                            hamilt::Velocity& velop, double* ct11, double* ct12, double* ct22)
{
    const int nbands = PARAM.inp.nbands;
    if (wg(ik, 0) - wg(ik, nbands - 1) < 1e-8 || nbands == 0) {
        return;
}
    const char transn = 'N';
    const char transc = 'C';
    const int ndim = 3;
    const int npwx = this->p_wfcpw->npwk_max;
    const double ef = this->p_elec->eferm.ef;
    const int npw = this->p_kv->ngk[ik];
    const int reducenb2 = (nbands - 1) * nbands / 2;
    const bool gamma_only = false; // ABACUS do not support gamma_only yet.
    std::complex<double>* levc = &(this->p_psi[0](ik, 0, 0));
    std::vector<std::complex<double>> prevc(ndim * npwx * nbands);
    std::vector<std::complex<double>> pij(nbands * nbands);
    std::vector<double> pij2(reducenb2, 0);
    // px|right>
    velop.act(this->p_psi, nbands * PARAM.globalv.npol, levc, prevc.data());
    for (int id = 0; id < ndim; ++id)
    {

        zgemm_(&transc, &transn, &nbands, &nbands, &npw, &ModuleBase::ONE, levc, &npwx,
               prevc.data() + id * npwx * nbands, &npwx, &ModuleBase::ZERO, pij.data(), &nbands);
#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, pij.data(), nbands * nbands, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
        if (!gamma_only) {
            for (int ib = 0, ijb = 0; ib < nbands; ++ib)
            {
                for (int jb = ib + 1; jb < nbands; ++jb, ++ijb)
                {
                    pij2[ijb] += norm(pij[ib * nbands + jb]);
                }
            }
}
    }

    if (GlobalV::RANK_IN_POOL == 0)
    {
        int nkstot = this->p_kv->get_nkstot();
        int ikglobal = K_Vectors::get_ik_global(ik, nkstot);
        std::stringstream ss;
        ss << PARAM.globalv.global_out_dir << "vmatrix" << ikglobal + 1 << ".dat";
        Binstream binpij(ss.str(), "w");
        binpij << 8 * reducenb2;
        binpij.write(pij2.data(), reducenb2);
        binpij << 8 * reducenb2;
    }

    int ntper = nt / GlobalV::NPROC_IN_POOL;
    int itstart = ntper * GlobalV::RANK_IN_POOL;
    if (nt % GlobalV::NPROC_IN_POOL > GlobalV::RANK_IN_POOL)
    {
        ntper++;
        itstart += GlobalV::RANK_IN_POOL;
    }
    else
    {
        itstart += nt % GlobalV::NPROC_IN_POOL;
    }

    for (int it = itstart; it < itstart + ntper; ++it)
    {
        double tmct11 = 0;
        double tmct12 = 0;
        double tmct22 = 0;
        double* enb = &(this->p_elec->ekb(ik, 0));
        for (int ib = 0, ijb = 0; ib < nbands; ++ib)
        {
            double ei = enb[ib];
            double fi = wg(ik, ib);
            for (int jb = ib + 1; jb < nbands; ++jb, ++ijb)
            {
                double ej = enb[jb];
                if (ej - ei > decut) {
                    continue;
}
                double fj = wg(ik, jb);
                double tmct = sin((ej - ei) * (it)*dt) * (fi - fj) * pij2[ijb];
                tmct11 += tmct;
                tmct12 += -tmct * ((ei + ej) / 2 - ef);
                tmct22 += tmct * pow((ei + ej) / 2 - ef, 2);
            }
        }
        ct11[it] += tmct11 / 2.0;
        ct12[it] += tmct12 / 2.0;
        ct22[it] += tmct22 / 2.0;
    }
    return;
}

void EleCond::calcondw(const int nt, const double dt, const int& smear_type, const double fwhmin, const double wcut,
                       const double dw_in, double* ct11, double* ct12, double* ct22)
{
    double factor = FACTOR;
    const int ndim = 3;
    int nw = ceil(wcut / dw_in);
    double dw = dw_in / ModuleBase::Ry_to_eV; // converge unit in eV to Ry
    const double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    const double gamma = fwhmin / 2.0 / ModuleBase::Ry_to_eV;
    std::vector<double> winfunc(nt);
    // 1: Gaussian, 2: Lorentzian
    if (smear_type == 1)
    {
        for (int it = 0; it < nt; ++it)
        {
            winfunc[it] = exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2));
        }
    }
    else if (smear_type == 2)
    {
        for (int it = 0; it < nt; ++it)
        {
            winfunc[it] = exp(-gamma * (it)*dt);
        }
    }

    std::ofstream ofscond("je-je.txt");
    ofscond << std::setw(8) << "#t(a.u.)" << std::setw(15) << "c11(t)" << std::setw(15) << "c12(t)" << std::setw(15)
            << "c22(t)" << std::setw(15) << "decay" << std::endl;
    for (int it = 0; it < nt; ++it)
    {
        ofscond << std::setw(8) << (it)*dt << std::setw(15) << -2 * ct11[it] << std::setw(15) << -2 * ct12[it]
                << std::setw(15) << -2 * ct22[it] << std::setw(15) << winfunc[it] << std::endl;
    }
    ofscond.close();
    std::vector<double> cw11(nw, 0);
    std::vector<double> cw12(nw, 0);
    std::vector<double> cw22(nw, 0);
    std::vector<double> kappa(nw);
    for (int iw = 0; iw < nw; ++iw)
    {
        for (int it = 0; it < nt; ++it)
        {
            cw11[iw] += -2 * ct11[it] * sin(-(iw + 0.5) * dw * it * dt) * winfunc[it] / (iw + 0.5) / dw * dt;
            cw12[iw] += -2 * ct12[it] * sin(-(iw + 0.5) * dw * it * dt) * winfunc[it] / (iw + 0.5) / dw * dt;
            cw22[iw] += -2 * ct22[it] * sin(-(iw + 0.5) * dw * it * dt) * winfunc[it] / (iw + 0.5) / dw * dt;
        }
    }
    ofscond.open("Onsager.txt");
    ofscond << std::setw(8) << "## w(eV) " << std::setw(20) << "sigma(Sm^-1)" << std::setw(20) << "kappa(W(mK)^-1)"
            << std::setw(20) << "L12/e(Am^-1)" << std::setw(20) << "L22/e^2(Wm^-1)" << std::endl;
    for (int iw = 0; iw < nw; ++iw)
    {
        cw11[iw] *= double(2) / ndim / this->p_ucell->omega * factor; // unit in Sm^-1
        cw12[iw]
            *= double(2) / ndim / this->p_ucell->omega * factor * 2.17987092759e-18 / 1.6021766208e-19; // unit in Am^-1
        cw22[iw] *= double(2) / ndim / this->p_ucell->omega * factor
                    * pow(2.17987092759e-18 / 1.6021766208e-19, 2); // unit in Wm^-1
        kappa[iw] = (cw22[iw] - pow(cw12[iw], 2) / cw11[iw]) / Occupy::gaussian_parameter / ModuleBase::Ry_to_eV
                    / 11604.518026;
        ofscond << std::setw(8) << (iw + 0.5) * dw * ModuleBase::Ry_to_eV << std::setw(20) << cw11[iw] << std::setw(20)
                << kappa[iw] << std::setw(20) << cw12[iw] << std::setw(20) << cw22[iw] << std::endl;
    }
    double sigma0 = cw11[0] - (cw11[1] - cw11[0]) * 0.5;
    double kappa0 = kappa[0] - (kappa[1] - kappa[0]) * 0.5;
    double Lorent0 = kappa0 / sigma0 / Occupy::gaussian_parameter / ModuleBase::Ry_to_eV / 11604.518026
                     * pow(1.6021766208e-19 / 1.3806505e-23, 2);
    std::cout << std::setprecision(6) << "DC electrical conductivity: " << sigma0 << " Sm^-1" << std::endl;
    std::cout << std::setprecision(6) << "Thermal conductivity: " << kappa0 << " W(mK)^-1" << std::endl;
    std::cout << std::setprecision(6) << "Lorenz number: " << Lorent0 << " k_B^2/e^2" << std::endl;
    ofscond.close();
}