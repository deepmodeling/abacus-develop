#include "esolver_ks_pw.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/binstream.h"
namespace ModuleESolver
{

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
template <typename FPTYPE, typename Device>
void ESolver_KS_PW<FPTYPE, Device>::KG(const double fwhmin,
                                       const double wcut,
                                       const double dw_in,
                                       const double dt_in,
                                       ModuleBase::matrix& wg)
{
    //-----------------------------------------------------------
    //               KS conductivity
    //-----------------------------------------------------------
    std::cout << "Calculating conductivity..." << std::endl;
    int nw = ceil(wcut / dw_in);
    double dw = dw_in / ModuleBase::Ry_to_eV; // converge unit in eV to Ry
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    double dt = dt_in;                    // unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    int nt = ceil(sqrt(46) / sigma / dt); // exp(-23) = 1e-10
    std::cout << "nw: " << nw << " ; dw: " << dw * ModuleBase::Ry_to_eV << " eV" << std::endl;
    std::cout << "nt: " << nt << " ; dt: " << dt << " a.u.(ry^-1)" << std::endl;
    assert(nw >= 1);
    assert(nt >= 1);
    const int nk = this->kv.nks;

    double* ct11 = new double[nt];
    double* ct12 = new double[nt];
    double* ct22 = new double[nt];
    ModuleBase::GlobalFunc::ZEROS(ct11, nt);
    ModuleBase::GlobalFunc::ZEROS(ct12, nt);
    ModuleBase::GlobalFunc::ZEROS(ct22, nt);

    hamilt::Velocity velop(this->pw_wfc, this->kv.isk.data(), &GlobalC::ppcell, &GlobalC::ucell, INPUT.cond_nonlocal);
    double decut = (wcut + 5*fwhmin)  / ModuleBase::Ry_to_eV;
    for (int ik = 0; ik < nk; ++ik)
    {
        velop.init(ik);
        jjcorr_ks(ik, nt, dt, decut, wg, velop, ct11, ct12, ct22);
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ct11, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct12, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ct22, nt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    //------------------------------------------------------------------
    //                    Output
    //------------------------------------------------------------------
    if (GlobalV::MY_RANK == 0)
    {
        calcondw(nt, dt, fwhmin, wcut, dw_in, ct11, ct12, ct22);
    }
    delete[] ct11;
    delete[] ct12;
    delete[] ct22;
}

template <typename FPTYPE, typename Device>
void ESolver_KS_PW<FPTYPE, Device>::jjcorr_ks(const int ik,
                                              const int nt,
                                              const double dt,
                                              const double decut,
                                              ModuleBase::matrix& wg,
                                              hamilt::Velocity& velop,
                                              double* ct11,
                                              double* ct12,
                                              double* ct22)
{
    char transn = 'N';
    char transc = 'C';
    const int ndim = 3;
    const int npwx = this->wf.npwx;
    const int nbands = GlobalV::NBANDS;
    const double ef = this->pelec->eferm.ef;
    const int npw = this->kv.ngk[ik];
    const int reducenb2 = (nbands - 1) * nbands / 2;
    bool gamma_only = false; // ABACUS do not support gamma_only yet.
    std::complex<double>* levc = &(this->psi[0](ik, 0, 0));
    std::complex<double>* prevc = new std::complex<double>[3 * npwx * nbands];
    std::complex<double>* pij = new std::complex<double>[nbands * nbands];
    double* pij2 = new double[reducenb2];
    ModuleBase::GlobalFunc::ZEROS(pij2, reducenb2);
    // px|right>
    velop.act(this->psi, nbands * GlobalV::NPOL, levc, prevc);
    for (int id = 0; id < ndim; ++id)
    {

        zgemm_(&transc,
               &transn,
               &nbands,
               &nbands,
               &npw,
               &ModuleBase::ONE,
               levc,
               &npwx,
               prevc + id * npwx * nbands,
               &npwx,
               &ModuleBase::ZERO,
               pij,
               &nbands);
#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, pij, nbands * nbands, MPI_DOUBLE_COMPLEX, MPI_SUM, POOL_WORLD);
#endif
        if (!gamma_only)
            for (int ib = 0, ijb = 0; ib < nbands; ++ib)
            {
                for (int jb = ib + 1; jb < nbands; ++jb, ++ijb)
                {
                    pij2[ijb] += norm(pij[ib * nbands + jb]);
                }
            }
    }

    if (GlobalV::RANK_IN_POOL == 0)
    {
        int nkstot = this->kv.nkstot;
        int ikglobal = this->kv.getik_global(ik);
        std::stringstream ss;
        ss << GlobalV::global_out_dir << "vmatrix" << ikglobal + 1 << ".dat";
        Binstream binpij(ss.str(), "w");
        binpij << 8 * reducenb2;
        binpij.write(pij2, reducenb2);
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
        double* enb = &(this->pelec->ekb(ik, 0));
        for (int ib = 0, ijb = 0; ib < nbands; ++ib)
        {
            double ei = enb[ib];
            double fi = wg(ik, ib);
            for (int jb = ib + 1; jb < nbands; ++jb, ++ijb)
            {
                double ej = enb[jb];
                if (ej - ei > decut )
                    continue;
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
    delete[] pij;
    delete[] prevc;
    delete[] pij2;
    return;
}

template <typename FPTYPE, typename Device>
void ESolver_KS_PW<FPTYPE, Device>::calcondw(const int nt,
                                             const double dt,
                                             const double fwhmin,
                                             const double wcut,
                                             const double dw_in,
                                             double* ct11,
                                             double* ct12,
                                             double* ct22)
{
    double factor = FACTOR;
    const int ndim = 3;
    int nw = ceil(wcut / dw_in);
    double dw = dw_in / ModuleBase::Ry_to_eV; // converge unit in eV to Ry
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    std::ofstream ofscond("je-je.txt");
    ofscond << std::setw(8) << "#t(a.u.)" << std::setw(15) << "c11(t)" << std::setw(15) << "c12(t)" << std::setw(15)
            << "c22(t)" << std::setw(15) << "decay" << std::endl;
    for (int it = 0; it < nt; ++it)
    {
        ofscond << std::setw(8) << (it)*dt << std::setw(15) << -2 * ct11[it] << std::setw(15) << -2 * ct12[it]
                << std::setw(15) << -2 * ct22[it] << std::setw(15)
                << exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) << std::endl;
    }
    ofscond.close();
    double* cw11 = new double[nw];
    double* cw12 = new double[nw];
    double* cw22 = new double[nw];
    double* kappa = new double[nw];
    ModuleBase::GlobalFunc::ZEROS(cw11, nw);
    ModuleBase::GlobalFunc::ZEROS(cw12, nw);
    ModuleBase::GlobalFunc::ZEROS(cw22, nw);
    for (int iw = 0; iw < nw; ++iw)
    {
        for (int it = 0; it < nt; ++it)
        {
            cw11[iw] += -2 * ct11[it] * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) / (iw + 0.5) / dw * dt;
            cw12[iw] += -2 * ct12[it] * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) / (iw + 0.5) / dw * dt;
            cw22[iw] += -2 * ct22[it] * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) / (iw + 0.5) / dw * dt;
        }
    }
    ofscond.open("Onsager.txt");
    ofscond << std::setw(8) << "## w(eV) " << std::setw(20) << "sigma(Sm^-1)" << std::setw(20) << "kappa(W(mK)^-1)"
            << std::setw(20) << "L12/e(Am^-1)" << std::setw(20) << "L22/e^2(Wm^-1)" << std::endl;
    for (int iw = 0; iw < nw; ++iw)
    {
        cw11[iw] *= double(2) / ndim / GlobalC::ucell.omega * factor; // unit in Sm^-1
        cw12[iw]
            *= double(2) / ndim / GlobalC::ucell.omega * factor * 2.17987092759e-18 / 1.6021766208e-19; // unit in Am^-1
        cw22[iw] *= double(2) / ndim / GlobalC::ucell.omega * factor
                    * pow(2.17987092759e-18 / 1.6021766208e-19, 2); // unit in Wm^-1
        kappa[iw] = (cw22[iw] - pow(cw12[iw], 2) / cw11[iw]) / Occupy::gaussian_parameter / ModuleBase::Ry_to_eV
                    / 11604.518026;
        ofscond << std::setw(8) << (iw + 0.5) * dw * ModuleBase::Ry_to_eV << std::setw(20) << cw11[iw] << std::setw(20)
                << kappa[iw] << std::setw(20) << cw12[iw] << std::setw(20) << cw22[iw] << std::endl;
    }
    std::cout << std::setprecision(6) << "DC electrical conductivity: " << cw11[0] - (cw11[1] - cw11[0]) * 0.5
              << " Sm^-1" << std::endl;
    std::cout << std::setprecision(6) << "Thermal conductivity: " << kappa[0] - (kappa[1] - kappa[0]) * 0.5
              << " W(mK)^-1" << std::endl;
    ;
    ofscond.close();

    delete[] cw11;
    delete[] cw12;
    delete[] cw22;
    delete[] kappa;
}

template class ESolver_KS_PW<float, psi::DEVICE_CPU>;
template class ESolver_KS_PW<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<float, psi::DEVICE_GPU>;
template class ESolver_KS_PW<double, psi::DEVICE_GPU>;
#endif
} // namespace ModuleESolver