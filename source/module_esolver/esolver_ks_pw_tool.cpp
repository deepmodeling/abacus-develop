#include "esolver_ks_pw.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "src_pw/global.h"
#include "src_pw/occupy.h"
#include "module_hamilt/ks_pw/velocity_pw.h"

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
#define FACTOR      1.839939223835727e7
template<typename FPTYPE, typename Device>
void ESolver_KS_PW<FPTYPE, Device>::KG(const int nche_KG, const FPTYPE fwhmin, const FPTYPE wcut, 
        const FPTYPE dw_in, const int times, ModuleBase::matrix& wg)
{
    //-----------------------------------------------------------
    //               KS conductivity
    //-----------------------------------------------------------
    cout << "Calculating conductivity..." << endl;
    char transn = 'N';
    char transc = 'C';
    int nw = ceil(wcut / dw_in);
    FPTYPE dw = dw_in / ModuleBase::Ry_to_eV; // converge unit in eV to Ry
    FPTYPE sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    FPTYPE dt = ModuleBase::PI / (dw * nw) / times; // unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    int nt = ceil(sqrt(20) / sigma / dt);
    cout << "nw: " << nw << " ; dw: " << dw * ModuleBase::Ry_to_eV << " eV" << endl;
    cout << "nt: " << nt << " ; dt: " << dt << " a.u.(ry^-1)" << endl;
    assert(nw >= 1);
    assert(nt >= 1);
    const int nk = GlobalC::kv.nks;
    const int ndim = 3;
    const int npwx = GlobalC::wf.npwx;
    const FPTYPE tpiba = GlobalC::ucell.tpiba;
    const int nbands = GlobalV::NBANDS;
    const FPTYPE ef = GlobalC::en.ef;

    FPTYPE *ct11 = new FPTYPE[nt];
    FPTYPE *ct12 = new FPTYPE[nt];
    FPTYPE *ct22 = new FPTYPE[nt];
    ModuleBase::GlobalFunc::ZEROS(ct11, nt);
    ModuleBase::GlobalFunc::ZEROS(ct12, nt);
    ModuleBase::GlobalFunc::ZEROS(ct22, nt);

    hamilt::Velocity velop(GlobalC::wfcpw, GlobalC::kv.isk.data(),&GlobalC::ppcell,&GlobalC::ucell, INPUT.cond_nonlocal);
    for (int ik = 0; ik < nk; ++ik)
    {
        velop.init(ik);
        const int npw = GlobalC::kv.ngk[ik];
        complex<FPTYPE> *levc = &(this->psi[0](ik, 0, 0));
        complex<FPTYPE> *prevc = new complex<FPTYPE>[3 * npwx * nbands];
        // px|right>
        velop.act(this->psi, nbands*GlobalV::NPOL, levc, prevc);
        for (int id = 0; id < ndim; ++id)
        {
            this->p_hamilt->updateHk(ik);
            complex<FPTYPE> *pij = new complex<FPTYPE>[nbands * nbands];
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
            MPI_Allreduce(MPI_IN_PLACE, pij, 2 * nbands * nbands, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
#endif
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
            // for(int it = 0 ; it < nt; ++it)
            {
                FPTYPE tmct11 = 0;
                FPTYPE tmct12 = 0;
                FPTYPE tmct22 = 0;
                FPTYPE *enb = &(this->pelec->ekb(ik, 0));
                for (int ib = 0; ib < nbands; ++ib)
                {
                    FPTYPE ei = enb[ib];
                    FPTYPE fi = wg(ik, ib);
                    for (int jb = ib + 1; jb < nbands; ++jb)
                    {
                        FPTYPE ej = enb[jb];
                        FPTYPE fj = wg(ik, jb);
                        FPTYPE tmct = sin((ej - ei) * (it)*dt) * (fi - fj) * norm(pij[ib * nbands + jb]);
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
        }
        delete[] prevc;
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
void ESolver_KS_PW<FPTYPE, Device>::calcondw(const int nt,
                             const FPTYPE dt,
                             const FPTYPE fwhmin,
                             const FPTYPE wcut,
                             const FPTYPE dw_in,
                             FPTYPE *ct11,
                             FPTYPE *ct12,
                             FPTYPE *ct22)
{
    FPTYPE factor = FACTOR;
    const int ndim = 3;
    int nw = ceil(wcut / dw_in);
    FPTYPE dw = dw_in / ModuleBase::Ry_to_eV; // converge unit in eV to Ry
    FPTYPE sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    ofstream ofscond("je-je.txt");
    ofscond << setw(8) << "#t(a.u.)" << setw(15) << "c11(t)" << setw(15) << "c12(t)" << setw(15) << "c22(t)" << setw(15)
            << "decay" << endl;
    for (int it = 0; it < nt; ++it)
    {
        ofscond << setw(8) << (it)*dt << setw(15) << -2 * ct11[it] << setw(15) << -2 * ct12[it] << setw(15)
                << -2 * ct22[it] << setw(15) << exp(-FPTYPE(1) / 2 * sigma * sigma * pow((it)*dt, 2)) << endl;
    }
    ofscond.close();
    FPTYPE *cw11 = new FPTYPE[nw];
    FPTYPE *cw12 = new FPTYPE[nw];
    FPTYPE *cw22 = new FPTYPE[nw];
    FPTYPE *kappa = new FPTYPE[nw];
    ModuleBase::GlobalFunc::ZEROS(cw11, nw);
    ModuleBase::GlobalFunc::ZEROS(cw12, nw);
    ModuleBase::GlobalFunc::ZEROS(cw22, nw);
    for (int iw = 0; iw < nw; ++iw)
    {
        for (int it = 0; it < nt; ++it)
        {
            cw11[iw] += -2 * ct11[it] * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-FPTYPE(1) / 2 * sigma * sigma * pow((it)*dt, 2)) / (iw + 0.5) / dw * dt;
            cw12[iw] += -2 * ct12[it] * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-FPTYPE(1) / 2 * sigma * sigma * pow((it)*dt, 2)) / (iw + 0.5) / dw * dt;
            cw22[iw] += -2 * ct22[it] * sin(-(iw + 0.5) * dw * it * dt)
                        * exp(-FPTYPE(1) / 2 * sigma * sigma * pow((it)*dt, 2)) / (iw + 0.5) / dw * dt;
        }
    }
    ofscond.open("Onsager.txt");
    ofscond << setw(8) << "## w(eV) " << setw(20) << "sigma(Sm^-1)" << setw(20) << "kappa(W(mK)^-1)" << setw(20)
            << "L12/e(Am^-1)" << setw(20) << "L22/e^2(Wm^-1)" << endl;
    for (int iw = 0; iw < nw; ++iw)
    {
        cw11[iw] *= FPTYPE(2) / ndim / GlobalC::ucell.omega * factor; // unit in Sm^-1
        cw12[iw]
            *= FPTYPE(2) / ndim / GlobalC::ucell.omega * factor * 2.17987092759e-18 / 1.6021766208e-19; // unit in Am^-1
        cw22[iw] *= FPTYPE(2) / ndim / GlobalC::ucell.omega * factor
                    * pow(2.17987092759e-18 / 1.6021766208e-19, 2); // unit in Wm^-1
        kappa[iw] = (cw22[iw] - pow(cw12[iw], 2) / cw11[iw]) / Occupy::gaussian_parameter / ModuleBase::Ry_to_eV
                    / 11604.518026;
        ofscond << setw(8) << (iw + 0.5) * dw * ModuleBase::Ry_to_eV << setw(20) << cw11[iw] << setw(20) << kappa[iw]
                << setw(20) << cw12[iw] << setw(20) << cw22[iw] << endl;
    }
    cout << setprecision(6) << "DC electrical conductivity: " << cw11[0] - (cw11[1] - cw11[0]) * 0.5 << " Sm^-1"
         << endl;
    cout << setprecision(6) << "Thermal conductivity: " << kappa[0] - (kappa[1] - kappa[0]) * 0.5 << " Wm^-1" << endl;
    ;
    ofscond.close();

    delete[] cw11;
    delete[] cw12;
    delete[] cw22;
    delete[] kappa;
}

template class ESolver_KS_PW<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<double, psi::DEVICE_GPU>;
#endif
} // namespace ModuleESolver