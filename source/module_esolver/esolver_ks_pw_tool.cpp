#include "esolver_ks_pw.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "src_pw/global.h"
#include "src_pw/occupy.h"

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
void ESolver_KS_PW::KG(const int nche_KG, const double fwhmin, const double wcut, const double dw_in, const int times)
{
    //-----------------------------------------------------------
    //               KS conductivity
    //-----------------------------------------------------------
    cout << "Calculating conductivity..." << endl;
    char transn = 'N';
    char transc = 'C';
    int nw = ceil(wcut / dw_in);
    double dw = dw_in / ModuleBase::Ry_to_eV; // converge unit in eV to Ry
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    double dt = ModuleBase::PI / (dw * nw) / times; // unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    int nt = ceil(sqrt(20) / sigma / dt);
    cout << "nw: " << nw << " ; dw: " << dw * ModuleBase::Ry_to_eV << " eV" << endl;
    cout << "nt: " << nt << " ; dt: " << dt << " a.u.(ry^-1)" << endl;
    assert(nw >= 1);
    assert(nt >= 1);
    const int nk = GlobalC::kv.nks;
    const int ndim = 3;
    const int npwx = GlobalC::wf.npwx;
    const double tpiba = GlobalC::ucell.tpiba;
    const int nbands = GlobalV::NBANDS;
    const double ef = GlobalC::en.ef;

    double *ct11 = new double[nt];
    double *ct12 = new double[nt];
    double *ct22 = new double[nt];
    ModuleBase::GlobalFunc::ZEROS(ct11, nt);
    ModuleBase::GlobalFunc::ZEROS(ct12, nt);
    ModuleBase::GlobalFunc::ZEROS(ct22, nt);

    for (int ik = 0; ik < nk; ++ik)
    {
        for (int id = 0; id < ndim; ++id)
        {
            this->phami->updateHk(ik);
            const int npw = GlobalC::kv.ngk[ik];

            complex<double> *pij = new complex<double>[nbands * nbands];
            complex<double> *prevc = new complex<double>[npw * nbands];
            complex<double> *levc = &(this->psi[0](ik, 0, 0));
            double *ga = new double[npw];
            for (int ig = 0; ig < npw; ig++)
            {
                ModuleBase::Vector3<double> v3 = GlobalC::wfcpw->getgpluskcar(ik, ig);
                ga[ig] = v3[id] * tpiba;
            }
            // px|right>
            for (int ib = 0; ib < nbands; ++ib)
            {
                for (int ig = 0; ig < npw; ++ig)
                {
                    prevc[ib * npw + ig] = ga[ig] * levc[ib * npwx + ig];
                }
            }
            zgemm_(&transc,
                   &transn,
                   &nbands,
                   &nbands,
                   &npw,
                   &ModuleBase::ONE,
                   levc,
                   &npwx,
                   prevc,
                   &npw,
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
                double tmct11 = 0;
                double tmct12 = 0;
                double tmct22 = 0;
                double *enb = &(this->pelec->ekb(ik, 0));
                for (int ib = 0; ib < nbands; ++ib)
                {
                    double ei = enb[ib];
                    double fi = GlobalC::wf.wg(ik, ib);
                    for (int jb = ib + 1; jb < nbands; ++jb)
                    {
                        double ej = enb[jb];
                        double fj = GlobalC::wf.wg(ik, jb);
                        double tmct = sin((ej - ei) * (it)*dt) * (fi - fj) * norm(pij[ib * nbands + jb]);
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
            delete[] ga;
        }
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

void ESolver_KS_PW::calcondw(const int nt,
                             const double dt,
                             const double fwhmin,
                             const double wcut,
                             const double dw_in,
                             double *ct11,
                             double *ct12,
                             double *ct22)
{
    double factor = FACTOR;
    const int ndim = 3;
    int nw = ceil(wcut / dw_in);
    double dw = dw_in / ModuleBase::Ry_to_eV; // converge unit in eV to Ry
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    ofstream ofscond("je-je.txt");
    ofscond << setw(8) << "#t(a.u.)" << setw(15) << "c11(t)" << setw(15) << "c12(t)" << setw(15) << "c22(t)" << setw(15)
            << "decay" << endl;
    for (int it = 0; it < nt; ++it)
    {
        ofscond << setw(8) << (it)*dt << setw(15) << -2 * ct11[it] << setw(15) << -2 * ct12[it] << setw(15)
                << -2 * ct22[it] << setw(15) << exp(-double(1) / 2 * sigma * sigma * pow((it)*dt, 2)) << endl;
    }
    ofscond.close();
    double *cw11 = new double[nw];
    double *cw12 = new double[nw];
    double *cw22 = new double[nw];
    double *kappa = new double[nw];
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
    ofscond << setw(8) << "## w(eV) " << setw(20) << "sigma(Sm^-1)" << setw(20) << "kappa(W(mK)^-1)" << setw(20)
            << "L12/e(Am^-1)" << setw(20) << "L22/e^2(Wm^-1)" << endl;
    for (int iw = 0; iw < nw; ++iw)
    {
        cw11[iw] *= double(2) / ndim / GlobalC::ucell.omega * factor; // unit in Sm^-1
        cw12[iw]
            *= double(2) / ndim / GlobalC::ucell.omega * factor * 2.17987092759e-18 / 1.6021766208e-19; // unit in Am^-1
        cw22[iw] *= double(2) / ndim / GlobalC::ucell.omega * factor
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
} // namespace ModuleESolver