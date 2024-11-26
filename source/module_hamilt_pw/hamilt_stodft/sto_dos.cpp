#include "sto_dos.h"

#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_parameter/parameter.h"
#include "sto_tool.h"
Sto_DOS::~Sto_DOS()
{
}

Sto_DOS::Sto_DOS(ModulePW::PW_Basis_K* p_wfcpw_in,
                 K_Vectors* p_kv_in,
                 elecstate::ElecState* p_elec_in,
                 psi::Psi<std::complex<double>>* p_psi_in,
                 hamilt::Hamilt<std::complex<double>>* p_hamilt_in,
                 StoChe<double>& stoche,
                 Stochastic_WF<std::complex<double>, base_device::DEVICE_CPU>* p_stowf_in)
{
    this->p_wfcpw = p_wfcpw_in;
    this->p_kv = p_kv_in;
    this->p_elec = p_elec_in;
    this->p_psi = p_psi_in;
    this->p_hamilt = p_hamilt_in;
    this->p_hamilt_sto = static_cast<hamilt::HamiltSdftPW<std::complex<double>>*>(p_hamilt_in);
    this->p_stowf = p_stowf_in;
    this->nbands_ks = p_psi_in->get_nbands();
    this->nbands_sto = p_stowf_in->nchi;
    this->method_sto = stoche.method_sto;
    this->stofunc.set_E_range(&stoche.emin_sto, &stoche.emax_sto);
}
void Sto_DOS::decide_param(const int& dos_nche,
                           const double& emin_sto,
                           const double& emax_sto,
                           const bool& dos_setemin,
                           const bool& dos_setemax,
                           const double& dos_emin_ev,
                           const double& dos_emax_ev,
                           const double& dos_scale)
{
    this->dos_nche = dos_nche;
    check_che(this->dos_nche, emin_sto, emax_sto, this->nbands_sto, this->p_kv, this->p_stowf, this->p_hamilt_sto);
    if (dos_setemax)
    {
        this->emax = dos_emax_ev;
    }
    else
    {
        this->emax = *p_hamilt_sto->emax * ModuleBase::Ry_to_eV;
    }
    if (dos_setemin)
    {
        this->emin = dos_emin_ev;
    }
    else
    {
        this->emin = *p_hamilt_sto->emin * ModuleBase::Ry_to_eV;
    }

    if (!dos_setemax && !dos_setemin)
    {
        double delta = (emax - emin) * dos_scale;
        this->emax = emax + delta / 2.0;
        this->emin = emin - delta / 2.0;
    }
}

void Sto_DOS::caldos(const double sigmain, const double de, const int npart)
{
    ModuleBase::TITLE("Sto_DOS", "caldos");
    ModuleBase::timer::tick("Sto_DOS", "caldos");
    std::cout << "=========================" << std::endl;
    std::cout << "###Calculating Dos....###" << std::endl;
    std::cout << "=========================" << std::endl;
    ModuleBase::Chebyshev<double> che(dos_nche);
    const int nk = p_kv->get_nks();
    const int npwx = p_wfcpw->npwk_max;

    std::vector<double> spolyv;
    std::vector<std::complex<double>> allorderchi;
    if (this->method_sto == 1)
    {
        spolyv.resize(dos_nche, 0);
    }
    else
    {
        spolyv.resize(dos_nche * dos_nche, 0);
        int nchip_new = ceil((double)this->p_stowf->nchip_max / npart);
        allorderchi.resize(nchip_new * npwx * dos_nche);
    }
    ModuleBase::timer::tick("Sto_DOS", "Tracepoly");
    std::cout << "1. TracepolyA:" << std::endl;
    for (int ik = 0; ik < nk; ik++)
    {
        std::cout << "ik: " << ik + 1 << std::endl;
        if (nk > 1)
        {
            this->p_hamilt->updateHk(ik);
        }
        const int npw = p_kv->ngk[ik];
        const int nchipk = this->p_stowf->nchip[ik];

        std::complex<double>* pchi;
        if (PARAM.inp.nbands > 0)
        {
            p_stowf->chiortho->fix_k(ik);
            pchi = p_stowf->chiortho->get_pointer();
        }
        else
        {
            p_stowf->chi0->fix_k(ik);
            pchi = p_stowf->chi0->get_pointer();
        }
        auto hchi_norm = std::bind(&hamilt::HamiltSdftPW<std::complex<double>>::hPsi_norm,
                                   p_hamilt_sto,
                                   std::placeholders::_1,
                                   std::placeholders::_2,
                                   std::placeholders::_3);
        if (this->method_sto == 1)
        {
            che.tracepolyA(hchi_norm, pchi, npw, npwx, nchipk);
            for (int i = 0; i < dos_nche; ++i)
            {
                spolyv[i] += che.polytrace[i] * p_kv->wk[ik] / 2;
            }
        }
        else
        {
            int N = dos_nche;
            double kweight = p_kv->wk[ik] / 2;
            char trans = 'T';
            char normal = 'N';
            double one = 1;
            for (int ipart = 0; ipart < npart; ++ipart)
            {
                int nchipk_new = nchipk / npart;
                int start_nchipk = ipart * nchipk_new + nchipk % npart;
                if (ipart < nchipk % npart)
                {
                    nchipk_new++;
                    start_nchipk = ipart * nchipk_new;
                }
                ModuleBase::GlobalFunc::ZEROS(allorderchi.data(), nchipk_new * npwx * dos_nche);
                std::complex<double>* tmpchi = pchi + start_nchipk * npwx;
                che.calpolyvec_complex(hchi_norm, tmpchi, allorderchi.data(), npw, npwx, nchipk_new);
                double* vec_all = (double*)allorderchi.data();
                int LDA = npwx * nchipk_new * 2;
                int M = npwx * nchipk_new * 2;
                dgemm_(&trans, &normal, &N, &N, &M, &kweight, vec_all, &LDA, vec_all, &LDA, &one, spolyv.data(), &N);
            }
        }
    }

    allorderchi.clear();
    allorderchi.shrink_to_fit();

    std::ofstream ofsdos;
    int ndos = int((emax - emin) / de) + 1;
    this->stofunc.sigma = sigmain / ModuleBase::Ry_to_eV;
    ModuleBase::timer::tick("Sto_DOS", "Tracepoly");

    std::cout << "2. Dos:" << std::endl;
    ModuleBase::timer::tick("Sto_DOS", "DOS Loop");
    int n10 = ndos / 10;
    int percent = 10;
    std::vector<double> sto_dos(ndos);
    std::vector<double> ks_dos(ndos);
    std::vector<double> error(ndos);
    for (int ie = 0; ie < ndos; ++ie)
    {
        double tmpks = 0;
        double tmpsto = 0;
        this->stofunc.targ_e = (emin + ie * de) / ModuleBase::Ry_to_eV;
        if (this->method_sto == 1)
        {
            auto ngauss = std::bind(&Sto_Func<double>::ngauss, &this->stofunc, std::placeholders::_1);
            che.calcoef_real(ngauss);
            tmpsto = BlasConnector::dot(dos_nche, che.coef_real, 1, spolyv.data(), 1);
        }
        else
        {
            auto nroot_gauss = std::bind(&Sto_Func<double>::nroot_gauss, &this->stofunc, std::placeholders::_1);
            che.calcoef_real(nroot_gauss);
            tmpsto = vTMv<double, base_device::DEVICE_CPU>(che.coef_real, spolyv.data(), dos_nche);
        }
        if (PARAM.inp.nbands > 0)
        {
            for (int ik = 0; ik < nk; ++ik)
            {
                double* en = &(this->p_elec->ekb(ik, 0));
                for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
                {
                    tmpks += this->stofunc.gauss(en[ib]) * p_kv->wk[ik] / 2;
                }
            }
        }
        tmpks /= GlobalV::NPROC_IN_POOL;

        double tmperror = 0;
        if (this->method_sto == 1)
        {
            tmperror = che.coef_real[dos_nche - 1] * spolyv[dos_nche - 1];
        }
        else
        {
            const int norder = dos_nche;
            double last_coef = che.coef_real[norder - 1];
            double last_spolyv = spolyv[norder * norder - 1];
            tmperror = last_coef
                       * (BlasConnector::dot(norder, che.coef_real, 1, spolyv.data() + norder * (norder - 1), 1)
                          + BlasConnector::dot(norder, che.coef_real, 1, spolyv.data() + norder - 1, norder)
                          - last_coef * last_spolyv);
        }

        if (ie % n10 == n10 - 1)
        {
            std::cout << percent << "%"
                      << " ";
            percent += 10;
        }
        sto_dos[ie] = tmpsto;
        ks_dos[ie] = tmpks;
        error[ie] = tmperror;
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ks_dos.data(), ndos, MPI_DOUBLE, MPI_SUM, STO_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, sto_dos.data(), ndos, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, error.data(), ndos, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (GlobalV::MY_RANK == 0)
    {
        std::string dosfile = PARAM.globalv.global_out_dir + "DOS1_smearing.dat";
        ofsdos.open(dosfile.c_str());
        double maxerror = 0;
        double sum = 0;
        ofsdos << std::setw(8) << "## E(eV) " << std::setw(20) << "dos(eV^-1)" << std::setw(20) << "sum"
               << std::setw(20) << "Error(eV^-1)" << std::endl;
        for (int ie = 0; ie < ndos; ++ie)
        {
            double tmperror = 2.0 * std::abs(error[ie]);
            if (maxerror < tmperror)
            {
                maxerror = tmperror;
            }
            double dos = 2.0 * (ks_dos[ie] + sto_dos[ie]) / ModuleBase::Ry_to_eV;
            sum += dos;
            ofsdos << std::setw(8) << emin + ie * de << std::setw(20) << dos << std::setw(20) << sum * de
                   << std::setw(20) << tmperror << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Finish DOS" << std::endl;
        std::cout << std::scientific << "DOS max absolute Chebyshev Error: " << maxerror << std::endl;
        ofsdos.close();
    }
    ModuleBase::timer::tick("Sto_DOS", "DOS Loop");
    ModuleBase::timer::tick("Sto_DOS", "caldos");
    return;
}
