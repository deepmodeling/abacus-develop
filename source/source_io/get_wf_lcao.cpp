#include "get_wf_lcao.h"

#include "source_io/cube_io.h"
#include "source_io/write_wfc_pw.h"
#include "source_base/memory.h"

#ifndef __OLD_GINT
#include "source_lcao/module_gint/temp_gint/gint_env_gamma.h"
#include "source_lcao/module_gint/temp_gint/gint_env_k.h"
#endif

Get_wf_lcao::Get_wf_lcao(const elecstate::ElecState* pes)
{
    pes_ = pes;
}

Get_wf_lcao::~Get_wf_lcao()
{
}

// For gamma_only
void Get_wf_lcao::begin(const UnitCell& ucell,
                        const psi::Psi<double>* psid,
                        const ModulePW::PW_Basis_K* pw_wfc,
                        const Parallel_Grid& pgrid,
                        const Parallel_Orbitals& para_orb,
                        Gint_Gamma& gg,
                        const int& out_wfc_pw,
                        const K_Vectors& kv,
                        const double nelec,
                        const std::vector<int>& out_wfc_norm,
                        const std::vector<int>& out_wfc_re_im,
                        const int nbands,
                        const int nspin,
                        const int nlocal,
                        const std::string& global_out_dir,
                        std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Get_wf_lcao", "begin");

    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    const int fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);

    prepare_get_wf(ofs_running);

#ifdef __OLD_GINT
    // allocate grid wave functions for gamma_only
    std::vector<double**> wfc_gamma_grid(nspin);
    for (int is = 0; is < nspin; ++is)
    {
        wfc_gamma_grid[is] = new double*[nbands];
        for (int ib = 0; ib < nbands; ++ib)
        {
            wfc_gamma_grid[is][ib] = new double[gg.gridt->lgd];
        }
    }
#endif

    // for pw_wfc in G space
    psi::Psi<std::complex<double>> psi_g;

    // if (out_wfc_pw || out_wfc_r)
    psi_g.resize(nspin, nbands, kv.ngk[0]);

#ifdef __OLD_GINT
    const double mem_size = sizeof(double) * double(gg.gridt->lgd) * double(nbands) * double(nspin) / 1024.0 / 1024.0;
    ModuleBase::Memory::record("Get_wf_lcao::begin", mem_size);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "On-the-fly memory consumption (MB)", mem_size);
#endif

    // Set this->bands_picked_
    this->select_bands(out_wfc_norm, nbands, fermi_band);

    // Calculate out_wfc_norm
    for (int is = 0; is < nspin; ++is)
    {
        psid->fix_k(is);
#ifdef __OLD_GINT
    #ifdef __MPI
        wfc_2d_to_grid(psid->get_pointer(), para_orb, wfc_gamma_grid[is], gg.gridt->trace_lo);
    #else
        // if not MPI enabled, it is the case psid holds a global matrix. 
        // use fix_k to switch between different spin channels (actually kpoints, 
        // because now the same kpoint in different spin channels are treated
        // as distinct kpoints)
        for (int i = 0; i < nbands; ++i)
        {
            for (int j = 0; j < nlocal; ++j)
            {
                wfc_gamma_grid[is][i][j] = psid[0](i, j);
            }
        }
    #endif
#else
        ModuleGint::Gint_env_gamma gint_env(psid->get_pointer(), &para_orb, nbands, nlocal, pes_->charge->rho[is]);
#endif
        for (int ib = 0; ib < nbands; ++ib)
        {
            if (bands_picked_[ib])
            {
            #ifdef __OLD_GINT
                ModuleBase::GlobalFunc::ZEROS(pes_->charge->rho[is], pw_wfc->nrxx);
                gg.cal_env(wfc_gamma_grid[is][ib], pes_->charge->rho[is], ucell);
            #else
                gint_env.cal_env_band(ib);
            #endif
                pes_->charge->save_rho_before_sum_band();

                // pint out information
                std::stringstream ss_file;
                ss_file << "wfi" << ib + 1 << "s" << is + 1 << ".cube";

                std::stringstream ss_out;
                ss_out << global_out_dir << ss_file.str();

                std::stringstream ss_info;
                ss_info << "Wave func. " << ib + 1 << " spin " << is + 1 << " saved in";

                ModuleBase::GlobalFunc::OUT(ofs_running, ss_info.str(), ss_file.str());

                const double ef_tmp = this->pes_->eferm.get_efval(is);
                ModuleIO::write_vdata_palgrid(pgrid,
                                              pes_->charge->rho_save[is],
                                              is,
                                              nspin,
                                              0,
                                              ss_out.str(),
                                              ef_tmp,
                                              &(ucell));
            }
        }
    }

    // Set this->bands_picked_
    this->select_bands(out_wfc_re_im, nbands, fermi_band);

    // Calculate out_wfc_re_im
    for (int is = 0; is < nspin; ++is)
    {
        psid->fix_k(is);
#ifdef __OLD_GINT
    #ifdef __MPI
        wfc_2d_to_grid(psid->get_pointer(), para_orb, wfc_gamma_grid[is], gg.gridt->trace_lo);
    #else
        // if not MPI enabled, it is the case psid holds a global matrix. use fix_k to switch between
        // different spin channels (actually kpoints, because now the same kpoint in different spin channels
        // are treated as distinct kpoints)
        for (int i = 0; i < nbands; ++i)
        {
            for (int j = 0; j < nlocal; ++j)
            {
                wfc_gamma_grid[is][i][j] = psid[0](i, j);
            }
        }
    #endif
#else
        ModuleGint::Gint_env_gamma gint_env(psid->get_pointer(), &para_orb, nbands, nlocal, pes_->charge->rho[is]);
#endif
        for (int ib = 0; ib < nbands; ++ib)
        {
            if (bands_picked_[ib])
            {
#ifdef __OLD_GINT
                ModuleBase::GlobalFunc::ZEROS(pes_->charge->rho[is], pw_wfc->nrxx);
                gg.cal_env(wfc_gamma_grid[is][ib], pes_->charge->rho[is], ucell);
#else
                gint_env.cal_env_band(ib);
#endif
                pes_->charge->save_rho_before_sum_band();

                const double ef_tmp = this->pes_->eferm.get_efval(is);

                // only for gamma_only now
                psi_g.fix_k(is);
                this->set_pw_wfc(pw_wfc, is, ib, nspin, pes_->charge->rho, psi_g);

                // Calculate real-space wave functions
                psi_g.fix_k(is);
                std::vector<std::complex<double>> wfc_r(pw_wfc->nrxx);
                pw_wfc->recip2real(&psi_g(ib, 0), wfc_r.data(), is);

                // Extract real and imaginary parts
                std::vector<double> wfc_real(pw_wfc->nrxx);
                std::vector<double> wfc_imag(pw_wfc->nrxx);
                for (int ir = 0; ir < pw_wfc->nrxx; ++ir)
                {
                    wfc_real[ir] = wfc_r[ir].real();
                    wfc_imag[ir] = wfc_r[ir].imag();
                }

                // Output real part
                std::stringstream ss_real;
                ss_real << global_out_dir << "wfi" << ib + 1 << "s" << is + 1 << "re.cube";
                ModuleIO::write_vdata_palgrid(pgrid, wfc_real.data(), is, nspin, 0, ss_real.str(), ef_tmp, &(ucell));

                // Output imaginary part
                std::stringstream ss_imag;
                ss_imag << global_out_dir << "wfi" << ib + 1 << "s" << is + 1 << "im.cube";
                ModuleIO::write_vdata_palgrid(pgrid, wfc_imag.data(), is, nspin, 0, ss_imag.str(), ef_tmp, &(ucell));
            }
        }
    }


    const int istep = -1; // -1 means ionic iteration number will not appear in file name
    const int iter = -1; // -1 means electronic iteration number will not appear in file name
    ModuleIO::write_wfc_pw(istep, iter, GlobalV::KPAR,
                           GlobalV::MY_POOL,
                           GlobalV::MY_RANK,
                           nbands,
                           nspin,
                           PARAM.globalv.npol,
                           GlobalV::RANK_IN_POOL,
                           GlobalV::NPROC_IN_POOL,
                           out_wfc_pw,
                           PARAM.inp.ecutwfc,
                           global_out_dir,
                           psi_g,
                           kv,
                           pw_wfc,
                           ofs_running);

#ifdef __OLD_GINT
    for (int is = 0; is < nspin; ++is)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            delete[] wfc_gamma_grid[is][ib];
        }
        delete[] wfc_gamma_grid[is];
    }
#endif
    return;
}

// For multi-k
void Get_wf_lcao::begin(const UnitCell& ucell,
                        const psi::Psi<std::complex<double>>* psi,
                        const ModulePW::PW_Basis_K* pw_wfc,
                        const Parallel_Grid& pgrid,
                        const Parallel_Orbitals& para_orb,
                        Gint_k& gk,
                        const int& out_wfc_pw,
                        const K_Vectors& kv,
                        const double nelec,
                        const std::vector<int>& out_wfc_norm,
                        const std::vector<int>& out_wfc_re_im,
                        const int nbands,
                        const int nspin,
                        const int nlocal,
                        const std::string& global_out_dir,
                        std::ofstream& ofs_running)
{
    ModuleBase::TITLE("Get_wf_lcao", "begin");

    const int fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);

    prepare_get_wf(ofs_running);

    // allocate grid wave functions for multi-k
    const int nks = kv.get_nks();
    std::vector<std::complex<double>**> wfc_k_grid(nks);
#ifdef __OLD_GINT
    for (int ik = 0; ik < nks; ++ik)
    {
        wfc_k_grid[ik] = new std::complex<double>*[nbands];
        for (int ib = 0; ib < nbands; ++ib)
        {
            wfc_k_grid[ik][ib] = new std::complex<double>[gk.gridt->lgd];
        }
    }

    const double mem_size
        = sizeof(std::complex<double>) * double(gk.gridt->lgd) * double(nbands) * double(nks) / 1024.0 / 1024.0;
    ModuleBase::Memory::record("Get_wf_lcao::begin", mem_size);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "On-the-fly memory consumption (MB)", mem_size);
#endif

    // for pw_wfc in G space
    psi::Psi<std::complex<double>> psi_g;

    // if (out_wfc_pw || out_wf_r)
    psi_g.resize(nks, nbands, pw_wfc->npwk_max);

    // Set this->bands_picked_
    this->select_bands(out_wfc_norm, nbands, fermi_band);

   // Calculate out_wfc_norm
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int ik = 0; ik < nks; ++ik) // the loop of nspin0 is included
    {
        const int ispin = kv.isk[ik];
        //  2d-to-grid conversion is unified into `wfc_2d_to_grid`.
        psi->fix_k(ik);

#ifdef __OLD_GINT
    #ifdef __MPI // need to deal with NSPIN=4 !!!!
        wfc_2d_to_grid(psi->get_pointer(), para_orb, wfc_k_grid[ik], gk.gridt->trace_lo);
    #else
        for (int i = 0; i < nbands; ++i)
        {
            for (int j = 0; j < nlocal; ++j)
            {
                wfc_k_grid[ik][i][j] = psi[0](i, j);
            }
        }
    #endif
#else
        ModuleGint::Gint_env_k gint_env(psi->get_pointer(), &para_orb, kv.kvec_c, kv.kvec_d,
                                        nbands, nlocal, ik, PARAM.inp.nspin, PARAM.globalv.npol, pes_->charge->rho[ispin]);
#endif
        
        for (int ib = 0; ib < nbands; ++ib)
        {
            if (bands_picked_[ib])
            {
#ifdef __OLD_GINT
                ModuleBase::GlobalFunc::ZEROS(pes_->charge->rho[ispin],
                                              pw_wfc->nrxx); // terrible, you make changes on another instance's data???

                // deal with NSPIN=4
                gk.cal_env_k(ik, wfc_k_grid[ik][ib], pes_->charge->rho[ispin], kv.kvec_c, kv.kvec_d, ucell);
#else
                gint_env.cal_env_band(ib);
#endif

                // ik0 is the real k-point index, starting from 0
                int ik0 = kv.ik2iktot[ik];
                if (nspin == 2)
                {
                    const int half_k = kv.get_nkstot() / 2;
                    if (ik0 >= half_k)
                    {
                        ik0 -= half_k;
                    }
                }

                // pint out information
                std::stringstream ss_file;
                ss_file << "wfi" << ib + 1 << "s" << ispin + 1 << "k" << ik0 + 1 << ".cube";

                std::stringstream ss_out;
                ss_out << global_out_dir << ss_file.str();

                std::stringstream ss_info;
                ss_info << "Wave func. " << ib + 1 << " spin " << ispin + 1 << " k-point " << ik0 + 1 << " saved in";

                ModuleBase::GlobalFunc::OUT(ofs_running, ss_info.str(), ss_file.str());

                const double ef_tmp = this->pes_->eferm.get_efval(ispin);

                ModuleIO::write_vdata_palgrid(pgrid,
                                              pes_->charge->rho[ispin],
                                              ispin,
                                              nspin,
                                              0,
                                              ss_out.str(),
                                              ef_tmp,
                                              &(ucell),
                                              3,
                                              1);

                // if (out_wfc_pw || out_wf_r)
                psi_g.fix_k(ik);
                this->set_pw_wfc(pw_wfc, ik, ib, nspin, pes_->charge->rho, psi_g);
            }
        }
    }

    const int istep = -1; // -1 means ionic iteration number will not appear in file name
    const int iter = -1; // -1 means electronic iteration number will not appear in file name
    ModuleIO::write_wfc_pw(istep, iter, GlobalV::KPAR,
                           GlobalV::MY_POOL,
                           GlobalV::MY_RANK,
                           nbands,
                           nspin,
                           PARAM.globalv.npol,
                           GlobalV::RANK_IN_POOL,
                           GlobalV::NPROC_IN_POOL,
                           out_wfc_pw,
                           PARAM.inp.ecutwfc,
                           global_out_dir,
                           psi_g,
                           kv,
                           pw_wfc,
                           ofs_running);

    // Set this->bands_picked_
    this->select_bands(out_wfc_re_im, nbands, fermi_band);

    // Calculate out_wfc_re_im
    for (int ib = 0; ib < nbands; ++ib)
    {
        if (bands_picked_[ib])
        {
            const int nspin0 = (nspin == 2) ? 2 : 1;
            for (int ik = 0; ik < nks; ++ik)
            {
                const int ispin = kv.isk[ik];

                psi_g.fix_k(ik);

                // Calculate real-space wave functions
                std::vector<std::complex<double>> wfc_r(pw_wfc->nrxx);
                pw_wfc->recip2real(&psi_g(ib, 0), wfc_r.data(), ik);

                // Extract real and imaginary parts
                std::vector<double> wfc_real(pw_wfc->nrxx);
                std::vector<double> wfc_imag(pw_wfc->nrxx);
                for (int ir = 0; ir < pw_wfc->nrxx; ++ir)
                {
                    wfc_real[ir] = wfc_r[ir].real();
                    wfc_imag[ir] = wfc_r[ir].imag();
                }

                // ik0 is the real k-point index, starting from 0
                int ik0 = kv.ik2iktot[ik];
                if (nspin == 2)
                {
                    const int half_k = kv.get_nkstot() / 2;
                    if (ik0 >= half_k)
                    {
                        ik0 -= half_k;
                    }
                }

                // Output real part
                std::stringstream ss_real;
                ss_real << global_out_dir << "wfi" << ib + 1 << "s" << ispin + 1 << "k" << ik0 + 1 << "re.cube";

                const double ef_tmp = this->pes_->eferm.get_efval(ispin);
                ModuleIO::write_vdata_palgrid(pgrid, wfc_real.data(), ispin, nspin, 0, ss_real.str(), ef_tmp, &(ucell));

                // Output imaginary part
                std::stringstream ss_imag;
                ss_imag << global_out_dir << "wfi" << ib + 1 << "s" << ispin + 1 << "k" << ik0 + 1 << "im.cube";
                ModuleIO::write_vdata_palgrid(pgrid, wfc_imag.data(), ispin, nspin, 0, ss_imag.str(), ef_tmp, &(ucell));
            }
        }
    }
#ifdef __OLD_GINT
    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib)
        {
            delete[] wfc_k_grid[ik][ib];
        }
        delete[] wfc_k_grid[ik];
    }
#endif
    return;
}

void Get_wf_lcao::select_bands(const std::vector<int>& out_wfc_kb, const int nbands, const int fermi_band)
{
    ModuleBase::TITLE("Get_wf_lcao", "select_bands");

    this->bands_picked_.resize(nbands);
    ModuleBase::GlobalFunc::ZEROS(bands_picked_.data(), nbands);

    // Select bands directly using parameter `out_wfc_norm` or `out_wfc_re_im`
    // Check if length of out_wfc_kb is valid
    if (static_cast<int>(out_wfc_kb.size()) > nbands)
    {
        ModuleBase::WARNING_QUIT("Get_wf_lcao::select_bands",
                                 "The number of bands specified by `out_wfc_norm` or `out_wfc_re_im` in the INPUT "
                                 "file exceeds `nbands`!");
    }
    // Check if all elements in out_wfc_kb are 0 or 1
    for (int value: out_wfc_kb)
    {
        if (value != 0 && value != 1)
        {
            ModuleBase::WARNING_QUIT(
                "Get_wf_lcao::select_bands",
                "The elements of `out_wfc_norm` or `out_wfc_re_im` must be either 0 or 1. Invalid values found!");
        }
    }
    // Fill bands_picked_ with values from out_wfc_kb
    // Remaining bands are already set to 0
    const int length = std::min(static_cast<int>(out_wfc_kb.size()), nbands);
    std::copy(out_wfc_kb.begin(), out_wfc_kb.begin() + length, bands_picked_.begin());

    // Check if there are selected bands below the Fermi surface
    bool has_below = false;
    for (int i = 0; i + 1 <= fermi_band; ++i)
    {
        if (bands_picked_[i] == 1)
        {
            has_below = true;
            break;
        }
    }
    if (has_below)
    {
        std::cout << " Plot wave functions below the Fermi surface: band ";
        for (int i = 0; i + 1 <= fermi_band; ++i)
        {
            if (bands_picked_[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
    }

    // Check if there are selected bands above the Fermi surface
    bool has_above = false;
    for (int i = fermi_band; i < nbands; ++i)
    {
        if (bands_picked_[i] == 1)
        {
            has_above = true;
            break;
        }
    }
    if (has_above)
    {
        std::cout << " Plot wave functions above the Fermi surface: band ";
        for (int i = fermi_band; i < nbands; ++i)
        {
            if (bands_picked_[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
    }
}

// for each band
void Get_wf_lcao::set_pw_wfc(const ModulePW::PW_Basis_K* pw_wfc,
                             const int& ik,
                             const int& ib,
                             const int& nspin,
                             const double* const* const rho,
                             psi::Psi<std::complex<double>>& wfc_g)
{
    if (ib == 0)
    {
        // once is enough
        ModuleBase::TITLE("Get_wf_lcao", "set_pw_wfc");
    }

    std::vector<std::complex<double>> Porter(pw_wfc->nrxx);
    // here I refer to v_hartree, but I don't know how to deal with NSPIN=4
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; ++is)
    {
        for (int ir = 0; ir < pw_wfc->nrxx; ++ir)
        {
            Porter[ir] += std::complex<double>(rho[is][ir], 0.0);
        }
    }

    // call FFT
    pw_wfc->real2recip(Porter.data(), &wfc_g(ib, 0), ik);
}

#ifdef __MPI
template <typename T>
int Get_wf_lcao::set_wfc_grid(const int naroc[2],
                              const int nb,
                              const int dim0,
                              const int dim1,
                              const int iprow,
                              const int ipcol,
                              const T* in,
                              T** out,
                              const std::vector<int>& trace_lo)
{
    ModuleBase::TITLE("Get_wf_lcao", "set_wfc_grid");
    if (!out)
    {
        return 0;
    }
    for (int j = 0; j < naroc[1]; ++j)
    {
        int igcol = globalIndex(j, nb, dim1, ipcol);
        if (igcol >= PARAM.inp.nbands)
        {
            continue;
        }
        for (int i = 0; i < naroc[0]; ++i)
        {
            int igrow = globalIndex(i, nb, dim0, iprow);
            int mu_local = trace_lo[igrow];
            if (out && mu_local >= 0)
            {
                out[igcol][mu_local] = in[j * naroc[0] + i];
            }
        }
    }
    return 0;
}

template int Get_wf_lcao::set_wfc_grid(const int naroc[2],
                                       const int nb,
                                       const int dim0,
                                       const int dim1,
                                       const int iprow,
                                       const int ipcol,
                                       const double* in,
                                       double** out,
                                       const std::vector<int>& trace_lo);
template int Get_wf_lcao::set_wfc_grid(const int naroc[2],
                                       const int nb,
                                       const int dim0,
                                       const int dim1,
                                       const int iprow,
                                       const int ipcol,
                                       const std::complex<double>* in,
                                       std::complex<double>** out,
                                       const std::vector<int>& trace_lo);

template <typename T>
void Get_wf_lcao::wfc_2d_to_grid(const T* lowf_2d,
                                 const Parallel_Orbitals& pv,
                                 T** lowf_grid,
                                 const std::vector<int>& trace_lo)
{
    ModuleBase::TITLE("Get_wf_lcao", "wfc_2d_to_grid");
    ModuleBase::timer::tick("Get_wf_lcao", "wfc_2d_to_grid");

    // dimension related
    const int nlocal = pv.desc_wfc[2];
    const int nbands = pv.desc_wfc[3];

    // MPI and memory related
    const int mem_stride = 1;
    int mpi_info = 0;

    // get the rank of the current process
    int rank = 0;
    MPI_Comm_rank(pv.comm(), &rank);

    // calculate the maximum number of nlocal over all processes in pv.comm() range
    long buf_size;
    mpi_info = MPI_Reduce(&pv.nloc_wfc, &buf_size, 1, MPI_LONG, MPI_MAX, 0, pv.comm());
    mpi_info = MPI_Bcast(&buf_size, 1, MPI_LONG, 0, pv.comm()); // get and then broadcast
    std::vector<T> lowf_block(buf_size);

    // this quantity seems to have the value returned by function numroc_ in ScaLAPACK?
    int naroc[2];

    // for BLACS broadcast
    char scope = 'A';
    char top = ' ';

    // loop over all processors
    for (int iprow = 0; iprow < pv.dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < pv.dim1; ++ipcol)
        {
            if (iprow == pv.coord[0] && ipcol == pv.coord[1])
            {
                BlasConnector::copy(pv.nloc_wfc, lowf_2d, mem_stride, lowf_block.data(), mem_stride);
                naroc[0] = pv.nrow;
                naroc[1] = pv.ncol_bands;
                Cxgebs2d(pv.blacs_ctxt, &scope, &top, 2, 1, naroc, 2);
                Cxgebs2d(pv.blacs_ctxt, &scope, &top, buf_size, 1, lowf_block.data(), buf_size);
            }
            else
            {
                Cxgebr2d(pv.blacs_ctxt, &scope, &top, 2, 1, naroc, 2, iprow, ipcol);
                Cxgebr2d(pv.blacs_ctxt, &scope, &top, buf_size, 1, lowf_block.data(), buf_size, iprow, ipcol);
            }

            // then use it to set the wfc_grid.
            mpi_info = this->set_wfc_grid(naroc,
                                          pv.nb,
                                          pv.dim0,
                                          pv.dim1,
                                          iprow,
                                          ipcol,
                                          lowf_block.data(),
                                          lowf_grid,
                                          trace_lo);
            // this operation will let all processors have the same wfc_grid
        }
    }
    ModuleBase::timer::tick("Get_wf_lcao", "wfc_2d_to_grid");
}

template void Get_wf_lcao::wfc_2d_to_grid(const double* lowf_2d,
                                          const Parallel_Orbitals& pv,
                                          double** lowf_grid,
                                          const std::vector<int>& trace_lo);
template void Get_wf_lcao::wfc_2d_to_grid(const std::complex<double>* lowf_2d,
                                          const Parallel_Orbitals& pv,
                                          std::complex<double>** lowf_grid,
                                          const std::vector<int>& trace_lo);
#endif

void Get_wf_lcao::prepare_get_wf(std::ofstream& ofs_running)
{
    ofs_running << "\n\n";
    ofs_running << " GET_WF CALCULATION BEGINS" << std::endl;

    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs_running << " |                                                                    |" << std::endl;
    ofs_running << " | Here we use real-space (r) grid integral technique to calculate    |" << std::endl;
    ofs_running << " | the electronic wave function psi(i,r) for each electronic state i. |" << std::endl;
    ofs_running << " | The |psi(i,r)|, Re[psi(i,r)], Im[psi(i,r)] are printed out using   |" << std::endl;
    ofs_running << " | numerical atomic orbitals as basis set.                            |" << std::endl;
    ofs_running << " |                                                                    |" << std::endl;
    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    ofs_running << "\n\n";

    ofs_running << std::setprecision(6);
}

int Get_wf_lcao::globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    const int iblock = localindex / nblk;
    const int gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

int Get_wf_lcao::localIndex(int globalindex, int nblk, int nprocs, int& myproc)
{
    myproc = int((globalindex % (nblk * nprocs)) / nblk);
    return int(globalindex / (nblk * nprocs)) * nblk + globalindex % nblk;
}
