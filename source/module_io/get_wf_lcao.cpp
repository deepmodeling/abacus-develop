#include "get_wf_lcao.h"

#include "module_parameter/parameter.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/cube_io.h"
#include "module_io/write_wfc_pw.h"
#include "module_io/write_wfc_r.h"
IState_Envelope::IState_Envelope(const elecstate::ElecState* pes)
{
    pes_ = pes;
}

IState_Envelope::~IState_Envelope()
{
}

void IState_Envelope::begin(const psi::Psi<double>* psid,
                            const ModulePW::PW_Basis* rhopw,
                            const ModulePW::PW_Basis_K* wfcpw,
                            const ModulePW::PW_Basis_Big* bigpw,
                            const Parallel_Orbitals& para_orb,
                            Gint_Gamma& gg,
                            const int& out_wfc_pw,
                            const int& out_wfc_r,
                            const K_Vectors& kv,
                            const double nelec,
                            const int nbands_istate,
                            const std::vector<int>& out_band_kb,
                            const int nbands,
                            const int nspin,
                            const int nlocal,
                            const std::string& global_out_dir)
{
    ModuleBase::TITLE("IState_Envelope", "begin");

    std::cout << " Perform |psi(band, r)| for selected bands." << std::endl;

    int mode = 0;
    if (nbands_istate > 0 && static_cast<int>(out_band_kb.size()) == 0)
    {
        mode = 1;
    }
    else if (static_cast<int>(out_band_kb.size()) > 0)
    {
        // If out_band_kb (bands_to_print) is not empty, set mode to 2
        mode = 2;
        std::cout << " Notice: INPUT parameter `nbands_istate` overwritten by `bands_to_print`!" << std::endl;
    }

    int fermi_band = 0;
    int bands_below = 0;
    int bands_above = 0;

    this->bands_picked_.resize(nbands);
    ModuleBase::GlobalFunc::ZEROS(bands_picked_.data(), nbands);

    // (1)
    // mohan update 2011-03-21
    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    std::cout << " number of electrons = " << nelec << std::endl;
    fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);
    std::cout << " number of occupied bands = " << fermi_band << std::endl;

    if (mode == 1)
    {
        bands_below = nbands_istate;
        bands_above = nbands_istate;

        std::cout << " Plot band decomposed charge density below Fermi surface with " << bands_below << " bands."
                  << std::endl;

        std::cout << " Plot band decomposed charge density above Fermi surface with " << bands_above << " bands."
                  << std::endl;

        for (int ib = 0; ib < nbands; ib++)
        {
            if (ib >= fermi_band - bands_below)
            {
                if (ib < fermi_band + bands_above)
                {
                    bands_picked_[ib] = 1;
                }
            }
        }
    }
    else if (mode == 2)
    {
        // Check if length of out_band_kb is valid
        if (static_cast<int>(out_band_kb.size()) > nbands)
        {
            ModuleBase::WARNING_QUIT(
                "IState_Envelope::begin",
                "The number of bands specified by `bands_to_print` in the INPUT file exceeds `nbands`!");
        }
        // Check if all elements in bands_picked_ are 0 or 1
        for (int value: out_band_kb)
        {
            if (value != 0 && value != 1)
            {
                ModuleBase::WARNING_QUIT(
                    "IState_Envelope::begin",
                    "The elements of `bands_to_print` must be either 0 or 1. Invalid values found!");
            }
        }
        // Fill bands_picked_ with values from out_band_kb
        // Remaining bands are already set to 0
        int length = std::min(static_cast<int>(out_band_kb.size()), nbands);
        for (int i = 0; i < length; ++i)
        {
            // out_band_kb rely on function parse_expression
            bands_picked_[i] = out_band_kb[i];
        }

        std::cout << " Plot band decomposed charge density below the Fermi surface: band ";
        for (int i = 0; i + 1 <= fermi_band; ++i)
        {
            if (bands_picked_[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
        std::cout << " Plot band decomposed charge density above the Fermi surface: band ";
        for (int i = fermi_band; i < nbands; ++i)
        {
            if (bands_picked_[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
    }
    else
    {
        ModuleBase::WARNING_QUIT("IState_Envelope::begin", "Invalid mode! Please check the code.");
    }

    // (2) cicle:

    // (2.1) calculate the selected density matrix
    // from wave functions.

    // (2.2) carry out the grid integration to
    // get the charge density.

    // (2.3) output the charge density in .cub format.

    // allocate grid wavefunction for gamma_only
    std::vector<double**> wfc_gamma_grid(nspin);
    for (int is = 0; is < nspin; ++is)
    {
        wfc_gamma_grid[is] = new double*[nbands];
        for (int ib = 0; ib < nbands; ++ib) {
            wfc_gamma_grid[is][ib] = new double[gg.gridt->lgd];
}
    }

    const double mem_size = sizeof(double) * double(gg.gridt->lgd) * double(nbands) * double(nspin) / 1024.0 / 1024.0;
    ModuleBase::Memory::record("IState_Envelope::begin::wfc_gamma_grid", mem_size);
    printf(" Estimated on-the-fly memory consuming by IState_Envelope::begin::wfc_gamma_grid: %f MB\n", mem_size);

    // for pw-wfc in G space
    psi::Psi<std::complex<double>> pw_wfc_g;

    if (out_wfc_pw || out_wfc_r)
    {
        pw_wfc_g.resize(1, nbands, kv.ngk[0]);
    }

    for (int ib = 0; ib < nbands; ib++)
    {
        if (bands_picked_[ib])
        {
            for (int is = 0; is < nspin; ++is) // loop over spin
            {
                std::cout << " Perform envelope function for band " << ib + 1 << std::endl;
                ModuleBase::GlobalFunc::ZEROS(pes_->charge->rho[is], wfcpw->nrxx);

                psid->fix_k(is);
#ifdef __MPI
                wfc_2d_to_grid(psid->get_pointer(), para_orb, wfc_gamma_grid[is], gg.gridt->trace_lo);
#else
                // if not MPI enabled, it is the case psid holds a global matrix. use fix_k to switch between different
                // spin channels (actually kpoints, because now the same kpoint in different spin channels are treated
                // as distinct kpoints)

                for (int i = 0; i < nbands; ++i)
                {
                    for (int j = 0; j < nlocal; ++j)
                        wfc_gamma_grid[is][i][j] = psid[0](i, j);
                }
#endif
                gg.cal_env(wfc_gamma_grid[is][ib], pes_->charge->rho[is], GlobalC::ucell);

                pes_->charge->save_rho_before_sum_band(); // xiaohui add 2014-12-09
                std::stringstream ss;
                ss << global_out_dir << "BAND" << ib + 1 << "_s_" << is + 1 << "_ENV.cube";
                const double ef_tmp = this->pes_->eferm.get_efval(is);
                ModuleIO::write_cube(
#ifdef __MPI
                    bigpw->bz,
                    bigpw->nbz,
                    rhopw->nplane,
                    rhopw->startz_current,
#endif
                    pes_->charge->rho_save[is],
                    is,
                    nspin,
                    0,
                    ss.str(),
                    rhopw->nx,
                    rhopw->ny,
                    rhopw->nz,
                    ef_tmp,
                    &(GlobalC::ucell),
                    3,
                    1);

                if (out_wfc_pw || out_wfc_r) { // only for gamma_only now
                    this->set_pw_wfc(wfcpw, 0, ib, nspin, pes_->charge->rho_save, pw_wfc_g);
}
            }
        }
    }

    if (out_wfc_pw)
    {
        std::stringstream ssw;
        ssw << global_out_dir << "WAVEFUNC";
        std::cout << " write G-space wavefunction into \"" << global_out_dir << "/" << ssw.str() << "\" files."
                  << std::endl;
        ModuleIO::write_wfc_pw(ssw.str(), pw_wfc_g, kv, wfcpw);
    }
    if (out_wfc_r)
    {
        ModuleIO::write_psi_r_1(pw_wfc_g, wfcpw, "wfc_realspace", false, kv);
    }

    for (int is = 0; is < nspin; ++is)
    {
        for (int ib = 0; ib < nbands; ++ib) {
            delete[] wfc_gamma_grid[is][ib];
}
        delete[] wfc_gamma_grid[is];
    }
    return;
}

void IState_Envelope::begin(const psi::Psi<std::complex<double>>* psi,
                            const ModulePW::PW_Basis* rhopw,
                            const ModulePW::PW_Basis_K* wfcpw,
                            const ModulePW::PW_Basis_Big* bigpw,
                            const Parallel_Orbitals& para_orb,
                            Gint_k& gk,
                            const int& out_wf,
                            const int& out_wf_r,
                            const K_Vectors& kv,
                            const double nelec,
                            const int nbands_istate,
                            const std::vector<int>& out_band_kb,
                            const int nbands,
                            const int nspin,
                            const int nlocal,
                            const std::string& global_out_dir)
{
    ModuleBase::TITLE("IState_Envelope", "begin");

    std::cout << " Perform |psi(band, r)| for selected bands." << std::endl;

    int mode = 0;
    if (nbands_istate > 0 && static_cast<int>(out_band_kb.size()) == 0)
    {
        mode = 1;
    }
    else if (static_cast<int>(out_band_kb.size()) > 0)
    {
        // If out_band_kb (bands_to_print) is not empty, set mode to 2
        mode = 2;
        std::cout << " Notice: INPUT parameter `nbands_istate` overwritten by `bands_to_print`!" << std::endl;
    }

    int fermi_band = 0;
    int bands_below = 0;
    int bands_above = 0;

    this->bands_picked_.resize(nbands);
    ModuleBase::GlobalFunc::ZEROS(bands_picked_.data(), nbands);

    // (1)
    // mohan update 2011-03-21
    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    // if NSPIN=4, each band only one electron, fermi_band should be nelec

    std::cout << " number of electrons = " << nelec << std::endl;
    fermi_band = nspin < 4 ? static_cast<int>((nelec + 1) / 2 + 1.0e-8) : nelec;
    std::cout << " number of occupied bands = " << fermi_band << std::endl;

    if (mode == 1)
    {
        bands_below = nbands_istate;
        bands_above = nbands_istate;

        std::cout << " Plot band decomposed charge density below Fermi surface with " << bands_below << " bands."
                  << std::endl;

        std::cout << " Plot band decomposed charge density above Fermi surface with " << bands_above << " bands."
                  << std::endl;

        for (int ib = 0; ib < nbands; ib++)
        {
            if (ib >= fermi_band - bands_below)
            {
                if (ib < fermi_band + bands_above)
                {
                    bands_picked_[ib] = 1;
                }
            }
        }
    }
    else if (mode == 2)
    {
        // Check if length of out_band_kb is valid
        if (static_cast<int>(out_band_kb.size()) > nbands)
        {
            ModuleBase::WARNING_QUIT(
                "IState_Envelope::begin",
                "The number of bands specified by `bands_to_print` in the INPUT file exceeds `nbands`!");
        }
        // Check if all elements in bands_picked_ are 0 or 1
        for (int value: out_band_kb)
        {
            if (value != 0 && value != 1)
            {
                ModuleBase::WARNING_QUIT(
                    "IState_Envelope::begin",
                    "The elements of `bands_to_print` must be either 0 or 1. Invalid values found!");
            }
        }
        // Fill bands_picked_ with values from out_band_kb
        // Remaining bands are already set to 0
        int length = std::min(static_cast<int>(out_band_kb.size()), nbands);
        for (int i = 0; i < length; ++i)
        {
            // out_band_kb rely on function parse_expression
            bands_picked_[i] = out_band_kb[i];
        }

        std::cout << " Plot band decomposed charge density below the Fermi surface: band ";
        for (int i = 0; i + 1 <= fermi_band; ++i)
        {
            if (bands_picked_[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
        std::cout << " Plot band decomposed charge density above the Fermi surface: band ";
        for (int i = fermi_band; i < nbands; ++i)
        {
            if (bands_picked_[i] == 1)
            {
                std::cout << i + 1 << " ";
            }
        }
        std::cout << std::endl;
    }
    else
    {
        ModuleBase::WARNING_QUIT("IState_Envelope::begin", "Invalid mode! Please check the code.");
    }

    // (2) cicle:

    // (2.1) calculate the selected density matrix
    // from wave functions.

    // (2.2) carry out the grid integration to
    // get the charge density.

    // (2.3) output the charge density in .cub format.

    // allocate grid wavefunction for gamma_only
    const int nks = kv.get_nks();
    std::vector<std::complex<double>**> wfc_k_grid(nks);
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
    ModuleBase::Memory::record("IState_Envelope::begin::wfc_k_grid", mem_size);
    printf(" Estimated on-the-fly memory consuming by IState_Envelope::begin::wfc_k_grid: %f MB\n", mem_size);

    // for pw-wfc in G space
    psi::Psi<std::complex<double>> pw_wfc_g(kv.ngk.data());

    if (out_wf || out_wf_r)
    {
        pw_wfc_g.resize(nks, nbands, wfcpw->npwk_max);
    }

    for (int ib = 0; ib < nbands; ib++)
    {
        if (bands_picked_[ib])
        {
            const int nspin0 = (nspin == 2) ? 2 : 1;
            for (int ik = 0; ik < nks; ++ik) // the loop of nspin0 is included
            {
                const int ispin = kv.isk[ik];
                ModuleBase::GlobalFunc::ZEROS(pes_->charge->rho[ispin],
                                              wfcpw->nrxx); // terrible, you make changes on another instance's data???
                std::cout << " Perform envelope function for kpoint " << ik << ",  band" << ib + 1 << std::endl;
                //  2d-to-grid conversion is unified into `wfc_2d_to_grid`.
                psi->fix_k(ik);
#ifdef __MPI // need to deal with NSPIN=4 !!!!
                wfc_2d_to_grid(psi->get_pointer(), para_orb, wfc_k_grid[ik], gk.gridt->trace_lo);
#else
                for (int i = 0; i < nbands; ++i)
                {
                    for (int j = 0; j < nlocal; ++j)
                        wfc_k_grid[ik][i][j] = psi[0](i, j);
                }
#endif
                // deal with NSPIN=4
                gk.cal_env_k(ik, wfc_k_grid[ik][ib], pes_->charge->rho[ispin], kv.kvec_c, kv.kvec_d, GlobalC::ucell);

                std::stringstream ss;
                ss << global_out_dir << "BAND" << ib + 1 << "_k_" << ik / nspin0 + 1 << "_s_" << ispin + 1
                   << "_ENV.cube";
                const double ef_tmp = this->pes_->eferm.get_efval(ispin);

                ModuleIO::write_cube(
#ifdef __MPI
                    bigpw->bz,
                    bigpw->nbz,
                    rhopw->nplane,
                    rhopw->startz_current,
#endif
                    pes_->charge->rho[ispin],
                    ispin,
                    nspin,
                    0,
                    ss.str(),
                    rhopw->nx,
                    rhopw->ny,
                    rhopw->nz,
                    ef_tmp,
                    &(GlobalC::ucell),
                    3,
                    1);

                if (out_wf || out_wf_r) // only for gamma_only now
                {
                    pw_wfc_g.fix_k(ik);
                    this->set_pw_wfc(wfcpw, ik, ib, nspin, pes_->charge->rho, pw_wfc_g);
                }
            }
        }
    }

    if (out_wf || out_wf_r)
    {
        if (out_wf)
        {
            std::stringstream ssw;
            ssw << global_out_dir << "WAVEFUNC";
            std::cout << " write G-space wavefunction into \"" << global_out_dir << "/" << ssw.str() << "\" files."
                      << std::endl;
            ModuleIO::write_wfc_pw(ssw.str(), pw_wfc_g, kv, wfcpw);
        }
        if (out_wf_r)
        {
            ModuleIO::write_psi_r_1(pw_wfc_g, wfcpw, "wfc_realspace", false, kv);
        }
    }

    for (int ik = 0; ik < nks; ++ik)
    {
        for (int ib = 0; ib < nbands; ++ib) {
            delete[] wfc_k_grid[ik][ib];
}
        delete[] wfc_k_grid[ik];
    }

    return;
}

// for each band
void IState_Envelope::set_pw_wfc(const ModulePW::PW_Basis_K* wfcpw,
                                 const int& ik,
                                 const int& ib,
                                 const int& nspin,
                                 const double* const* const rho,
                                 psi::Psi<std::complex<double>>& wfc_g)
{
    if (ib == 0) { // once is enough
        ModuleBase::TITLE("IState_Envelope", "set_pw_wfc");
}

    std::vector<std::complex<double>> Porter(wfcpw->nrxx);
    // here I refer to v_hartree, but I don't know how to deal with NSPIN=4
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; is++) {
        for (int ir = 0; ir < wfcpw->nrxx; ir++) {
            Porter[ir] += std::complex<double>(rho[is][ir], 0.0);
}
}

    // call FFT
    wfcpw->real2recip(Porter.data(), &wfc_g(ib, 0), ik);
}

#ifdef __MPI
template <typename T>
int IState_Envelope::set_wfc_grid(const int naroc[2],
                                  const int nb,
                                  const int dim0,
                                  const int dim1,
                                  const int iprow,
                                  const int ipcol,
                                  const T* in,
                                  T** out,
                                  const std::vector<int>& trace_lo)
{
    ModuleBase::TITLE(" Local_Orbital_wfc", "set_wfc_grid");
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

template int IState_Envelope::set_wfc_grid(const int naroc[2],
                                           const int nb,
                                           const int dim0,
                                           const int dim1,
                                           const int iprow,
                                           const int ipcol,
                                           const double* in,
                                           double** out,
                                           const std::vector<int>& trace_lo);
template int IState_Envelope::set_wfc_grid(const int naroc[2],
                                           const int nb,
                                           const int dim0,
                                           const int dim1,
                                           const int iprow,
                                           const int ipcol,
                                           const std::complex<double>* in,
                                           std::complex<double>** out,
                                           const std::vector<int>& trace_lo);

template <typename T>
void IState_Envelope::wfc_2d_to_grid(const T* lowf_2d,
                                     const Parallel_Orbitals& pv,
                                     T** lowf_grid,
                                     const std::vector<int>& trace_lo)
{
    ModuleBase::TITLE(" Local_Orbital_wfc", "wfc_2d_to_grid");
    ModuleBase::timer::tick("Local_Orbital_wfc", "wfc_2d_to_grid");

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
    ModuleBase::timer::tick("Local_Orbital_wfc", "wfc_2d_to_grid");
}

template void IState_Envelope::wfc_2d_to_grid(const double* lowf_2d,
                                              const Parallel_Orbitals& pv,
                                              double** lowf_grid,
                                              const std::vector<int>& trace_lo);
template void IState_Envelope::wfc_2d_to_grid(const std::complex<double>* lowf_2d,
                                              const Parallel_Orbitals& pv,
                                              std::complex<double>** lowf_grid,
                                              const std::vector<int>& trace_lo);
#endif

int IState_Envelope::globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

int IState_Envelope::localIndex(int globalindex, int nblk, int nprocs, int& myproc)
{
    myproc = int((globalindex % (nblk * nprocs)) / nblk);
    return int(globalindex / (nblk * nprocs)) * nblk + globalindex % nblk;
}
