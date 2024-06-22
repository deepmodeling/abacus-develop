#include "istate_envelope.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/rho_io.h"
#include "module_io/write_wfc_pw.h"
#include "module_io/write_wfc_r.h"
IState_Envelope::IState_Envelope(const elecstate::ElecState* pes_in)
{
    pes = pes_in;
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
                            int& out_wfc_pw,
                            int& out_wfc_r,
                            const K_Vectors& kv,
                            const double nelec,
                            const int nbands_istate,
                            const int nbands,
                            const int nspin,
                            const int nlocal,
                            const std::string& global_out_dir)
{
    ModuleBase::TITLE("IState_Envelope", "begin");

    std::cout << " perform |psi(band, r)| for selected bands." << std::endl;

    // (1)
    // mohan update 2011-03-21
    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    int fermi_band = static_cast<int>((nelec + 1) / 2 + 1.0e-8);
    int bands_below = nbands_istate;
    int bands_above = nbands_istate;

    std::cout << " number of electrons = " << nelec << std::endl;
    std::cout << " number of occupied bands = " << fermi_band << std::endl;
    std::cout << " plot band decomposed charge density below fermi surface with " << bands_below << " bands."
              << std::endl;

    std::cout << " plot band decomposed charge density above fermi surface with " << bands_above << " bands."
              << std::endl;

    // (2) cicle:

    // (2.1) calculate the selected density matrix
    // from wave functions.

    // (2.2) carry out the grid integration to
    // get the charge density.

    // (2.3) output the charge density in .cub format.
    this->bands_picked = new bool[nbands];
    ModuleBase::GlobalFunc::ZEROS(bands_picked, nbands);
    for (int ib = 0; ib < nbands; ib++)
    {
        if (ib >= fermi_band - bands_below)
        {
            if (ib < fermi_band + bands_above)
            {
                bands_picked[ib] = true;
            }
        }
    }

    // allocate grid wavefunction for gamma_only
    std::vector<double**> wfc_gamma_grid(nspin);
    for (int is = 0; is < nspin; ++is)
    {
        wfc_gamma_grid[is] = new double*[nbands];
        for (int ib = 0; ib < nbands; ++ib)
            wfc_gamma_grid[is][ib] = new double[gg.gridt->lgd];
    }
    const size_t mem_size = sizeof(double) * gg.gridt->lgd * nbands * nspin / 1024 / 1024;
    ModuleBase::Memory::record("IState_Envelope::begin::wfc_gamma_grid", mem_size);
    printf("Estimated on-the-fly memory consuming by IState_Envelope::begin::wfc_gamma_grid: %ld MB\n", mem_size);
    // for pw-wfc in G space
    psi::Psi<std::complex<double>> pw_wfc_g;

    if (out_wfc_pw || out_wfc_r)
    {
        pw_wfc_g.resize(1, nbands, kv.ngk[0]);
    }

    for (int ib = 0; ib < nbands; ib++)
    {
        if (bands_picked[ib])
        {
            for (int is = 0; is < nspin; ++is) // loop over spin
            {
                std::cout << " Perform envelope function for band " << ib + 1 << std::endl;
                ModuleBase::GlobalFunc::ZEROS(pes->charge->rho[is], wfcpw->nrxx);

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
                gg.cal_env(wfc_gamma_grid[is][ib], pes->charge->rho[is], GlobalC::ucell);

                pes->charge->save_rho_before_sum_band(); // xiaohui add 2014-12-09
                std::stringstream ss;
                ss << global_out_dir << "BAND" << ib + 1 << "_s_" << is + 1 << "_ENV.cube";
                const double ef_tmp = this->pes->eferm.get_efval(is);
                ModuleIO::write_rho(
#ifdef __MPI
                    bigpw->bz,
                    bigpw->nbz,
                    rhopw->nplane,
                    rhopw->startz_current,
#endif
                    pes->charge->rho_save[is],
                    is,
                    nspin,
                    0,
                    ss.str(),
                    rhopw->nx,
                    rhopw->ny,
                    rhopw->nz,
                    ef_tmp,
                    &(GlobalC::ucell),
                    3);

                if (out_wfc_pw || out_wfc_r) // only for gamma_only now
                    this->set_pw_wfc(wfcpw, 0, ib, nspin, pes->charge->rho_save, pw_wfc_g);
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

    delete[] bands_picked;
    for (int is = 0; is < nspin; ++is)
    {
        for (int ib = 0; ib < nbands; ++ib)
            delete[] wfc_gamma_grid[is][ib];
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
                            int& out_wf,
                            int& out_wf_r,
                            const K_Vectors& kv,
                            const double nelec,
                            const int nbands_istate,
                            const int nbands,
                            const int nspin,
                            const int nlocal,
                            const std::string& global_out_dir)
{
    ModuleBase::TITLE("IState_Envelope", "begin");

    std::cout << " perform |psi(band, r)| for selected bands." << std::endl;

    // (1)
    // mohan update 2011-03-21
    // if ucell is odd, it's correct,
    // if ucell is even, it's also correct.
    // +1.0e-8 in case like (2.999999999+1)/2
    // if NSPIN=4, each band only one electron, fermi_band should be nelec
    int fermi_band = nspin < 4 ? static_cast<int>((nelec + 1) / 2 + 1.0e-8) : nelec;
    int bands_below = nbands_istate;
    int bands_above = nbands_istate;

    std::cout << " number of electrons = " << nelec << std::endl;
    std::cout << " number of occupied bands = " << fermi_band << std::endl;
    std::cout << " plot band decomposed charge density below fermi surface with " << bands_below << " bands."
              << std::endl;

    // (2) cicle:

    // (2.1) calculate the selected density matrix
    // from wave functions.

    // (2.2) carry out the grid integration to
    // get the charge density.

    // (2.3) output the charge density in .cub format.
    this->bands_picked = new bool[nbands];
    ModuleBase::GlobalFunc::ZEROS(bands_picked, nbands);
    for (int ib = 0; ib < nbands; ib++)
    {
        if (ib >= fermi_band - bands_below)
        {
            if (ib < fermi_band + bands_above)
            {
                bands_picked[ib] = true;
            }
        }
    }

    // allocate grid wavefunction for gamma_only
    std::vector<std::complex<double>**> wfc_k_grid(nspin);
    for (int ik = 0; ik < kv.get_nks(); ++ik)
    {
        wfc_k_grid[ik] = new std::complex<double>*[nbands];
        for (int ib = 0; ib < nbands; ++ib)
            wfc_k_grid[ik][ib] = new std::complex<double>[gk.gridt->lgd];
    }
    const double mem_size
        = sizeof(std::complex<double>) * double(gk.gridt->lgd) * double(nbands) * double(nspin) / 1024.0 / 1024.0;
    ModuleBase::Memory::record("IState_Envelope::begin::wfc_k_grid", mem_size);
    printf(" Estimated on-the-fly memory consuming by IState_Envelope::begin::wfc_k_grid: %f MB\n", mem_size);
    assert(mem_size > 0);
    // for pw-wfc in G space
    psi::Psi<std::complex<double>> pw_wfc_g(kv.ngk.data());

    if (out_wf || out_wf_r)
    {
        pw_wfc_g.resize(kv.get_nks(), nbands, wfcpw->npwk_max);
    }

    for (int ib = 0; ib < nbands; ib++)
    {
        if (bands_picked[ib])
        {
            const int nspin0 = (nspin == 2) ? 2 : 1;
            for (int ik = 0; ik < kv.get_nks(); ++ik) // the loop of nspin0 is included
            {
                const int ispin = kv.isk[ik];
                ModuleBase::GlobalFunc::ZEROS(pes->charge->rho[ispin],
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
                gk.cal_env_k(ik, wfc_k_grid[ik][ib], pes->charge->rho[ispin], kv.kvec_c, kv.kvec_d, GlobalC::ucell);

                std::stringstream ss;
                ss << global_out_dir << "BAND" << ib + 1 << "_k_" << ik / nspin0 + 1 << "_s_" << ispin + 1
                   << "_ENV.cube";
                const double ef_tmp = this->pes->eferm.get_efval(ispin);
                ModuleIO::write_rho(
#ifdef __MPI
                    bigpw->bz,
                    bigpw->nbz,
                    rhopw->nplane,
                    rhopw->startz_current,
#endif
                    pes->charge->rho[ispin],
                    ispin,
                    nspin,
                    0,
                    ss.str(),
                    rhopw->nx,
                    rhopw->ny,
                    rhopw->nz,
                    ef_tmp,
                    &(GlobalC::ucell),
                    3);

                if (out_wf || out_wf_r) // only for gamma_only now
                {
                    pw_wfc_g.fix_k(ik);
                    this->set_pw_wfc(wfcpw, ik, ib, nspin, pes->charge->rho, pw_wfc_g);
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

    delete[] bands_picked;
    for (int is = 0; is < nspin; ++is)
    {
        for (int ib = 0; ib < nbands; ++ib)
            delete[] wfc_k_grid[is][ib];
        delete[] wfc_k_grid[is];
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
    if (ib == 0) // once is enough
        ModuleBase::TITLE("IState_Envelope", "set_pw_wfc");

    std::vector<std::complex<double>> Porter(wfcpw->nrxx);
    // here I refer to v_hartree, but I don't know how to deal with NSPIN=4
    const int nspin0 = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin0; is++)
        for (int ir = 0; ir < wfcpw->nrxx; ir++)
            Porter[ir] += std::complex<double>(rho[is][ir], 0.0);

    // call FFT
    wfcpw->real2recip(Porter.data(), &wfc_g(ib, 0), ik);
}

#ifdef __MPI
template <typename T>
int IState_Envelope::set_wfc_grid(const Parallel_2D& p2d,
                                  const int nbands,
                                  const std::vector<int>& trace_lo,
                                  const T* in,
                                  T** out)
{
    ModuleBase::TITLE(" Local_Orbital_wfc", "set_wfc_grid");
    if (!out)
    {
        return 0;
    }
    const int naroc[2] = {p2d.nrow, p2d.ncol};
    for (int j = 0; j < naroc[1]; ++j)
    {
        const int igcol = p2d.local2global_col(j);
        if (igcol >= nbands)
        {
            continue;
        }
        for (int i = 0; i < naroc[0]; ++i)
        {
            const int igrow = p2d.local2global_row(i);
            const int mu_local = trace_lo[igrow];
            out[igcol][mu_local] = in[j * naroc[0] + i];
        }
    }
    return 0;
}

template int IState_Envelope::set_wfc_grid(const Parallel_2D& p2d,
                                           const int nbands,
                                           const std::vector<int>& trace_lo,
                                           const double* in,
                                           double** out);
template int IState_Envelope::set_wfc_grid(const Parallel_2D& p2d,
                                           const int nbands,
                                           const std::vector<int>& trace_lo,
                                           const std::complex<double>* in,
                                           std::complex<double>** out);

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
    auto mpi_dtype = std::is_same<T, double>::value ? MPI_DOUBLE : MPI_DOUBLE_COMPLEX;

    // get the rank of the current process
    int rank = 0;
    MPI_Comm_rank(pv.comm_2D, &rank);

    // calculate the maximum number of nlocal over all processes in pv.comm_2D range
    long buf_size;
    mpi_info = MPI_Reduce(&pv.nloc_wfc, &buf_size, 1, MPI_LONG, MPI_MAX, 0, pv.comm_2D);
    mpi_info = MPI_Bcast(&buf_size, 1, MPI_LONG, 0, pv.comm_2D); // get and then broadcast
    std::vector<T> lowf_block(buf_size);

    // this quantity seems to have the value returned by function numroc_ in ScaLAPACK?
    int naroc[2];

    // loop over all processors
    for (int iprow = 0; iprow < pv.dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < pv.dim1; ++ipcol)
        {
            // get the rank of the processor at the given coordinate
            int rank_at_coord;
            const int mpi_cart_coord[2] = {iprow, ipcol};
            mpi_info = MPI_Cart_rank(pv.comm_2D, mpi_cart_coord, &rank_at_coord); // get the MPI rank

            // keep in mind present function is concurrently called by all processors, thus
            // the following code block will only be executed once for each processor, which means
            // for each processor, get its MPI rank and MPI coord, then assign the naroc[0] and naroc[1]
            // with the value which should have been calculated automatically by ScaLAPACK function
            // numroc_.
            if (rank == rank_at_coord)
            {
                BlasConnector::copy(pv.nloc_wfc, lowf_2d, mem_stride, lowf_block.data(), mem_stride);
                naroc[0] = pv.nrow;
                naroc[1] = pv.ncol_bands;
            }

            // broadcast the number of row and column
            mpi_info = MPI_Bcast(naroc, 2, MPI_INT, rank_at_coord, pv.comm_2D);

            // broadcast the data, this means the data owned by one processor is broadcast
            // to all other processors in the communicator.
            mpi_info = MPI_Bcast(lowf_block.data(), buf_size, mpi_dtype, rank_at_coord, pv.comm_2D);

            // then use it to set the wfc_grid.
            Parallel_2D p2d;
            p2d.init(nlocal, nbands, pv.nb, pv.comm_2D);
            mpi_info = this->set_wfc_grid(p2d, nbands, trace_lo, lowf_block.data(), lowf_grid);
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