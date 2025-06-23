#ifndef STO_TOOL_H
#define STO_TOOL_H
#include "source_cell/klist.h"
#include "module_hamilt_pw/hamilt_stodft/hamilt_sdft_pw.h"
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"
#include "source_base/module_device/memory_op.h"
#include "module_psi/psi.h"

/**
 * @brief Check if Emin and Emax are converged
 *
 * @param nche_in N order of Chebyshev expansion
 * @param try_emin trial Emin
 * @param try_emax trial Emax
 * @param nbands_sto number of stochastic bands
 */
template <typename FPTYPE, typename Device>
struct check_che_op
{
    using syncmem_complex_h2d_op = base_device::memory::synchronize_memory_op<std::complex<FPTYPE>, Device, base_device::DEVICE_CPU>;
    void operator()(const int& nche_in,
                    const double& try_emin,
                    const double& try_emax,
                    const int& nbands_sto,
                    K_Vectors* p_kv,
                    Stochastic_WF<std::complex<FPTYPE>, Device>* p_stowf,
                    hamilt::HamiltSdftPW<std::complex<FPTYPE>, Device>* p_hamilt_sto);
};

/**
 * @brief structure to distribute calculation among processors
 *
 * @param start start index of this processor
 * @param num_per number of elements for this processor
 *
 */
struct parallel_distribution
{
    parallel_distribution(const int& num_all, const int& np, const int myrank)
    {
        int num_per = num_all / np;
        int st_per = num_per * myrank;
        int re = num_all % np;
        if (myrank < re)
        {
            ++num_per;
            st_per += myrank;
        }
        else
        {
            st_per += re;
        }
        this->start = st_per;
        this->num_per = num_per;
    }
    int start;
    int num_per;
};

#ifdef __MPI
/**
 * @brief gather information from all processors
 *
 */
struct info_gatherv
{
    info_gatherv(const int& ngroup_per, const int& np, const int& num_in_group, MPI_Comm comm_world)
    {
        nrecv = new int[np];
        displs = new int[np];
        MPI_Allgather(&ngroup_per, 1, MPI_INT, nrecv, 1, MPI_INT, comm_world);
        displs[0] = 0;
        for (int i = 1; i < np; ++i)
        {
            displs[i] = displs[i - 1] + nrecv[i - 1];
        }
        for (int i = 0; i < np; ++i)
        {
            nrecv[i] *= num_in_group;
            displs[i] *= num_in_group;
        }
    }
    ~info_gatherv()
    {
        delete[] nrecv;
        delete[] displs;
    }
    int* nrecv = nullptr;
    int* displs = nullptr;
};
#endif

/**
 * @brief convert psi from double to float
 *
 * @param psi_in input psi of double
 * @param psi_out output psi of float
 */
template <typename FPTYPE_IN, typename FPTYPE_OUT, typename Device>
struct convert_psi_op
{
    using castmem_complex_op
        = base_device::memory::cast_memory_op<std::complex<FPTYPE_OUT>, std::complex<FPTYPE_IN>, Device, Device>;
    void operator()(const psi::Psi<std::complex<FPTYPE_IN>, Device>& psi_in,
                    psi::Psi<std::complex<FPTYPE_OUT>, Device>& psi_out)
    {
        psi_in.fix_k(0);
        psi_out.fix_k(0);
        castmem_complex_op()(psi_out.get_pointer(), psi_in.get_pointer(), psi_in.size());
    }
};

/**
 * @brief gather chi from all processors
 *
 * @param chi stochasitc wave function of this processor
 * @param chi_all gathered stochastic wave function
 * @param npwx maximum number of plane waves on all processors
 * @param nrecv_sto number of stochastic orbitals on each processor
 * @param displs_sto displacement of stochastic orbitals on each processor
 * @param perbands_sto number of stochastic bands of this processor
 * @return psi::Psi<std::complex<float>> pointer to gathered stochastic wave function
 *
 */
template <typename FPTYPE, typename Device>
struct gatherchi_op
{
    psi::Psi<std::complex<FPTYPE>, Device>* operator()(psi::Psi<std::complex<FPTYPE>, Device>& chi,
                                                       psi::Psi<std::complex<FPTYPE>, Device>& chi_all,
                                                       const int& npwx,
                                                       int* nrecv_sto,
                                                       int* displs_sto,
                                                       const int perbands_sto);
};

#endif