#include "gint_common.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_parameter/parameter.h"

#ifdef __MPI
#include "module_base/blacs_connector.h"
#include <mpi.h>
#endif

namespace ModuleGint
{

void compose_hr_gint(std::shared_ptr<HContainer<double>> hr_gint)
{
    for (int iap = 0; iap < hr_gint->size_atom_pairs(); iap++)
    {
        auto& ap = hr_gint->get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        if (iat1 > iat2)
        {
            // fill lower triangle matrix with upper triangle matrix
            // the upper <IJR> is <iat2, iat1>
            const hamilt::AtomPair<double>* upper_ap = hr_gint->find_pair(iat2, iat1);
            const hamilt::AtomPair<double>* lower_ap = hr_gint->find_pair(iat1, iat2);
#ifdef __DEBUG
            assert(upper_ap != nullptr);
#endif
            for (int ir = 0; ir < ap.get_R_size(); ir++)
            {
                auto R_index = ap.get_R_index(ir);
                auto upper_mat = upper_ap->find_matrix(-R_index);
                auto lower_mat = lower_ap->find_matrix(R_index);
                for (int irow = 0; irow < upper_mat->get_row_size(); ++irow)
                {
                    for (int icol = 0; icol < upper_mat->get_col_size(); ++icol)
                    {
                        lower_mat->get_value(icol, irow) = upper_ap->get_value(irow, icol);
                    }
                }
            }
        }
    }
}

void compose_hr_gint(std::vector<std::shared_ptr<HContainer<double>>> hr_gint_part,
        std::shared_ptr<HContainer<std::complex<double>>> hr_gint_full)
{
    for (int iap = 0; iap < hr_gint_full->size_atom_pairs(); iap++)
    {
        auto* ap = &hr_gint_full->get_atom_pair(iap);
        const int iat1 = ap->get_atom_i();
        const int iat2 = ap->get_atom_j();
        if (iat1 <= iat2)
        {
            hamilt::AtomPair<std::complex<double>>* upper_ap = ap;
            hamilt::AtomPair<std::complex<double>>* lower_ap = hr_gint_full->find_pair(iat2, iat1);
            const hamilt::AtomPair<double>* ap_nspin_0 = hr_gint_part[0]->find_pair(iat1, iat2);
            const hamilt::AtomPair<double>* ap_nspin_3 = hr_gint_part[3]->find_pair(iat1, iat2);
            for (int ir = 0; ir < upper_ap->get_R_size(); ir++)
            {
                const auto R_index = upper_ap->get_R_index(ir);
                auto upper_mat = upper_ap->find_matrix(R_index);
                auto mat_nspin_0 = ap_nspin_0->find_matrix(R_index);
                auto mat_nspin_3 = ap_nspin_3->find_matrix(R_index);

                // The row size and the col size of upper_matrix is double that of matrix_nspin_0
                for (int irow = 0; irow < mat_nspin_0->get_row_size(); ++irow)
                {
                    for (int icol = 0; icol < mat_nspin_0->get_col_size(); ++icol)
                    {
                        upper_mat->get_value(2*irow, 2*icol) = mat_nspin_0->get_value(irow, icol) + mat_nspin_3->get_value(irow, icol);
                        upper_mat->get_value(2*irow+1, 2*icol+1) = mat_nspin_0->get_value(irow, icol) - mat_nspin_3->get_value(irow, icol);
                    }
                }

                if (PARAM.globalv.domag)
                {
                    const hamilt::AtomPair<double>* ap_nspin_1 = hr_gint_part[1]->find_pair(iat1, iat2);
                    const hamilt::AtomPair<double>* ap_nspin_2 = hr_gint_part[2]->find_pair(iat1, iat2);
                    const auto mat_nspin_1 = ap_nspin_1->find_matrix(R_index);
                    const auto mat_nspin_2 = ap_nspin_2->find_matrix(R_index);
                    for (int irow = 0; irow < mat_nspin_1->get_row_size(); ++irow)
                    {
                        for (int icol = 0; icol < mat_nspin_1->get_col_size(); ++icol)
                        {
                            upper_mat->get_value(2*irow, 2*icol+1) = mat_nspin_1->get_value(irow, icol) +  std::complex<double>(0.0, 1.0) * mat_nspin_2->get_value(irow, icol);
                            upper_mat->get_value(2*irow+1, 2*icol) = mat_nspin_1->get_value(irow, icol) -  std::complex<double>(0.0, 1.0) * mat_nspin_2->get_value(irow, icol);
                        }
                    }
                }

                // fill the lower triangle matrix
                if (iat1 < iat2)
                {
                    auto lower_mat = lower_ap->find_matrix(-R_index);
                    for (int irow = 0; irow < upper_mat->get_row_size(); ++irow)
                    {
                        for (int icol = 0; icol < upper_mat->get_col_size(); ++icol)
                        {
                            lower_mat->get_value(icol, irow) = conj(upper_mat->get_value(irow, icol));
                        }
                    }
                }
            }
        }
    }
}

template <typename T>
void transfer_hr_gint_to_hR(std::shared_ptr<const HContainer<T>> hr_gint, HContainer<T>* hR)
{
#ifdef __MPI
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1)
    {
        hR->add(*hr_gint);
    }
    else
    {
        hamilt::transferSerials2Parallels(*hr_gint, hR);
    }
#else
    hR->add(*hr_gint);
#endif
}

// gint_info should not have been a parameter, but it was added to initialize dm_gint_full
// In the future, we might try to remove the gint_info parameter
template<typename T>
void transfer_dm_2d_to_gint(
    std::shared_ptr<const GintInfo> gint_info,
    std::vector<HContainer<T>*> dm,
    std::vector<std::shared_ptr<HContainer<T>>> dm_gint)
{
    // To check whether input parameter dm_2d has been initialized
    assert(PARAM.inp.nspin == dm_gint.size()
           && "The size of dm should be equal to the number of spins!");
    if(PARAM.inp.nspin != 4)
    {
        assert(dm.size() == PARAM.inp.nspin);
    } else
    {
        assert(dm.size() == 1);
    }

    if (PARAM.inp.nspin != 4)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
#ifdef __MPI
            hamilt::transferParallels2Serials(*dm[is], dm_gint[is].get());
#else
            dm_gint[is]->set_zero();
            dm_gint[is]->add(*dm[is]);
#endif
        }
    } else  // NSPIN=4 case
    {
#ifdef __MPI
        const int npol = 2;
        std::shared_ptr<HContainer<T>> dm_full = gint_info->get_hr<T>(npol);
        hamilt::transferParallels2Serials(*dm[0], dm_full.get());
#else
        HContainer<T>* dm_full = dm[0];
#endif
        std::vector<T*> tmp_pointer(4, nullptr);
        for (int iap = 0; iap < dm_full->size_atom_pairs(); iap++)
        {
            auto& ap = dm_full->get_atom_pair(iap);
            const int iat1 = ap.get_atom_i();
            const int iat2 = ap.get_atom_j();
            for (int ir = 0; ir < ap.get_R_size(); ir++)
            {
                const ModuleBase::Vector3<int> r_index = ap.get_R_index(ir);
                for (int is = 0; is < 4; is++)
                {
                    tmp_pointer[is] =
                        dm_gint[is]->find_matrix(iat1, iat2, r_index)->get_pointer();
                }
                T* data_full = ap.get_pointer(ir);
                for (int irow = 0; irow < ap.get_row_size(); irow += 2)
                {
                    for (int icol = 0; icol < ap.get_col_size(); icol += 2)
                    {
                        *(tmp_pointer[0])++ = data_full[icol];
                        *(tmp_pointer[1])++ = data_full[icol + 1];
                    }
                    data_full += ap.get_col_size();
                    for (int icol = 0; icol < ap.get_col_size(); icol += 2)
                    {
                        *(tmp_pointer[2])++ = data_full[icol];
                        *(tmp_pointer[3])++ = data_full[icol + 1];
                    }
                    data_full += ap.get_col_size();
                }
            }
        }
    }
}

int globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    const int iblock = localindex / nblk;
    const int gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

int localIndex(int globalindex, int nblk, int nprocs, int& myproc)
{
    myproc = int((globalindex % (nblk * nprocs)) / nblk);
    return int(globalindex / (nblk * nprocs)) * nblk + globalindex % nblk;
}

template <typename T>
void wfc_2d_to_gint(const T* wfc_2d,
                    const Parallel_Orbitals& pv,
                    T* wfc_gint,
                    std::shared_ptr<const GintInfo> gint_info)
{
    ModuleBase::TITLE("Module_gint", "wfc_2d_to_gint");
    ModuleBase::timer::tick("Module_gint", "wfc_2d_to_gint");

    // dimension related
    const int nlocal = pv.desc_wfc[2];
    const int nbands = pv.desc_wfc[3];

#ifdef __MPI
    const std::vector<int>& trace_lo = gint_info->get_trace_lo();

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
    std::vector<T> wfc_block(buf_size);

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
                BlasConnector::copy(pv.nloc_wfc, wfc_2d, mem_stride, wfc_block.data(), mem_stride);
                naroc[0] = pv.nrow;
                naroc[1] = pv.ncol_bands;
                Cxgebs2d(pv.blacs_ctxt, &scope, &top, 2, 1, naroc, 2);
                Cxgebs2d(pv.blacs_ctxt, &scope, &top, buf_size, 1, wfc_block.data(), buf_size);
            }
            else
            {
                Cxgebr2d(pv.blacs_ctxt, &scope, &top, 2, 1, naroc, 2, iprow, ipcol);
                Cxgebr2d(pv.blacs_ctxt, &scope, &top, buf_size, 1, wfc_block.data(), buf_size, iprow, ipcol);
            }

            // then use it to set the wfc_grid.
            const int nb = pv.nb;
            const int dim0 = pv.dim0;
            const int dim1 = pv.dim1;
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
                    if (wfc_gint && mu_local >= 0)
                    {
                        wfc_gint[igcol * nlocal + mu_local] = wfc_block[j * naroc[0] + i];
                    }
                }
            }
            // this operation will let all processors have the same wfc_grid
        }
    }
#else
    for (int i = 0; i < nbands; ++i)
    {
        for (int j = 0; j < nlocal; ++j)
        {
            wfc_k_grid[i * nlocal + j] = psi[0](i, j);
        }
    }
#endif
    ModuleBase::timer::tick("Module_gint", "wfc_2d_to_gint");
}

template void transfer_hr_gint_to_hR(
    std::shared_ptr<const HContainer<double>> hr_gint,
    HContainer<double>* hR);
template void transfer_hr_gint_to_hR(
    std::shared_ptr<const HContainer<std::complex<double>>> hr_gint,
    HContainer<std::complex<double>>* hR);
template void transfer_dm_2d_to_gint(
    std::shared_ptr<const GintInfo> gint_info,
    std::vector<HContainer<double>*> dm,
    std::vector<std::shared_ptr<HContainer<double>>> dm_gint);
template void transfer_dm_2d_to_gint(
    std::shared_ptr<const GintInfo> gint_info,
    std::vector<HContainer<std::complex<double>>*> dm,
    std::vector<std::shared_ptr<HContainer<std::complex<double>>>> dm_gint);
template void wfc_2d_to_gint(
    const double* wfc_2d,
    const Parallel_Orbitals& pv,
    double* wfc_grid,
    std::shared_ptr<const GintInfo> gint_info);
template void wfc_2d_to_gint(
    const std::complex<double>* wfc_2d,
    const Parallel_Orbitals& pv,
    std::complex<double>* wfc_grid,
    std::shared_ptr<const GintInfo> gint_info);
}