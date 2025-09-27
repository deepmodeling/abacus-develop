#include "gint_common.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_io/module_parameter/parameter.h"

#ifdef __MPI
#include "source_base/module_external/blacs_connector.h"
#include <mpi.h>
#endif

namespace ModuleGint
{

void compose_hr_gint(HContainer<double>& hr_gint)
{
    ModuleBase::TITLE("Gint", "compose_hr_gint");
    ModuleBase::timer::tick("Gint", "compose_hr_gint");
    for (int iap = 0; iap < hr_gint.size_atom_pairs(); iap++)
    {
        auto& ap = hr_gint.get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        if (iat1 > iat2)
        {
            // fill lower triangle matrix with upper triangle matrix
            // the upper <IJR> is <iat2, iat1>
            const hamilt::AtomPair<double>* upper_ap = hr_gint.find_pair(iat2, iat1);
            const hamilt::AtomPair<double>* lower_ap = hr_gint.find_pair(iat1, iat2);
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
    ModuleBase::timer::tick("Gint", "compose_hr_gint");
}

template <typename T>
void transfer_hr_gint_to_hR(const HContainer<T>& hr_gint, HContainer<T>& hR)
{
    ModuleBase::TITLE("Gint", "transfer_hr_gint_to_hR");
    ModuleBase::timer::tick("Gint", "transfer_hr_gint_to_hR");
#ifdef __MPI
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1)
    {
        hR.add(hr_gint);
    }
    else
    {
        hamilt::transferSerials2Parallels(hr_gint, &hR);
    }
#else
    hR.add(hr_gint);
#endif
    ModuleBase::timer::tick("Gint", "transfer_hr_gint_to_hR");
}


void merge_hr_part_to_hR(const std::vector<hamilt::HContainer<double>>& hr_gint_tmp ,
                         hamilt::HContainer<std::complex<double>>* hR,
                         const GintInfo& gint_info){
    ModuleBase::TITLE("Gint_k", "transfer_pvpR");
    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");

    const UnitCell* ucell_in = gint_info.get_ucell();
    int mg = hR->get_paraV()->get_global_row_size()/2;
    int ng = hR->get_paraV()->get_global_col_size()/2;
    int nb = hR->get_paraV()->get_block_size()/2;
    hamilt::HContainer<std::complex<double>>* hR_tmp;


#ifdef __MPI
    int blacs_ctxt = hR->get_paraV()->blacs_ctxt;
    std::vector<int> iat2iwt(ucell_in->nat);
    for (int iat = 0; iat < ucell_in->nat; iat++) {
        iat2iwt[iat] = ucell_in->get_iat2iwt()[iat]/2;
    }
    Parallel_Orbitals *pv = new Parallel_Orbitals();
    pv->set(mg, ng, nb, blacs_ctxt);
    pv->set_atomic_trace(iat2iwt.data(), ucell_in->nat, mg);
    auto ijr_info = hR->get_ijr_info();
    hR_tmp = new hamilt::HContainer<std::complex<double>>(pv, nullptr, &ijr_info);
#endif

    //select hr_gint_tmp 
    std::vector<int> first = {0, 1, 1, 0};
    std::vector<int> second= {3, 2, 2, 3};
    //select position in the big matrix
    std::vector<int> row_set = {0, 0, 1, 1};
    std::vector<int> col_set = {0, 1, 0, 1};
    //construct complex matrix
    std::vector<int> clx_i = {1, 0, 0, -1};
    std::vector<int> clx_j = {0, 1, -1, 0};
    for (int is = 0; is < 4; is++){
        if(!PARAM.globalv.domag && (is==1 || is==2)) continue;
        hR_tmp->set_zero();
        hamilt::HContainer<std::complex<double>>* hRGint_tmpCd = new hamilt::HContainer<std::complex<double>>(ucell_in->nat);
        hRGint_tmpCd->insert_ijrs( &(gint_info.get_ijr_info()), *(ucell_in));
        hRGint_tmpCd->allocate(nullptr, true);
        hRGint_tmpCd->set_zero();
        for (int iap = 0; iap < hRGint_tmpCd->size_atom_pairs(); iap++)
        {
            auto* ap = &hRGint_tmpCd->get_atom_pair(iap);
            const int iat1 = ap->get_atom_i();
            const int iat2 = ap->get_atom_j();
            if (iat1 <= iat2)
            {
                hamilt::AtomPair<std::complex<double>>* upper_ap = ap;
                hamilt::AtomPair<std::complex<double>>* lower_ap = hRGint_tmpCd->find_pair(iat2, iat1);
                const hamilt::AtomPair<double>* ap_nspin1 = hr_gint_tmp [first[is]].find_pair(iat1, iat2);
                const hamilt::AtomPair<double>* ap_nspin2 = hr_gint_tmp [second[is]].find_pair(iat1, iat2);
                for (int ir = 0; ir < upper_ap->get_R_size(); ir++)
                {   
                    const auto R_index = upper_ap->get_R_index(ir);
                    auto upper_mat = upper_ap->find_matrix(R_index);
                    auto mat_nspin1 = ap_nspin1->find_matrix(R_index);
                    auto mat_nspin2 = ap_nspin2->find_matrix(R_index);
                    // The row size and the col size of upper_matrix is double that of matrix_nspin_0
                    for (int irow = 0; irow < mat_nspin1->get_row_size(); ++irow)
                    {
                        for (int icol = 0; icol < mat_nspin1->get_col_size(); ++icol)
                        {
                            upper_mat->get_value(irow, icol) = mat_nspin1->get_value(irow, icol) 
                            + std::complex<double>(clx_i[is], clx_j[is]) * mat_nspin2->get_value(irow, icol);
                        }
                    }
                    //fill the lower triangle matrix
                    //When is=0 or 3, the real part does not need conjugation; 
                    //when is=1 or 2, the small matrix is not Hermitian, so conjugation is not needed
                    if (iat1 < iat2)
                    {
                        auto lower_mat = lower_ap->find_matrix(-R_index);
                        for (int irow = 0; irow < upper_mat->get_row_size(); ++irow)
                        {
                            for (int icol = 0; icol < upper_mat->get_col_size(); ++icol)
                            {
                                lower_mat->get_value(icol, irow) = upper_mat->get_value(irow, icol);
                            }
                        }
                    }

                }
            }
        }
        // transfer hRGint_tmpCd to parallel hR_tmp
#ifdef __MPI
        hamilt::transferSerials2Parallels( *hRGint_tmpCd, hR_tmp);
#else
        hR_tmp = hRGint_tmpCd;
#endif
        // merge hR_tmp to hR
        for (int iap = 0; iap < hR->size_atom_pairs(); iap++)
        {
            auto* ap = &hR->get_atom_pair(iap);
            const int iat1 = ap->get_atom_i();
            const int iat2 = ap->get_atom_j();
            auto* ap_nspin = hR_tmp ->find_pair(iat1, iat2);
            for (int ir = 0; ir < ap->get_R_size(); ir++)
            {   
                const auto R_index = ap->get_R_index(ir);
                auto upper_mat = ap->find_matrix(R_index);
                auto mat_nspin = ap_nspin->find_matrix(R_index);
                // The row size and the col size of upper_matrix is double that of matrix_nspin_0
                for (int irow = 0; irow < mat_nspin->get_row_size(); ++irow)
                {
                    for (int icol = 0; icol < mat_nspin->get_col_size(); ++icol)
                    {
                        upper_mat->get_value(2*irow+row_set[is], 2*icol+col_set[is]) = 
                        mat_nspin->get_value(irow, icol);
                    }
                }
            }
        }
        delete hRGint_tmpCd;
    }
    ModuleBase::timer::tick("Gint_k", "transfer_pvpR");
    return;
}



// gint_info should not have been a parameter, but it was added to initialize dm_gint_full
// In the future, we might try to remove the gint_info parameter
template<typename T>
void transfer_dm_2d_to_gint(
    const GintInfo& gint_info,
    std::vector<HContainer<T>*> dm,
    std::vector<HContainer<T>>& dm_gint)
{
    ModuleBase::TITLE("Gint", "transfer_dm_2d_to_gint");
    ModuleBase::timer::tick("Gint", "transfer_dm_2d_to_gint");

    if (PARAM.inp.nspin != 4)
    {
        // dm_gint.size() usually equals to PARAM.inp.nspin,
        // but there is exception within source_lcao/module_lr
        for (int is = 0; is < dm_gint.size(); is++)
        {
#ifdef __MPI
            hamilt::transferParallels2Serials(*dm[is], &dm_gint[is]);
#else
            dm_gint[is].set_zero();
            dm_gint[is].add(*dm[is]);
#endif
        }
    } else  // NSPIN=4 case
    {

        // is=0:↑↑, 1:↑↓, 2:↓↑, 3:↓↓
        const int row_set[4] = {0, 0, 1, 1};
        const int col_set[4] = {0, 1, 0, 1};
        int mg = dm[0]->get_paraV()->get_global_row_size()/2;
        int ng = dm[0]->get_paraV()->get_global_col_size()/2;
        int nb = dm[0]->get_paraV()->get_block_size()/2;
        const UnitCell* ucell = gint_info.get_ucell();
        auto ijr_info = dm[0]->get_ijr_info();
#ifdef __MPI
        int blacs_ctxt = dm[0]->get_paraV()->blacs_ctxt;
        std::vector<int> iat2iwt(ucell->nat);
        for (int iat = 0; iat < ucell->nat; iat++) {
            iat2iwt[iat] = ucell->get_iat2iwt()[iat]/2;
        }
        Parallel_Orbitals pv{};
        pv.set(mg, ng, nb, blacs_ctxt);
        pv.set_atomic_trace(iat2iwt.data(), ucell->nat, mg);
        HContainer<T> dm2d_tmp(&pv, nullptr, &ijr_info);
#else
        auto* dm2d_tmp = new hamilt::HContainer<T>(ucell->nat);
        dm2d_tmp -> insert_ijrs(&ijr_info, *ucell);
        dm2d_tmp -> allocate(nullptr, true);
#endif
         for (int is = 0; is < 4; is++){
            for (int iap = 0; iap < dm[0]->size_atom_pairs(); ++iap) {
                auto& ap = dm[0]->get_atom_pair(iap);
                int iat1 = ap.get_atom_i();
                int iat2 = ap.get_atom_j();
                for (int ir = 0; ir < ap.get_R_size(); ++ir) {
                    const ModuleBase::Vector3<int> r_index = ap.get_R_index(ir);
#ifdef __MPI
                    T* matrix_out = dm2d_tmp.find_matrix(iat1, iat2, r_index)->get_pointer();
#else
                    T* matrix_out = dm2d_tmp->find_matrix(iat1, iat2, r_index)->get_pointer();
#endif
                    T* matrix_in = ap.get_pointer(ir);
                    for (int irow = 0; irow < ap.get_row_size()/2; irow ++) {
                        for (int icol = 0; icol < ap.get_col_size()/2; icol ++) {
                            int index_i = irow* ap.get_col_size()/2 + icol;
                            int index_j = (irow*2+row_set[is]) * ap.get_col_size() + icol*2+col_set[is];
                            matrix_out[index_i] = matrix_in[index_j];
                        }
                    }
                }
            }
#ifdef __MPI         
            hamilt::transferParallels2Serials(dm2d_tmp, &dm_gint[is]);
#else
            dm_gint[is].set_zero();
            dm_gint[is].add(*dm2d_tmp);
#endif
        }//is=4
    }
    ModuleBase::timer::tick("Gint", "transfer_dm_2d_to_gint");
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
                    int nbands,  // needed if MPI is disabled
                    int nlocal,  // needed if MPI is disabled
                    const Parallel_Orbitals& pv,
                    T* wfc_gint,
                    const GintInfo& gint_info)
{
    ModuleBase::TITLE("Gint", "wfc_2d_to_gint");
    ModuleBase::timer::tick("Gint", "wfc_2d_to_gint");

#ifdef __MPI
    // dimension related
    nlocal = pv.desc_wfc[2];
    nbands = pv.desc_wfc[3];

    const std::vector<int>& trace_lo = gint_info.get_trace_lo();

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
            wfc_gint[i * nlocal + j] = wfc_2d[i * nlocal + j];
        }
    }
#endif
    ModuleBase::timer::tick("Gint", "wfc_2d_to_gint");
}

template void transfer_hr_gint_to_hR(
    const HContainer<double>& hr_gint,
    HContainer<double>& hR);
template void transfer_hr_gint_to_hR(
    const HContainer<std::complex<double>>& hr_gint,
    HContainer<std::complex<double>>& hR);
template void transfer_dm_2d_to_gint(
    const GintInfo& gint_info,
    std::vector<HContainer<double>*> dm,
    std::vector<HContainer<double>>& dm_gint);
template void transfer_dm_2d_to_gint(
    const GintInfo& gint_info,
    std::vector<HContainer<std::complex<double>>*> dm,
    std::vector<HContainer<std::complex<double>>>& dm_gint);
template void wfc_2d_to_gint(
    const double* wfc_2d,
    int nbands,
    int nlocal,
    const Parallel_Orbitals& pv,
    double* wfc_grid,
    const GintInfo& gint_info);
template void wfc_2d_to_gint(
    const std::complex<double>* wfc_2d,
    int nbands,
    int nlocal,
    const Parallel_Orbitals& pv,
    std::complex<double>* wfc_grid,
    const GintInfo& gint_info);
}