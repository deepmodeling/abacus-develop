#include "overlap.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/module_rt/td_folding.h"
#include "source_lcao/module_rt/td_info.h"

#include <functional>
#include <vector>

namespace
{

/// @brief Populate an HContainer with atom pairs based on neighbor search
/// @details Searches for adjacent atoms within orbital cutoff distance and inserts
/// valid atom pairs into the container. This is the common logic shared between
/// SR (overlap matrix) and SR_async (asynchronous overlap for Hefei-NAMD).
/// @tparam TR Data type for real-space matrix elements (double or complex<double>)
/// @param container Target HContainer to populate with atom pairs
/// @param ucell Unit cell containing atomic structure
/// @param grid_driver Grid driver for neighbor searching
/// @param orb_cutoff Orbital cutoff radii for each atom type
/// @param dtau_modifier Optional function to modify dtau before cutoff check (for velocity shifts in NAMD)
///                      Signature: void(int iat1, int I1, int T1, ModuleBase::Vector3<double>& dtau)
template <typename TR>
void populate_atom_pairs(hamilt::HContainer<TR>* container,
                         const UnitCell* ucell,
                         const Grid_Driver* grid_driver,
                         const std::vector<double>& orb_cutoff,
                         std::function<void(int, int, int, ModuleBase::Vector3<double>&)> dtau_modifier = nullptr)
{
    const Parallel_Orbitals* paraV = container->get_paraV();

    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        const auto tau1 = ucell->get_tau(iat1);
        int T1 = 0;
        int I1 = 0;
        ucell->iat2iait(iat1, &I1, &T1);

        AdjacentAtomInfo adjs;
        grid_driver->Find_atom(*ucell, tau1, T1, I1, &adjs);

        // Loop over adjacent atoms (including self: ad=0 corresponds to iat1 itself)
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = ucell->itia2iat(T2, I2);

            // Skip if atom pair has no local orbitals in parallel distribution
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }

            const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            ModuleBase::Vector3<double> dtau = ucell->cal_dtau(iat1, iat2, R_index);

            // Apply optional dtau modification (e.g., velocity shift for NAMD)
            if (dtau_modifier)
            {
                dtau_modifier(iat1, I1, T1, dtau);
            }

            // Skip atoms beyond orbital cutoff distance
            // Note: Use strict inequality (<) because at exactly the cutoff distance,
            // the theoretical matrix element is zero, but numerical errors can produce
            // non-zero values that affect reproducibility.
            const double distance = dtau.norm() * ucell->lat0;
            const double cutoff_sum = orb_cutoff[T1] + orb_cutoff[T2];
            if (distance >= cutoff_sum)
            {
                continue;
            }

            hamilt::AtomPair<TR> atom_pair(iat1, iat2, R_index, paraV);
            container->insert_pair(atom_pair);
        }
    }

    // Allocate memory for matrix elements and initialize to zero
    container->allocate(nullptr, true);
}

} // anonymous namespace

template <typename TK, typename TR>
hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>::Overlap(HS_Matrix_K<TK>* hsk_in,
                                                             const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                                             hamilt::HContainer<TR>* hR_in,
                                                             hamilt::HContainer<TR>* SR_in,
                                                             const UnitCell* ucell_in,
                                                             const std::vector<double>& orb_cutoff,
                                                             const Grid_Driver* GridD_in,
                                                             const TwoCenterIntegrator* intor)
    : hamilt::OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in), orb_cutoff_(orb_cutoff), intor_(intor), gridD(GridD_in)
{
    this->cal_type = calculation_type::lcao_overlap;
    this->ucell = ucell_in;
    this->SR = SR_in;
#ifdef __DEBUG
    assert(this->ucell != nullptr);
    assert(this->SR != nullptr);
#endif
    // Initialize SR to allocate sparse overlap matrix memory.
    // Only initialize if SR_in is not nullptr (for force calculation, SR_in can be nullptr).
    if (SR_in != nullptr)
    {
        this->initialize_SR(GridD_in);
    }
}

template <typename TK, typename TR>
hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>::~Overlap()
{
}

template <typename TK, typename TR>
void hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>::initialize_SR(const Grid_Driver* GridD)
{
    ModuleBase::TITLE("OverlapNew", "initialize_SR");
    ModuleBase::timer::start("OverlapNew", "initialize_SR");

    populate_atom_pairs(this->SR, this->ucell, GridD, this->orb_cutoff_);

    ModuleBase::timer::end("OverlapNew", "initialize_SR");
}

template <typename TK, typename TR>
void hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>::calculate_SR()
{
    ModuleBase::TITLE("Overlap", "calculate_SR");
    ModuleBase::timer::start("Overlap", "calculate_SR");
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iap = 0; iap < this->SR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<TR>& tmp = this->SR->get_atom_pair(iap);
        const int iat1 = tmp.get_atom_i();
        const int iat2 = tmp.get_atom_j();
        const Parallel_Orbitals* paraV = tmp.get_paraV();

        for (int iR = 0; iR < tmp.get_R_size(); ++iR)
        {
            const ModuleBase::Vector3<int> R_index = tmp.get_R_index(iR);
            auto dtau = ucell->cal_dtau(iat1, iat2, R_index);
            TR* data_pointer = tmp.get_pointer(iR);
            this->cal_SR_IJR(iat1, iat2, paraV, dtau, data_pointer);
        }
    }
    // if TK == double, then SR should be fixed to gamma case
    // judge type of TK equal to double
    if (std::is_same<TK, double>::value)
    {
        this->SR->fix_gamma();
    }
    ModuleBase::timer::end("Overlap", "calculate_SR");
}

// cal_SR_IJR()
template <typename TK, typename TR>
void hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>::cal_SR_IJR(const int& iat1,
                                                                  const int& iat2,
                                                                  const Parallel_Orbitals* paraV,
                                                                  const ModuleBase::Vector3<double>& dtau,
                                                                  TR* data_pointer)
{
    // ---------------------------------------------
    // get info of orbitals of atom1 and atom2 from ucell
    // ---------------------------------------------
    int T1=0;
    int I1=0;
    this->ucell->iat2iait(iat1, &I1, &T1);
    int T2=0;
    int I2=0;
    this->ucell->iat2iait(iat2, &I2, &T2);
    Atom& atom1 = this->ucell->atoms[T1];
    Atom& atom2 = this->ucell->atoms[T2];

    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();

    const int* iw2l1 = atom1.iw2l.data();
    const int* iw2n1 = atom1.iw2n.data();
    const int* iw2m1 = atom1.iw2m.data();
    const int* iw2l2 = atom2.iw2l.data();
    const int* iw2n2 = atom2.iw2n.data();
    const int* iw2m2 = atom2.iw2m.data();

    // ---------------------------------------------
    // calculate the overlap matrix for each pair of orbitals
    // ---------------------------------------------
    double olm[3] = {0, 0, 0};
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    const int step_trace = col_indexes.size() + 1;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const int iw1 = row_indexes[iw1l] / npol;
        const int L1 = iw2l1[iw1];
        const int N1 = iw2n1[iw1];
        const int m1 = iw2m1[iw1];

        // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
        int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const int iw2 = col_indexes[iw2l] / npol;
            const int L2 = iw2l2[iw2];
            const int N2 = iw2n2[iw2];
            const int m2 = iw2m2[iw2];

            // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
            int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;
            intor_->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * this->ucell->lat0, olm);
            for (int ipol = 0; ipol < npol; ipol++)
            {
                data_pointer[ipol * step_trace] += olm[0];
            }
            data_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
    }
}

// contributeHR()
template <typename TK, typename TR>
void hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
    if (this->SR_fixed_done)
    {
        return;
    }
    this->calculate_SR();
    this->SR_fixed_done = true;
}

// contributeHk()
template <>
void hamilt::Overlap<hamilt::OperatorLCAO<double, double>>::contributeHk(int ik)
{
    //! if k vector is not changed, then do nothing and return, only for gamma_only case
    if (this->kvec_d[ik] == this->kvec_d_old)
    {
        return;
    }
    ModuleBase::TITLE("Overlap", "contributeHk");
    ModuleBase::timer::start("Overlap", "contributeHk");
    
    //! set SK to zero and then calculate SK for each k vector
    this->hsk->set_zero_sk();
    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
    {
        const int nrow = this->SR->get_atom_pair(0).get_paraV()->get_row_size();
        hamilt::folding_HR(*this->SR, this->hsk->get_sk(), this->kvec_d[ik], nrow, 1);
    }
    else
    {
        const int ncol = this->SR->get_atom_pair(0).get_paraV()->get_col_size();
        hamilt::folding_HR(*this->SR, this->hsk->get_sk(), this->kvec_d[ik], ncol, 0);
    }
    
    // update kvec_d_old
    this->kvec_d_old = this->kvec_d[ik];

    ModuleBase::timer::end("Overlap", "contributeHk");
}
template <typename TK, typename TR>
void hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    ModuleBase::TITLE("Overlap", "contributeHk");
    ModuleBase::timer::start("Overlap", "contributeHk");
    
    //! set SK to zero and then calculate SK for each k vector
    this->hsk->set_zero_sk();
    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
    {
        const int nrow = this->SR->get_atom_pair(0).get_paraV()->get_row_size();
        if(PARAM.inp.td_stype == 2)
        {
            module_rt::folding_HR_td(*ucell, *this->SR, this->hsk->get_sk(), this->kvec_d[ik], TD_info::cart_At, nrow, 1);
        }
        else
        {
            hamilt::folding_HR(*this->SR, this->hsk->get_sk(), this->kvec_d[ik], nrow, 1);
        }
    }
    else
    {
        const int ncol = this->SR->get_atom_pair(0).get_paraV()->get_col_size();
        if(PARAM.inp.td_stype == 2)
        {
            module_rt::folding_HR_td(*ucell, *this->SR, this->hsk->get_sk(), this->kvec_d[ik], TD_info::cart_At, ncol, 0);
        }
        else
        {
            hamilt::folding_HR(*this->SR, this->hsk->get_sk(), this->kvec_d[ik], ncol, 0);
        }
    }
    
    // update kvec_d_old
    this->kvec_d_old = this->kvec_d[ik];

    ModuleBase::timer::end("Overlap", "contributeHk");
}
template <typename TK, typename TR>
TK* hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>::getSk()
{
    if (this->hsk != nullptr)
    {
        return this->hsk->get_sk();
    }
    return nullptr;
}

//==============================================================================
// Asynchronous overlap matrix methods for Hefei-NAMD interface
// These methods calculate <phi(t-dt)|phi(t)> by shifting atom positions backward
// using atomic velocities, enabling non-adiabatic molecular dynamics calculations.
//==============================================================================

template <typename TK, typename TR>
hamilt::HContainer<TR>* hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>::calculate_SR_async(const UnitCell& ucell_in,
                                                                                            const double md_dt,
                                                                                            const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("OverlapNew", "calculate_SR_async");
    ModuleBase::timer::start("OverlapNew", "calculate_SR_async");

    // Initialize SR_async for Hefei-NAMD asynchronous overlap calculation
    // This is done here to use the exact dtau with velocity shifts
    hamilt::HContainer<TR>* SR_async = new hamilt::HContainer<TR>(paraV);

    // Define velocity shift modifier for dtau
    // This shifts atom1 backward to its position at (t - dt),
    // giving the overlap <phi(t-dt)|phi(t)> needed for NAMD.
    const double vel_to_dtau = md_dt / ModuleBase::AU_to_FS / ucell_in.lat0;
    auto velocity_shift = [&ucell_in, vel_to_dtau](int iat1, int I1, int T1, ModuleBase::Vector3<double>& dtau) {
        const Atom* atom1 = &ucell_in.atoms[T1];
        for (int dim = 0; dim < 3; dim++)
        {
            dtau[dim] += atom1->vel[I1][dim] * vel_to_dtau;
        }
    };

    // Populate atom pairs with velocity-shifted dtau for cutoff checking
    populate_atom_pairs(SR_async, &ucell_in, this->gridD, this->orb_cutoff_, velocity_shift);

    // Calculate overlap matrix elements with velocity-shifted dtau
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iap = 0; iap < SR_async->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<TR>& atom_pair = SR_async->get_atom_pair(iap);
        const int iat1 = atom_pair.get_atom_i();
        const int iat2 = atom_pair.get_atom_j();
        const Parallel_Orbitals* paraV_local = atom_pair.get_paraV();

        // Get atom type and index for velocity lookup
        int T1 = 0;
        int I1 = 0;
        ucell_in.iat2iait(iat1, &I1, &T1);
        const Atom* atom1 = &ucell_in.atoms[T1];

        for (int iR = 0; iR < atom_pair.get_R_size(); ++iR)
        {
            const ModuleBase::Vector3<int> R_index = atom_pair.get_R_index(iR);
            ModuleBase::Vector3<double> dtau = ucell_in.cal_dtau(iat1, iat2, R_index);

            // Apply velocity shift to dtau
            for (int dim = 0; dim < 3; dim++)
            {
                dtau[dim] += atom1->vel[I1][dim] * vel_to_dtau;
            }

            TR* data_pointer = atom_pair.get_pointer(iR);
            this->cal_SR_IJR(iat1, iat2, paraV_local, dtau, data_pointer);
        }
    }

    // For gamma-only calculations, apply symmetry constraint
    if (std::is_same<TK, double>::value)
    {
        SR_async->fix_gamma();
    }

    ModuleBase::timer::end("OverlapNew", "calculate_SR_async");
    return SR_async;
}

template <typename TK, typename TR>
void hamilt::Overlap<hamilt::OperatorLCAO<TK, TR>>::output_SR_async_csr(const int istep,
                                                                         hamilt::HContainer<TR>* SR_async,
                                                                         const int precision)
{
    if (SR_async == nullptr)
    {
        return;
    }

    ModuleBase::TITLE("OverlapNew", "output_SR_async_csr");
    ModuleBase::timer::start("OverlapNew", "output_SR_async_csr");

#ifdef __MPI
    // Gather distributed SR_async to rank 0 for serial output
    const Parallel_Orbitals* paraV = SR_async->get_paraV();
    const int nbasis = SR_async->get_nbasis();

    Parallel_Orbitals serial_paraV;
    serial_paraV.init(nbasis, nbasis, nbasis, paraV->comm());
    serial_paraV.set_serial(nbasis, nbasis);
    serial_paraV.set_atomic_trace(this->ucell->get_iat2iwt(), this->ucell->nat, nbasis);

    hamilt::HContainer<TR> SR_async_serial(&serial_paraV);
    hamilt::gatherParallels(*SR_async, &SR_async_serial, 0);
#else
    hamilt::HContainer<TR>& SR_async_serial = *SR_async;
#endif

    // Only rank 0 writes the output file
    if (GlobalV::MY_RANK == 0)
    {
        const std::string filename = PARAM.globalv.global_out_dir + "syns_nao.csr";
        std::ofstream ofs;

        // First step creates new file, subsequent steps append
        if (istep <= 0)
        {
            ofs.open(filename);
        }
        else
        {
            ofs.open(filename, std::ios::app);
        }

        // Write header information
        ofs << "IONIC_STEP: " << istep + 1 << std::endl;
        ofs << "Matrix Dimension of S_async(R): " << SR_async_serial.get_nbasis() << std::endl;
        ofs << "Matrix number of S_async(R): " << SR_async_serial.size_R_loop() << std::endl;

        // Write matrix data in CSR format
        const double sparse_threshold = 1e-10;
        hamilt::Output_HContainer<TR> output_handler(&SR_async_serial, ofs, sparse_threshold, precision);
        output_handler.write();

        ofs.close();
    }

    ModuleBase::timer::end("OverlapNew", "output_SR_async_csr");
}

// Include force/stress implementation
#include "overlap_force_stress.hpp"

template class hamilt::Overlap<hamilt::OperatorLCAO<double, double>>;
template class hamilt::Overlap<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::Overlap<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;
