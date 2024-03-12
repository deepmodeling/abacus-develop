#include "module_io/to_qo.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_base/ylm.h"
#include "module_base/parallel_common.h"

toQO::toQO(std::string qo_basis, std::vector<std::string> strategies)
{
    qo_basis_ = qo_basis;
    strategies_ = strategies;
}

toQO::~toQO()
{
}

void toQO::initialize(UnitCell* p_ucell,
                      const std::vector<ModuleBase::Vector3<double>>& kvecs_d)
{
    // BEGIN: "Two-center bundle build"
    read_abacus_variables(p_ucell, kvecs_d, GlobalV::MY_RANK, GlobalV::NPROC);
    // build two-center overlap calculator
    overlap_calculator_ = std::unique_ptr<TwoCenterIntegrator>(new TwoCenterIntegrator);
    // build orbitals
    /* 
       orbitals building operation is to (for example) read all orbital files and save to
       radial collection data structure. These kinds of information is irrelevant to the
       structure simulated. Once two atoms' one certain overlap integral is needed, get
       orbitals from radial collection by type-l-izeta, so that get the radial function.
       Then specify the angular part Ylm l and m, also l' and m', and the correct distance.
       Then the overlap integral is calculated.
     */
    // build the numerical atomic orbital basis
    // PARALLELIZATION STRATEGY: use RANK-0 to read in the files, then broadcast
    build_nao(ntype_, 
              GlobalV::global_orbital_dir,
              p_ucell_->orbital_fn,
              GlobalV::MY_RANK);
    // build another atomic orbital
    build_ao(ntype_, 
             GlobalV::global_pseudo_dir,
             p_ucell_->pseudo_fn, 
             GlobalV::qo_screening_coeff, 
             GlobalV::qo_thr,
             GlobalV::ofs_running,
             GlobalV::MY_RANK);
    // neighbor list search
    scan_supercell(GlobalV::MY_RANK, GlobalV::NPROC);
    // build grids, for all processes
    double rcut_max = std::max(nao_->rcut_max(), ao_->rcut_max());
    int ngrid = int(rcut_max / 0.01) + 1;
    double cutoff = 2.0*rcut_max;
    nao_->to_file("qo-mpi-test-nao");
    ao_->to_file("qo-mpi-test-ao");
    nao_->set_uniform_grid(true, ngrid, cutoff, 'i', true);
    ao_->set_uniform_grid(true, ngrid, cutoff, 'i', true);
    overlap_calculator_->tabulate(*ao_, *nao_, 'S', ngrid, cutoff);
    // prepare for Ylm
    ModuleBase::Ylm::set_coefficients();
    // END: "Two-center bundle build"
    // allocate memory for ovlp_ao_nao_R_ and ovlp_ao_nao_k_
    allocate_ovlp(true); allocate_ovlp(false);
}

void toQO::build_nao(const int ntype, 
                     const std::string orbital_dir,
                     const std::string* const orbital_fn,
                     const int rank)
{
    // build the numerical atomic orbital basis
    ModuleBase::SphericalBesselTransformer sbt;
    nao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    // add GlobalV::global_orbital_dir ahead of orbital_fn
    int ntype_ = ntype;
#ifdef __MPI
    Parallel_Common::bcast_int(ntype_);
#endif
    std::string* orbital_fn_ = new std::string[ntype_];
    if(rank == 0)
    {
        for(int it = 0; it < ntype_; it++)
        {
            orbital_fn_[it] = orbital_dir + orbital_fn[it];
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_string(orbital_fn_, ntype_);
#endif

    nao_->build(ntype_, orbital_fn_, 'o');
    nao_->set_transformer(sbt);
    // indexing
    radialcollection_indexing(*nao_, na_, index_nao_, rindex_nao_);
    nphi_ = index_nao_.size();

    delete[] orbital_fn_;
}

bool toQO::orbital_filter(const int l, const std::string spec)
{
    std::vector<std::string> l2symbol = {"s", "p", "d", "f", "g"}; // seems enough
    if(spec == "all") return true;
    else if(spec.find_first_of(l2symbol[l]) != std::string::npos) return true;
    else return false;
}

void toQO::build_hydrogen(const int ntype, 
                          const double* const charges, 
                          const bool slater_screening,
                          const int* const nmax,
                          const double qo_thr,
                          const int rank)
{
    ModuleBase::SphericalBesselTransformer sbt;
    ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    // for this method, all processes CAN do together
    ao_->build(ntype, 
               charges, 
               slater_screening, 
               nmax, 
               symbols_.data(), 
               qo_thr, 
               strategies_.data(),
               rank);
    ao_->set_transformer(sbt);
    // indexing
    radialcollection_indexing(*ao_, na_, index_ao_, rindex_ao_);
    nchi_ = index_ao_.size();
}

void toQO::build_pswfc(const int ntype, 
                       const std::string pseudo_dir,
                       const std::string* const pspot_fn, 
                       const double* const screening_coeffs,
                       const double qo_thr,
                       const int rank)
{
    ModuleBase::SphericalBesselTransformer sbt;
    ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    int ntype_ = ntype;
#ifdef __MPI
    Parallel_Common::bcast_int(ntype_);
#endif
    std::string* pspot_fn_ = new std::string[ntype_];
    for(int it = 0; it < ntype; it++)
    {
        pspot_fn_[it] = pseudo_dir + pspot_fn[it];
    }
#ifdef __MPI
    Parallel_Common::bcast_string(pspot_fn_, ntype_);
#endif
    // for this method, all processes MIGHT NOT do together, because of possible conflict of reading files
    // in the following build function, the file reading is done by RANK-0, then broadcast to other processes
    ao_->build(ntype, pspot_fn_, screening_coeffs, qo_thr, rank);
    ao_->set_transformer(sbt);
    // indexing
    radialcollection_indexing(*ao_, na_, index_ao_, rindex_ao_);
    nchi_ = index_ao_.size();
    delete[] pspot_fn_;
}
/*
void toQO::build_szv(const int ntype)
{
    // build the numerical atomic orbital basis
    ModuleBase::SphericalBesselTransformer sbt;
    ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    // add GlobalV::global_orbital_dir ahead of orbital_fn
    int ntype_ = ntype;
    std::vector<std::string> orbital_fn_(ntype_, "szv.orb");
    ao_->build(ntype_, orbital_fn_.data(), 'o');
    ao_->set_transformer(sbt);
    for(int itype = 0; itype < ntype; itype++)
    {
        int _nchi_it = 0;
        for(int l = 0; l <= ao_->lmax(itype); l++)
        {
            _nchi_it += (2*l+1)*ao_->nzeta(itype, l);
        }
        nchi_ += _nchi_it * na_[itype];
    }
}
*/
void toQO::build_ao(const int ntype, 
                    const std::string pseudo_dir,
                    const std::string* const pspot_fn,
                    const std::vector<double> screening_coeffs,
                    const double qo_thr,
                    const std::ofstream& ofs_running,
                    const int rank)
{
    if(qo_basis_ == "hydrogen")
    {
        bool with_slater_screening = std::find_if(screening_coeffs.begin(), screening_coeffs.end(), 
            [](double sc) { return sc > 1e-10; }) != screening_coeffs.end();
        build_hydrogen(ntype_,                  /// ntype
                       charges_.data(),         /// charges
                       with_slater_screening,   /// slater_screening
                       nmax_.data(),            /// nmax
                       qo_thr,                  /// qo_thr
                       rank);                   /// rank
    }
    else if(qo_basis_ == "pswfc")
    {
        build_pswfc(ntype_,                     /// ntype
                    pseudo_dir,                 /// pseudo_dir
                    pspot_fn,                   /// pspot_fn
                    screening_coeffs.data(),    /// screening_coeffs
                    qo_thr,                     /// qo_thr
                    rank);                      /// rank
    }
}

void toQO::calculate_ovlpR(const int iR)
{
    assert (rindex_ao_.size() == nchi_);
    assert (rindex_nao_.size() == nphi_);
    for(int irow = 0; irow < nchi_; irow++)
    {
        //         it,  ia,  li,  izeta, mi
        std::tuple<int, int, int, int, int> orb1 = rindex_ao_[irow];
        int it = std::get<0>(orb1);
        int ia = std::get<1>(orb1);
        int li = std::get<2>(orb1);
        int izeta = std::get<3>(orb1);
        int mi = std::get<4>(orb1);
        for(int icol = 0; icol < nphi_; icol++)
        {
            //         jt,  ja,  lj,  jzeta, mj
            std::tuple<int, int, int, int, int> orb2 = rindex_nao_[icol];
            int jt = std::get<0>(orb2);
            int ja = std::get<1>(orb2);
            int lj = std::get<2>(orb2);
            int jzeta = std::get<3>(orb2);
            int mj = std::get<4>(orb2);
            ModuleBase::Vector3<double> rij = p_ucell_->atoms[jt].tau[ja] - p_ucell_->atoms[it].tau[ia];
            // there is waste here, but for easy to understand, I don't optimize it.
            ModuleBase::Vector3<int> R = supercells_[iR];
            ModuleBase::Vector3<double> Rij;
            Rij.x = rij.x + double(R.x) * p_ucell_->a1.x 
                          + double(R.y) * p_ucell_->a2.x 
                          + double(R.z) * p_ucell_->a3.x;
            Rij.y = rij.y + double(R.x) * p_ucell_->a1.y 
                          + double(R.y) * p_ucell_->a2.y 
                          + double(R.z) * p_ucell_->a3.y;
            Rij.z = rij.z + double(R.x) * p_ucell_->a1.z 
                          + double(R.y) * p_ucell_->a2.z 
                          + double(R.z) * p_ucell_->a3.z;
            Rij *= p_ucell_->lat0; // convert to Bohr
            // print two-center integrals related variables in one line
            printf("it: %d, ia: %d, li: %d, iz: %d, mi: %d, jt: %d, ja: %d, lj: %d, jz: %d, mj: %d, iR: %d, Rij: %f %f %f\n", 
                it, ia, li, izeta, mi, jt, ja, lj, jzeta, mj, iR, Rij.x, Rij.y, Rij.z);
            overlap_calculator_->calculate(
                it, li, izeta, mi,
                jt, lj, jzeta, mj,
                Rij, &ovlpR_[irow*nphi_+icol]
            );
        }
    }
}

void toQO::calculate_ovlpk(int ik)
{
    // first calculate all Rs corresponding two-center integrals on own process
    if(ik == iks_[0])
    {
        // std::string info = "Calculate two-center integrals in realspace for given pair of orbitals\n";
        // info += "Current process: " + std::to_string(GlobalV::MY_RANK) + "\n";
        // info += "Total number of Rs: " + std::to_string(nR_tot_) + "\n";
        // info += "Local number of Rs: " + std::to_string(iRs_.size()) + "\n";
        // info += "Rs indices: ";
        // for(auto iR: iRs_) info += std::to_string(iR) + " ";
        // info += "\n";
        // printf("%s", info.c_str());
        for(auto iR: iRs_)
        {
            calculate_ovlpR(iR);
            write_ovlp<double>(GlobalV::global_out_dir, ovlpR_, nchi_, nphi_, true, iR);
        }
#ifdef __MPI
        // wait for all processes to finish the calculation
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    // then read and calculate all ks corresponding two-center integrals on own process
    for(int iR = 0; iR < nR_tot_; iR++)
    {
        int barrier_iR = (iR + GlobalV::MY_RANK) % GlobalV::NPROC;
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        read_ovlp(GlobalV::global_out_dir, nchi_, nphi_, true, barrier_iR);
#ifdef __MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        append_ovlpR_eiRk(ik, iR);
    }
}
void toQO::calculate()
{
    // std::string info = "Calculate two-center integrals in kspace for given pair of orbitals\n";
    // info += "Current process: " + std::to_string(GlobalV::MY_RANK) + "\n";
    // info += "Total number of kpoints: " + std::to_string(nks_tot_) + "\n";
    // info += "Local number of kpoints: " + std::to_string(iks_.size()) + "\n";
    // info += "kpoints indices: ";
    // for(auto ik: iks_) info += std::to_string(ik) + " ";
    // info += "\n";
    // printf("%s", info.c_str());
    for(auto ik: iks_)
    {
        zero_out_ovlps(false);
        calculate_ovlpk(ik);
        write_ovlp<std::complex<double>>(GlobalV::global_out_dir, ovlpk_, nchi_, nphi_, false, ik);
    }
#ifdef __MPI
    // wait for all processes to finish the calculation
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    // delete all QO_ovlpR_* files
    // if(GlobalV::MY_RANK == 0)
    // {
    //     for(int iR = 0; iR < nR_tot_; iR++)
    //     {
    //         std::string filename = GlobalV::global_out_dir + "/QO_ovlpR_" + std::to_string(iR) + ".dat";
    //         std::remove(filename.c_str());
    //     }
    // }
}
