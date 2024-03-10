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
    #ifdef __MPI
    if(GlobalV::MY_RANK == 0)
    {
    #endif
    printf("\n---- Quasiatomic Orbital (QO) Analysis Initialization ----\n");
    #ifdef __MPI
    }
    #endif
    kvecs_d_ = kvecs_d;
    nkpts_ = kvecs_d.size();

    // BEGIN: "Two-center bundle build"
    unwrap_unitcell(p_ucell);
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
    build_nao(p_ucell_->ntype, 
              GlobalV::global_orbital_dir,
              p_ucell_->orbital_fn,
              GlobalV::MY_RANK);
    // build another atomic orbital
    // PARALLELIZATION STRATEGY: only RANK-0 works
    #ifdef __MPI
    if(GlobalV::MY_RANK == 0)
    {
    #endif
    build_ao(ntype_, 
             GlobalV::global_pseudo_dir,
             p_ucell_->pseudo_fn, 
             GlobalV::qo_screening_coeff, 
             GlobalV::qo_thr,
             GlobalV::ofs_running,
             GlobalV::MY_RANK);
    // neighbor list search
    scan_supercell();
    // build grids
    double rcut_max = std::max(nao_->rcut_max(), ao_->rcut_max());
    int ngrid = int(rcut_max / 0.01) + 1;
    double cutoff = 2.0*rcut_max;
    nao_->set_uniform_grid(true, ngrid, cutoff, 'i', true);
    ao_->set_uniform_grid(true, ngrid, cutoff, 'i', true);
    overlap_calculator_->tabulate(*ao_, *nao_, 'S', ngrid, cutoff);
    // prepare for Ylm
    ModuleBase::Ylm::set_coefficients();
    // END: "Two-center bundle build"
    // allocate memory for ovlp_ao_nao_R_ and ovlp_ao_nao_k_
    allocate_ovlp(true); allocate_ovlp(false);
    printf("---- Quasiatomic Orbital (QO) Analysis Initialization Done ----\n");
    #ifdef __MPI
    }
    #endif
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
    #ifdef __MPI
    if(rank == 0)
    {
    #endif
    printf("Build numerical atomic orbital basis done.\n");
    #ifdef __MPI
    }
    #endif
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
    ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    ao_->build(ntype, 
               charges, 
               slater_screening, 
               nmax, 
               symbols_.data(), 
               qo_thr, 
               strategies_.data());
    ModuleBase::SphericalBesselTransformer sbt;
    ao_->set_transformer(sbt);
}

void toQO::build_pswfc(const int ntype, 
                       const std::string pseudo_dir,
                       const std::string* const pspot_fn, 
                       const double* const screening_coeffs,
                       const double qo_thr,
                       const int rank)
{
    ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    std::string* pspot_fn_ = new std::string[ntype_];
    for(int it = 0; it < ntype; it++)
    {
        pspot_fn_[it] = pseudo_dir + pspot_fn[it];
    }
    ao_->build(ntype, pspot_fn_, screening_coeffs, qo_thr);
    ModuleBase::SphericalBesselTransformer sbt;
    ao_->set_transformer(sbt);
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
        build_hydrogen(ntype_, 
                       charges_.data(),
                       with_slater_screening, 
                       nmax_.data(),
                       qo_thr,
                       rank);
    }
    else if(qo_basis_ == "pswfc")
    {
        build_pswfc(ntype_, 
                    pseudo_dir,
                    pspot_fn, 
                    screening_coeffs.data(),
                    qo_thr,
                    rank);
    }
    /*
    else if(qo_basis_ == "szv")
    {
        build_szv(ntype_);
    }
    */
    else
    {
        #ifdef __MPI
        if(rank == 0)
        {
        #endif
        // Not implemented error
        ModuleBase::WARNING_QUIT("toQO::initialize", "Error: " + qo_basis_ + " is not implemented yet.");
        #ifdef __MPI
        }
        #endif
    }
    // indexing
    radialcollection_indexing(*ao_, na_, index_ao_, rindex_ao_);
    nchi_ = index_ao_.size();
    // radial functions generation completed
    #ifdef __MPI
    if(rank == 0)
    {
    #endif
    printf("Build atom-centered orbital basis done.\n");
    #ifdef __MPI
    }
    #endif
}

void toQO::calculate_ovlpR(const int iR)
{
    for(int irow = 0; irow < nchi_; irow++)
    {
        //         it,  ia,  li,  izeta, mi
        std::tuple<int, int, int, int, int> orb1 = rindex_ao_[irow];
        for(int icol = 0; icol < nphi_; icol++)
        {
            //         jt,  ja,  lj,  jzeta, mj
            std::tuple<int, int, int, int, int> orb2 = rindex_nao_[icol];
            int it = std::get<0>(orb1);
            int ia = std::get<1>(orb1);
            int jt = std::get<0>(orb2);
            int ja = std::get<1>(orb2);
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
            overlap_calculator_->calculate(
                it, std::get<2>(orb1), std::get<3>(orb1), std::get<4>(orb1),
                jt, std::get<2>(orb2), std::get<3>(orb2), std::get<4>(orb2),
                Rij, &ovlpR_[irow*nphi_+icol]
            );
        }
    }
}

void toQO::calculate_ovlpk(int ik)
{
    for(int iR = 0; iR < nR_; iR++)
    {
        calculate_ovlpR(iR); // calculate S(R) for each R, save to ovlp_ao_nao_R_
        append_ovlpR_eiRk(ik, iR);
    }
}
void toQO::calculate()
{
    #ifdef __MPI
    if(GlobalV::MY_RANK == 0)
    {
    #endif
    printf("Calculating overlap integrals for kpoints.\n");
    if(nkpts_ < nR_)
    {
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl
                  << "! Warning: number of kpoints is less than number of supercells, " << std::endl
                  << "! this will cause information loss when transform matrix R -> k. " << std::endl
                  << "! The further conversion k -> R cannot recover full information." << std::endl
                  << "! Number of kpoints: " << nkpts_ << std::endl
                  << "! Number of supercells: " << nR_ << std::endl
                  << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    }
    write_supercells();
    for(int ik = 0; ik < nkpts_; ik++)
    {
        zero_out_ovlps(false);
        calculate_ovlpk(ik);
        write_ovlp<std::complex<double>>(GlobalV::global_out_dir, /// dir
                                         ovlpk_,                  /// ovlp
                                         nchi_,                   /// nrows  
                                         nphi_,                   /// ncols
                                         false,                   /// is_R
                                         ik);                     /// ik or iR
    }
    printf("Calculating overlap integrals for kpoints done.\n\n");
    #ifdef __MPI
    }
    #endif
}
