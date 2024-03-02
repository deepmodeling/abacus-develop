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
    build_nao(p_ucell_->ntype, p_ucell_->orbital_fn);
    // build another atomic orbital
    // PARALLELIZATION STRATEGY: only RANK-0 works
    #ifdef __MPI
    if(GlobalV::MY_RANK == 0)
    {
    #endif
    build_ao(ntype_, p_ucell_->pseudo_fn);
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

void toQO::build_nao(const int ntype, const std::string* const orbital_fn)
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
    if(GlobalV::MY_RANK == 0)
    {
        for(int it = 0; it < ntype_; it++)
        {
            orbital_fn_[it] = GlobalV::global_orbital_dir + orbital_fn[it];
        }
    }
#ifdef __MPI
    Parallel_Common::bcast_string(orbital_fn_, ntype_);
#endif

    nao_->build(ntype_, orbital_fn_, 'o');
    nao_->set_transformer(sbt);
    for(int it = 0; it < ntype_; it++)
    {
        int _nphi_it = 0;
        for(int l = 0; l <= nao_->lmax(it); l++)
        {
            for(int izeta = 0; izeta < nao_->nzeta(it, l); izeta++)
            {
                _nphi_it += 2*l + 1;
            }
        }
        nphi_ += _nphi_it*na_[it];
    }
    #ifdef __MPI
    if(GlobalV::MY_RANK == 0)
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

void toQO::build_hydrogen(const int ntype, const double* const charges, const int* const nmax)
{
    ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    ao_->build(ntype, charges, nmax, symbols_.data(), GlobalV::qo_thr, strategies_.data());
    ModuleBase::SphericalBesselTransformer sbt;
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

    #ifdef __MPI
    if(GlobalV::MY_RANK == 0)
    {
    #endif
    printf("Build arbitrary atomic orbital basis done.\n");
    #ifdef __MPI
    }
    #endif
}

void toQO::build_pswfc(const int ntype, const std::string* const pspot_fn, const double* const screening_coeffs)
{
    ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    std::string* pspot_fn_ = new std::string[ntype_];
    for(int it = 0; it < ntype; it++)
    {
        pspot_fn_[it] = GlobalV::global_pseudo_dir + pspot_fn[it];
    }
    ao_->build(ntype, pspot_fn_, screening_coeffs, GlobalV::qo_thr);
    ModuleBase::SphericalBesselTransformer sbt;
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

    #ifdef __MPI
    if(GlobalV::MY_RANK == 0)
    {
    #endif
    printf("Build arbitrary atomic orbital basis done.\n");
    #ifdef __MPI
    }
    #endif
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
void toQO::build_ao(const int ntype, const std::string* const pspot_fn)
{
    if(qo_basis_ == "hydrogen")
    {
        build_hydrogen(ntype_, charges_.data(), nmax_.data());
    }
    else if(qo_basis_ == "pswfc")
    {
        build_pswfc(ntype_, pspot_fn, GlobalV::qo_screening_coeff.data());
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
        if(GlobalV::MY_RANK == 0)
        {
        #endif
        // Not implemented error
        GlobalV::ofs_running << "Error: " << qo_basis_ << " is not implemented yet." << std::endl;
        ModuleBase::WARNING_QUIT("toQO::initialize", "Error: " + qo_basis_ + " is not implemented yet.");
        #ifdef __MPI
        }
        #endif
    } // radial functions generation completed
}

void toQO::calculate_ovlp_R(const int iR)
{
    // save memory mode: only write to ovlp_ao_nao_R_[0]
    int iR_save = save_mem_? 0 : iR;

    int irow = 0; // row and column index of ovlp_ao_nao_R_
    for(int it = 0; it < p_ucell_->ntype; it++)
    {
    // FOR EACH TYPE it, GET THE MAXIMUM l
        int lmaxi = atom_database_.principle_quantum_number[p_ucell_->atoms[it].ncpp.psd] - 1;
        for(int ia = 0; ia < p_ucell_->atoms[it].na; ia++)
        {
    // FOR EACH ATOM ia OF PRESENT TYPE it, SPECIFIES AN ATOM itia
    // BUT SPECIFYING AN ATOM HERE IS NOT NECESSARY, THE ONLY REASON IS THE ARRANGEMENT OF ovlp_ao_nao_R_
            for(int li = 0; li <= lmaxi; li++)
            {
                // orbitals arrange in the way stated: https://abacus.deepmodeling.com/en/latest/advanced/pp_orb.html#basis-set
                // generate the magnetic quantum number mi list
                std::vector<int> mis;
                for(int mi_abs = 0; mi_abs <= li; mi_abs++)
                {
                    mis.push_back(mi_abs);
                    if(mi_abs != 0) mis.push_back(-mi_abs);
                }
                if((!orbital_filter(li, strategies_[it]))&&(qo_basis_ == "pswfc")) continue;
    // RADIAL FUNCTIONS ARE ORGANIZED BY (l, zeta), SO FOR EACH l, GET THE MAXIMUM zeta
                int nzetai = ao_->nzeta(it, li);
    // FOR (l, zeta) OF ATOM itia, SPECIFY A RADIAL ATOMIC ORBITAL
                for(int izetai = 0; izetai < nzetai; izetai++)
                {
    // FOR EACH RADIAL ATOMIC ORBITAL, SPECIFY A SPHERICAL HARMONIC
                    //for(int mi = -li; mi <= li; mi++) // natural but it's not how ABACUS arrange orbitals
                    for(int mi : mis)
                    {
    // HERE WE GET flzeta(r)*Ylm(theta, phi),
    // THEN ANOTHER ORBITAL...(jt, ja, lj, izetaj, mj)
                        int icol = 0; 
                        for(int jt = 0; jt < p_ucell_->ntype; jt++)
                        {
                            for(int ja = 0; ja < p_ucell_->atoms[jt].na; ja++)
                            {
                                int lmaxj = p_ucell_->atoms[jt].nwl;
                                for(int lj = 0; lj <= lmaxj; lj++)
                                {
                                    // orbitals arrange in the way stated: https://abacus.deepmodeling.com/en/latest/advanced/pp_orb.html#basis-set
                                    // generate the magnetic quantum number mj list
                                    std::vector<int> mjs;
                                    for(int mj_abs = 0; mj_abs <= lj; mj_abs++)
                                    {
                                        mjs.push_back(mj_abs);
                                        if(mj_abs != 0) mjs.push_back(-mj_abs);
                                    }
                                    int nzetaj = nao_->nzeta(jt, lj);
                                    for(int izetaj = 0; izetaj < nzetaj; izetaj++)
                                    {
                                        //for(int mj = -lj; mj <= lj; mj++) // natural but it's not how ABACUS arrange orbitals
                                        for(int mj : mjs)
                                        {
    // TWO ATOMIC ORBITALS ARE SPECIFIED, THEN WE NEED TO CALCULATE THE OVERLAP IN SUPERCELL
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
                                                it, li, izetai, mi,
                                                jt, lj, izetaj, mj,
                                                Rij, &ovlp_R_[iR_save][irow][icol]
                                            );
                                            icol++; // CHARNGE ORBITAL2: (jt, ja, lj, izetaj, mj)
                                        }
                                    }
                                }
                            }
                        }
                        irow++; // CHARNGE ORBITAL1: (it, ia, li, izetai, mi)
                    }
                }
            }
        }
    }
}

void toQO::calculate_ovlp_k(int ik)
{
    for(int iR = 0; iR < nR_; iR++)
    {
        calculate_ovlp_R(iR); // calculate S(R) for each R, save to ovlp_ao_nao_R_
        if(save_mem_) append_ovlp_R_eiRk(ik, iR);
    }
    if(!save_mem_) fold_ovlp_R(ik);
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
        calculate_ovlp_k(ik);
        write_ovlp<std::complex<double>>(ovlp_k_, ik);
    }
    printf("Calculating overlap integrals for kpoints done.\n\n");
    #ifdef __MPI
    }
    #endif
}
