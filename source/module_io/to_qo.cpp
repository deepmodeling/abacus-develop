#include "module_io/to_qo.h"
#include "module_basis/module_nao/two_center_integrator.h"


toQO::toQO(std::string qo_basis, std::string strategy)
{
    qo_basis_ = qo_basis;
    strategy_ = strategy;
}

toQO::~toQO()
{
}

void toQO::initialize(UnitCell *p_ucell,
                      const std::vector<ModuleBase::Vector3<double>> kvecs_c)
{
    kvecs_c_ = kvecs_c;
    nkpts_ = kvecs_c.size();
    unwrap_unitcell(p_ucell);
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
    build_nao(p_ucell_->ntype, p_ucell_->orbital_fn);

    // build another atomic orbital
    std::vector<int> nmax = std::vector<int>(p_ucell_->ntype);
    for(int itype = 0; itype < ntype_; itype++)
    {
        nmax[itype] = atom_database_.principle_quantum_number[symbols_[itype]];
    }
    build_ao(p_ucell_->ntype, charges_.data(), nmax.data());

    // neighbor list search
    scan_supercell();
    
    // allocate memory for ovlp_ao_nao_R_ and ovlp_ao_nao_k_
    clean_up();
}

void toQO::build_nao(const int ntype, const std::string* const orbital_fn)
{
    // build the numerical atomic orbital basis
    ModuleBase::SphericalBesselTransformer sbt;
    nao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    nao_->build(ntype, orbital_fn, 'o');
    nao_->set_transformer(sbt);
    nphi_ = nao_->nchi();
}

void toQO::build_ao(const int ntype, const double* const charges, const int* const nmax)
{
    // build another atomic orbital
    if(qo_basis_ == "hydrogen")
    {
        ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
        ao_->build(ntype, charges, nmax, strategy_);
        ModuleBase::SphericalBesselTransformer sbt;
        ao_->set_transformer(sbt);
        nchi_ = ao_->nchi();
    }
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
    double rcut_max = std::max(nao_->rcut_max(), ao_->rcut_max());
    int ngrid = int(rcut_max / 0.01) + 1;
    double cutoff = 2*rcut_max;
    nao_->set_uniform_grid(true, ngrid, cutoff, 'i', true);
    ao_->set_uniform_grid(true, ngrid, cutoff, 'i', true);

    overlap_calculator_->tabulate(*ao_, *nao_, 'S', ngrid, rcut_max);

    // save memory mode: only write to ovlp_ao_nao_R_[0]
    int iR_save = save_mem_? 0 : iR;

    int irow = 0; int icol = 0; // row and column index of ovlp_ao_nao_R_
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
    // RADIAL FUNCTIONS ARE ORGANIZED BY (l, zeta), SO FOR EACH l, GET THE MAXIMUM zeta
                int nzetai = ao_->nzeta(it, li);
    // FOR (l, zeta) OF ATOM itia, SPECIFY A RADIAL ATOMIC ORBITAL
                for(int izetai = 0; izetai < nzetai; izetai++)
                {
    // FOR EACH RADIAL ATOMIC ORBITAL, SPECIFY A SPHERICAL HARMONIC
                    for(int mi = -li; mi <= li; mi++)
                    {
    // HERE WE GET flzeta(r)*Ylm(theta, phi),
    // THEN ANOTHER ORBITAL...(jt, ja, lj, izetaj, mj)
                        for(int jt = 0; jt < p_ucell_->ntype; jt++)
                        {
                            for(int ja = 0; ja < p_ucell_->atoms[jt].na; ja++)
                            {
                                int lmaxj = p_ucell_->atoms[jt].nwl;
                                for(int lj = 0; lj <= lmaxj; lj++)
                                {        
                                    int nzetaj = nao_->nzeta(jt, lj);
                                    for(int izetaj = 0; izetaj < nzetaj; izetaj++)
                                    {
                                        for(int mj = -lj; mj <= lj; mj++)
                                        {
    // TWO ATOMIC ORBITALS ARE SPECIFIED, THEN WE NEED TO CALCULATE THE OVERLAP IN SUPERCELL
                                            ModuleBase::Vector3<double> rij = p_ucell_->atoms[it].tau[ia] - p_ucell_->atoms[jt].tau[ja];
                                            // there is waste here, but for easy to understand, I don't optimize it.
                                            ModuleBase::Vector3<int> R = supercells_[iR];
                                            ModuleBase::Vector3<double> Rij;
                                            Rij.x = rij.x + R.x * p_ucell_->a1.x 
                                                          + R.y * p_ucell_->a2.x 
                                                          + R.z * p_ucell_->a3.x;
                                            Rij.y = rij.y + R.x * p_ucell_->a1.y 
                                                          + R.y * p_ucell_->a2.y 
                                                          + R.z * p_ucell_->a3.y;
                                            Rij.z = rij.z + R.x * p_ucell_->a1.z 
                                                          + R.y * p_ucell_->a2.z 
                                                          + R.z * p_ucell_->a3.z;
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

void toQO::calculate_ovlp_k(ModuleBase::Vector3<double> kvec_c)
{
    for(int iR = 0; iR < nR_; iR++)
    {
        calculate_ovlp_R(iR); // calculate S(R) for each R, save to ovlp_ao_nao_R_
        if(save_mem_) append_ovlp_R_eiRk(kvec_c, iR);
    }
    if(!save_mem_) fold_ovlp_R(kvec_c);
}

void toQO::calculate(std::vector<ModuleBase::Vector3<double>> kvecs_c)
{
    for(int ik = 0; ik < nkpts_; ik++)
    {
        calculate_ovlp_k(kvecs_c[ik]);
        write_ovlp(ovlp_k_, "QO_OVLP_K_" + std::to_string(ik) + ".dat"  );
    }
}

using CplxMatrix = std::vector<std::vector<std::complex<double>>>;
using RealMatrix = std::vector<std::vector<double>>;