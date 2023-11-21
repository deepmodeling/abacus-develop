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
                      int nkpts)
{
    p_ucell_ = p_ucell;
    nkpts_ = nkpts;
    // build orbitals
    ModuleBase::SphericalBesselTransformer sbt;

    // build the numerical atomic orbital basis
    nao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
    nao_->build(p_ucell_->ntype, p_ucell_->orbital_fn, 'o');
    nao_->set_transformer(sbt);
    
    // build another atomic orbital
    if(qo_basis_ == "hydrogen")
    {
        int ntype_ = p_ucell_->ntype;
        double* charges = new double[ntype_];
        std::string* symbols = new std::string[ntype_];
        int* nmax = new int[ntype_];
        for(int itype = 0; itype < ntype_; itype++)
        {
            charges[itype] = p_ucell_->atoms[itype].ncpp.zv;
            symbols[itype] = p_ucell_->atoms[itype].ncpp.psd;
            nmax[itype] = atom_database_.principle_quantum_number[symbols[itype]];
        }
        ao_ = std::unique_ptr<RadialCollection>(new RadialCollection);
        ao_->build(p_ucell_->ntype, charges, nmax, strategy_);
        ao_->set_transformer(sbt);
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
    }
    // neighbor list search
    scan_supercell();
    // allocate memory for ovlp_ao_nao_R_ and ovlp_ao_nao_k_
    nchi_ = ao_->nchi();
    nphi_ = nao_->nchi();
    nR_ = supercells_.size();
    allocate_ovlps();
}

void toQO::cal_ovlp_ao_nao_R(const int iR)
{
    double rcut_max = std::max(nao_->rcut_max(), ao_->rcut_max());
    int ngrid = int(rcut_max / 0.01) + 1;
    double cutoff = 2*rcut_max;
    nao_->set_uniform_grid(true, ngrid, cutoff, 'i', true);
    ao_->set_uniform_grid(true, ngrid, cutoff, 'i', true);

    overlap_calculator_->tabulate(*ao_, *nao_, 'S', ngrid, rcut_max);
    
    int irow = 0; int icol = 0;
    for(int it = 0; it < p_ucell_->ntype; it++)
    {
    // FOR EACH TYPE it
        for(int ia = 0; ia < p_ucell_->atoms[it].na; ia++)
        {
    // FOR EACH ATOM ia OF PRESENT TYPE it
            for(int jt = 0; jt < p_ucell_->ntype; jt++)
            {
    // FOR ANOTHER TYPE jt
                for(int ja = 0; ja < p_ucell_->atoms[jt].na; ja++)
                {
    // FOR ANOTHER ATOM ja OF ANOTHER TYPE jt
    // calculate displacement vector here
                    ModuleBase::Vector3<double> rij = p_ucell_->atoms[it].tau[ia] - p_ucell_->atoms[jt].tau[ja];
    // EXTRACTING RADIAL ORBITAL INFORMATION, resolution (l, izeta)
                    int lmaxi = atom_database_.principle_quantum_number[p_ucell_->atoms[it].ncpp.psd] - 1;
                    for(int li = 0; li <= lmaxi; li++)
                    {
                        int nzetai = ao_->nzeta(it, li);
                        for(int izetai = 0; izetai < nzetai; izetai++)
                        {
    // FOR EACH RADIAL ORBITAL OF ATOM itia (li, izetai)
                            int lmaxj = p_ucell_->atoms[jt].nwl;
                            for(int lj = 0; lj <= lmaxj; lj++)
                            {
                                int nzetaj = nao_->nzeta(jt, lj);
                                for(int izetaj = 0; izetaj < nzetaj; izetaj++)
                                {
    // FOR EACH RADIAL ORBITAL OF ATOM jtja (lj, izetaj)
    // GET ORBITAL-PAIR (li, izetai)-(lj, izetaj)
    // LOOP OVER Ylm...
                                    for(int mi = -li; mi <= li; mi++)
                                    {
                                        for(int mj = -lj; mj <= lj; mj++)
                                        {
    // Ylm OF ATOM itia, jtja DETERMINED HERE
                                            ModuleBase::Vector3<int> R = supercells_[iR];

                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}