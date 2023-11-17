#include "module_io/to_qo.h"
#include "module_basis/module_nao/two_center_integrator.h"


toQO::toQO(std::string qo_basis)
{
    qo_basis_ = qo_basis;
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
        ao_->build(p_ucell_->ntype, charges, nmax);
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
}

void toQO::cal_ovlp_ao_nao()
{
    double rcut_max = std::max(nao_->rcut_max(), ao_->rcut_max());
    int ngrid = int(rcut_max / 0.01) + 1;
    nao_->set_uniform_grid(true, ngrid, rcut_max, 'i', true);
    ao_->set_uniform_grid(true, ngrid, rcut_max, 'i', true);
    for(int ikpt = 0; ikpt < nkpts_; ikpt++)
    {
        // calculate overlap matrix between atomic orbitals and numerical atomic orbitals for each kpoint

    }
}