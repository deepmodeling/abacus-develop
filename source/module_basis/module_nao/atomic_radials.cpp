#include "module_basis/module_nao/atomic_radials.h"

#include "module_base/math_integral.h"
#include "module_base/parallel_common.h"
#include "module_base/tool_quit.h"
#include "module_io/orb_io.h"
#include "projgen.h"

#include <fstream>
#include <iostream>
#include <string>
#include <numeric>

AtomicRadials& AtomicRadials::operator=(const AtomicRadials& rhs)
{
    RadialSet::operator=(rhs);
    orb_ecut_ = rhs.orb_ecut_;
    return *this;
}

void AtomicRadials::build(const std::string& file, const int itype, std::ofstream* ptr_log, const int rank)
{
    // deallocates all arrays and reset variables (excluding sbt_)
    cleanup();

    std::ifstream ifs;
    bool is_open = false;

    if (rank == 0)
    {
        ifs.open(file);
        is_open = ifs.is_open();
    }

#ifdef __MPI
    Parallel_Common::bcast_bool(is_open);
#endif

    if (!is_open)
    {
        ModuleBase::WARNING_QUIT("AtomicRadials::build", "Couldn't open orbital file: " + file);
    }

    if (ptr_log)
    {
        (*ptr_log) << "\n\n\n\n";
        (*ptr_log) << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        (*ptr_log) << " |                                                                   |" << std::endl;
        (*ptr_log) << " |               SETUP NUMERICAL ATOMIC ORBITALS                     |" << std::endl;
        (*ptr_log) << " |                                                                   |" << std::endl;
        (*ptr_log) << " | Orbital information includes the cutoff radius, angular momentum, |" << std::endl;
        (*ptr_log) << " | zeta number and numerical values on a radial grid.                |" << std::endl;
        (*ptr_log) << " |                                                                   |" << std::endl;
        (*ptr_log) << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        (*ptr_log) << "\n\n\n\n";
    }

    itype_ = itype;
    read_abacus_orb(ifs, ptr_log, rank);
    set_rcut_max();

    if (rank == 0)
    {
        ifs.close();
    }
}

void AtomicRadials::build(RadialSet* const other, const int itype, const double rcut)
{
    this->symbol_ = other->symbol();
    this->lmax_ = other->lmax();
    this->nchi_ = other->nchi();
    this->nzeta_max_ = other->nzeta_max();
    this->itype_ = itype;
    this->symbol_ = other->symbol();
    this->nzeta_ = new int[this->lmax_ + 1];
    for (int l = 0; l <= this->lmax_; ++l)
    {
        this->nzeta_[l] = other->nzeta(l);
    }
    this->indexing();
    this->chi_ = new NumericalRadial[nchi_];
    for (int ichi = 0; ichi < this->nchi_; ichi++)
    {
        const int l = other->cbegin()[ichi].l();
        int ngrid = other->cbegin()[ichi].nr();
        const double* rgrid = other->cbegin()[ichi].rgrid();
        const double* rvalue = other->cbegin()[ichi].rvalue();
        const int izeta = other->cbegin()[ichi].izeta();
        // if the cutoff radius is larger than the original one, just copy the orbitals
        if (rcut >= other->cbegin()[ichi].rcut())
        {
            this->chi_[ichi].build(l, true, ngrid, rgrid, rvalue, 0, izeta, symbol_, itype, false);
        }
        else
        {
            // call smoothgen to modify the orbitals to the local projections
            std::vector<double> rvalue_new;
            smoothgen(ngrid, rgrid, rvalue, rcut, rvalue_new);
            ngrid = rvalue_new.size();
            // projgen(l, ngrid, rgrid, rvalue, rcut, 20, rvalue_new);
            // build the new on-site orbitals
            this->chi_[ichi].build(l, true, ngrid, rgrid, rvalue_new.data(), 0, izeta, symbol_, itype, false);
        }
    }
}

void AtomicRadials::read_abacus_orb(std::ifstream& ifs, std::ofstream* ptr_log, const int rank)
{
    /*
     * Read the orbital file.
     *
     * For orbital file format, see
     * (new) abacus-develop/tools/SIAB/PyTorchGradient/source/IO/print_orbital.py
     * (old) abacus-develop/tools/SIAB/SimulatedAnnealing/source/src_spillage/Plot_Psi.cpp
     *                                                                                  */
    int ngrid = 0; // number of grid points
    double dr = 0; // grid spacing
    std::string tmp;

    int nr;
    std::vector<int> nzeta;
    std::vector<std::vector<double>> radials;

    ModuleIO::read_abacus_orb(ifs, symbol_, orb_ecut_, nr, dr, nzeta, radials, rank);

    lmax_ = nzeta.size() - 1;
    nzeta_ = new int[lmax_ + 1];
    std::copy(nzeta.begin(), nzeta.end(), nzeta_);

    nchi_ = std::accumulate(nzeta.begin(), nzeta.end(), 0);
    nzeta_max_ = *std::max_element(nzeta.begin(), nzeta.end());

    indexing();
    
    std::vector<double> rgrid(nr);
    std::iota(rgrid.begin(), rgrid.end(), 0);
    std::for_each(rgrid.begin(), rgrid.end(), [dr](double& r) { r *= dr; });
    chi_ = new NumericalRadial[nchi_];
    int ichi = 0;
    for (int l = 0; l <= lmax_; ++l)
    {
        for (int izeta = 0; izeta < nzeta[l]; ++izeta)
        {
            chi_[index(l, izeta)].build(l, true, nr, rgrid.data(), radials[ichi].data(), 0, izeta, symbol_, itype_, false);
            chi_[index(l, izeta)].normalize();
            ++ichi;
        }
    }
}
