#include "module_basis/module_nao/hydrogen_radials.h"
#include "module_base/global_variable.h"
#include <iostream>
#include <algorithm>

HydrogenRadials& HydrogenRadials::operator=(const HydrogenRadials& rhs)
{
    RadialSet::operator=(rhs);
    return *this;
}

void HydrogenRadials::build(const int itype,
                            const double charge,
                            const int nmax,
                            const double rcut,
                            const double dr,
                            const int rank,
                            std::ofstream* ptr_log)
{
    cleanup();
    itype_ = itype;
    if(rcut != GlobalV::qo_rcut)
    {
        #ifdef __MPI
        if(GlobalV::MY_RANK == 0)
        {
        #endif
            GlobalV::ofs_running << "Quasiatomic orbital analysis: use rcut = " 
                                 << GlobalV::qo_rcut 
                                 << " a.u., which is different from the default setting 10 a.u."
                                 << std::endl;
        #ifdef __MPI
        }
        #endif
    }
    #ifdef __MPI
    if(GlobalV::MY_RANK == 0)
    {
    #endif
        std::cout << "--- ENABLE QUASIATOMIC ORBITAL (QO) ANALYSIS ---" << std::endl;
    #ifdef __MPI
    }
    #endif
    generate_hydrogen_radials(charge, nmax, rcut, dr, rank, ptr_log);
}

void HydrogenRadials::generate_hydrogen_radials(const double charge,
                                                const int nmax,
                                                const double rcut,
                                                const double dr,
                                                const int rank,
                                                std::ofstream* ptr_log)
{
    // because all of these are parallelized, no need to bcast
    double a0 = 1.0; // Bohr radius

    lmax_ = nmax - 1;
    // count orbitals
    nchi_ = nmax * (nmax + 1) / 2;
    /*
    1 1s
    2 2s 2p
    3 3s 3p 3d
    */
    nzeta_ = new int[nmax];
    for(int i = 0; i != nmax; ++i)
    {
        nzeta_[i] = nmax - i;
    }
    nzeta_max_ = nmax;
    indexing();

    int ngrid = static_cast<int>(rcut / dr) + 1;
    double* rvalue = new double[ngrid];
    double* rgrid = new double[ngrid];
    for(int ir = 0; ir != ngrid; ++ir)
    {
        rgrid[ir] = ir * dr;
    }

    chi_ = new NumericalRadial[nchi_];

    int ichi = 0;
    for(int l_ = 0; l_ != nmax; ++l_)
    {
        for(int n_ = 1; n_ <= nmax; ++n_)
        {
            if(n_ <= l_)
            {
                // the loop sequence is in this way is because the multiplicity of l
                // is actually nmax - l, so if n < l, it is unphysical
                continue;
            }
            double norm_factor = sqrt(
                4.0*std::pow(charge, 3)*
                static_cast<double>(this->assoc_laguerre_.factorial(n_ - l_ - 1)) /
                std::pow(n_, 4)*
                static_cast<double>(this->assoc_laguerre_.factorial(n_ + l_)) /
                std::pow(a0, 3)
            );
            for(int ir = 0; ir != ngrid; ++ir)
            {
                // Bohr radius is 1.0
                double rho = 2.0 * rgrid[ir] * charge / n_ / a0;
                rvalue[ir] = norm_factor * pow(rho, l_) * exp(-rho/2.0) * this->assoc_laguerre_.value(
                    n_, l_, rho
                );
            }
            // what izeta, symbol and itype should be?
            chi_[ichi].build(l_, true, ngrid, rgrid, rvalue, 0, n_ - 1 - l_, "", itype_, false);
            chi_[ichi].normalize();
            ++ichi;
        }
    }

    delete[] rvalue;
    delete[] rgrid;
}