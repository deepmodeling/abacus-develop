#include "module_basis/module_nao/atomic_radials.h"

#include "module_base/math_integral.h"
#include "module_base/parallel_common.h"
#include "module_base/tool_quit.h"

// FIXME: should update with pyabacus
// #include "module_io/orb_io.h"

#include "module_base/projgen.h"

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

    if (rank == 0)
    {
        /*
         * read the header & grid information, including
         *
         * 1. element symbol --> symbol_
         * 2. energy cutoff --> orb_ecut_
         * 3. maximum angular momentum --> lmax_
         * 4. number of radial functions for each angular momentum --> nzeta_
         * 5. number of grid points --> ngrid
         * 6. grid spacing --> dr
         *                                                                              */
        while (ifs >> tmp)
        {
            if (tmp == "Element")
            {
                ifs >> symbol_;
            }
            else if (tmp == "Cutoff(Ry)")
            {
                ifs >> orb_ecut_;
            }
            else if (tmp == "Lmax")
            {
                ifs >> lmax_;
#ifdef __DEBUG
                assert(lmax_ >= 0);
#endif
                nzeta_ = new int[lmax_ + 1];
                for (int l = 0; l <= lmax_; ++l)
                {
                    ifs >> tmp >> tmp >> tmp >> nzeta_[l]; // skip "Number" "of" "Xorbital-->"
                }
            }
            else if (tmp == "Mesh")
            {
                ifs >> ngrid;
                continue;
            }
            else if (tmp == "dr")
            {
                ifs >> dr;
                break;
            }
        }

        /*
         * calculate:
         *
         * 1. the total number of radial functions --> nchi_
         * 2. maximum number of radial functions for each angular momentum --> nzeta_max_
         * 3. a map from (l, izeta) to 1-d array index in chi_
         *                                                                              */
        nchi_ = 0;
        for (int l = 0; l <= lmax_; ++l)
        {
            nchi_ += nzeta_[l];
        }
        nzeta_max_ = *std::max_element(nzeta_, nzeta_ + lmax_ + 1);
        indexing(); // build index_map_
    }

#ifdef __MPI
    Parallel_Common::bcast_string(symbol_);
    Parallel_Common::bcast_double(orb_ecut_);
    Parallel_Common::bcast_int(lmax_);

    Parallel_Common::bcast_int(nchi_);
    Parallel_Common::bcast_int(nzeta_max_);

    Parallel_Common::bcast_int(ngrid);
    Parallel_Common::bcast_double(dr);
#endif

    if (rank != 0)
    {
        nzeta_ = new int[lmax_ + 1];
        index_map_ = new int[(lmax_ + 1) * nzeta_max_];
    }

#ifdef __MPI
    Parallel_Common::bcast_int(nzeta_, lmax_ + 1);
    Parallel_Common::bcast_int(index_map_, (lmax_ + 1) * nzeta_max_);
#endif

    double* rvalue = new double[ngrid];
    double* rgrid = new double[ngrid];
    for (int ir = 0; ir != ngrid; ++ir)
    {
        rgrid[ir] = ir * dr;
    }
    chi_ = new NumericalRadial[nchi_];

    // record whether an orbital has been read or not
    bool* is_read = new bool[nchi_];
    for (int i = 0; i != nchi_; ++i)
    {
        is_read[i] = false;
    }

    for (int l = 0; l <= lmax_; ++l)
    {
        for (int izeta = 0; izeta < nzeta_[l]; ++izeta)
        {
            if (rank == 0)
            {
                /*
                 * read the orbital information, including
                 *
                 * 1. angular momentum
                 * 2. zeta number
                 * 3. values on the grid
                 *                                                                              */
                while (ifs.good())
                {
                    while (ifs >> tmp)
                    {
                        if (tmp == "N")
                        {
                            break;
                        }
                    }
                    int read_l, read_izeta;
                    ifs >> tmp >> read_l >> read_izeta;
                    if (l == read_l && izeta == read_izeta)
                    {
                        break;
                    }
                }

                for (int ir = 0; ir != ngrid; ++ir)
                {
                    ifs >> rvalue[ir];
                }
            }

#ifdef __MPI
            Parallel_Common::bcast_double(rvalue, ngrid);
#endif
#ifdef __DEBUG
            assert(index(l, izeta) >= 0 && index(l, izeta) < nchi_);
            assert(!is_read[index(l, izeta)]);
#endif
            is_read[index(l, izeta)] = true;

            // skip the initialization of sbt_ in this stage
            chi_[index(l, izeta)].build(l, true, ngrid, rgrid, rvalue, 0, izeta, symbol_, itype_, false);
            chi_[index(l, izeta)].normalize();
        }
    }

    delete[] is_read;
    delete[] rvalue;
    delete[] rgrid;
}

// FIXME: should update with pyabacus
// void AtomicRadials::read_abacus_orb(std::ifstream& ifs, std::ofstream* ptr_log, const int rank)
// {
//     /*
//      * Read the orbital file.
//      *
//      * For orbital file format, see
//      * (new) abacus-develop/tools/SIAB/PyTorchGradient/source/IO/print_orbital.py
//      * (old) abacus-develop/tools/SIAB/SimulatedAnnealing/source/src_spillage/Plot_Psi.cpp
//      *                                                                                  */
//     int ngrid = 0; // number of grid points
//     double dr = 0; // grid spacing
//     std::string tmp;

//     int nr;
//     std::vector<int> nzeta;
//     std::vector<std::vector<double>> radials;

//     ModuleIO::read_abacus_orb(ifs, symbol_, orb_ecut_, nr, dr, nzeta, radials, rank);

//     lmax_ = nzeta.size() - 1;
//     nzeta_ = new int[lmax_ + 1];
//     std::copy(nzeta.begin(), nzeta.end(), nzeta_);

//     nchi_ = std::accumulate(nzeta.begin(), nzeta.end(), 0);
//     nzeta_max_ = *std::max_element(nzeta.begin(), nzeta.end());

//     indexing();
    
//     std::vector<double> rgrid(nr);
//     std::iota(rgrid.begin(), rgrid.end(), 0);
//     std::for_each(rgrid.begin(), rgrid.end(), [dr](double& r) { r *= dr; });
//     chi_ = new NumericalRadial[nchi_];
//     int ichi = 0;
//     for (int l = 0; l <= lmax_; ++l)
//     {
//         for (int izeta = 0; izeta < nzeta[l]; ++izeta)
//         {
//             chi_[index(l, izeta)].build(l, true, nr, rgrid.data(), radials[ichi].data(), 0, izeta, symbol_, itype_, false);
//             chi_[index(l, izeta)].normalize();
//             ++ichi;
//         }
//     }
// }
