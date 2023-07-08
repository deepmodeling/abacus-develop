#include "module_basis/module_nao/atomic_radials.h"

#include <fstream>
#include <iostream>
#include <string>

#include "module_base/parallel_common.h"
#include "module_base/tool_quit.h"

AtomicRadials& AtomicRadials::operator=(const AtomicRadials& rhs)
{
    RadialSet::operator=(rhs);
    orb_ecut_ = rhs.orb_ecut_;
    return *this;
}

void AtomicRadials::build(const std::string& file, const int itype, std::ofstream* ptr_log, const int rank)
{
    // deallocates all arrays and reset variables
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
        ModuleBase::WARNING_QUIT("AtomicRadials::read", "Couldn't open orbital file: " + file);
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

    if (rank == 0)
    {
        ifs.close();
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
                assert(lmax_ >= 0);
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
        indexing(); // calculate nzeta_max_ and build index_map_
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

    int l = 0;
    int izeta = 0;
    for (int i = 0; i != nchi_; ++i)
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
            // ifs >> tmp >> tmp >> tmp; // skip "Type" "L" "N"
            ifs >> tmp >> tmp >> tmp;
            assert(tmp == "N");

            ifs >> tmp >> l >> izeta;
            assert(l >= 0 && l <= lmax_);
            assert(izeta >= 0 && izeta < nzeta_[l]);

            for (int ir = 0; ir != ngrid; ++ir)
            {
                ifs >> rvalue[ir];
            }
        }

#ifdef __MPI
        Parallel_Common::bcast_int(l);
        Parallel_Common::bcast_int(izeta);
        Parallel_Common::bcast_double(rvalue, ngrid);
#endif
        assert(index(l, izeta) >= 0 && index(l, izeta) < nchi_);
        assert(!is_read[index(l, izeta)]);
        is_read[index(l, izeta)] = true;

        chi_[index(l, izeta)].build(l, true, ngrid, rgrid, rvalue, 0, izeta, symbol_, itype_);
        chi_[index(l, izeta)].normalize();
    }

    delete[] is_read;
    delete[] rvalue;
    delete[] rgrid;
}
