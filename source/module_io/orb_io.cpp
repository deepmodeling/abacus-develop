#include "module_io/orb_io.h"
#include "module_base/tool_quit.h"
#ifdef __MPI
#include "module_base/parallel_common.h"
#endif

void ModuleIO::read_abacus_orb(std::ifstream& ifs,
                               std::string& elem,
                               double& ecut,
                               int& nr,
                               double& dr,
                               std::vector<int>& nzeta,
                               std::vector<std::vector<double>>& radials,
                               const int rank)
{
    nr = 0; // number of grid points
    dr = 0; // grid spacing
    int lmax = 0, nchi = 0; // number of radial functions
    std::vector<std::vector<int>> radial_map_; // build a map from [l][izeta] to 1-d array index
    std::string tmp;
    // first read the header
    if (rank == 0)
    {
        if (!ifs.is_open())
        {
            ModuleBase::WARNING_QUIT("AtomicRadials::read_abacus_orb", "Couldn't open orbital file.");
        }
        while (ifs >> tmp)
        {
            if (tmp == "Element")
            {
                ifs >> elem;
            }
            else if (tmp == "Cutoff(Ry)")
            {
                ifs >> ecut;
            }
            else if (tmp == "Lmax")
            {
                ifs >> lmax;
                nzeta.resize(lmax + 1);
                for (int l = 0; l <= lmax; ++l)
                {
                    ifs >> tmp >> tmp >> tmp >> nzeta[l];
                }
            }
            else if (tmp == "Mesh")
            {
                ifs >> nr;
                continue;
            }
            else if (tmp == "dr")
            {
                ifs >> dr;
                break;
            }
        }
        radial_map_.resize(lmax + 1);
        for (int l = 0; l <= lmax; ++l)
        {
            radial_map_[l].resize(nzeta[l]);
        }
        int ichi = 0;
        for (int l = 0; l <= lmax; ++l)
        {
            for (int iz = 0; iz < nzeta[l]; ++iz)
            {
                radial_map_[l][iz] = ichi++; // return the value of ichi, then increment
            }
        }
        nchi = ichi; // total number of radial functions
        radials.resize(nchi);
        std::for_each(radials.begin(), radials.end(), [nr](std::vector<double>& v) { v.resize(nr); });
    }

    // broadcast the header information
#ifdef __MPI
    Parallel_Common::bcast_string(elem);
    Parallel_Common::bcast_double(ecut);
    Parallel_Common::bcast_int(lmax);
    Parallel_Common::bcast_int(nchi);
    Parallel_Common::bcast_int(nr);
    Parallel_Common::bcast_double(dr);
#endif

    // then adjust the size of the vectors
    if (rank != 0)
    {
        nzeta.resize(lmax + 1);
        radials.resize(nchi);
        std::for_each(radials.begin(), radials.end(), [nr](std::vector<double>& v) { v.resize(nr); });
    }
    // broadcast the number of zeta functions for each angular momentum
#ifdef __MPI
    Parallel_Common::bcast_int(nzeta.data(), lmax + 1);
#endif

    // read the radial functions by rank0
    int ichi = 0;
    for (int i = 0; i != nchi; ++i)
    {
        if (rank == 0)
        {
            int l, izeta;
            ifs >> tmp >> tmp >> tmp;
            ifs >> tmp >> l >> izeta;
            ichi = radial_map_[l][izeta];
            for (int ir = 0; ir != nr; ++ir)
            {
                ifs >> radials[ichi][ir];
            }
        }
    // broadcast the radial functions
#ifdef __MPI
        Parallel_Common::bcast_int(ichi); // let other ranks know where to store the radial function
        Parallel_Common::bcast_double(radials[ichi].data(), nr);
#endif
    }
}

void ModuleIO::write_abacus_orb(std::ofstream& ofs,
                                const std::string& elem,
                                const double& ecut,
                                const int nr,
                                const double dr,
                                const std::vector<int>& nzeta,
                                const std::vector<std::vector<double>>& radials,
                                const int rank)
{
    const std::vector<std::string> spec = {"S", "P", "D", "F", "G", "H", "I", "J", "K"};

    if (rank == 0)
    {
        if (!ofs.is_open())
        {
            ModuleBase::WARNING_QUIT("AtomicRadials::write_abacus_orb", "Couldn't open orbital file.");
        }
        const int lmax = nzeta.size() - 1;

        for (int i = 0; i < 75; ++i)
        {
            ofs << "-";
        }
        ofs << std::endl;
        // left aligned
        ofs << std::left << std::setw(28) << "Element" << elem << std::endl;
        ofs << std::left << std::setw(28) << "Energy Cutoff(Ry)" << ecut << std::endl;
        // rcut .1f, not scientific
        ofs << std::left << std::setw(28) << "Radius Cutoff(a.u.)" 
            << std::fixed << std::setprecision(1) << dr * (nr - 1) << std::endl;
        ofs << std::left << std::setw(28) << "Lmax" << lmax << std::endl;
        for (int l = 0; l != nzeta.size(); ++l)
        {
            std::string title = "Number of " + spec[l] + "orbital-->";
            ofs << std::left << std::setw(28) << title << nzeta[l] << std::endl;
        }
        for (int i = 0; i < 75; ++i)
        {
            ofs << "-";
        }
        ofs << std::endl;
        ofs << "SUMMARY  END\n\n";
        ofs << std::left << std::setw(28) << "Mesh" << nr << std::endl;
        ofs << std::left << std::setw(28) << "dr" << std::fixed << std::setprecision(2) << dr << std::endl;

        int ichi = 0;
        for (int l = 0; l <= lmax; l++)
        {
            for (int izeta = 0; izeta < nzeta[l]; izeta++)
            {
                ofs << std::right << std::setw(20) << "Type"
                    << std::right << std::setw(20) << "L"
                    << std::right << std::setw(20) << "N" << std::endl;
                ofs << std::right << std::setw(20) << 0
                    << std::right << std::setw(20) << l
                    << std::right << std::setw(20) << izeta;
                for (int i = 0; i < nr; i++)
                {
                    if (i % 4 == 0)
                    {
                        ofs << std::endl;
                    }
                    ofs << std::left << std::setw(22) << std::setprecision(14) 
                        << std::scientific << radials[ichi][i];
                }
                ofs << std::endl;
                ichi++;
            }
        }
    }
    // ofs.close(); // like read_abacus_orb, who opens it, who closes it
}
