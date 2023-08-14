#include "module_basis/module_nao/beta_radials.h"

#include <regex>

#include "module_base/parallel_common.h"
#include "module_base/tool_quit.h"

void BetaRadials::build(const std::string& file, const int itype, std::ofstream* ptr_log, const int rank)
{
    /*
     * Build a BetaRadials object from beta functions read from a pseudopotential file
     *
     * NOTE: only the rank == 0 process reads the file
     *                                                                                      */
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
        ModuleBase::WARNING_QUIT("BetaRadials::read", "Couldn't open pseudopotential file: " + file);
    }

    if (ptr_log)
    {
        (*ptr_log) << "\n\n\n\n";
        (*ptr_log) << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        (*ptr_log) << " |                                                                   |" << std::endl;
        (*ptr_log) << " |               SETUP BETA PROJECTORS                               |" << std::endl;
        (*ptr_log) << " |                                                                   |" << std::endl;
        (*ptr_log) << " | Projector information includes the cutoff radius, angular         |" << std::endl;
        (*ptr_log) << " | momentum, zeta number and numerical values on a radial grid.      |" << std::endl;
        (*ptr_log) << " |                                                                   |" << std::endl;
        (*ptr_log) << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        (*ptr_log) << "\n\n\n\n";
    }

    itype_ = itype;

    // identify the version of UPF file
    // UPF 2.0.1 format starts with <UPF version="2.0.1">
    int upf_version = 100;
    if (rank == 0)
    {
        std::string first_line;
        std::getline(ifs, first_line);
        if (first_line.find("2.0.1") != std::string::npos)
        {
            upf_version = 201;
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_int(upf_version);
#endif

    switch (upf_version)
    {
    case 100:
        read_beta_upf100(ifs, ptr_log, rank);
        break;
    case 201:
        read_beta_upf201(ifs, ptr_log, rank);
        break;
    default: /* not supposed to happen */;
    }

    if (rank == 0)
    {
        ifs.close();
    }

    for (int i = 0; i < nchi_; i++)
    {
        chi_[i].set_transformer(sbt_, 0);
    }
}

void BetaRadials::read_beta_upf100(std::ifstream& ifs, std::ofstream* ptr_log, const int rank)
{
    /*
     * Read the nonlocal beta functions from the pseudopotential file (in old UPF format)
     *                                                                                      */
    std::stringstream ss;
    std::string tmp;
    std::string line;
    int ngrid_max = 0;
    bool is_good = false;

    if (rank == 0)
    {
        /*
         * Read the header, including
         *
         * 1. Element symbol
         * 2. Maximum angular momentum
         * 3. Number of radial grid points
         * 4. Number of beta functions
         *                                                                                  */
        while (std::getline(ifs, line))
        {
            if (line.find("Element") != std::string::npos)
            {
                ss.str("");
                ss << line;
                ss >> symbol_;
            }
            else if (line.find("Max angular momentum component") != std::string::npos)
            {
                ss.str("");
                ss << line;
                ss >> lmax_;
            }
            else if (line.find("Number of points in mesh") != std::string::npos)
            {
                ss.str("");
                ss << line;
                ss >> ngrid_max;
            }
            else if (line.find("Number of Projectors") != std::string::npos)
            {
                ss.str("");
                ss << line;
                ss >> tmp >> nchi_; // nchi_ is the number of beta functions
            }

            // should obtain valid nchi_, symbol_ & ngrid_max upon reaching the end of header
            if (line.find("</PP_HEADER>") != std::string::npos)
            {
                // lmax could be -1 when there is no beta projectors, see, e.g., H.pz-vbc.UPF
                is_good = (nchi_ >= 0) && symbol_.size() && (ngrid_max > 0);
                break;
            }
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_bool(is_good);
    Parallel_Common::bcast_string(symbol_);
    Parallel_Common::bcast_int(lmax_);
    Parallel_Common::bcast_int(nchi_);
    Parallel_Common::bcast_int(ngrid_max);
#endif

    if (!is_good)
    {
        ModuleBase::WARNING_QUIT("BetaRadials::read_beta_upf100", "PP_HEADER error");
    }

    // In case some pseudopotential file does not have any beta function
    if (nchi_ == 0)
    {
        return;
    }

    double* rgrid = new double[ngrid_max];
    if (rank == 0)
    {
        /*
         * Read the radial grid
         *                                                                                  */
        while (ifs >> tmp)
        {
            if (tmp == "<PP_R>")
            {
                break;
            }
        }
        assert(!ifs.eof());

        for (int ir = 0; ir != ngrid_max; ++ir)
        {
            ifs >> rgrid[ir];
        }

        ifs >> tmp;
        assert(tmp == "</PP_R>");
    }

#ifdef __MPI
    Parallel_Common::bcast_double(rgrid, ngrid_max);
#endif

    assert(lmax_ >= 0);
    nzeta_ = new int[lmax_ + 1];
    for (int l = 0; l <= lmax_; ++l)
    {
        nzeta_[l] = 0;
    }

    // rbeta, l & ngrid are directly read from file
    double* rbeta = new double[ngrid_max]; // beta function is given as r*beta(r)
    int l = 0;
    int ngrid = 0;

    int l_last = -1;
    int izeta = 0;

    chi_ = new NumericalRadial[nchi_];
    for (int i = 0; i != nchi_; ++i)
    {
        if (rank == 0)
        {
            /*
             * Read the beta functions, including
             *
             * 1. Angular momentum
             * 2. Number of actual radial grid points (should be smaller than ngrid_max)
             * 3. Numerical values (r*beta(r)) on the radial grid
             *
             * and record the zeta number (izeta).
             *                                                                              */
            while (std::getline(ifs, line))
            {
                if (line.find("<PP_BETA>") != std::string::npos)
                {
                    break;
                }
            }

            ifs >> tmp >> l;
            ifs >> tmp >> tmp; // skip "Beta" "L"
            ifs >> ngrid;

            if (l == l_last)
            {
                izeta += 1;
            }
            else
            {
                izeta = 0;
                l_last = l;
            }

            for (int ir = 0; ir != ngrid; ++ir)
            {
                ifs >> rbeta[ir];
            }

            nzeta_[l] += 1;

            ifs >> tmp;
            assert(tmp == "</PP_BETA>");
        } // rank == 0

#ifdef __MPI
        Parallel_Common::bcast_int(l);
        Parallel_Common::bcast_int(ngrid);
        Parallel_Common::bcast_int(izeta);
        Parallel_Common::bcast_double(rbeta, ngrid);
#endif
        chi_[i].build(l, true, ngrid, rgrid, rbeta, 1, izeta, symbol_, itype_);
    }

#ifdef __MPI
    Parallel_Common::bcast_int(nzeta_, lmax_ + 1);
#endif

    indexing();

    delete[] rgrid;
    delete[] rbeta;
}

void BetaRadials::read_beta_upf201(std::ifstream& ifs, std::ofstream* ptr_log, const int rank)
{
    /*
     * Read the nonlocal beta functions from the pseudopotential file (in UPF 2.0.1 format)
     *                                                                                      */
    std::string tmp;
    std::string line;
    int ngrid_max = 0;
    bool is_good = false;

    if (rank == 0)
    {
        /*
         * Read the header, including
         *
         * 1. Element symbol
         * 2. Maximum angular momentum
         * 3. Number of radial grid points
         * 4. Number of beta functions
         *                                                                                  */
        while (std::getline(ifs, line))
        {
            if (line.find("element=") != std::string::npos)
            {
                symbol_ = trim201(line);
            }
            else if (line.find("l_max=") != std::string::npos)
            {
                lmax_ = std::stoi(trim201(line));
            }
            else if (line.find("mesh_size=") != std::string::npos)
            {
                ngrid_max = std::stoi(trim201(line));
            }
            else if (line.find("number_of_proj=") != std::string::npos)
            {
                nchi_ = std::stoi(trim201(line));
            }

            if (line.find("/>") != std::string::npos)
            {
                is_good = (nchi_ >= 0) && (ngrid_max >= 0) && symbol_.size();
                break;
            }
        }
    }

#ifdef __MPI
    Parallel_Common::bcast_bool(is_good);
    Parallel_Common::bcast_string(symbol_);
    Parallel_Common::bcast_int(lmax_);
    Parallel_Common::bcast_int(nchi_);
    Parallel_Common::bcast_int(ngrid_max);
#endif

    if (!is_good)
    {
        ModuleBase::WARNING_QUIT("BetaRadials::read_beta_upf201", "PP_HEADER error");
    }

    // In case some pseudopotential file does not have any beta function
    if (nchi_ == 0)
    {
        return;
    }

    double* rgrid = new double[ngrid_max];
    if (rank == 0)
    {
        /*
         * Read the radial grid
         *                                                                                  */
        while (std::getline(ifs, line))
        {
            if (line.find("<PP_R") != std::string::npos)
            {
                break;
            }
        }
        assert(!ifs.eof());

        for (int ir = 0; ir != ngrid_max; ++ir)
        {
            ifs >> rgrid[ir];
        }

        ifs >> tmp;
        assert(tmp == "</PP_R>");
    }

#ifdef __MPI
    Parallel_Common::bcast_double(rgrid, ngrid_max);
#endif

    nzeta_ = new int[lmax_ + 1];
    for (int l = 0; l <= lmax_; ++l)
    {
        nzeta_[l] = 0;
    }

    // rbeta, l & ngrid are directly read from file
    double* rbeta = new double[ngrid_max];
    int l = 0;
    int ngrid = 0;

    int l_last = -1;
    int izeta = 0;

    chi_ = new NumericalRadial[nchi_];
    for (int i = 0; i != nchi_; ++i)
    {
        if (rank == 0)
        {
            /*
             * Read the beta functions, including
             *
             * 1. Angular momentum
             * 2. Number of actual radial grid points (should be smaller than ngrid_max)
             * 3. Numerical values on the radial grid
             *
             * and record the zeta number (izeta).
             *                                                                              */
            while (std::getline(ifs, line))
            {
                if (line.find("<PP_BETA") != std::string::npos)
                {
                    break;
                }
            }
            assert(!ifs.eof());

            // multiple key-value pairs may appear in the same line
            // and they may appear in the same line as "<PP_BETA"
            // see, e.g., Si.pz-n-nc.UPF
            do
            {
                if (line.find("angular_momentum") != std::string::npos)
                {
                    l = std::stoi(extract201(line, "angular_momentum"));
                }

                if (line.find("size") != std::string::npos)
                {
                    ngrid = std::stoi(extract201(line, "size"));
                }
                // neither "cutoff_radius_index" nor "cutoff_radius" is reliable!
                // the code will read all the values first and then reverse scan to determine the grid size

                if (line.find(">") != std::string::npos)
                {
                    is_good &= (l >= 0) && (l <= lmax_) && (ngrid > 0) && (ngrid <= ngrid_max);
                    break;
                }
            } while (std::getline(ifs, line));

            if (l == l_last)
            {
                izeta += 1;
            }
            else
            {
                izeta = 0;
                l_last = l;
            }

            for (int ir = 0; ir != ngrid; ++ir)
            {
                ifs >> rbeta[ir];
            }

            // reverse scan to determine the grid size
            for (; ngrid > 0; --ngrid)
            {
                if (std::abs(rbeta[ngrid - 1]) > 1e-12)
                {
                    break;
                }
            }

            nzeta_[l] += 1;

            ifs >> tmp;
            assert(tmp.find("</PP_BETA") != std::string::npos);
            is_good &= ifs.good() && (ngrid > 0);
        } // rank == 0

#ifdef __MPI
        Parallel_Common::bcast_bool(is_good);
        Parallel_Common::bcast_int(l);
        Parallel_Common::bcast_int(ngrid);
        Parallel_Common::bcast_int(izeta);
        Parallel_Common::bcast_double(rbeta, ngrid);
#endif

        if (!is_good)
        {
            ModuleBase::WARNING_QUIT("BetaRadials::read_beta_upf201", "PP_BETA error");
        }

        chi_[i].build(l, true, ngrid, rgrid, rbeta, 1, izeta, symbol_, itype_);
    }

#ifdef __MPI
    Parallel_Common::bcast_int(nzeta_, lmax_ + 1);
#endif

    indexing();

    delete[] rgrid;
    delete[] rbeta;
}

std::string BetaRadials::trim201(std::string const& str)
{

    // extract the substring between quotation marks (with whitespace trimmed)
    // str MUST contain exactly a pair of quotation marks

    std::string::size_type start = str.find('"');
    std::string::size_type end = str.find_last_of('"');
    std::string tmp = str.substr(start + 1, end - start - 1);

    if (tmp.length() == 0)
    {
        return tmp;
    }

    start = tmp.find_first_not_of(" \t");
    end = tmp.find_last_not_of(" \t");
    return tmp.substr(start, end + 1 - start);
}

std::string BetaRadials::extract201(std::string const& str, std::string const& keyword) {
    std::smatch match;
    std::string regex_string = ".*" + keyword + "=\" *(.+) *\".*";
    std::regex re(regex_string);
    std::regex_match(str, match, re);
    return match[1].str();
}
