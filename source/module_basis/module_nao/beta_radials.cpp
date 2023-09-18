#include "module_basis/module_nao/beta_radials.h"

#include <regex>

#include "module_base/parallel_common.h"
#include "module_base/tool_quit.h"
#include "module_base/global_variable.h"

void BetaRadials::build(const Numerical_Nonlocal& nl, const int itype, std::ofstream* const ptr_log)
{
    cleanup();

    itype_ = itype;
    assert(itype_ == nl.getType());

    symbol_ = nl.getLabel();
    lmax_ = nl.getLmax();
    nchi_ = nl.get_nproj();

    chi_ = new NumericalRadial[nchi_];
    nzeta_ = new int[lmax_ + 1];
    std::fill(nzeta_, nzeta_ + lmax_ + 1, 0);

    for (int ichi = 0; ichi != nchi_; ++ichi)
    {
        Numerical_Nonlocal_Lm& beta = nl.Proj[ichi];
        int l = beta.getL();
        chi_[ichi].build(l, true, beta.getNr(), beta.getRadial(), beta.getBeta_r(), 1, nzeta_[l], symbol_, itype_);
        nzeta_[l] += 1;
    }
    nzeta_max_ = *std::max_element(nzeta_, nzeta_ + lmax_ + 1);

    indexing();
    set_rcut_max();
}

//void BetaRadials::build(const std::string& file, const int itype, std::ofstream* ptr_log, const int rank)
//{
//    /*
//     * Build a BetaRadials object from beta functions read from a pseudopotential file
//     *
//     * NOTE: only the rank == 0 process reads the file
//     *                                                                                      */
//    cleanup();
//
//    std::ifstream ifs;
//    bool is_open = false;
//
//    if (rank == 0)
//    {
//        ifs.open(file);
//        is_open = ifs.is_open();
//    }
//
//#ifdef __MPI
//    Parallel_Common::bcast_bool(is_open);
//#endif
//
//    if (!is_open)
//    {
//        ModuleBase::WARNING_QUIT("BetaRadials::read", "Couldn't open pseudopotential file: " + file);
//    }
//
//    if (ptr_log)
//    {
//        (*ptr_log) << "\n\n\n\n";
//        (*ptr_log) << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
//        (*ptr_log) << " |                                                                   |" << std::endl;
//        (*ptr_log) << " |       SETUP BETA PROJECTORS FOR TWO-CENTER INTEGRALS              |" << std::endl;
//        (*ptr_log) << " |                                                                   |" << std::endl;
//        (*ptr_log) << " | Projector information includes the cutoff radius, angular         |" << std::endl;
//        (*ptr_log) << " | momentum, zeta number and numerical values on a radial grid.      |" << std::endl;
//        (*ptr_log) << " |                                                                   |" << std::endl;
//        (*ptr_log) << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
//        (*ptr_log) << "\n\n\n\n";
//    }
//
//    itype_ = itype;
//
//    // identify the version of UPF file
//    // UPF 2.0.1 format starts with <UPF version="2.0.1">
//    int upf_version = 100;
//    if (rank == 0)
//    {
//        std::string first_line;
//        std::getline(ifs, first_line);
//        if (first_line.find("2.0.1") != std::string::npos)
//        {
//            upf_version = 201;
//        }
//    }
//
//#ifdef __MPI
//    Parallel_Common::bcast_int(upf_version);
//#endif
//
//    switch (upf_version)
//    {
//    case 100:
//        read_beta_upf100(ifs, ptr_log, rank);
//        break;
//    case 201:
//        read_beta_upf201(ifs, ptr_log, rank);
//        break;
//    default: /* not supposed to happen */;
//    }
//
//    if (rank == 0)
//    {
//        ifs.close();
//    }
//
//    for (int i = 0; i < nchi_; i++)
//    {
//        chi_[i].set_transformer(sbt_, 0);
//    }
//}
//
//void BetaRadials::read_beta_upf100(std::ifstream& ifs, std::ofstream* ptr_log, const int rank)
//{
//    /*
//     * Read the nonlocal beta functions from the pseudopotential file (in old UPF format)
//     *                                                                                      */
//    std::stringstream ss;
//    std::string tmp;
//    std::string line;
//    int ngrid_max = 0;
//    bool is_good = false;
//
//    if (rank == 0)
//    {
//        /*
//         * Read the header, including
//         *
//         * 1. Element symbol
//         * 2. Maximum angular momentum
//         * 3. Number of radial grid points
//         * 4. Number of beta functions
//         *                                                                                  */
//        while (std::getline(ifs, line))
//        {
//            if (line.find("Element") != std::string::npos)
//            {
//                ss.str("");
//                ss << line;
//                ss >> symbol_;
//            }
//            else if (line.find("Max angular momentum component") != std::string::npos)
//            {
//                ss.str("");
//                ss << line;
//                ss >> lmax_;
//            }
//            else if (line.find("Number of points in mesh") != std::string::npos)
//            {
//                ss.str("");
//                ss << line;
//                ss >> ngrid_max;
//            }
//            else if (line.find("Number of Projectors") != std::string::npos)
//            {
//                ss.str("");
//                ss << line;
//                ss >> tmp >> nchi_; // nchi_ is the number of beta functions
//            }
//
//            // should obtain valid nchi_, symbol_ & ngrid_max upon reaching the end of header
//            if (line.find("</PP_HEADER>") != std::string::npos)
//            {
//                // lmax could be -1 when there is no beta projectors, see, e.g., H.pz-vbc.UPF
//                is_good = (nchi_ >= 0) && symbol_.size() && (ngrid_max > 0);
//                break;
//            }
//        }
//    }
//
//#ifdef __MPI
//    Parallel_Common::bcast_bool(is_good);
//    Parallel_Common::bcast_string(symbol_);
//    Parallel_Common::bcast_int(lmax_);
//    Parallel_Common::bcast_int(nchi_);
//    Parallel_Common::bcast_int(ngrid_max);
//#endif
//
//    if (!is_good)
//    {
//        ModuleBase::WARNING_QUIT("BetaRadials::read_beta_upf100", "PP_HEADER error");
//    }
//
//    // In case some pseudopotential file does not have any beta function
//    if (nchi_ == 0)
//    {
//        return;
//    }
//
//    double* rgrid = new double[ngrid_max];
//    if (rank == 0)
//    {
//        /*
//         * Read the radial grid
//         *                                                                                  */
//        while (ifs >> tmp)
//        {
//            if (tmp == "<PP_R>")
//            {
//                break;
//            }
//        }
//        assert(!ifs.eof());
//
//        for (int ir = 0; ir != ngrid_max; ++ir)
//        {
//            ifs >> rgrid[ir];
//        }
//
//        ifs >> tmp;
//        assert(tmp == "</PP_R>");
//    }
//
//#ifdef __MPI
//    Parallel_Common::bcast_double(rgrid, ngrid_max);
//#endif
//
//    assert(lmax_ >= 0);
//    nzeta_ = new int[lmax_ + 1];
//    for (int l = 0; l <= lmax_; ++l)
//    {
//        nzeta_[l] = 0;
//    }
//
//    // rbeta, l & ngrid are directly read from file
//    double* rbeta = new double[ngrid_max]; // beta function is given as r*beta(r)
//    int l = 0;
//    int ngrid = 0;
//
//    int l_last = -1;
//    int izeta = 0;
//
//    chi_ = new NumericalRadial[nchi_];
//    for (int i = 0; i != nchi_; ++i)
//    {
//        if (rank == 0)
//        {
//            /*
//             * Read the beta functions, including
//             *
//             * 1. Angular momentum
//             * 2. Number of actual radial grid points (should be smaller than ngrid_max)
//             * 3. Numerical values (r*beta(r)) on the radial grid
//             *
//             * and record the zeta number (izeta).
//             *                                                                              */
//            while (std::getline(ifs, line))
//            {
//                if (line.find("<PP_BETA>") != std::string::npos)
//                {
//                    break;
//                }
//            }
//
//            ifs >> tmp >> l;
//            ifs >> tmp >> tmp; // skip "Beta" "L"
//            ifs >> ngrid;
//
//            if (l == l_last)
//            {
//                izeta += 1;
//            }
//            else
//            {
//                izeta = 0;
//                l_last = l;
//            }
//
//            for (int ir = 0; ir != ngrid; ++ir)
//            {
//                ifs >> rbeta[ir];
//            }
//
//            nzeta_[l] += 1;
//
//            ifs >> tmp;
//            assert(tmp == "</PP_BETA>");
//        } // rank == 0
//
//#ifdef __MPI
//        Parallel_Common::bcast_int(l);
//        Parallel_Common::bcast_int(ngrid);
//        Parallel_Common::bcast_int(izeta);
//        Parallel_Common::bcast_double(rbeta, ngrid);
//#endif
//        chi_[i].build(l, true, ngrid, rgrid, rbeta, 1, izeta, symbol_, itype_);
//    }
//
//#ifdef __MPI
//    Parallel_Common::bcast_int(nzeta_, lmax_ + 1);
//#endif
//
//    indexing();
//
//    delete[] rgrid;
//    delete[] rbeta;
//}
//
///*======================================================================
// *
// *  Read beta functions (UPF 2.0.1) for two-center integrals
// *
// *======================================================================*/
//void BetaRadials::read_beta_upf201(std::ifstream& ifs, std::ofstream* ptr_log, const int rank)
//{
//    /*
//     * Read the nonlocal beta functions from the pseudopotential file (in UPF 2.0.1 format)
//     *                                                                                      */
//
//    // temporary variables
//    std::string tmpword;
//    std::string tmpline;
//
//    // status flag
//    bool is_good = true;
//
//
//    /*===========================================================
//     *
//     *                  Read the header
//     *
//     * Necessary information:
//     * 1. Element symbol
//     * 2. Maximum angular momentum
//     * 3. Number of radial grid points
//     * 4. Number of beta functions
//     * 5. Whether spin-orbit coupling is included
//     *
//     *===========================================================*/
//
//    // variables needed to be extracted from the header
//    // element & lmax read from file are directly used to update member variables symbol_ & lmax_
//    int ngrid_max = 0;
//    bool has_so = false;
//    int nbeta = 0; 
//    // NOTE: The number of PP_BETA from a pseudopotential file would not necessarily 
//    // be equal to the final number of beta projectors used in two-center integrals (nchi_).
//    // This would happen if a pseudopotential file contains spin-orbit information but
//    // lspinorb is not switched on, in which case SOC-related beta functions would be averaged.
//
//    if (rank == 0)
//    {
//        while (std::getline(ifs, tmpline))
//        {
//            if (tmpline.find("element=") != std::string::npos)
//            {
//                symbol_ = trim201(tmpline);
//            }
//            else if (tmpline.find("l_max=") != std::string::npos)
//            {
//                lmax_ = std::stoi(trim201(tmpline));
//            }
//            else if (tmpline.find("mesh_size=") != std::string::npos)
//            {
//                ngrid_max = std::stoi(trim201(tmpline));
//            }
//            else if (tmpline.find("number_of_proj=") != std::string::npos)
//            {
//                nbeta = std::stoi(trim201(tmpline));
//            }
//            else if (tmpline.find("has_so") != std::string::npos)
//            {
//                std::string has_so_str = trim201(tmpline);
//                if ( has_so_str == "T" || has_so_str == "TRUE" || has_so_str == "True" || has_so_str == "true" )
//                    has_so = true;
//            }
//
//            if (tmpline.find("/>") != std::string::npos)
//            {
//                is_good = (nchi_ >= 0) && (ngrid_max >= 0) && symbol_.size();
//                break;
//            }
//        }
//        printf("element = %s, lmax = %i, ngrid_max = %i, nchi = %i, has_so = %d, is_good = %d\n",
//                symbol_.c_str(), lmax_, ngrid_max, nbeta, has_so, is_good);
//    }
//
//#ifdef __MPI
//    Parallel_Common::bcast_bool(is_good);
//#endif
//    if (!is_good)
//    {
//        ModuleBase::WARNING_QUIT("BetaRadials::read_beta_upf201", "PP_HEADER error");
//    }
//
//#ifdef __MPI
//    Parallel_Common::bcast_string(symbol_);
//    Parallel_Common::bcast_int(lmax_);
//    Parallel_Common::bcast_int(ngrid_max);
//    Parallel_Common::bcast_int(nbeta);
//    Parallel_Common::bcast_bool(has_so);
//#endif
//
//    // It is an error if lspinorb is set to true but the pseudopotential file does not contain spin-orbit information
//    if (!has_so && GlobalV::LSPINORB)
//    {
//        ModuleBase::WARNING_QUIT("BetaRadials::read_beta_upf201",
//                "lspinorb is set to true but the pseudopotential file does not contain spin-orbit information");
//    }
//
//    // In case some pseudopotential file does not have any beta function
//    if (nbeta == 0)
//    {
//        return;
//    }
//
//    /*===========================================================
//     *
//     *              Read the real-space grid
//     *
//     *===========================================================*/
//    double* rgrid = new double[ngrid_max];
//    if (rank == 0)
//    {
//        /*
//         * Read the radial grid
//         *                                                                                  */
//        while (std::getline(ifs, tmpline))
//        {
//            if (tmpline.find("<PP_R") != std::string::npos)
//            {
//                break;
//            }
//        }
//        assert(!ifs.eof());
//
//        for (int ir = 0; ir != ngrid_max; ++ir)
//        {
//            ifs >> rgrid[ir];
//        }
//
//        ifs >> tmpword;
//        assert(tmpword == "</PP_R>");
//    }
//
//#ifdef __MPI
//    Parallel_Common::bcast_double(rgrid, ngrid_max);
//#endif
//
//    /*===========================================================
//     *
//     *          Read the beta functions (rank-0 only)
//     * 
//     * Information to read from file:
//     * 1. Angular momentum
//     * 2. r*beta(r) values on the radial grid
//     *
//     * Information to find out:
//     * 1. Number of radial grid points with non-zero r*beta(r) values
//     *
//     * Since an averaging over SOC-related beta functions might
//     * be needed, the broadcast is performed at this stage.
//     *
//     *===========================================================*/
//    double* rbeta = new double[nbeta * ngrid_max];
//    int* l = new int[nbeta];
//    int* ngrid = new int[nbeta]; // number of grid points with non-zero r*beta(r) values
//
//    if (rank == 0)
//    {
//        for (int ibeta = 0; ibeta != nbeta; ++ibeta)
//        {
//            while (std::getline(ifs, tmpline))
//            {
//                if (tmpline.find("<PP_BETA") != std::string::npos)
//                {
//                    break;
//                }
//            }
//            assert(!ifs.eof());
//
//            // extract the angular momentum and the number of grid points
//            // from the little header of each PP_BETA block.
//            //
//            // NOTE: multiple key-value pairs may appear in the same line
//            // as "<PP_BETA", see, e.g., Si.pz-n-nc.UPF
//            do
//            {
//                if (tmpline.find("angular_momentum") != std::string::npos)
//                {
//                    l[ibeta] = std::stoi(extract201(tmpline, "angular_momentum"));
//                }
//
//                if (tmpline.find("size") != std::string::npos)
//                {
//                    ngrid[ibeta] = std::stoi(extract201(tmpline, "size"));
//                }
//
//                // neither "cutoff_radius_index" nor "cutoff_radius" is reliable!
//                // the code will read all the values first and then reverse scan to determine the grid size
//
//                if (tmpline.find(">") != std::string::npos)
//                {
//                    is_good &= (l[ibeta] >= 0) && (l[ibeta] <= lmax_) && 
//                               (ngrid[ibeta] > 0) && (ngrid[ibeta] <= ngrid_max);
//                    break;
//                }
//            } while (std::getline(ifs, tmpline));
//
//            // read in r*beta(r)
//            for (int ir = 0; ir != ngrid[ibeta]; ++ir)
//            {
//                ifs >> rbeta[ir];
//            }
//
//            // reverse scan to determine the non-zero grid size
//            for (; ngrid[ibeta] > 0; --ngrid[ibeta])
//            {
//                if (std::abs(rbeta[ngrid[ibeta] - 1]) > 1e-12)
//                {
//                    break;
//                }
//            }
//
//            ifs >> tmpword;
//            assert(tmpword.find("</PP_BETA") != std::string::npos);
//            is_good &= ifs.good() && (ngrid[ibeta] > 0);
//
//            printf("l = %i, ngrid = %i, is_good = %d\n", l[ibeta], ngrid[ibeta], is_good);
//        } // for loop over beta functions
//    }// rank == 0
//
//#ifdef __MPI
//    Parallel_Common::bcast_bool(is_good);
//#endif
//    if (!is_good)
//    {
//        ModuleBase::WARNING_QUIT("BetaRadials::read_beta_upf201", "PP_BETA error");
//    }
//
//    /*===========================================================
//     *
//     *          Finalize beta functions (rank-0 only)
//     * 
//     * If lspinorb is set to false but the pseudopotential file
//     * contains spin-orbit information, an averaging over SOC-
//     * related beta functions is performed. Otherwise, the "final"
//     * variables are simply the raw ones.
//     *
//     * Necessary information for averaging:
//     * 1. PP_DIJ block
//     * 2. PP_SPIN_ORB block
//     *
//     *===========================================================*/
//    int nbeta_final = 0;
//    int* l_final = nullptr;
//    int* ngrid_final = nullptr;
//    double* rbeta_final = nullptr;
//    if (rank == 0)
//    {
//        if (!GlobalV::LSPINORB && has_so)
//        {
//            /*
//             * read the PP_DIJ block for averaging
//             *                                                                                 */
//            double* dij = new double[nbeta * nbeta];
//            while (std::getline(ifs, tmpline))
//            {
//                if (tmpline.find("<PP_DIJ") != std::string::npos)
//                {
//                    break;
//                }
//            }
//            assert(!ifs.eof());
//
//            for (int i = 0; i != nbeta * nbeta; ++i)
//            {
//                ifs >> dij[i];
//            }
//            ifs >> tmpword;
//            assert(tmpword == "</PP_DIJ>");
//
//            /*
//             * read the PP_SPIN_ORB block for averaging
//             *                                                                                 */
//            int* lll = new int[nbeta];
//            int* jjj = new int[nbeta];
//            while (std::getline(ifs, tmpline))
//            {
//                if (tmpline.find("<PP_SPIN_ORB") != std::string::npos)
//                {
//                    break;
//                }
//            }
//            assert(!ifs.eof());
//
//            for (int i = 0; i != nbeta; ++i)
//            {
//                std::getline(ifs, tmpline);
//                lll[i] = std::stoi(extract201(tmpline, "lll"));
//                jjj[i] = std::stod(extract201(tmpline, "jjj"));
//            }
//            ifs >> tmpword;
//            assert(tmpword == "</PP_SPIN_ORB>");
//
//            /*
//             * Find out the final number of beta functions
//             *
//             * For each non-zero l, there are two associated beta functions (j = l +/- 0.5)
//             * to be averaged.
//             *                                                                                 */
//            int num_0 = std::count(lll, lll + nbeta, 0);
//            nbeta_final = (nbeta - num_0) / 2 + num_0;
//
//            /*
//             * Average over beta functions
//             *                                                                                 */
//            l_final = new int[nbeta_final];
//            ngrid_final = new int[nbeta_final];
//            rbeta_final = new double[nbeta_final * ngrid_max];
//
//            int ibeta = 0;
//            for (int i = 0; i < nbeta_final; ++i)
//            {
//                int l = lll[ibeta];
//                l_final[i] = l;
//
//                if (l == 0)
//                {
//                    // no averaging for l = 0
//                    ngrid_final[i] = ngrid[ibeta];
//                    std::memcpy(&rbeta_final[i * ngrid_max], &rbeta[ibeta * ngrid_max], ngrid[ibeta] * sizeof(double));
//                    ibeta += 1;
//                }
//                else
//                {
//                    // check that beta functions do come in pairs
//                    assert( lll[ibeta] == lll[ibeta + 1] && 
//                            std::abs(std::abs(jjj[ibeta] - jjj[ibeta+1]) - 1.0) < 1e-6 &&
//                            std::abs(std::abs(jjj[ibeta] - lll[ibeta]) - 0.5) < 1e-6 );
//
//                    // m/p stands for minus/plus 1/2
//                    int im = ibeta;
//                    int ip = ibeta + 1;
//                    if (jjj[ibeta] > lll[ibeta]) std::swap(im, ip);
//
//                    double v = ( (l + 1.0) * dij[ip * nbeta + ip] + l * dij[im * nbeta + im] )
//                                / (2.0 * l + 1.0);
//                    if (std::abs(v) < 1e-8) v = 0.1;
//
//                    // weight of averaging
//                    double wm = 1.0 / (2.0 * l + 1.0) * std::sqrt(std::abs(dij[im * nbeta + im] / v)) * l;
//                    double wp = 1.0 / (2.0 * l + 1.0) * std::sqrt(std::abs(dij[ip * nbeta + ip] / v)) * (l + 1.0);
//                    
//                    // beta(final) = wm * beta(j=l-1/2) + wp * beta(j=l+1/2)
//                    ngrid_final[i] = std::max(ngrid[im], ngrid[ip]);
//                    std::transform(&rbeta[im * ngrid_max], &rbeta[im * ngrid_max] + ngrid_final[i],
//                                   &rbeta[ip * ngrid_max], &rbeta_final[i * ngrid_max],
//                                   [wm, wp](double rbeta_m, double rbeta_p){ return wm * rbeta_m + wp * rbeta_p; });
//
//                    ibeta += 2;
//                }
//            }
//
//            // clean up
//            delete[] dij;
//            delete[] lll;
//            delete[] jjj;
//        }
//        else
//        {
//            nbeta_final = nbeta;
//            l_final = l;
//            ngrid_final = ngrid;
//            rbeta_final = rbeta;
//        }
//    }
//#ifdef __MPI
//    Parallel_Common::bcast_bool(is_good);
//#endif
//    if (!is_good)
//    {
//        ModuleBase::WARNING_QUIT("BetaRadials::read_beta_upf201", "PP_BETA error");
//    }
//
//    /*===========================================================
//     *
//     *          Broadcast final beta functions
//     * 
//     *===========================================================*/
//#ifdef __MPI
//    Parallel_Common::bcast_int(nbeta_final); 
//#endif
//
//    if (rank != 0)
//    {
//        rbeta_final = new double[nbeta_final * ngrid_max];
//        l_final = new int[nbeta_final];
//        ngrid_final = new int[nbeta_final];
//    }
//
//#ifdef __MPI
//    Parallel_Common::bcast_int(l_final, nbeta_final);
//    Parallel_Common::bcast_int(ngrid_final, nbeta_final);
//    Parallel_Common::bcast_double(rbeta_final, nbeta_final * ngrid_max);
//#endif
//
//    /*===========================================================
//     *
//     *              Build BetaRadials object
//     * 
//     *===========================================================*/
//    nchi_ = nbeta_final;
//    chi_ = new NumericalRadial[nchi_];
//    int izeta = 0;
//    for (int i = 0; i != nchi_; ++i)
//    {
//        izeta = (i == 0 || l_final[i] != l_final[i - 1]) ? 0 : izeta + 1;
//        chi_[i].build(l_final[i], true, ngrid_final[i], rgrid, &rbeta_final[i * ngrid_max], 1, izeta, symbol_, itype_);
//    }
//
//    indexing();
//
//    // clean up
//    delete[] l_final;
//    delete[] ngrid_final;
//    delete[] rbeta_final;
//
//    if (rank == 0 && !GlobalV::LSPINORB && has_so)
//    {
//        delete[] l;
//        delete[] ngrid;
//        delete[] rbeta;
//    }
//}
//
//std::string BetaRadials::trim201(std::string const& str)
//{
//    // extract the substring between quotation marks (with whitespace trimmed)
//    // str MUST contain exactly a pair of quotation marks
//
//    std::string::size_type start = str.find('"');
//    std::string::size_type end = str.find_last_of('"');
//    std::string tmp = str.substr(start + 1, end - start - 1);
//
//    if (tmp.length() == 0)
//    {
//        return tmp;
//    }
//
//    start = tmp.find_first_not_of(" \t");
//    end = tmp.find_last_not_of(" \t");
//    return tmp.substr(start, end + 1 - start);
//}
//
//std::string BetaRadials::extract201(std::string const& str, std::string const& keyword) {
//    std::smatch match;
//    std::string regex_string = ".*" + keyword + "=\" *([^=]+) *\".*";
//    std::regex re(regex_string);
//    std::regex_match(str, match, re);
//    return match[1].str();
//}
