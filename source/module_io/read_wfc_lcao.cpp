#include "module_io/read_wfc_lcao.h"
#include "module_base/formatter.h"
#include "module_base/tool_quit.h"
#include <fstream>
#include <regex>
#include <cassert>

template<typename T>
void cast_to_stdcomplex(const std::string& src, std::vector<T>& dest, const std::string& delimiter = " ")
{
    // there may be parentheses in the string
    // so we need to remove them first
    std::string str = src;
    str.erase(std::remove(str.begin(), str.end(), '('), str.end());
    str.erase(std::remove(str.begin(), str.end(), ')'), str.end());
    const std::vector<std::string> tokens= FmtCore::split(str, delimiter);
    for (const auto& token : tokens)
    {
        dest.push_back(std::stod(token));
    }
}

void ModuleIO::read_abacus_lowf(const std::string& flowf, 
                                int& ik,
                                ModuleBase::Vector3<double>& kvec_c, 
                                int& nbands, 
                                int& nbasis, 
                                std::vector<std::complex<double>>& lowf, 
                                std::vector<double>& ekb, 
                                std::vector<double>& occ,
                                double& wk)  //<[out] wavefunction coefficients
{
    std::ifstream ifs(flowf.c_str());
    if(!ifs) ModuleBase::WARNING_QUIT("read_abacus_lowf", "open file failed: " + flowf);
    // will use line-by-line parse
    std::string line;
    bool read_kvec = false;
    int iband = 0;
    int ilocal = 0;
    while (std::getline(ifs, line))
    {
        // remove leading and trailing whitespaces
        line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");
        if(FmtCore::endswith(line, "(index of k points)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            ik = std::stoi(result[0]);
            read_kvec = true;
            continue;
        }
        if(read_kvec)
        {
            const std::vector<std::string> result = FmtCore::split(line);
            kvec_c.x = std::stod(result[0]);
            kvec_c.y = std::stod(result[1]);
            kvec_c.z = std::stod(result[2]);
            read_kvec = false;
            continue;
        }
        if(FmtCore::endswith(line, "(number of bands)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            nbands = std::stoi(result[0]);
            ekb.resize(nbands, 0.0); // initialize ekb
            occ.resize(nbands, 0.0); // initialize occ
        }
        else if(FmtCore::endswith(line, "(number of orbitals)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            nbasis = std::stoi(result[0]);
            lowf.resize(nbands * nbasis, std::complex<double>(0.0, 0.0)); // initialize lowf
        }
        else if(FmtCore::endswith(line, "(band)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            #ifdef __DEBUG
            assert (ilocal == 0)||(ilocal == nlocal);
            #endif
            iband = std::stoi(result[0]) - 1;
            ilocal = 0; // reset ilocal
        }
        else if(FmtCore::endswith(line, "(Ry)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            ekb[iband] = std::stod(result[0]);
        }
        else if(FmtCore::endswith(line, "(Occupations)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            occ[iband] = std::stod(result[0]);
            #ifdef __DEBUG
            assert (ilocal == 0);
            #endif
        }
        else// read wavefunction coefficients
        {
            const std::vector<std::string> result = FmtCore::split(line);
            // for the case the complex number is written as a b
            for(int i = 0; i < result.size(); i += 2)
            {
                lowf[iband * nbasis + ilocal] = std::complex<double>(std::stod(result[i]), std::stod(result[i + 1]));
                ilocal += 1;
            }
            // // for the case the complex number is written as (a, b)
            // for (const auto& token : result)
            // {
            //     std::vector<double> temp;
            //     cast_to_stdcomplex(token, temp);
            //     for (int i = 0; i < temp.size(); i += 2)
            //     {
            //         lowf.push_back(std::complex<double>(temp[i], temp[i + 1]));
            //     }
            //     ilocal += temp.size() / 2;
            // }
        }
    }
    #ifdef __DEBUG
    assert (lowf.size() == nbands * nlocal);
    assert (iband == nbands);
    assert (ilocal == nlocal);
    #endif
}

void ModuleIO::read_abacus_lowf(const std::string& flowf, 
                                int& ik,
                                ModuleBase::Vector3<double>& kvec_c, 
                                int& nbands, 
                                int& nbasis, 
                                std::vector<double>& lowf, 
                                std::vector<double>& ekb, 
                                std::vector<double>& occ,
                                double& wk)
{
    std::ifstream ifs(flowf.c_str());
    if(!ifs) ModuleBase::WARNING_QUIT("read_abacus_lowf", "open file failed: " + flowf);
    // will use line-by-line parse
    std::string line;
    bool read_kvec = false;
    int iband = 0;
    int ilocal = 0;
    while (std::getline(ifs, line))
    {
        // remove leading and trailing whitespaces
        line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");
        if(FmtCore::endswith(line, "(index of k points)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            ik = std::stoi(result[0]);
            read_kvec = true;
            continue;
        }
        if(read_kvec)
        {
            const std::vector<std::string> result = FmtCore::split(line);
            kvec_c.x = std::stod(result[0]);
            kvec_c.y = std::stod(result[1]);
            kvec_c.z = std::stod(result[2]);
            read_kvec = false;
            continue;
        }
        if(FmtCore::endswith(line, "(number of bands)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            nbands = std::stoi(result[0]);
        }
        else if(FmtCore::endswith(line, "(number of orbitals)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            nbasis = std::stoi(result[0]);
        }
        else if(FmtCore::endswith(line, "(band)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            #ifdef __DEBUG
            assert (ilocal == 0)||(ilocal == nlocal);
            #endif
            iband = std::stoi(result[0]);
            ilocal = 0; // reset ilocal
        }
        else if(FmtCore::endswith(line, "(Ry)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            ekb.push_back(std::stod(result[0]));
        }
        else if(FmtCore::endswith(line, "(Occupations)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            occ.push_back(std::stod(result[0]));
            #ifdef __DEBUG
            assert (ilocal == 0);
            #endif
        }
        else// read wavefunction coefficients
        {
            const std::vector<std::string> result = FmtCore::split(line);
            for (const auto& token : result)
            {
                lowf.push_back(std::stod(token));
                ilocal += 1;
            }
        }
    }
    #ifdef __DEBUG
    assert (lowf.size() == nbands * nlocal);
    assert (iband == nbands);
    assert (ilocal == nlocal);
    #endif
}

#ifdef __MPI
void ModuleIO::scatter_lowf(const int& nbands, 
                            const int& nbasis, 
                            std::vector<std::complex<double>>& lowf_glb,
                            const MPI_Comm& comm,  // MPI comm world for 2d bcd
                            const int desc[9],     // the descriptor of the matrix
                            const int& blacs_ctxt, // the one of orb_con.ParaV
                            std::vector<std::complex<double>>& lowf_loc)
{
    // assert(lowf_glb.size() == nbands * nbasis); // this is only valid for source rank
    // while the para2d_loc has been determined externally, it stores in Parallel_Orbitals, but
    // that class is not well-designed, therefore, we directly get information of interest in
    // param list of this function: COMM_WORLD, nb, desc, blacs_ctxt
    Parallel_2D para2d_loc, para2d_glb;

    // get nbands_loc and nbasis_loc from desc. They should be both controlled by a param nb2d
    const int nb = desc[4]; // the block size
    const int mb = desc[5]; // the block size
    para2d_loc.set(nbasis, nbands, std::max(nb, mb), comm, blacs_ctxt);
    // set means: "para2d_loc is the thing that indicates a 2dbcd from nbands*nbasis to
    // at maximum the data of size nb*nb, in comm world comm and blacs_ctxt blacs_ctxt"
    para2d_glb.init(nbasis, nbands, std::max(nbands, nbasis), comm);

    // the other things is this is just scatter of one single k-point. Loop over all kpoints
    // to get the full lowf to construct psi

    lowf_loc.resize(para2d_loc.nrow * para2d_loc.ncol, std::complex<double>(0.0, 0.0));
    const int one = 1;
    pzgemr2d_(&nbasis, &nbands, lowf_glb.data(), &one, &one, para2d_glb.desc, 
                                lowf_loc.data(), &one, &one, desc, 
              &(para2d_glb.blacs_ctxt));
    // alternative would be Cpxgemr2d, which means, C-style, complex, general matrix, rowwise distribution
    // Cpzgemr2d(nbasis, nbands, lowf_glb.data(), 1, 1, para2d_glb.desc, 
    //                           lowf_loc.data(), 1, 1, desc, para2d_glb.blacs_ctxt);
    // RESULTING PROCESSOR GRID IS COLUNM-MAJOR, why???
}
#endif

#ifdef __MPI
void ModuleIO::restart_from_file(const std::string& out_dir, // hard-code the file name to be LOWF_K_*.txt?
                                 const MPI_Comm& comm,  // MPI comm world for 2d bcd
                                 const int desc[9],
                                 const int& blacs_ctxt, // the one of orb_con.ParaV
                                 const int& nks,
                                 int& nbands,
                                 int& nbasis,
                                 std::vector<std::complex<double>>& lowf_loc,
                                 std::vector<double>& ekb,
                                 std::vector<double>& occ,
                                 std::vector<ModuleBase::Vector3<double>>& kvec_c,
                                 std::vector<double>& wk)
{
    int nbands_ = -1, nbasis_ = -1;
    // get nbands_loc and nbasis_loc from desc. They should be both controlled by a param nb2d
    const int nb = desc[4]; // the block size
    const int mb = desc[5]; // the block size
    for(int ik = 0; ik < nks; ik++)
    {
        // check existence of file
        const std::string flowf = "LOWF_K_" + std::to_string(ik) + ".txt";
        std::ifstream ifs(out_dir + "/" + flowf);
        if(!ifs) ModuleBase::WARNING_QUIT("restart_from_file", "open file failed: " + flowf);
        
        int nbands, nbasis;
        std::vector<std::complex<double>> lowf_glb;
        std::vector<std::complex<double>> loc_;
        std::vector<double> ekb_, occ_;
        ModuleBase::Vector3<double> kvec;
        double wk_;
        read_abacus_lowf(flowf, ik, kvec, nbands, nbasis, lowf_glb, ekb_, occ_, wk_);
        assert(nbands == nbands_ || nbands_ == -1); // check the consistency of nbands
        assert(nbasis == nbasis_ || nbasis_ == -1); // check the consistency of nbasis
        nbands_ = nbands_ == -1 ? nbands : nbands_;
        nbasis_ = nbasis_ == -1 ? nbasis : nbasis_;
        // scatter the lowf_glb to lowf_loc
        scatter_lowf(nbands, nbasis, lowf_glb, comm, desc, blacs_ctxt, loc_);
        // append to the global lowf_loc
        lowf_loc.insert(lowf_loc.end(), loc_.begin(), loc_.end());
        ekb.insert(ekb.end(), ekb_.begin(), ekb_.end());
        occ.insert(occ.end(), occ_.begin(), occ_.end());
        wk.push_back(wk_);
        kvec_c.push_back(kvec);
    }
    assert(lowf_loc.size() == nks * nb * nb);
}
#endif