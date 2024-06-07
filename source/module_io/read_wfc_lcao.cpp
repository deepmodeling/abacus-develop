/**
 * @file read_wfc_lcao.cpp
 * @author kirk0830
 * @brief this file is part of the project of removal of class Local_Orbital_Wfc (LOWF)
 * @version 0.1
 * @date 2024-06-07
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "module_io/read_wfc_lcao.h"
#include "module_base/formatter.h"
#include "module_base/tool_quit.h"
#include <fstream>
#include <regex>
#include <cassert>

// only when you know why you need the T, you can write function with template parameter, 
// otherwise, you should overload the function for different types
// For example in this case, ONLY lowf support to be std::complex<double> and std::complex<float>,
// not ekb, occ, wk and kvec_c.
/**
 * @brief Read the wavefunction coefficients from the file (for complex wavefunction coefficients)
 * 
 * @tparam T 
 * @param flowf [in] file name like "LOWF_K_*.txt", dumped from ABACUS INPUT out_wfc_lcao 1 
 * @param ik [out] the index of k points
 * @param kvec_c [out] the k vector in Cartesian coordinates
 * @param nbands [out] the number of bands
 * @param nbasis [out] the number of orbitals
 * @param lowf [out] wavefunction coefficients
 * @param ekb [out] eigenvalues
 * @param occ [out] occupations
 * @param wk [out] weight of k points
 */
template<typename T>
void ModuleIO::read_abacus_lowf(const std::string& flowf, 
                                int& ik,
                                ModuleBase::Vector3<double>& kvec_c, 
                                int& nbands, 
                                int& nbasis, 
                                std::vector<std::complex<T>>& lowf, 
                                std::vector<double>& ekb, 
                                std::vector<double>& occ,
                                double& wk)  //<[out] wavefunction coefficients
{
    // assert the T must be double or float
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
            lowf.resize(nbands * nbasis, std::complex<T>(0.0, 0.0)); // initialize lowf
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
                lowf[iband * nbasis + ilocal] = std::complex<T>(std::stod(result[i]), std::stod(result[i + 1]));
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
// instantiate the template function
template void ModuleIO::read_abacus_lowf(const std::string& flowf, int& ik, ModuleBase::Vector3<double>& kvec_c, 
int& nbands, int& nbasis, std::vector<std::complex<double>>& lowf, std::vector<double>& ekb, 
std::vector<double>& occ, double& wk);
template void ModuleIO::read_abacus_lowf(const std::string& flowf, int& ik, ModuleBase::Vector3<double>& kvec_c, 
int& nbands, int& nbasis, std::vector<std::complex<float>>& lowf, std::vector<double>& ekb, 
std::vector<double>& occ, double& wk);
/**
 * @brief Read the wavefunction coefficients from the file (for real wavefunction coefficients)
 * 
 * @tparam T 
 * @param flowf [in] file name like "LOWF_K_*.txt", dumped from ABACUS INPUT out_wfc_lcao 1 
 * @param ik [out] the index of k points
 * @param kvec_c [out] the k vector in Cartesian coordinates
 * @param nbands [out] the number of bands
 * @param nbasis [out] the number of orbitals
 * @param lowf [out] wavefunction coefficients
 * @param ekb [out] eigenvalues
 * @param occ [out] occupations
 * @param wk [out] weight of k points
 */
template<typename T>
void ModuleIO::read_abacus_lowf(const std::string& flowf, 
                                int& ik,
                                ModuleBase::Vector3<double>& kvec_c, 
                                int& nbands, 
                                int& nbasis, 
                                std::vector<T>& lowf, 
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

    ik = 0;
    kvec_c = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    wk = 1.0;
    while (std::getline(ifs, line))
    {
        // remove leading and trailing whitespaces
        line = std::regex_replace(line, std::regex("^ +| +$|( ) +"), "$1");
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
            lowf.resize(nbands * nbasis, 0.0); // initialize lowf
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
            for (const auto& token : result)
            {
                lowf[iband * nbasis + ilocal] = static_cast<T>(std::stod(token));
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
// instantiate the template function
template void ModuleIO::read_abacus_lowf(const std::string& flowf, int& ik, ModuleBase::Vector3<double>& kvec_c,
int& nbands, int& nbasis, std::vector<double>& lowf, std::vector<double>& ekb, std::vector<double>& occ, double& wk);
template void ModuleIO::read_abacus_lowf(const std::string& flowf, int& ik, ModuleBase::Vector3<double>& kvec_c,
int& nbands, int& nbasis, std::vector<float>& lowf, std::vector<double>& ekb, std::vector<double>& occ, double& wk);

#ifdef __MPI
/**
 * @brief Restart the wavefunction coefficients from the file (MPI 2D-BCD version)
 * 
 * @tparam T: datatype of the wavefunction coefficients, can be double, float, std::complex<double> or std::complex<float>
 * @param out_dir [in] the directory where the wavefunction coefficients are stored
 * @param p2d [in] the 2D parallel distribution
 * @param nks [in] the number of k points
 * @param nbands [out] the number of bands
 * @param nbasis [out] the number of orbitals
 * @param lowf_loc [out] the local wavefunction coefficients, can be used to construct psi, see constructor No.8
 * @param ekb [out] the eigenvalues
 * @param occ [out] the occupations
 * @param kvec_c [out] the k vectors in Cartesian coordinates
 * @param wk [out] the weight of k points
 * 
 * @warning Cpxgemr2d not implemented yet
 */
template<typename T>
void ModuleIO::restart_from_file(const std::string& out_dir, // hard-code the file name to be LOWF_K_*.txt?
                                 const Parallel_2D& p2d,
                                 const int& nks,
                                 int& nbands,
                                 int& nbasis,
                                 std::vector<T>& lowf_loc,
                                 std::vector<double>& ekb,
                                 std::vector<double>& occ,
                                 std::vector<ModuleBase::Vector3<double>>& kvec_c,
                                 std::vector<double>& wk)
{
    int nbands_ = -1, nbasis_ = -1;
    for(int ik = 0; ik < nks; ik++)
    {
        // check existence of file
        const std::string flowf = out_dir + "/LOWF_K_" + std::to_string(ik + 1) + ".txt";
        std::ifstream ifs(flowf);
        if(!ifs) ModuleBase::WARNING_QUIT("restart_from_file", "open file failed: " + flowf);

        std::vector<T> lowf_glb;
        std::vector<T> lowf_loc_k;
        std::vector<double> ekb_, occ_;
        ModuleBase::Vector3<double> kvec;
        double wk_;
        read_abacus_lowf(flowf, ik, kvec, nbands, nbasis, lowf_glb, ekb_, occ_, wk_);
        assert(nbands == nbands_ || nbands_ == -1); // check the consistency of nbands
        assert(nbasis == nbasis_ || nbasis_ == -1); // check the consistency of nbasis
        nbands_ = (nbands_ == -1) ? nbands : nbands_;
        nbasis_ = (nbasis_ == -1) ? nbasis : nbasis_;
        // scatter the lowf_glb to lowf_loc
        Parallel_2D p2d_glb;
        p2d_glb.init(nbasis, nbands, std::max(nbasis, nbands), p2d.comm_2D);
        lowf_loc_k.resize(p2d.nrow * p2d.ncol);
        // after PR#4268, Cpxgemr2d can be called to scatter from lowf_glb to lowf_loc_k
        // Cpxgemr2d(nbasis, nbands,
        //           lowf_glb.data(), 1, 1, p2d_glb.desc,
        //           lowf_loc_k.data(), 1, 1, p2d.desc,
        //           &(p2d_glb.blacs_ctxt));
        // now raise not implemented error directly
        ModuleBase::WARNING_QUIT("restart_from_file", "Cpxgemr2d not implemented yet");
        // append to the global lowf_loc
        lowf_loc.insert(lowf_loc.end(), lowf_loc_k.begin(), lowf_loc_k.end());
        ekb.insert(ekb.end(), ekb_.begin(), ekb_.end());
        occ.insert(occ.end(), occ_.begin(), occ_.end());
        wk.push_back(wk_);
        kvec_c.push_back(kvec);
    }
    assert(lowf_loc.size() == nks * p2d.nrow * p2d.ncol);
}
// instantiate the template function
template void ModuleIO::restart_from_file(const std::string& out_dir, const Parallel_2D& p2d, const int& nks, int& nbands, int& nbasis,
std::vector<double>& lowf_loc, std::vector<double>& ekb, std::vector<double>& occ, std::vector<ModuleBase::Vector3<double>>& kvec_c, std::vector<double>& wk);
template void ModuleIO::restart_from_file(const std::string& out_dir, const Parallel_2D& p2d, const int& nks, int& nbands, int& nbasis,
std::vector<float>& lowf_loc, std::vector<double>& ekb, std::vector<double>& occ, std::vector<ModuleBase::Vector3<double>>& kvec_c, std::vector<double>& wk);
template void ModuleIO::restart_from_file(const std::string& out_dir, const Parallel_2D& p2d, const int& nks, int& nbands, int& nbasis,
std::vector<std::complex<double>>& lowf_loc, std::vector<double>& ekb, std::vector<double>& occ, std::vector<ModuleBase::Vector3<double>>& kvec_c, std::vector<double>& wk);
template void ModuleIO::restart_from_file(const std::string& out_dir, const Parallel_2D& p2d, const int& nks, int& nbands, int& nbasis,
std::vector<std::complex<float>>& lowf_loc, std::vector<double>& ekb, std::vector<double>& occ, std::vector<ModuleBase::Vector3<double>>& kvec_c, std::vector<double>& wk);
#endif
/**
 * @brief Restart the wavefunction coefficients from the file (serial version)
 * 
 * @tparam T: datatype of the wavefunction coefficients, can be double, float, std::complex<double> or std::complex<float>
 * @param out_dir [in] the directory where the wavefunction coefficients are stored
 * @param nks [in] the number of k points
 * @param nbands [out] the number of bands
 * @param nbasis [out] the number of orbitals
 * @param lowf_loc [out] the local wavefunction coefficients, can be used to construct psi, see constructor No.8
 * @param ekb [out] the eigenvalues
 * @param occ [out] the occupations
 * @param kvec_c [out] the k vectors in Cartesian coordinates
 * @param wk [out] the weight of k points
 */
template <typename T>
void ModuleIO::restart_from_file(const std::string& out_dir, // hard-code the file name to be LOWF_K_*.txt?
                                 const int& nks,
                                 int& nbands,
                                 int& nbasis,
                                 std::vector<T>& lowf,
                                 std::vector<double>& ekb,
                                 std::vector<double>& occ,
                                 std::vector<ModuleBase::Vector3<double>>& kvec_c,
                                 std::vector<double>& wk)
{
    int nbands_ = -1, nbasis_ = -1;
    for(int ik = 0; ik < nks; ik++)
    {
        // check existence of file
        const std::string flowf = out_dir + "/LOWF_K_" + std::to_string(ik + 1) + ".txt";
        const std::ifstream ifs(flowf);
        if(!ifs) ModuleBase::WARNING_QUIT("restart_from_file", "open file failed: " + flowf);

        std::vector<T> lowf_;
        std::vector<double> ekb_, occ_;
        ModuleBase::Vector3<double> kvec_;
        double wk_;
        int ik_;
        read_abacus_lowf(flowf, ik_, kvec_, nbands, nbasis, lowf_, ekb_, occ_, wk_);
        assert(nbands == nbands_ || nbands_ == -1); // check the consistency of nbands
        assert(nbasis == nbasis_ || nbasis_ == -1); // check the consistency of nbasis
        nbands_ = (nbands_ == -1) ? nbands : nbands_;
        nbasis_ = (nbasis_ == -1) ? nbasis : nbasis_;
        assert(ik_ == ik + 1); // check the consistency of ik
        // append to the global lowf_loc
        lowf.insert(lowf.end(), lowf_.begin(), lowf_.end());
        ekb.insert(ekb.end(), ekb_.begin(), ekb_.end());
        occ.insert(occ.end(), occ_.begin(), occ_.end());
        wk.push_back(wk_);
        kvec_c.push_back(kvec_);
    }
    assert(ekb.size() == nks * nbands);
    assert(occ.size() == nks * nbands);
    assert(lowf.size() == nks * nbands * nbasis);
}
// instantiate the template function
template void ModuleIO::restart_from_file(const std::string& out_dir, const int& nks, int& nbands, int& nbasis,
std::vector<double>& lowf, std::vector<double>& ekb, std::vector<double>& occ, std::vector<ModuleBase::Vector3<double>>& kvec_c, std::vector<double>& wk);
template void ModuleIO::restart_from_file(const std::string& out_dir, const int& nks, int& nbands, int& nbasis,
std::vector<float>& lowf, std::vector<double>& ekb, std::vector<double>& occ, std::vector<ModuleBase::Vector3<double>>& kvec_c, std::vector<double>& wk);
template void ModuleIO::restart_from_file(const std::string& out_dir, const int& nks, int& nbands, int& nbasis,
std::vector<std::complex<double>>& lowf, std::vector<double>& ekb, std::vector<double>& occ, std::vector<ModuleBase::Vector3<double>>& kvec_c, std::vector<double>& wk);
template void ModuleIO::restart_from_file(const std::string& out_dir, const int& nks, int& nbands, int& nbasis,
std::vector<std::complex<float>>& lowf, std::vector<double>& ekb, std::vector<double>& occ, std::vector<ModuleBase::Vector3<double>>& kvec_c, std::vector<double>& wk);
