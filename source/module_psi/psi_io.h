/**
 * @file psi_io.cpp
 * @author kirk0830
 * @brief This file collects both I/O of wavefunction for planewave and numerical atomic orbitals basis sets
 * @version 0.1
 * @date 2024-05-30
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include <regex>
#include <cassert>
#include "module_base/formatter.h"
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "module_base/tool_quit.h"

namespace psi
{
namespace io
{
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

// read file named as LOWF_K_*.txt, the name should not be hard-coded
// structure of LOWF_K_*.txt:
// [ik] (index of k points)
// [xk] [yk] [zk]
// [nb] (number of bands)
// [nlocal] (number of orbitals)
// [ib] (band)
// [energy] (energy of the band)
// [occ] (Occupations)
// [real] [imag] [real] [imag] ... (wavefunction coefficients)
// ...
// [ib] (band)
// [energy] (energy of the band)
// [occ] (Occupations)
// [real] [imag] [real] [imag] ... (wavefunction coefficients)
// ...
/**
 * @brief read ABACUS output wavefunction file named as LOWF_K_*.txt for std::complex<double> wavefunction coefficients
 * 
 * @param flowf file name under convention of LOWF_K_*.txt
 * @param ik index of k points, will be returned
 * @param kvec_c k vector in Cartesian coordinates, will be returned
 * @param nbands number of bands, will be returned
 * @param nlocal number of orbitals, will be returned
 * @param energies energies of bands, will be returned
 * @param occs occupations, will be returned
 * @param lowf wavefunction coefficients, will be returned. Note! 1D array of complex numbers
 */
void read_abacus_lowf(const std::string& flowf,                 //<[in]  LOWF_K_*.txt
                      int& ik,                                  //<[out] index of k points
                      std::vector<int>& kvec_c,                 //<[out] k vector
                      int& nbands,                              //<[out] number of bands
                      int& nlocal,                              //<[out] number of orbitals
                      std::vector<double>& energies,            //<[out] energies of bands
                      std::vector<double>& occs,                //<[out] occupations
                      std::vector<std::complex<double>>& lowf)  //<[out] wavefunction coefficients
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
            for (const auto& token : result)
            {
                kvec_c.push_back(std::stod(token));
            }
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
            nlocal = std::stoi(result[0]);
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
            energies.push_back(std::stod(result[0]));
        }
        else if(FmtCore::endswith(line, "(Occupations)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            occs.push_back(std::stod(result[0]));
            #ifdef __DEBUG
            assert (ilocal == 0);
            #endif
        }
        else// read wavefunction coefficients
        {
            const std::vector<std::string> result = FmtCore::split(line);
            for (const auto& token : result)
            {
                std::vector<double> temp;
                cast_to_stdcomplex(token, temp);
                for (int i = 0; i < temp.size(); i += 2)
                {
                    lowf.push_back(std::complex<double>(temp[i], temp[i + 1]));
                }
                ilocal += temp.size() / 2;
            }
        }
    }
    #ifdef __DEBUG
    assert (lowf.size() == nbands * nlocal);
    assert (iband == nbands);
    assert (ilocal == nlocal);
    #endif
}
// overload for the case of double wavefunction coefficients
/**
 * @brief double overloaded version of read_abacus_lowf, read wavefunction coefficients as double
 * 
 * @param flowf file name under convention of LOWF_K_*.txt
 * @param ik index of k points, will be returned
 * @param kvec_c k vector in Cartesian coordinates, will be returned
 * @param nbands number of bands, will be returned
 * @param nlocal number of orbitals, will be returned
 * @param energies energies of bands, will be returned
 * @param occs occupations, will be returned
 * @param lowf wavefunction coefficients, will be returned. Note! 1D array of complex numbers 
 */
void read_abacus_lowf(const std::string& flowf,                 //<[in]  LOWF_K_*.txt
                      int& ik,                                  //<[out] index of k points
                      std::vector<int>& kvec_c,                 //<[out] k vector
                      int& nbands,                              //<[out] number of bands
                      int& nlocal,                              //<[out] number of orbitals
                      std::vector<double>& energies,            //<[out] energies of bands
                      std::vector<double>& occs,                //<[out] occupations
                      std::vector<double>& lowf)                //<[out] wavefunction coefficients
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
            for (const auto& token : result)
            {
                kvec_c.push_back(std::stod(token));
            }
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
            nlocal = std::stoi(result[0]);
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
            energies.push_back(std::stod(result[0]));
        }
        else if(FmtCore::endswith(line, "(Occupations)"))
        {
            std::vector<std::string> result = FmtCore::split(line);
            occs.push_back(std::stod(result[0]));
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

} // namespace io
} // namespace psi