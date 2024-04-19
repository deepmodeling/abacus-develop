#ifndef LOCAL_ORBITAL_WFC
#define LOCAL_ORBITAL_WFC

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/ORB_control.h" // mohan add 2021-05-24
#include "module_elecstate/elecstate.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_psi/psi.h"

class Local_Orbital_wfc
{
public:

	Local_Orbital_wfc();    
	~Local_Orbital_wfc();
    // refactor new implementation: RAII
    // a more look-looking name would be LocalOrbitalWfc, I suppose...
    Local_Orbital_wfc(const int& nspin,
                      const int& nks,
                      const int& nbands,
                      const int& nlocal,
                      const int& gamma_only,
                      const int& nb2d,
                      const std::string& ks_solver,
                      const std::string& readin_dir);
    //
    void initialize();
    ///=========================================
    /// grid wfc
    /// used to generate density matrix: LOC.DM_R,
    /// which is used to calculate the charge density.
    /// which is got after the diagonalization of
    /// std::complex Hamiltonian matrix.
    ///=========================================
    std::complex<double>*** wfc_k_grid; // [NK, GlobalV::NBANDS, GlobalV::NLOCAL]
    std::complex<double>* wfc_k_grid2;  // [NK*GlobalV::NBANDS*GlobalV::NLOCAL]

    // pointer to const Parallel_Orbitals object
    const Parallel_Orbitals* ParaV;
    // pointer to const Grid_Technique object, although have no idea about what it is...
    // the name of Grid_Technique should be changed to be more informative
    const Grid_Technique* gridt;

    /// read wavefunction coefficients: LOWF_*.txt
    void gamma_file(psi::Psi<double>* psid, elecstate::ElecState* pelec);
    void allocate_k(const int& lgd,
        psi::Psi<std::complex<double>>* psi,
        elecstate::ElecState* pelec,
        const int& nks,
        const int& nkstot,
        const std::vector<ModuleBase::Vector3<double>>& kvec_c,
        const int& istep);

    //=========================================
    // Init Cij, make it satisfy 2 conditions:
    // (1) Unit
    // (2) Orthogonal <i|S|j>= \delta{ij}
    //=========================================
    // void init_Cij(const bool change_c = 1);

    ///=========================================
    /// Parallel: map of index in 2D distribution: global<->local
    ///=========================================
    static int globalIndex(int localindex, int nblk, int nprocs, int myproc);
    static int localIndex(int globalindex, int nblk, int nprocs, int& myproc);

#ifdef __MPI
    ///=========================================
    /// Parallel: convert the distribution of wavefunction from 2D to grid
    ///=========================================
    /// For gamma_only, T = double;
    /// For multi-k, T = complex<double>;
    /// Set myid and ctot when output is needed;
    /// Set wfc as nullptr when 2d-to-grid convertion is not needed.

    // Notice: here I reload this function rather than use template
    // (in which the implementation should be put in header file )
    // because sub-function `write_wfc_nao_complex`contains GlobalC declared in `global.h`
    // which will cause lots of "not defined" if included in a header file.
    void wfc_2d_to_grid(const int istep,
                        const int out_wfc_lcao,
                        const double* wfc_2d,
                        double** wfc_grid,
                        const ModuleBase::matrix& ekb,
                        const ModuleBase::matrix& wg);
    void wfc_2d_to_grid(const int istep,
                        const int out_wfc_lcao,
                        const std::complex<double>* wfc_2d,
                        std::complex<double>** wfc_grid,
                        int ik,
                        const ModuleBase::matrix& ekb,
                        const ModuleBase::matrix& wg,
                        const std::vector<ModuleBase::Vector3<double>>& kvec_c);
#endif

    int error = 0;

  private:

    template <typename T>
    int set_wfc_grid(int naroc[2],
                     int nb,
                     int dim0,
                     int dim1,
                     int iprow,
                     int ipcol,
                     T* work,
                     T** wfc,
                     int myid = -1,
                     T** ctot = nullptr);

    bool wfck_flag;
    bool complex_flag;
    bool allocate_flag;
    int nks;
};


// the function should not be defined here!! mohan 2024-03-28
template <typename T>
int Local_Orbital_wfc::set_wfc_grid(int naroc[2],
                                    int nb,
                                    int dim0,
                                    int dim1,
                                    int iprow,
                                    int ipcol,
                                    T* work,
                                    T** wfc,
                                    int myid,
                                    T** ctot)
{
#ifdef __MPI
    ModuleBase::TITLE(" Local_Orbital_wfc", "set_wfc_grid");
    if (!wfc && !ctot)
    {
        return 0;
    }
    for (int j = 0; j < naroc[1]; ++j)
    {
        int igcol = globalIndex(j, nb, dim1, ipcol);
		if (igcol >= GlobalV::NBANDS)
		{
			continue;
		}
		for (int i = 0; i < naroc[0]; ++i)
		{
			int igrow = globalIndex(i, nb, dim0, iprow);
			int mu_local = this->gridt->trace_lo[igrow];
			if (wfc && mu_local >= 0)
			{
				wfc[igcol][mu_local] = work[j * naroc[0] + i];
			}
			if (ctot && myid == 0)
			{
				ctot[igcol][igrow] = work[j * naroc[0] + i];
			}
        }
    }
#endif
    return 0;
}
#endif

namespace LOWF{
// standalone python style string split function
void split(const std::string& src, const char& delimiter, std::vector<std::string>& dest)
{
    std::string str = src;
    std::string substring;
    std::string::size_type start = 0, index;
    do
    {
        index = str.find_first_of(delimiter, start);
        if (index != std::string::npos)
        {
            substring = str.substr(start, index - start);
            dest.push_back(substring);
            start = str.find_first_not_of(delimiter, index);
            if (start == std::string::npos) return;
        }
    } while (index != std::string::npos);
    // the last token
    substring = str.substr(start);
    dest.push_back(substring);
}
// standalone python style string startswith function
bool startswith(const std::string& src, const std::string& prefix)
{
    return src.find(prefix) == 0;
}
// standalone python style string endswith function
bool endswith(const std::string& src, const std::string& suffix)
{
    return src.rfind(suffix) == src.size() - suffix.size();
}

template<typename T>
void cast_to_stdcomplex(const std::string& src, std::vector<T>& dest, const char& delimiter = ' ')
{
    // there may be parentheses in the string
    // so we need to remove them first
    std::string str = src;
    str.erase(std::remove(str.begin(), str.end(), '('), str.end());
    str.erase(std::remove(str.begin(), str.end(), ')'), str.end());
    std::vector<std::string> tokens;
    split(str, delimiter, tokens);
    for (const auto& token : tokens)
    {
        dest.push_back(std::stod(token));
    }
}

#include <regex>
#include <cassert>
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
        if(endswith(line, "(index of k points)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            ik = std::stoi(result[0]);
            read_kvec = true;
            continue;
        }
        if(read_kvec)
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            for (const auto& token : result)
            {
                kvec_c.push_back(std::stod(token));
            }
            read_kvec = false;
            continue;
        }
        if(endswith(line, "(number of bands)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            nbands = std::stoi(result[0]);
        }
        else if(endswith(line, "(number of orbitals)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            nlocal = std::stoi(result[0]);
        }
        else if(endswith(line, "(band)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            #ifdef __DEBUG
            assert (ilocal == 0)||(ilocal == nlocal);
            #endif
            iband = std::stoi(result[0]);
            ilocal = 0; // reset ilocal
        }
        else if(endswith(line, "(Ry)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            energies.push_back(std::stod(result[0]));
        }
        else if(endswith(line, "(Occupations)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            occs.push_back(std::stod(result[0]));
            #ifdef __DEBUG
            assert (ilocal == 0);
            #endif
        }
        else// read wavefunction coefficients
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            for (const auto& token : result)
            {
                std::vector<double> temp;
                cast_to_stdcomplex(token, temp, ' ');
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
        if(endswith(line, "(index of k points)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            ik = std::stoi(result[0]);
            read_kvec = true;
            continue;
        }
        if(read_kvec)
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            for (const auto& token : result)
            {
                kvec_c.push_back(std::stod(token));
            }
            read_kvec = false;
            continue;
        }
        if(endswith(line, "(number of bands)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            nbands = std::stoi(result[0]);
        }
        else if(endswith(line, "(number of orbitals)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            nlocal = std::stoi(result[0]);
        }
        else if(endswith(line, "(band)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            #ifdef __DEBUG
            assert (ilocal == 0)||(ilocal == nlocal);
            #endif
            iband = std::stoi(result[0]);
            ilocal = 0; // reset ilocal
        }
        else if(endswith(line, "(Ry)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            energies.push_back(std::stod(result[0]));
        }
        else if(endswith(line, "(Occupations)"))
        {
            std::vector<std::string> result;
            split(line, ' ', result);
            occs.push_back(std::stod(result[0]));
            #ifdef __DEBUG
            assert (ilocal == 0);
            #endif
        }
        else// read wavefunction coefficients
        {
            std::vector<std::string> result;
            split(line, ' ', result);
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

// function convert local index to global index
int to_ixx_global(const int& ixx_local, //<[in] local index of something
                  const int& nxx_block, //<[in] number of something per block
                  const int& iproc,     //<[in] processor id
                  const int& nprocs)    //<[in] number of processors
{
    int iblock = ixx_local / nxx_block;
    // locali seems to be index of something but more than the number of blocks,
    // then the index of block it belongs to can be calculated.

    // then there seems to be a (iblock, iproc) 2-dimensional coordinate
    // and flattened to a one-dimensional index
    int temp = iblock * nprocs + iproc;
    // means the 2d mapping would be:
    //         proc1, proc2, proc3, proc4, ...
    // block1
    // block2
    // block3
    // block4
    // ...

    // then return the index of a two-dimensional mapping like:
    // irow = (iblock*nprocs + iproc)
    // icol = locali % nblock
    // (irow, icol) in a two-dimensional mapping
    //         block1, block2, block3, block4, ...
    // ?1
    // ?2                    (irow, icol)
    // ...

}
} // end of namespace LOWF
