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

/**
 * @brief Given there is a two-dimensional strategy (iblock, jproc), which means one block can
 * distribute its data to all processors. This function convert local index (ixx_local, iproc) pair 
 * to global 1d flattened index.
 * in which: 
 * ixx_local is local index of something on one processor but across all blocks, 
 * iproc is the index of processor. The need of nxx_block and nprocs will determine the 
 * stepsize of "pointer" to access the memory.
 * THIS FUNCTION PARSES THE LAYOUT OF THE MEMORY ACCORDING TO NPROC AND NXX_BLOCK
 * THIS FUNCTION INDICATES A HARD-CORDED MEMORY LAYOUT/DISTRIBUTION STRATEGY
 * 
 * @param ixx_local [in] local index of something on one processor but across all blocks
 * @param iproc [in] processor id
 * @param ixx_global [out] global 1d flattened index
 * @param nxx_block_local [in] number of something per block on one processor
 * @param nprocs [in] number of processors
 */
void to_ixx_global(const int& ixx_local,
                   const int& iproc,
                   int& ixx_global,
                   const int& nxx_block_local,
                   const int& nprocs)
{
    // ixx_local seems to be index of something but more than the number of blocks,
    // then the index of block it belongs to can be calculated.
    // for example the ixx_local is 10 and nxx_block is 3, then iblock = 10 / 3 = 3
    // which means the ixx_local indiced element on the 4th block
    int iblock = ixx_local / nxx_block_local;          // block index
    int ixx_block_local = ixx_local % nxx_block_local; // index within a block on one processor

    // then there seems to be a (iblock, iproc) 2-dimensional coordinate
    // and flattened to a one-dimensional index
    // means the 2d mapping would be:
    //          proc1,   proc2,   proc3,  proc4, ...
    // block1    0         1        2       3    ...
    // block2  nprocs+0 nprocs+1 nprocs+2 nprocs+3 ...
    // block3 2*nprocs+0 2*nprocs+1 2*nprocs+2 2*nprocs+3 ...
    // block4
    // ...
    // it means each block distributes its data to all processors.
    // (iblock, jproc) -> ij_blockproc
    int temp = iblock * nprocs + iproc; // block is the leading dimension

    // ixx_local % nxx_block is the local index with in a block
    // temp * nxx_block is then the shift, then there seems to be a two-dimensional
    // coordinate like (temp, ixx_block_local).
    //              xx1_inblock, xx2_inblock, xx3_inblock, ...
    // block1-proc1
    // block1-proc2
    // block1-proc3                RETURN                       <- temp
    // ...
    // block2-proc1
    // block2-proc2
    // ...
    // block3-proc1
    // ...
    
    // In summary it is a mapping from (ixx_local, iproc) to a one-dimensional index
    // named ixx_global
    ixx_global = temp * nxx_block_local + ixx_block_local;
}

/**
 * @brief Reverse function of to_ixx_global, convert global 1d flattened index to local index (ixx_local, iproc) pair
 * 
 * @param ixx_local [out] local index of something on one processor but across all blocks
 * @param iproc [out] processor id
 * @param ixx_global [in] global 1d flattened index
 * @param nxx_block_local [in] number of something per block, or say capacity of something per block
 * @param nprocs [in] number of processors
 */
void to_ixx_local(int& ixx_local,
                  int& iproc,
                  const int& ixx_global,
                  const int& nxx_block_local,
                  const int& nprocs)
{
    int nxx_block = nxx_block_local * nprocs; // number of something collected over all processors
    int iblock = ixx_global / nxx_block;      // block index
    int ixx_block = ixx_global % nxx_block;   // index within a block across all processors
    iproc = ixx_block / nxx_block_local;      // processor index

    int shift = iblock * nxx_block_local;     // shift: how many "something" in front of the "thing" of interest
    int ixx_block_local = ixx_global % nxx_block_local; // index within a block on one processor
    ixx_local = shift + ixx_block_local;      // local index
}

/**
 * @brief refactored version of function ModuleIO::CTOT2q. Slicing a part of data 
 * specified by (iproc_dim1, iproc_dim2) tuple from `src` and each element is indiced 
 * by `global indices` (`ndim1_src` and `ndim2_src` are used to control the range)
 * to `dst` on 2-dimensional memory space (whose size is ndim1_dst x ndim2_dst).
 * With this function, any thing can be distributed to all processors, if want the
 * present distribution strategy.
 * 
 * (iproc_dim1, iproc_dim2) identifies one processor, then copy corresponding block
 * of data from src to dst.
 * 
 * @param iproc index of present processor
 * @param nxx_block_local number of "something" per block per processor
 * @param dst destination where data copied to
 * @param ndim1_dst first dimension of destination, always to be nrow
 * @param ndim2_dst second dimension of destination, always to be ncol
 * @param src source storing data
 * @param iproc_dim1 (on MPI 2D Parallel Grid) coordinate i, index of processor for dim1
 * @param iproc_dim2 (on MPI 2D Parallel Grid) coordinate j, index of processor for dim2
 * @param nprocs_dim1 (on MPI 2D Parallel Grid) bound of iproc_dim1
 * @param nprocs_dim2 (on MPI 2D Parallel Grid) bound of iproc_dim2
 * @param ndim1_src first dimension of source, similar with ndim1_dst but will be larger
 * @param ndim2_src second dimension of source, similar with ndim2_dst but will be larger
 */
void distribute_2d(const int& iproc,             //< [in] (myid)
                   const int& nxx_block_local,   //< [in] (nb)
                   double* dst,                  //< [out] (work)
                   const int& ndim1_dst,         //< [in] (naroc[0])
                   const int& ndim2_dst,         //< [in] (naroc[1])
                   double** src,                 //< [in] (CTOT)
                   const int& iproc_dim1,        //< [in] (iprow)
                   const int& iproc_dim2,        //< [in] (ipcol)
                   const int& nprocs_dim1,       //< [in] (dim0)
                   const int& nprocs_dim2,       //< [in] (dim1)
                   const int& ndim1_src = -1,    //< [in] (nbands)
                   const int& ndim2_src = -1)    //< [in] new parameter
{
    // FIRST DEFINE THE PROCESSOR INFORMATION
    // workflow consistency judgement: only processor 0 can output the data
    // but all processors have the same workflow
    if(iproc != 0) return;
    for(int i = 0; i < ndim1_dst; i++)
    {
        int ixx_global;
        to_ixx_global(i, iproc_dim1, ixx_global, nxx_block_local, nprocs_dim1);
        // check if the global index is out of range
        if((ixx_global >= ndim1_src)||(ndim1_src < 0)) continue;
        // if not
        for(int j = 0; j < ndim2_dst; j++)
        {
            int jxx_global;
            to_ixx_global(j, iproc_dim2, jxx_global, nxx_block_local, nprocs_dim2);
            // check if the global index is out of range
            if((jxx_global >= ndim2_src)||(ndim2_src < 0)) continue;
            // if not, access the memory
            int flatten_index = i * ndim2_dst + j;
            dst[flatten_index] = src[ixx_global][jxx_global];
        }
    }
}
/**
 * @brief overloaded version of `distribute_2d` function, which is used to distribute complex data
 * for more information and annotation, see the above function.
 * 
 * @param iproc index of present processor
 * @param nxx_block_local number of "something" per block per processor
 * @param dst destination where data copied to
 * @param ndim1_dst first dimension of destination, always to be nrow
 * @param ndim2_dst second dimension of destination, always to be ncol
 * @param src source storing data
 * @param iproc_dim1 index of processor for dim1
 * @param iproc_dim2 index of processor for dim2
 * @param nprocs_dim1 number of processors for dim1
 * @param nprocs_dim2 number of processors for dim2
 * @param ndim1_src maximal number of the first dimension of source
 * @param ndim2_src maximal number of the second dimension of source
 */
void distribute_2d(const int& iproc,             //< [in] (myid)
                   const int& nxx_block_local,   //< [in] (nb)
                   std::complex<double>* dst,    //< [out] (work)
                   const int& ndim1_dst,         //< [in] (naroc[0])
                   const int& ndim2_dst,         //< [in] (naroc[1])
                   std::complex<double>** src,   //< [in] (CTOT)
                   const int& iproc_dim1,        //< [in] (iprow)
                   const int& iproc_dim2,        //< [in] (ipcol)
                   const int& nprocs_dim1,       //< [in] (dim0)
                   const int& nprocs_dim2,       //< [in] (dim1)
                   const int& ndim1_src = -1,    //< [in] (nbands)
                   const int& ndim2_src = -1)    //< [in] new parameter
{
    if(iproc != 0) return;
    for(int i = 0; i < ndim1_dst; i++)
    {
        int ixx_global;
        to_ixx_global(i, iproc_dim1, ixx_global, nxx_block_local, nprocs_dim1);
        if((ixx_global >= ndim1_src)||(ndim1_src < 0)) continue;
        for(int j = 0; j < ndim2_dst; j++)
        {
            int jxx_global;
            to_ixx_global(j, iproc_dim2, jxx_global, nxx_block_local, nprocs_dim2);
            if((jxx_global >= ndim2_src)||(ndim2_src < 0)) continue;
            int flatten_index = i * ndim2_dst + j;
            dst[flatten_index] = src[ixx_global][jxx_global];
        }
    }
}

/**
 * @brief map data distributed from CTOT onto 2D parallel grid to exact wavefunction container
 * 
 * @param nxx_block_local number of basis function coefficients per block per processor
 * @param orb_indices mapping from global index in CTOT to index of numerical atomic orbital
 * @param src source storing data, which is distributed from CTOT
 * @param ndim1_src maximal number of the first dimension of source
 * @param ndim2_src maximal number of the second dimension of source
 * @param wfc destination storing wavefunction coefficients, is nbands x nlocal
 * @param iproc_dim1 index of processor for dim1
 * @param iproc_dim2 index of processor for dim2
 * @param nprocs_dim1 number of processors for dim1
 * @param nprocs_dim2 number of processors for dim2
 */
void map_onto_wfc(const int& nxx_block_local,
                  const std::vector<int>& orb_indices,
                  double* src,
                  const int& ndim1_src,
                  const int& ndim2_src,
                  double** wfc,
                  const int& iproc_dim1,
                  const int& iproc_dim2,
                  const int& nprocs_dim1 = -1,
                  const int& nprocs_dim2 = -1)
{
    if(wfc == nullptr) return;
    for(int i = 0; i < ndim1_src; i++)
    {
        int ixx_global;
        to_ixx_global(i, iproc_dim1, ixx_global, nxx_block_local, nprocs_dim1);
        for(int j = 0; j < ndim2_src; j++)
        {
            int jxx_global;
            to_ixx_global(j, iproc_dim2, jxx_global, nxx_block_local, nprocs_dim2);
            int iorb = orb_indices[jxx_global];
            #ifdef __DEBUG
            assert (iorb >= 0);
            assert (jxx_global >= 0);
            #endif
            wfc[ixx_global][iorb] = src[i * ndim2_src + j];
        }
    }
}
// overloaded version for complex data, for more information and annotation, see the above function
void map_onto_wfc(const int& nxx_block_local,
                  const std::vector<int>& orb_indices,
                  std::complex<double>* src,
                  const int& ndim1_src,
                  const int& ndim2_src,
                  std::complex<double>** wfc,
                  const int& iproc_dim1,
                  const int& iproc_dim2,
                  const int& nprocs_dim1 = -1,
                  const int& nprocs_dim2 = -1)
{
    if(wfc == nullptr) return;
    for(int i = 0; i < ndim1_src; i++)
    {
        int ixx_global;
        to_ixx_global(i, iproc_dim1, ixx_global, nxx_block_local, nprocs_dim1);
        for(int j = 0; j < ndim2_src; j++)
        {
            int jxx_global;
            to_ixx_global(j, iproc_dim2, jxx_global, nxx_block_local, nprocs_dim2);
            int iorb = orb_indices[jxx_global];
            #ifdef __DEBUG
            assert (iorb >= 0);
            assert (jxx_global >= 0);
            #endif
            wfc[ixx_global][iorb] = src[i * ndim2_src + j];
        }
    }
}

} // end of namespace LOWF
