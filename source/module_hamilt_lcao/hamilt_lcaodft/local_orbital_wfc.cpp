#include "local_orbital_wfc.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/write_wfc_nao.h"
#include "module_io/read_wfc_nao.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

Local_Orbital_wfc::Local_Orbital_wfc()
{
	allocate_flag = false;
	wfck_flag = false;	
	complex_flag = false;
    nks = 0;
}

Local_Orbital_wfc::~Local_Orbital_wfc()
{

	// used for k-points.
	if(this->complex_flag)
	{
        delete[] this->wfc_k_grid2;
    }
    if(this->wfck_flag)
    {
		for(int i=0; i<nks; i++)
		{
			delete[] this->wfc_k_grid[i];
		}
		delete[] this->wfc_k_grid;
	}
}

void Local_Orbital_wfc::gamma_file(psi::Psi<double>* psid, elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "gamma_file");
    std::cout << " Read in gamma point wave function files " << std::endl;

    double** ctot;

    //allocate psi
    int ncol = this->ParaV->ncol_bands;

    if (GlobalV::KS_SOLVER == "genelpa" 
     || GlobalV::KS_SOLVER == "lapack_gvx" 
     || GlobalV::KS_SOLVER == "scalapack_gvx" 
     || GlobalV::KS_SOLVER == "cg_in_lcao"
#ifdef __CUDA
        || GlobalV::KS_SOLVER == "cusolver"
#endif
        )
    {
        ncol = this->ParaV->ncol;
    }

    if (psid == nullptr)
    {
        ModuleBase::WARNING_QUIT("gamma_file", "psid should be allocated first!");
    }
    else
    {
        psid->resize(GlobalV::NSPIN, ncol, this->ParaV->nrow);
    }
    ModuleBase::GlobalFunc::ZEROS(psid->get_pointer(), psid->size());

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->error = ModuleIO::read_wfc_nao(ctot, is, GlobalV::GAMMA_ONLY_LOCAL, GlobalV::NB2D, GlobalV::NBANDS, 
                                             GlobalV::NLOCAL, GlobalV::global_readin_dir, this->ParaV, psid, pelec);
#ifdef __MPI
        Parallel_Common::bcast_int(this->error);
#endif
        switch (this->error)
        {
        case 1:
            std::cout << "Can't find the wave function file: LOWF_GAMMA_S" << is + 1 << ".txt" << std::endl;
            break;
        case 2:
            std::cout << "In wave function file, band number doesn't match" << std::endl;
            break;
        case 3:
            std::cout << "In wave function file, nlocal doesn't match" << std::endl;
            break;
        case 4:
            std::cout << "In k-dependent wave function file, k point is not correct" << std::endl;
            break;
        default:
            std::cout << " Successfully read in wave functions " << is << std::endl;
        }
        if (this->error)
        {
            std::cout << "WARNING: Failed to read in wavefunction, use default initialization instead." << std::endl;
            break;
        }
    }//loop ispin
}

void Local_Orbital_wfc::allocate_k(const int& lgd,
    psi::Psi<std::complex<double>>* psi,
    elecstate::ElecState* pelec,
    const int& nks,
    const int& nkstot,
    const std::vector<ModuleBase::Vector3<double>>& kvec_c,
    const int& istep)
{
    this->nks = nks;

    ModuleBase::TITLE("Local_Orbital_wfc", "allocate_k");
    if(GlobalV::NLOCAL < GlobalV::NBANDS)
	{
		ModuleBase::WARNING_QUIT("Local_Orbital_wfc::allocate","NLOCAL<NBANDS");
	}

	// mohan add the flag 2011-03-02
	// allocate the first part (only once!).
	if(this->wfck_flag == false)
	{
		this->wfc_k_grid = new std::complex<double>**[nks];
		for(int ik=0; ik<nks; ik++)
		{
			this->wfc_k_grid[ik] = new std::complex<double>*[GlobalV::NBANDS];
		}
		this->wfck_flag = true;
	}
	
	if(this->complex_flag)
	{
		delete[] this->wfc_k_grid2;
		this->complex_flag = false;
	}
	// allocate the second part.
	//if(lgd != 0) xiaohui modify 2015-02-04, fixed memory bug
	if(lgd != 0)
	{
		const int page=GlobalV::NBANDS*lgd;
		this->wfc_k_grid2=new std::complex<double> [nks*page];
		ModuleBase::GlobalFunc::ZEROS(wfc_k_grid2, nks*page);
		for(int ik=0; ik<nks; ik++)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ib++)
			{
				this->wfc_k_grid[ik][ib] = &wfc_k_grid2[ik*page+ib*lgd];
			}
			ModuleBase::Memory::record("LOWF::wfc_k_grid", sizeof(std::complex<double>) * GlobalV::NBANDS*GlobalV::NLOCAL);
			this->complex_flag = true;
		}
	}

    if (INPUT.init_wfc == "atomic")
    {
    }
    else if (INPUT.init_wfc == "file")
    {
		if (istep > 0)
		{
			return;
		}
        std::cout << " Read in wave functions files: " << nkstot << std::endl;
        if (psi == nullptr)
        {
            ModuleBase::WARNING_QUIT("allocate_k", "psi should be allocated first!");
        }
        else
        {
            psi->resize(nkstot, this->ParaV->ncol_bands, this->ParaV->nrow);
        }
        for (int ik = 0; ik < nkstot; ++ik)
        {
            std::complex<double>** ctot;
            this->error = ModuleIO::read_wfc_nao_complex(ctot, ik, GlobalV::NB2D, GlobalV::NBANDS, GlobalV::NLOCAL, 
                                    GlobalV::global_readin_dir, kvec_c[ik], this->ParaV, psi, pelec);
#ifdef __MPI
            Parallel_Common::bcast_int(this->error);
#endif
            switch (this->error)
            {
            case 1:
                std::cout << "Can't find the wave function file: LOWF_K_" << ik + 1 << ".txt" << std::endl;
                break;
            case 2:
                std::cout << "In wave function file, band number doesn't match" << std::endl;
                break;
            case 3:
                std::cout << "In wave function file, nlocal doesn't match" << std::endl;
                break;
            case 4:
                std::cout << "In k-dependent wave function file, k point is not correct" << std::endl;
                break;
            default:
                std::cout << " Successfully read in wave functions " << ik + 1 << std::endl;
            }
            if (this->error)
            {
                std::cout << "WARNING: Failed to read in wavefunction, use default initialization instead." << std::endl;
                break;
            }
        }
    }
    else
    {
		ModuleBase::WARNING_QUIT("Local_Orbital_wfc","check the parameter: init_wfc");
	}

	return;
}


int Local_Orbital_wfc::globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock=localindex/nblk;
    gIndex=(iblock*nprocs+myproc)*nblk+localindex%nblk;
    return gIndex;
}


int Local_Orbital_wfc::localIndex(int globalindex, int nblk, int nprocs, int& myproc)
{
    myproc=int((globalindex%(nblk*nprocs))/nblk);
    return int(globalindex/(nblk*nprocs))*nblk+globalindex%nblk;
}

#ifdef __MPI
void Local_Orbital_wfc::wfc_2d_to_grid(const int istep,
                                       const int out_wfc_lcao,
                                       const double* wfc_2d,
                                       double** wfc_grid,
                                       const ModuleBase::matrix& ekb,
                                       const ModuleBase::matrix& wg)
{
    ModuleBase::TITLE(" Local_Orbital_wfc", "wfc_2d_to_grid");
    ModuleBase::timer::tick("Local_Orbital_wfc","wfc_2d_to_grid");

    const Parallel_Orbitals* pv = this->ParaV;
    const int inc = 1;
    int myid=0;
    MPI_Comm_rank(pv->comm_2D, &myid);
    int info=0;
    
    //calculate maxnloc for bcasting 2d-wfc
    long maxnloc; // maximum number of elements in local matrix
    info=MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
    info=MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
    double *work=new double[maxnloc]; // work/buffer matrix

    double** ctot;
    if (out_wfc_lcao && myid == 0)
    {
        ctot = new double* [GlobalV::NBANDS];
        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            ctot[i] = new double[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
        }
        ModuleBase::Memory::record("LOWF::ctot", sizeof(double) * GlobalV::NBANDS * GlobalV::NLOCAL);
    }

    int naroc[2]; // maximum number of row or column
    for(int iprow=0; iprow<pv->dim0; ++iprow)
    {
        for(int ipcol=0; ipcol<pv->dim1; ++ipcol)
        {
            const int coord[2]={iprow, ipcol};
            int src_rank;
            info=MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
            if(myid==src_rank)
            {
                BlasConnector::copy(pv->nloc_wfc, wfc_2d, inc, work, inc);
                naroc[0]=pv->nrow;
                naroc[1]=pv->ncol_bands;
            }
            info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
            info=MPI_Bcast(work, maxnloc, MPI_DOUBLE, src_rank, pv->comm_2D);

			if (out_wfc_lcao)
			{
				info = this->set_wfc_grid(naroc, pv->nb,
						pv->dim0, pv->dim1, iprow, ipcol,
						work, wfc_grid, myid, ctot);
			}
			else
			{
				info = this->set_wfc_grid(naroc, pv->nb,
						pv->dim0, pv->dim1, iprow, ipcol,
						work, wfc_grid);
			}

        }//loop ipcol
    }//loop iprow
    if(out_wfc_lcao && myid == 0)
    {
        std::stringstream ss;
        if (GlobalV::out_app_flag)
        {
            if (out_wfc_lcao == 1)
            {
                ss << GlobalV::global_out_dir << "LOWF_GAMMA_S" << GlobalV::CURRENT_SPIN + 1 << ".txt";
            }
            else if (out_wfc_lcao == 2)
            {
                ss << GlobalV::global_out_dir << "LOWF_GAMMA_S" << GlobalV::CURRENT_SPIN + 1 << ".dat";
            }
        }
        else
        {
            if (out_wfc_lcao == 1)
            {
                ss << GlobalV::global_out_dir << istep << "_"
                    << "LOWF_GAMMA_S" << GlobalV::CURRENT_SPIN + 1 << ".txt";
            }
            else if (out_wfc_lcao == 2)
            {
                ss << GlobalV::global_out_dir << istep << "_"
                    << "LOWF_GAMMA_S" << GlobalV::CURRENT_SPIN + 1 << ".dat";
            }
        }
        if (out_wfc_lcao == 1)
        {
            ModuleIO::write_wfc_nao(ss.str(), ctot, ekb, wg);
        }
        else if (out_wfc_lcao == 2)
        {
            ModuleIO::write_wfc_nao(ss.str(), ctot, ekb, wg, true);
        }
        for (int i = 0; i < GlobalV::NBANDS; i++)
        {
            delete[] ctot[i];
        }
        delete[] ctot;
    }
    delete[] work;
    ModuleBase::timer::tick("Local_Orbital_wfc","wfc_2d_to_grid");
}

void Local_Orbital_wfc::wfc_2d_to_grid(const int istep,
                                       const int out_wfc_lcao,
                                       const std::complex<double>* wfc_2d,
                                       std::complex<double>** wfc_grid,
                                       int ik,
                                       const ModuleBase::matrix& ekb,
                                       const ModuleBase::matrix& wg,
                                       const std::vector<ModuleBase::Vector3<double>>& kvec_c)
{
    ModuleBase::TITLE("Local_Orbital_wfc", "wfc_2d_to_grid");
    ModuleBase::timer::tick("Local_Orbital_wfc","wfc_2d_to_grid");

    const Parallel_Orbitals* pv = this->ParaV;
    const int inc = 1;
    int myid=0;
    MPI_Comm_rank(pv->comm_2D, &myid);
    int info=0;
    
    //calculate maxnloc for bcasting 2d-wfc
    long maxnloc=0; // maximum number of elements in local matrix
    info=MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
    info=MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
    std::complex<double> *work=new std::complex<double>[maxnloc]; // work/buffer matrix

    std::complex<double> **ctot;
    if (out_wfc_lcao && myid == 0)
    {
        ctot = new std::complex<double>*[GlobalV::NBANDS];
        for (int i=0; i<GlobalV::NBANDS; i++)
        {
            ctot[i] = new std::complex<double>[GlobalV::NLOCAL];
            ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
        }
        ModuleBase::Memory::record("LOWF::ctot", sizeof(std::complex<double>) * GlobalV::NBANDS * GlobalV::NLOCAL);
    }
    
    int naroc[2] = {0}; // maximum number of row or column
    for(int iprow=0; iprow<pv->dim0; ++iprow)
    {
        for(int ipcol=0; ipcol<pv->dim1; ++ipcol)
        {
            const int coord[2]={iprow, ipcol};
            int src_rank;
            info=MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
            if(myid==src_rank)
            {
                BlasConnector::copy(pv->nloc_wfc, wfc_2d, inc, work, inc);
                naroc[0]=pv->nrow;
                naroc[1]=pv->ncol_bands;
            }
            info=MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
            info = MPI_Bcast(work, maxnloc, MPI_DOUBLE_COMPLEX, src_rank, pv->comm_2D);
            
			if (out_wfc_lcao)
			{
				info = this->set_wfc_grid(naroc, pv->nb,
						pv->dim0, pv->dim1, iprow, ipcol,
						work, wfc_grid, myid, ctot);
			}
			else
			{
				// mohan update 2021-02-12, delte BFIELD option
				info = this->set_wfc_grid(naroc, pv->nb,
						pv->dim0, pv->dim1, iprow, ipcol,
						work, wfc_grid);
			}
        }//loop ipcol
    }//loop iprow

    if (out_wfc_lcao && myid == 0)
    {
        std::stringstream ss;
        if (GlobalV::out_app_flag)
        {
            if (out_wfc_lcao == 1)
            {
                ss << GlobalV::global_out_dir << "LOWF_K_" << ik + 1 << ".txt";
            }
            else if (out_wfc_lcao == 2)
            {
                ss << GlobalV::global_out_dir << "LOWF_K_" << ik + 1 << ".dat";
            }
        }
        else
        {
            if (out_wfc_lcao == 1)
            {
                ss << GlobalV::global_out_dir << istep << "_"
                   << "LOWF_K_" << ik + 1 << ".txt";
            }
            else if (out_wfc_lcao == 2)
            {
                ss << GlobalV::global_out_dir << istep << "_"
                   << "LOWF_K_" << ik + 1 << ".dat";
                   std::cout << __LINE__ << " " << ss.str() << std::endl;
            }
        }
        if (out_wfc_lcao == 1)
        {
            ModuleIO::write_wfc_nao_complex(ss.str(), ctot, ik, kvec_c[ik], ekb, wg);
        }
        else if (out_wfc_lcao == 2)
        {
            ModuleIO::write_wfc_nao_complex(ss.str(), ctot, ik, kvec_c[ik], ekb, wg, true);
        }
        for (int i = 0; i < GlobalV::NBANDS; i++)
        {
            delete[] ctot[i];
        }
        delete[] ctot;
    }
    delete[] work;
    ModuleBase::timer::tick("Local_Orbital_wfc","wfc_2d_to_grid");
}
#endif

/**
 * @brief Temporary namespace for Scalapack 2D Block-Cyclic-Distribution (2D-BCD) strategy, storing detached functions from Local_Orbital_wfc and related classes. (The most important) parts of them will be moved to interface-to-ScaLAPACK class in the future, and the rest will be moved elsewhere.
 * 
 */
namespace ScaLAPACK2DBCD{
// standalone python style string split function
/**
 * @brief split a string by a delimiter
 * 
 * @param src string to be split
 * @param delimiter delimiter
 * @param dest vector to store the splitted strings
 */
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
/**
 * @brief check if a string starts with a prefix
 * 
 * @param src string to be checked
 * @param prefix prefix
 * @return true 
 * @return false 
 */
bool startswith(const std::string& src, const std::string& prefix)
{
    return src.find(prefix) == 0;
}
// standalone python style string endswith function
/**
 * @brief check if a string ends with a suffix
 * 
 * @param src string to be checked
 * @param suffix suffix
 * @return true 
 * @return false 
 */
bool endswith(const std::string& src, const std::string& suffix)
{
    return src.rfind(suffix) == src.size() - suffix.size();
}
/**
 * @brief cast a string to a vector of T, the string should be in the format of "real imag real imag ..."
 * 
 * @tparam T 
 * @param src string to be cast
 * @param dest vector to store the casted values
 * @param delimiter delimiter between real and imag parts
 */
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
 * @brief Developer notes on Scalapack 2D Block-Cyclic-Distribution (2D-BCD) strategy
 * One matrix should first be divided into blocks with size MB*NB, then blocks will be distributed to different processors (onto processor grid).
 * A schematic view of 2D-BCD is shown below. Firstly a matrix is like:
 * 
 *    0   1   2   3   4   5   6   7   8  ...
 * 0 a11 a12 a13 a14 a15 a16 a17 a18 a19 ...
 * 1 a21 a22 a23 a24 a25 a26 a27 a28 a29 ...
 * 2 a31 a32 a33 a34 a35 a36 a37 a38 a39 ...
 * 3 a41 a42 a43 a44 a45 a46 a47 a48 a49 ...
 * 4 a51 a52 a53 a54 a55 a56 a57 a58 a59 ...
 * 5 a61 a62 a63 a64 a65 a66 a67 a68 a69 ...
 * 6 a71 a72 a73 a74 a75 a76 a77 a78 a79 ...
 * ...
 * But remember it is stored in a 1D array in memory, so all elements are stored in a continuous manner. This means a11 -> a0, ...:
 *    0   1    2    3    4    5    6    7   8  ...
 * 0 a0  a1   a2   a3   a4   a5   a6   a7  a8  ...
 * 1 an  an+1 ...  ...  ...  ...  ...  ... ...  
 * 
 * if MB=2, NB=3, then the matrix is divided into 2*3 blocks:
 *  |  0   1   2  |  3   4   5  |  6   7   8  | ...
 * -+-------------+-------------+-------------+ ...
 * 0| a11 a12 a13 | a14 a15 a16 | a17 a18 a19 | ...
 * 1| a21 a22 a23 | a24 a25 a26 | a27 a28 a29 | ...
 * -+-------------+-------------+-------------+ ...
 * 2| a31 a32 a33 | a34 a35 a36 | a37 a38 a39 | ...
 * 3| a41 a42 a43 | a44 a45 a46 | a47 a48 a49 | ...
 * -+-------------+-------------+-------------+ ...
 * 4| a51 a52 a53 | a54 a55 a56 | a57 a58 a59 | ...
 * 5| a61 a62 a63 | a64 a65 a66 | a67 a68 a69 | ...
 * -+-------------+-------------+-------------+ ...
 * 6| a71 a72 a73 | a74 a75 a76 | a77 a78 a79 | ...
 * ...
 * 
 * if there are 4 processors, then the blocks will be distributed to the processor grid:
 * 
 *  +-------------+-------------+-------------+ ...
 *  | Processor 0 | Processor 1 | Processor 0 | ...
 *  |   (0, 0)    |   (0, 1)    |   (0, 0)    | ...
 *  +-------------+-------------+-------------+ ...
 *  | Processor 2 | Processor 3 | Processor 2 | ...
 *  |   (1, 0)    |   (1, 1)    |   (1, 0)    | ...
 *  +-------------+-------------+-------------+ ...
 *  | Processor 0 | Processor 1 | Processor 0 | ...
 *  |   (0, 0)    |   (0, 1)    |   (0, 0)    | ...
 *  +-------------+-------------+-------------+ ...
 *  | Processor 2 | Processor 3 | Processor 2 | ...
 * ...
 * 
 * The processor grid is a 2D grid with size Pr*Pc, if there are 20 blocks in one line and there are only 4 processors, then from 0th to 19th block, they will correspond to processor 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3.
 * Additionaly, present Scalapack requires continuous memory allocation of blocks, therefore, the blocks should be allocated in a continuous manner.
 * 
 * For example on Processor 0, the distributed matrix should be like:
 * 
 *  a11 a12 a13 | a17 a18 a19 ...
 *  a21 a22 a23 | a27 a28 a29 ...
 *  ------------+------------ ...
 *  a51 a52 a53 | a57 a58 a59 ...
 *  a61 a62 a63 | a67 a68 a69 ...
 * ...
 * 
 * or say                               it will be locally indiced as:
 *  +-------------+-------------+...    +-------------+-------------+...
 *  |   Block 0   |   Block 2   |...    |   Block 0   |   Block 1   |...
 *  |             |             |...    |             |             |...
 *  +-------------+-------------+...    +-------------+-------------+...
 *  |   Block 2n  |  Block 2n+2 |...    |   Block n   |  Block n+1  |...
 *  |             |             |...    |             |             |...
 *  +-------------+-------------+...    +-------------+-------------+...
 *  ...
 * 
 * Moreover it should be like (1D flatten, row-major, continuous memory allocation, leading dimension specified to let it be 2D matrix):
 * a11 a12 a13 a17 a18 a19 ... a21 a22 a23 a27 a28 a29 ... a51 a52 a53 a57 a58 a59 ...
 * they will be locally indiced as:
 * a1  a2  a3  a4  a5  a6  ...
 * 
 * The following two functions are for converting between global index and local index.
 */

/**
 * @brief Find global index of a matrix element (index before 2DBCD) from local index (index after 2DBCD) and processor index. (2DBCD: 2D Block-Cyclic-Distribution)
 * 
 * @param ixx_local [in] index of matrix element on one processor (across blocks distributed on it)
 * @param iproc [in] index of processor
 * @param ixx_global [out] global index of matrix element
 * @param nxx_block_local [in] number of "something" per block per processor
 * @param nprocs [in] number of processors
 */
void index_before2DBCD(const int& ixx_local,
                       const int& iproc,
                       int& ixx_global,
                       const int& nxx_block_local,
                       const int& nprocs)
{
    // calculate the index of block the element belongs to.
    // note that it is just the index within the domain of processor
    int iblock_local = ixx_local / nxx_block_local;
    // calculate the index within a block
    int ixx_block_local = ixx_local % nxx_block_local;

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

    // it is convenient to assume each processor has same number of blocks
    // but in general it is not necessary
    // however if our assumption is true, the `temp` would be GLOBAL block index
    int temp = iblock_local * nprocs + iproc; // block is the leading dimension

    // ixx_local % nxx_block is the local index with in a block
    // temp * nxx_block is then the shift, then there seems to be a two-dimensional
    // coordinate like (temp, ixx_block_local).
    //                 xx1_inblock, xx2_inblock, xx3_inblock, ...
    // block1 of proc1
    // block1 of proc2
    // block1 of proc3                RETURN                       <- temp
    // ...
    // block2 of proc1
    // block2 of proc2
    // ...
    // block3 of proc1
    // ...
    
    // In summary it is a mapping from (ixx_local, iproc) to a one-dimensional index
    // named ixx_global
    ixx_global = temp * nxx_block_local + ixx_block_local;
}

/**
 * @brief Reverse function of index_before2DBCD, convert global 1d flattened index to local index (ixx_local, iproc) pair
 * 
 * @param ixx_local [out] local index of something on one processor but across all blocks
 * @param iproc [out] processor id
 * @param ixx_global [in] global 1d flattened index
 * @param nxx_block_local [in] number of something per block, or say capacity of something per block
 * @param nprocs [in] number of processors
 */
void index_after2DBCD(int& ixx_local,
                      int& iproc,
                      const int& ixx_global,
                      const int& nxx_block_local,
                      const int& nprocs)
{
    // calculate the total number of elements in blocks with the same index across all processors
    int nxx_block = nxx_block_local * nprocs;
    // therefore, the index of block is (the index of block on each processor)
    int iblock = ixx_global / nxx_block;
    // and index within the concatenation of all blocks is
    int ixx_block = ixx_global % nxx_block;
    // then the processor index is
    iproc = ixx_block / nxx_block_local;

    // so on each processor, there are nxx_block_local elements in each "local block",
    // the index iblock also indicates there are "iblock" of "local blocks" ahead of 
    // the element of interest, therefore, the total number of elements ahead will be
    // the shift on index
    int shift = iblock * nxx_block_local;
    // calculate the index within a "local block"
    int ixx_block_local = ixx_global % nxx_block_local;
    // add together
    ixx_local = shift + ixx_block_local;
}

/**
 * @brief distribute manually the data from a 2D array to a 1D array, the data is distributed according to the 2D Block-Cyclic-Distribution (2D-BCD) strategy
 * 
 * @param iproc present processor index
 * @param nxx_block_local number of element in a block per processor
 * @param dst destination where data copied to
 * @param nrows_dst (ScaLAPACK: MB), number of rows of block per processor
 * @param ncols_dst (ScaLAPACK: NB), number of columns of block per processor
 * @param src source storing data
 * @param iproc_row (ScaLAPACK: MYPROW), index of processor for row in processor grid
 * @param iproc_col (ScaLAPACK: MYPCOL), index of processor for column in processor grid
 * @param nprocs_row (ScaLAPACK: NPROW), number of processors for row in processor grid
 * @param nprocs_col (ScaLAPACK: NPCOL), number of processors for column in processor grid
 * @param nrows_src maximal number of rows of source
 * @param ncols_src maximal number of columns of source
 * 
 * @attention FORTRAN and C++ style data layout is different, the data layout in C++ is row-major, but in FORTRAN it is column-major.
 * 
 * @note find more information on BLACS quick reference: https://netlib.org/blacs/BLACS/QRef.html
 */
void distribute(const int& iproc,             //< [in] (myid)
                const int& nxx_block_local,   //< [in] (nb)
                double* dst,                  //< [out] (work)
                const int& nrows_dst,         //< [in] (naroc[0])
                const int& ncols_dst,         //< [in] (naroc[1])
                double** src,                 //< [in] (CTOT)
                const int& iproc_row,         //< [in] (iprow)
                const int& iproc_col,         //< [in] (ipcol)
                const int& nprocs_row,        //< [in] (dim0)
                const int& nprocs_col,        //< [in] (dim1)
                const int& nrows_src = -1,    //< [in] (nbands)
                const int& ncols_src = -1)    //< [in] new parameter
{
    // FIRST DEFINE THE PROCESSOR INFORMATION
    // workflow consistency judgement: only processor 0 can output the data
    // but all processors have the same workflow
    if(iproc != 0) return;
    for(int i = 0; i < nrows_dst; i++) // loop over the leading dimension
    {
        int irow_global;
        // reversely find what is the element correspond to the element (i, j) on
        // processor (iproc_row, iproc_col)
        // so calculate the global index, or say index of the large matrix to 
        // distribute
        index_before2DBCD(i,                // local index
                          iproc_row,        // processor index for row
                          irow_global,      // global row index
                          nxx_block_local,  // number of elements in a block
                          nprocs_row);      // number of processors for row
        // check if the global index is out of range
        // WHY PRETEND LIKE NOTHING HAPPENED?
        if(irow_global >= nrows_src) continue;
        // if valid...
        for(int j = 0; j < ncols_dst; j++) // loop over the second dimension
        {
            int icol_global;
            index_before2DBCD(j, iproc_col, icol_global, nxx_block_local, nprocs_col);
            // check if the global index is out of range
            // WHY PRETEND LIKE NOTHING HAPPENED?
            if(icol_global >= ncols_src) continue;
            // if not, access the memory. 
            // HERE IS THE PLACE CONTROLLING THE DIFFERENT DATA LAYOUT BETWEEN C++ AND FORTRAN
            // C++ style: row-major, memory of elements in a row is continuous
            // int flatten_index = i * ncols_dst + j;
            // FORTRAN style: column-major, memory of elements in a column is continuous
            int flatten_index = j * nrows_dst + i;
            dst[flatten_index] = src[irow_global][icol_global];
        }
    }
}
/**
 * @brief overloaded version of `distribute` for complex data
 * 
 * @param iproc present processor index
 * @param nxx_block_local number of element in a block per processor
 * @param dst destination where data copied to
 * @param nrows_dst (ScaLAPACK: MB), number of rows of block per processor
 * @param ncols_dst (ScaLAPACK: NB), number of columns of block per processor
 * @param src source storing data
 * @param iproc_row (ScaLAPACK: MYPROW), index of processor for row in processor grid
 * @param iproc_col (ScaLAPACK: MYPCOL), index of processor for column in processor grid
 * @param nprocs_row (ScaLAPACK: NPROW), number of processors for row in processor grid
 * @param nprocs_col (ScaLAPACK: NPCOL), number of processors for column in processor grid
 * @param nrows_src maximal number of rows of source
 * @param ncols_src maximal number of columns of source
 * 
 * @attention FORTRAN and C++ style data layout is different, the data layout in C++ is row-major, but in FORTRAN it is column-major.
 */
void distribute(const int& iproc,             //< [in] (myid)
                const int& nxx_block_local,   //< [in] (nb)
                std::complex<double>* dst,    //< [out] (work)
                const int& nrows_dst,         //< [in] (naroc[0])
                const int& ncols_dst,         //< [in] (naroc[1])
                std::complex<double>** src,   //< [in] (CTOT)
                const int& iproc_row,         //< [in] (iprow)
                const int& iproc_col,         //< [in] (ipcol)
                const int& nprocs_row,        //< [in] (dim0)
                const int& nprocs_col,        //< [in] (dim1)
                const int& nrows_src = -1,    //< [in] (nbands)
                const int& ncols_src = -1)    //< [in] new parameter
{
    if(iproc != 0) return;
    for(int i = 0; i < nrows_dst; i++)
    {
        int irow_global;
        index_before2DBCD(i, iproc_row, irow_global, nxx_block_local, nprocs_row);
        if(irow_global >= nrows_src) continue;
        for(int j = 0; j < ncols_dst; j++)
        {
            int icol_global;
            index_before2DBCD(j, iproc_col, icol_global, nxx_block_local, nprocs_col);
            if(icol_global >= ncols_src) continue;
            int flatten_index = j * nrows_dst + i;
            dst[flatten_index] = src[irow_global][icol_global];
        }
    }
}

/**
 * @brief Reverse function of distribute, collect the data from a 1D array to a 2D array, the data is collected according to the 2D Block-Cyclic-Distribution (2D-BCD) strategy
 * 
 * @param nxx_block_local number of element in a block per processor
 * @param src 
 * @param nrows_src 
 * @param ncols_src 
 * @param dst 
 * @param iproc_row 
 * @param iproc_col 
 * @param nprocs_row 
 * @param nprocs_col 
 * @param iorbs 
 */
void collect(const int& nxx_block_local,
             double* src,
             const int& nrows_src,
             const int& ncols_src,
             double** dst,
             const int& iproc_row,
             const int& iproc_col,
             const int& nprocs_row,
             const int& nprocs_col,
             const int* iorbs = nullptr)
{
    if(dst == nullptr) return;
    for(int i = 0; i < nrows_src; i++) // do you remember the difference between the FORTRAN and C++ matrix layouts?
    {
        int irow_global;
        index_before2DBCD(i, iproc_row, irow_global, nxx_block_local, nprocs_row);
        for(int j = 0; j < ncols_src; j++)
        {
            int icol_global;
            index_before2DBCD(j, iproc_col, icol_global, nxx_block_local, nprocs_col);
            const int icol = iorbs == nullptr ? icol_global : iorbs[icol_global];
            const int flatten_index = j * nrows_src + i;
            dst[irow_global][icol] = src[flatten_index];
        }
    }
}

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
        index_before2DBCD(i, iproc_dim1, ixx_global, nxx_block_local, nprocs_dim1);
        for(int j = 0; j < ndim2_src; j++)
        {
            int jxx_global;
            index_before2DBCD(j, iproc_dim2, jxx_global, nxx_block_local, nprocs_dim2);
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
