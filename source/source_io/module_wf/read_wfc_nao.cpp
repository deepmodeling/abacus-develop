#include "read_wfc_nao.h"

#include "source_base/parallel_common.h"
#include "source_base/timer.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/module_out/filename.h"
#include "source_base/tool_title.h" // use title
#include "source_base/global_function.h" // use READ_VALUE

#include <type_traits>

namespace
{

template <typename T>
bool read_binary_value(std::ifstream& ifs, T& value)
{
    ifs.read(reinterpret_cast<char*>(&value), sizeof(T));
    return static_cast<bool>(ifs);
}

template <typename T>
bool read_record_value(std::ifstream& ifs, T& value, const bool binary)
{
    if (binary)
    {
        return read_binary_value(ifs, value);
    }
    ModuleBase::GlobalFunc::READ_VALUE(ifs, value);
    return static_cast<bool>(ifs);
}

bool read_binary_wfc_data(std::ifstream& ifs, double& data)
{
    return read_binary_value(ifs, data);
}

bool read_binary_wfc_data(std::ifstream& ifs, float& data)
{
    double value = 0.0;
    if (!read_binary_value(ifs, value))
    {
        return false;
    }
    data = static_cast<float>(value);
    return true;
}

bool read_binary_wfc_data(std::ifstream& ifs, std::complex<double>& data)
{
    double real = 0.0;
    double imag = 0.0;
    if (!read_binary_value(ifs, real) || !read_binary_value(ifs, imag))
    {
        return false;
    }
    data = std::complex<double>(real, imag);
    return true;
}

bool read_binary_wfc_data(std::ifstream& ifs, std::complex<float>& data)
{
    double real = 0.0;
    double imag = 0.0;
    if (!read_binary_value(ifs, real) || !read_binary_value(ifs, imag))
    {
        return false;
    }
    data = std::complex<float>(static_cast<float>(real), static_cast<float>(imag));
    return true;
}

} // namespace

// mohan add 2025-10-19
void ModuleIO::read_wfc_nao_one_data(std::ifstream& ifs, float& data)
{
    ifs >> data;
}

void ModuleIO::read_wfc_nao_one_data(std::ifstream& ifs, double& data)
{
    ifs >> data;
}

void ModuleIO::read_wfc_nao_one_data(std::ifstream& ifs, std::complex<double>& data)
{
    double a = 0.0;
    double b = 0.0;
    ifs >> a >> b;
    data = std::complex<double>(a, b);
}

void ModuleIO::read_wfc_nao_one_data(std::ifstream& ifs, std::complex<float>& data)
{
    float a = 0.0;
    float b = 0.0;
    ifs >> a >> b;
    data = std::complex<float>(a, b);
}

template <typename T>
bool ModuleIO::read_wfc_nao(
    const std::string& global_readin_dir,
    const Parallel_Orbitals& ParaV,
    psi::Psi<T>& psid,
	ModuleBase::matrix& ekb,
    ModuleBase::matrix& wg,
	const std::vector<int> &ik2iktot,
	const int nkstot,
	const int nspin,
    const int read_type,
    const int skip_band,
    const int istep)
{
    ModuleBase::TITLE("ModuleIO", "read_wfc_nao");
    ModuleBase::timer::start("ModuleIO", "read_wfc_nao");

    const int nk = ekb.nr;

    const bool gamma_only = std::is_same<T, double>::value || std::is_same<T, float>::value;
    const bool binary = (read_type == 2);
    bool read_success = true;
    int myrank = 0;
#ifdef __MPI
    MPI_Comm_rank(ParaV.comm(), &myrank);
#endif
    if (read_type != 1 && read_type != 2)
    {
        if (myrank == 0)
        {
            std::cout << " Error in reading wave function files!\n"
                      << " read_wfc_lcao should be 1 or 2, but got " << read_type << std::endl;
        }
        ModuleBase::timer::end("ModuleIO", "read_wfc_nao");
        return false;
    }
    if (skip_band < 0)
    {
        if (myrank == 0)
        {
            std::cout << " Error in reading wave function files!\n"
                      << " skip_band should not be negative, but got " << skip_band << std::endl;
        }
        ModuleBase::timer::end("ModuleIO", "read_wfc_nao");
        return false;
    }

    int nbands = ParaV.get_wfc_global_nbands(); // the global number of bands
    int nlocal = ParaV.get_wfc_global_nbasis(); // the global number of basis functions
    int nbands_local = ParaV.ncol_bands; // the number of bands in the local process
    int nlocal_local = ParaV.nrow_bands; // the number of basis functions in the local process

    if (gamma_only)
    {
        // I don't know why, but in orther places, the init of psi is using ParaV.ncol for gamma_only case
        // It seems that the diagonalization of double case need a full matrix of psi.
        nbands_local = ParaV.ncol;
    }
    psid.resize(nk, nbands_local, nlocal_local);

    // lambda function to read one file
	auto read_one_file = [&](const std::string& ss, 
			std::stringstream& error_message, 
			const int ik, 
			std::vector<T>& ctot)
    {
        std::ifstream ifs;
        const std::ios_base::openmode mode
            = binary ? (std::ios::in | std::ios::binary) : std::ios::in;
        ifs.open(ss.c_str(), mode);
        if (!ifs)
        {
            error_message << " Can't open file:" << ss << std::endl;
            return false;
        }
        else
		{
            std::cout << " Read NAO wave functions from " << ss << std::endl;
		}

        const auto incomplete_file = [&](const std::string& field) {
            error_message << "The wave function file is incomplete or corrupted while reading "
                          << field << ": " << ss << std::endl;
            ifs.close();
            return false;
        };

        if (!gamma_only)
        {
            int ik_file = 0;
			double kx = 0.0;
			double ky = 0.0;
			double kz = 0.0;
            if (!read_record_value(ifs, ik_file, binary))
            {
                return incomplete_file("the k-point index");
            }
            if (binary)
            {
                if (!read_binary_value(ifs, kx) || !read_binary_value(ifs, ky)
                    || !read_binary_value(ifs, kz))
                {
                    return incomplete_file("the k-point vector");
                }
            }
            else
            {
                ifs >> kx >> ky >> kz;
                if (!ifs)
                {
                    return incomplete_file("the k-point vector");
                }
            }
            if (ik_file != ik + 1)
            {
                error_message << "The k index read in from file do not match the k index generated by ABACUS!\n";
                error_message << " read in k index=" << ik_file;
                error_message << " ABACUS k index =" << ik+1 << std::endl;
                ifs.close();
                return false;
            }
        }
        int nbands_file = 0, nlocal_file = 0;
        if (!read_record_value(ifs, nbands_file, binary)
            || !read_record_value(ifs, nlocal_file, binary))
        {
            return incomplete_file("the dimensions");
        }
        if (nbands_file < 0 || skip_band > nbands_file || nbands > nbands_file - skip_band)
        {
            error_message << "The number of bands to be read exceeds the number of bands in the file generated by ABACUS!\n";
            error_message << " nbands in the existing file=" << nbands_file;
            error_message << " skip_band=" << skip_band;
            error_message << " nbands to be read into ABACUS=" << nbands << std::endl;
            ifs.close();
            return false;
        }
        if (nlocal != nlocal_file)
        {
            error_message << "The nlocal read in from file do not match the nlocal generated by ABACUS!\n";
            error_message << " read in nlocal=" << nlocal_file;
            error_message << " ABACUS nlocal =" << nlocal << std::endl;
            ifs.close();
            return false;
        }
        for (int i = 0; i < skip_band + nbands; i++)
        {
            // the first skip_bands useless bands are read into 0th band to be overwritten
            const int ib_read = std::max(i - skip_band, 0);
            int ib = 0;
            if (!read_record_value(ifs, ib, binary)
                || !read_record_value(ifs, ekb(ik, ib_read), binary)
                || !read_record_value(ifs, wg(ik, ib_read), binary))
            {
                return incomplete_file("band " + std::to_string(i + 1) + " metadata");
            }
            if (i+1 != ib)
            {
                error_message << "The band index read in from file do not match the global parameter band index!\n";
                error_message << " read in band index=" << ib;
                error_message << " band index=" << i+1 << std::endl;
                ifs.close();
                return false;
            }
            for (int j = 0; j < nlocal; j++)
            {
                bool data_read = false;
                if (binary)
                {
                    data_read = read_binary_wfc_data(ifs, ctot[ib_read * nlocal + j]);
                }
                else
                {
                    read_wfc_nao_one_data(ifs, ctot[ib_read * nlocal + j]);
                    data_read = static_cast<bool>(ifs);
                }
                if (!data_read)
                {
                    return incomplete_file("band " + std::to_string(i + 1) + " coefficients");
                }
            }
        }
        ifs.close();
        return true;
    }; // end read one file
        

    std::string errors;

	std::vector<T> ctot;
	if (myrank == 0) 
	{
		ctot.resize(nbands * nlocal);
	}
	else
	{
		ctot.resize(0);
	}

    for(int ik=0;ik<nk;ik++)
    {
        if (myrank == 0)
        {
            const bool out_app_flag = false;
            std::stringstream error_message;
            std::string readin_dir = global_readin_dir;
            if(istep >= 0)
            {
                readin_dir = readin_dir + "WFC/";
            }
            std::string ss = ModuleIO::filename_output(readin_dir,"wf","nao",
                    ik,ik2iktot,nspin,nkstot,read_type,out_app_flag,gamma_only,istep);

            read_success = read_one_file(ss, error_message, ik, ctot);
            errors = error_message.str();
        }   
#ifdef __MPI
        Parallel_Common::bcast_bool(read_success);
        Parallel_Common::bcast_string(errors);
#endif 
        if (!read_success)
        {
            std::cout << " Error in reading wave function files!\n";
            std::cout << errors << std::endl;
            ModuleBase::timer::end("ModuleIO", "read_wfc_nao");
            return false;
        }

        psid.fix_k(ik);
#ifdef __MPI
        Parallel_2D pv_glb;
        pv_glb.set(nlocal, nbands, std::max(nlocal, nbands), ParaV.blacs_ctxt);
        Cpxgemr2d(nlocal,
                  nbands,
                  ctot.data(),
                  1,
                  1,
                  pv_glb.desc,
                  psid.get_pointer(),
                  1,
                  1,
                  const_cast<int*>(ParaV.desc_wfc),
                  pv_glb.blacs_ctxt);
        Parallel_Common::bcast_double(&ekb(ik, 0), nbands);
        Parallel_Common::bcast_double(&wg(ik, 0), nbands);
#else
        BlasConnector::copy(nbands*nlocal, ctot.data(), 1, psid.get_pointer(), 1);
#endif
    }// end of loop over k-points
    
    ModuleBase::timer::end("ModuleIO", "read_wfc_nao");
    return true;    
};

template bool ModuleIO::read_wfc_nao<double>(const std::string& global_readin_dir,
    const Parallel_Orbitals& ParaV,
    psi::Psi<double>& psid,
	ModuleBase::matrix& ekb,
    ModuleBase::matrix& wg,
	const std::vector<int> &ik2iktot,
	const int nkstot,
	const int nspin,
    const int read_type,
    const int skip_band,
    const int istep);

// mohan add 2025-10-19
template bool ModuleIO::read_wfc_nao<float>(const std::string& global_readin_dir,
    const Parallel_Orbitals& ParaV,
    psi::Psi<float>& psid,
	ModuleBase::matrix& ekb,
    ModuleBase::matrix& wg,
	const std::vector<int> &ik2iktot,
	const int nkstot,
	const int nspin,
    const int read_type,
    const int skip_band,
    const int istep);

template bool ModuleIO::read_wfc_nao<std::complex<double>>(const std::string& global_readin_dir,
    const Parallel_Orbitals& ParaV,
	psi::Psi<std::complex<double>>& psid,
	ModuleBase::matrix& ekb,
    ModuleBase::matrix& wg,
	const std::vector<int> &ik2iktot,
	const int nkstot,
	const int nspin,
    const int read_type,
	const int skip_band,
    const int istep);

// mohan add 2025-10-19
template bool ModuleIO::read_wfc_nao<std::complex<float>>(const std::string& global_readin_dir,
    const Parallel_Orbitals& ParaV,
	psi::Psi<std::complex<float>>& psid,
	ModuleBase::matrix& ekb,
    ModuleBase::matrix& wg,
	const std::vector<int> &ik2iktot,
	const int nkstot,
	const int nspin,
    const int read_type,
	const int skip_band,
    const int istep);
