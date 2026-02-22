#include "source_io/module_dm/write_dmk.h"

#include "source_base/parallel_common.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_output/ucell_io.h"

std::string ModuleIO::dmk_gen_fname(const bool gamma_only, const int ispin, const int nspin, const int ik, const int istep)
{
    // set istep = -1 if you don't want the 'g' index appears in the file name
    assert(istep>=-1);

    // ik should be the correct one

    std::string fname = "dm";

    if (!gamma_only)
    {
        fname += "k" + std::to_string(ik + 1);
    }

    if (nspin == 2)
    {
        fname += "s" + std::to_string(ispin + 1);
    }

    if( istep >= 0 )
    {
        fname += "g" + std::to_string(istep + 1);
    }

    fname += "_nao.txt";

    return fname;
}



void ModuleIO::dmk_readData(std::ifstream& ifs, double& data)
{
    ifs >> data;
}

void ModuleIO::dmk_readData(std::ifstream& ifs, std::complex<double>& data)
{
    std::string complex_str;
    ifs >> complex_str;

    size_t comma_pos = complex_str.find(',');
    if (complex_str.front() == '(' && complex_str.back() == ')' && comma_pos != std::string::npos)
    {
        double real = std::stod(complex_str.substr(1, comma_pos - 1));
        double imag = std::stod(complex_str.substr(comma_pos + 1, complex_str.size() - comma_pos - 2));
        data = std::complex<double>(real, imag);
    }
    else
    {
        ModuleBase::WARNING_QUIT("ModuleIO::dmk_readData",
        "Invalid complex number format in dmk: " + complex_str);
    }
}

template <typename T>
bool ModuleIO::read_dmk(const int nspin,
		const int nk,
		const K_Vectors &kv,
		const Parallel_2D& pv,
		const std::string& dmk_dir,
		std::vector<std::vector<T>>& dmk,
		std::ofstream &ofs_running)
{
    ModuleBase::TITLE("ModuleIO", "read_dmk");
    ModuleBase::timer::tick("ModuleIO", "read_dmk");

    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(pv.comm(), &my_rank);
#endif

    int nlocal = pv.get_global_row_size();
    bool gamma_only = std::is_same<double, T>::value;
    std::vector<std::vector<T>> dmk_global(nspin * nk, std::vector<T>(nlocal * nlocal, 0));

    bool read_success = true;
    std::string tmp;
    if (my_rank == 0)
    {
        for (int ispin = 0; ispin < nspin; ispin++)
        {
            for (int ik = 0; ik < nk; ik++)
            {
                // to read density matrix in k space, remember to delete the step information 'g'
                // set istep = -1 if you don't want the 'g' index appears in the file name
                const int istep = -1;
                std::string fn = dmk_dir + dmk_gen_fname(gamma_only, ispin, nspin, ik, istep);
                std::ifstream ifs(fn.c_str());

                if (!ifs)
                {
                    ofs_running << " Cannot find density matrix file " << fn << " for k-point " << ik+1 << std::endl;
                    ModuleBase::WARNING("ModuleIO::read_dmk", "Can't open density matrix (k) file < " + fn + " >.");
                    read_success = false;
                    break;
                }
                else
                {
                    ofs_running << " Read density matrix file " << fn << " for k-point " << ik+1 << std::endl;
                }

                int spin_tmp = 0;
                ModuleBase::GlobalFunc::READ_VALUE(ifs, spin_tmp);

                double fermi_tmp = 0.0;
                ModuleBase::GlobalFunc::READ_VALUE(ifs, fermi_tmp);

                int nlocal_tmp = 0;
                ModuleBase::GlobalFunc::READ_VALUE(ifs, nlocal_tmp);

                if(nlocal_tmp==nlocal)
                {
                    ofs_running << " number of basis (nlocal) is correct: " << nlocal << std::endl;
                }
                else
                {
                    ModuleBase::WARNING_QUIT("ModuleIO::read_dmk","nlocal does not match!");
                }

                // read the DMK data
                const size_t index_k = ik + nk * ispin;
                for (int i = 0; i < nlocal; ++i)
                {
                    const size_t index_i = i * nlocal;
                    for (int j = 0; j < nlocal; ++j)
                    {
                        dmk_readData(ifs, dmk_global[index_k][index_i + j]);
                    }
                }
                ifs.close();
            } // ik
            if (!read_success)
            {
                break;
            }
        } // ispin
    }     // rank0

#ifdef __MPI
    MPI_Bcast(&read_success, 1, MPI_C_BOOL, 0, pv.comm());
#endif

    if (read_success)
    {
#ifdef __MPI
        // seperate dmk data to each processor with 2D block distribution
        dmk.resize(nspin * nk, std::vector<T>(pv.get_row_size() * pv.get_col_size()));
        Parallel_2D pv_glb;
        pv_glb.set(nlocal, nlocal, nlocal, pv.blacs_ctxt);
        for (int ik = 0; ik < nspin * nk; ik++)
        {
            Cpxgemr2d(
            nlocal,
            nlocal,
            dmk_global[ik].data(),
            1,
            1,
            pv_glb.desc,
            dmk[ik].data(),
            1,
            1,
            const_cast<int*>(pv.desc),
            pv_glb.blacs_ctxt);
        }
#else
        dmk = dmk_global;
#endif
    }
    ModuleBase::timer::tick("ModuleIO", "read_dmk");
    return read_success;
}

template <typename T>
void ModuleIO::write_dmk(const std::vector<std::vector<T>>& dmk,
		const K_Vectors &kv,
		const int precision,
		const std::vector<double>& efs,
		const UnitCell* ucell,
		const Parallel_2D& pv,
		const int istep)
{
    ModuleBase::TITLE("ModuleIO", "write_dmk");
    ModuleBase::timer::tick("ModuleIO", "write_dmk");

    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(pv.comm(), &my_rank);
#endif

    bool gamma_only = std::is_same<double, T>::value;
    const int nlocal = pv.get_global_row_size();
    const int nspin = efs.size();
    assert(nspin > 0);
    const int nk = dmk.size() / nspin;
    const double dm_thr = 1.0e-16; // mohan set 2025-09-02

    if (nk * nspin != dmk.size())
    {
        ModuleBase::WARNING_QUIT("ModuleIO::write_dmk", "The size of dmk is not consistent with nspin and nk.");
    }

    Parallel_2D pv_glb;

    for (int ispin = 0; ispin < nspin; ispin++)
    {
        for (int ik = 0; ik < nk; ik++)
        {
            // gather dmk[ik] to dmk_global
            std::vector<T> dmk_global(my_rank == 0 ? nlocal * nlocal : 0);
#ifdef __MPI
            pv_glb.set(nlocal, nlocal, nlocal, pv.blacs_ctxt);
            Cpxgemr2d(
            nlocal,
            nlocal,
            const_cast<T*>(dmk[ik + nk * ispin].data()),
            1,
            1,
            const_cast<int*>(pv.desc),
            dmk_global.data(),
            1,
            1,
            pv_glb.desc,
            pv_glb.blacs_ctxt);
#else
            dmk_global = dmk[ik + nk * ispin];
#endif

            if (my_rank == 0)
            {
                std::string fn = PARAM.globalv.global_out_dir 
			+ dmk_gen_fname(gamma_only, ispin, nspin, ik, istep);

                std::ofstream ofs(fn.c_str());

                ofs << " --- Ionic Step " << istep+1 << " ---" << std::endl;

                if (!ofs)
                {
                    ModuleBase::WARNING("ModuleIO::write_dmk", "Can't create DENSITY MATRIX File < " + fn + " >.");
                    continue;
                }
                else
                {
                    //std::cout << " Write the density matrix to file " << fn << std::endl;
                }


		// information about density matrix at this k-point
		ofs << " " << nspin << " # number of spin directions" << std::endl;
                ofs << " " << ispin+1 << " # spin index" << std::endl; 
		ofs << " " << kv.get_nkstot_full() << " # total k points " << std::endl;
		ofs << " " << kv.get_nkstot() << " # total k points after symmetrized (if open) " << std::endl;
		ofs << " " << ik+1 << " # k-point index " << std::endl;
                ofs << " " << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z  
			<< " # k point coordinate (Cartesian) " << std::endl;
                ofs << " " << kv.kvec_d[ik].x << " " << kv.kvec_d[ik].y << " " << kv.kvec_d[ik].z  
			<< " # k point coordinate (direct) " << std::endl;
		ofs << " " << kv.wk[ik] << " # weight of this k point" << std::endl;
                ofs << " " << efs[ispin] << " # Fermi energy in Ry " << std::endl;
		ofs << " " << nlocal << " # number of localized basis " << std::endl;
		ofs << " " << nlocal << " " << nlocal << " # size of this matrix " << std::endl; 
		ofs << std::endl;

		// write ucell
		ModuleIO::UcellIO::write_ucell(ofs, ucell);

                ofs << std::fixed;
                ofs << std::scientific;
                ofs << std::setprecision(precision);
                ofs << std::right;
                for (int i = 0; i < nlocal; ++i)
                {
		    const size_t ii = i * nlocal;
                    for (int j = 0; j < nlocal; ++j)
                    {
                        if (std::is_same<double, T>::value)
                        {
                            if (j % 8 == 0)
                            {
                                ofs << "\n";
                            }
                            ofs << " " << dmk_global[ii + j];
                        }
                        else if (std::is_same<std::complex<double>, T>::value)
                        {
                            if (j % 4 == 0)
                            {
                                ofs << "\n";
                            }

                            double real_v = std::real(dmk_global[ii + j]);
                            if(std::abs(real_v) < dm_thr)
                            {
                                real_v = 0.0;
                            }
                            double imag_v = std::imag(dmk_global[ii + j]);
                            if(std::abs(imag_v) < dm_thr)
                            {
                                imag_v = 0.0;
                            }

                            ofs << " (" << real_v << "," << imag_v << ")";
                        }
                    }
                }
                ofs.close();
            } // rank0
        }     // ik
    }         // ispin

    ModuleBase::timer::tick("ModuleIO", "write_dmk");
}

template bool ModuleIO::read_dmk<double>(const int nspin,
		const int nk,
		const K_Vectors &kv,
		const Parallel_2D& pv,
		const std::string& dmk_dir,
		std::vector<std::vector<double>>& dmk,
		std::ofstream &ofs);

template bool ModuleIO::read_dmk<std::complex<double>>(const int nspin,
		const int nk,
		const K_Vectors &kv,
		const Parallel_2D& pv,
		const std::string& dmk_dir,
		std::vector<std::vector<std::complex<double>>>& dmk,
		std::ofstream &ofs);

template void ModuleIO::write_dmk<double>(const std::vector<std::vector<double>>& dmk,
		const K_Vectors &kv,
		const int precision,
		const std::vector<double>& efs,
		const UnitCell* ucell,
		const Parallel_2D& pv,
		const int istep);

template void ModuleIO::write_dmk<std::complex<double>>(const std::vector<std::vector<std::complex<double>>>& dmk,
		const K_Vectors &kv,
		const int precision,
		const std::vector<double>& efs,
		const UnitCell* ucell,
		const Parallel_2D& pv,
		const int istep);

