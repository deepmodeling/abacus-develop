#include "hsk_writer.h"
#include "hs_dense_io.h"
#include "source_base/module_out/filename.h" // use filename_output function
#include <complex>

template <typename T>
void ModuleIO::write_hsk(
        const std::string &global_out_dir,
        const int nspin,
        const int nks,
        const int nkstot,
        const std::vector<int> &ik2iktot,
        const std::vector<int> &isk,
        hamilt::Hamilt<T>* p_hamilt,
        const Parallel_Orbitals &pv,
        const bool gamma_only,
        const bool out_app_flag,
        const int istep,
        const int out_type,
        const int precision,
        const int nlocal,
        const std::string &ks_solver,
        const int drank,
        std::ofstream &ofs_running)
{

    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        ">>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs_running << " |                                            "
        "                        |" << std::endl;
    ofs_running << " | Write Hamiltonian matrix H(k) or overlap matrix S(k) in numerical  |" << std::endl;
    ofs_running << " | atomic orbitals at each k-point.                                   |" << std::endl;
    ofs_running << " |                                            "
        "                        |" << std::endl;
    ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        ">>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs_running << "\n WRITE H(k) OR S(k)" << std::endl;

    for (int ik = 0; ik < nks; ++ik)
    {
        p_hamilt->updateHk(ik);
        const bool binary = (out_type == 2);

        hamilt::MatrixBlock<T> h_mat;
        hamilt::MatrixBlock<T> s_mat;

        p_hamilt->matrix(h_mat, s_mat);

        std::string h_fn = ModuleIO::filename_output(global_out_dir,
                "hk","nao",ik,ik2iktot,nspin,nkstot,
                out_type,out_app_flag,gamma_only,istep);

        ModuleIO::save_mat(istep,
                h_mat.p,
                nlocal,
                binary,
                precision,
                1,
                out_app_flag,
                h_fn,
                pv,
                drank,
                ks_solver);

        // mohan note 2025-06-02
        // for overlap matrix, the two spin channels yield the same matrix
        // so we only need to print matrix from one spin channel.
        const int current_spin = isk[ik];
        if(current_spin == 1)
        {
            continue;
        }

        std::string s_fn = ModuleIO::filename_output(global_out_dir,
                "sk","nao",ik,ik2iktot,nspin,nkstot,
                out_type,out_app_flag,gamma_only,istep);

        ofs_running << " The output filename is " << s_fn << std::endl;

        ModuleIO::save_mat(istep,
                s_mat.p,
                nlocal,
                binary,
                precision,
                1,
                out_app_flag,
                s_fn,
                pv,
                drank,
                ks_solver);
    } // end ik
}

// Explicit instantiations
template void ModuleIO::write_hsk<double>(
    const std::string&, const int, const int, const int,
    const std::vector<int>&, const std::vector<int>&,
    hamilt::Hamilt<double>*, const Parallel_Orbitals&,
    const bool, const bool, const int, const int, const int,
    const int, const std::string&, const int, std::ofstream&);
template void ModuleIO::write_hsk<std::complex<double>>(
    const std::string&, const int, const int, const int,
    const std::vector<int>&, const std::vector<int>&,
    hamilt::Hamilt<std::complex<double>>*, const Parallel_Orbitals&,
    const bool, const bool, const int, const int, const int,
    const int, const std::string&, const int, std::ofstream&);
