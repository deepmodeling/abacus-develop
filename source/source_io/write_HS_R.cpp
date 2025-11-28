#include "write_HS_R.h"

#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/LCAO_HS_arrays.hpp"
#include "source_lcao/spar_dh.h"
#include "source_lcao/spar_hsr.h"
#include "source_lcao/spar_st.h"
#include "write_HS_sparse.h"

// if 'binary=true', output binary file.
// The 'sparse_thr' is the accuracy of the sparse matrix.
// If the absolute value of the matrix element is less than or equal to the
// 'sparse_thr', it will be ignored.
template <typename TK>
void ModuleIO::output_HSR(const UnitCell& ucell,
                          const int& istep,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
						  const K_Vectors& kv,
						  Plus_U &dftu, // mohan add 20251107
						  hamilt::Hamilt<TK>* p_ham,
#ifdef __EXX
                          const std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd,
                          const std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc,
#endif
                          const std::string& SR_filename,
                          const std::string& HR_filename_up,
                          const std::string HR_filename_down,
                          const bool& binary,
                          const double& sparse_thr)
{

    ModuleBase::TITLE("ModuleIO", "output_HSR");
    ModuleBase::timer::tick("ModuleIO", "output_HSR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |     #Print out Hamiltonian matrix H(R) or overlap matrix S(R)#     |" << std::endl;
    GlobalV::ofs_running << " | Use numerical atomic orbitals basis. Here R is the Bravis lattice  |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    const int nspin = PARAM.inp.nspin;

    if (nspin == 1 || nspin == 4)
    {
        const int spin_now = 0;
        // jingan add 2021-6-4, modify 2021-12-2
		sparse_format::cal_HSR(ucell,
				dftu, // mohan add 20251107
				pv,
				HS_Arrays,
				grid,
				spin_now,
				sparse_thr,
				kv.nmp,
				p_ham
#ifdef __EXX
				,
				Hexxd,
				Hexxc
#endif
        );
    }
    else if (nspin == 2)
    {
        int spin_now = 1;

        // save HR of spin down first (the current spin always be down)
		sparse_format::cal_HSR(ucell,
				dftu,
				pv,
				HS_Arrays,
				grid,
				spin_now,
				sparse_thr,
				kv.nmp,
				p_ham
#ifdef __EXX
				,
				Hexxd,
				Hexxc
#endif
				);

        // cal HR of the spin up
        if (PARAM.inp.vl_in_h)
        {
            const int ik = 0;
            p_ham->refresh();
            p_ham->updateHk(ik);
            spin_now = 0;
        }

		sparse_format::cal_HSR(ucell,
				dftu,
				pv,
				HS_Arrays,
				grid,
				spin_now,
				sparse_thr,
				kv.nmp,
				p_ham
#ifdef __EXX
				,
				Hexxd,
				Hexxc
#endif
        );
    }

    ModuleIO::save_HSR_sparse(istep, pv, HS_Arrays, sparse_thr, binary, SR_filename, HR_filename_up, HR_filename_down);

    sparse_format::destroy_HS_R_sparse(HS_Arrays);

    ModuleBase::timer::tick("ModuleIO", "output_HSR");
    return;
}

void ModuleIO::output_dSR(const int& istep,
                          const UnitCell& ucell,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
                          const TwoCenterBundle& two_center_bundle,
                          const LCAO_Orbitals& orb,
                          const K_Vectors& kv,
                          const bool& binary,
                          const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_dSR");
    ModuleBase::timer::tick("ModuleIO", "output_dSR");

    sparse_format::cal_dS(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, sparse_thr);

    // mohan update 2024-04-01
    ModuleIO::save_dH_sparse(istep, pv, HS_Arrays, sparse_thr, binary, "s");

    sparse_format::destroy_dH_R_sparse(HS_Arrays);

    ModuleBase::timer::tick("ModuleIO", "output_dSR");
    return;
}

void ModuleIO::output_dHR(const int& istep,
                          const ModuleBase::matrix& v_eff,
                          const UnitCell& ucell,
                          const Parallel_Orbitals& pv,
                          LCAO_HS_Arrays& HS_Arrays,
                          const Grid_Driver& grid, // mohan add 2024-04-06
                          const TwoCenterBundle& two_center_bundle,
                          const LCAO_Orbitals& orb,
                          const K_Vectors& kv,
                          const bool& binary,
                          const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_dHR");
    ModuleBase::timer::tick("ModuleIO", "output_dHR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |                         #Print out dH/dR#                          |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    const int nspin = PARAM.inp.nspin;

    if (nspin == 1 || nspin == 4)
    {
        // mohan add 2024-04-01
        const int cspin = 0;

        sparse_format::cal_dH(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, cspin, sparse_thr, v_eff);
    }
    else if (nspin == 2)
    {
        for (int cspin = 0; cspin < 2; cspin++)
        {
            sparse_format::cal_dH(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, cspin, sparse_thr, v_eff);
        }
    }
    // mohan update 2024-04-01
    ModuleIO::save_dH_sparse(istep, pv, HS_Arrays, sparse_thr, binary);

    sparse_format::destroy_dH_R_sparse(HS_Arrays);

    ModuleBase::timer::tick("ModuleIO", "output_dHR");
    return;
}

template <typename TK>
void ModuleIO::output_SR(Parallel_Orbitals& pv,
                         const Grid_Driver& grid,
                         hamilt::Hamilt<TK>* p_ham,
                         const std::string& SR_filename,
                         const bool& binary,
                         const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_SR");
    ModuleBase::timer::tick("ModuleIO", "output_SR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |                 #Print out overlap matrix S(R)#                    |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    std::cout << " Overlap matrix file is in " << SR_filename << std::endl;
    GlobalV::ofs_running << " Overlap matrix file is in " << SR_filename << std::endl;

    LCAO_HS_Arrays HS_Arrays;

    sparse_format::cal_SR(pv,
                          HS_Arrays.all_R_coor,
                          HS_Arrays.SR_sparse,
                          HS_Arrays.SR_soc_sparse,
                          grid,
                          sparse_thr,
                          p_ham);

    const int istep = 0;

    if (PARAM.inp.nspin == 4)
    {
        ModuleIO::save_sparse(HS_Arrays.SR_soc_sparse,
                              HS_Arrays.all_R_coor,
                              sparse_thr,
                              binary,
                              SR_filename,
                              pv,
                              "S",
                              istep);
    }
    else
    {
        ModuleIO::save_sparse(HS_Arrays.SR_sparse,
                              HS_Arrays.all_R_coor,
                              sparse_thr,
                              binary,
                              SR_filename,
                              pv,
                              "S",
                              istep);
    }

    sparse_format::destroy_HS_R_sparse(HS_Arrays);

    ModuleBase::timer::tick("ModuleIO", "output_SR");
    return;
}

void ModuleIO::output_TR(const int istep,
                         const UnitCell& ucell,
                         const Parallel_Orbitals& pv,
                         LCAO_HS_Arrays& HS_Arrays,
                         const Grid_Driver& grid,
                         const TwoCenterBundle& two_center_bundle,
                         const LCAO_Orbitals& orb,
                         const std::string& TR_filename,
                         const bool& binary,
                         const double& sparse_thr)
{
    ModuleBase::TITLE("ModuleIO", "output_TR");
    ModuleBase::timer::tick("ModuleIO", "output_TR");

    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " |           #Print out kinetic energy term matrix T(R)#              |" << std::endl;
    GlobalV::ofs_running << " |                                                                    |" << std::endl;
    GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;

    std::stringstream sst;
    if (PARAM.inp.calculation == "md" && !PARAM.inp.out_app_flag)
    {
        sst << PARAM.globalv.global_matrix_dir << TR_filename << "g" << istep;
        GlobalV::ofs_running << " T(R) data are in file: " << sst.str() << std::endl;
    }
    else
    {
        sst << PARAM.globalv.global_out_dir << TR_filename;
        GlobalV::ofs_running << " T(R) data are in file: " << sst.str() << std::endl;
    }

    sparse_format::cal_TR(ucell, pv, HS_Arrays, grid, two_center_bundle, orb, sparse_thr);

    ModuleIO::save_sparse(HS_Arrays.TR_sparse,
                          HS_Arrays.all_R_coor,
                          sparse_thr,
                          binary,
                          sst.str().c_str(),
                          pv,
                          "T",
                          istep);

    sparse_format::destroy_T_R_sparse(HS_Arrays);

    ModuleBase::timer::tick("ModuleIO", "output_TR");
    return;
}

template void ModuleIO::output_HSR<double>(
    const UnitCell& ucell,
    const int& istep,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& HS_Arrays,
    const Grid_Driver& grid,
    const K_Vectors& kv,
	Plus_U &dftu, // mohan add 20251107
    hamilt::Hamilt<double>* p_ham,
#ifdef __EXX
    const std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd,
    const std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc,
#endif
    const std::string& SR_filename,
    const std::string& HR_filename_up,
    const std::string HR_filename_down,
    const bool& binary,
    const double& sparse_thr);

template void ModuleIO::output_HSR<std::complex<double>>(
    const UnitCell& ucell,
    const int& istep,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& HS_Arrays,
    const Grid_Driver& grid,
    const K_Vectors& kv,
	Plus_U &dftu, // mohan add 20251107
    hamilt::Hamilt<std::complex<double>>* p_ham,
#ifdef __EXX
    const std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd,
    const std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc,
#endif
    const std::string& SR_filename,
    const std::string& HR_filename_up,
    const std::string HR_filename_down,
    const bool& binary,
    const double& sparse_thr);

template void ModuleIO::output_SR<double>(Parallel_Orbitals& pv,
                                          const Grid_Driver& grid,
                                          hamilt::Hamilt<double>* p_ham,
                                          const std::string& SR_filename,
                                          const bool& binary,
                                          const double& sparse_thr);
template void ModuleIO::output_SR<std::complex<double>>(Parallel_Orbitals& pv,
                                                        const Grid_Driver& grid,
                                                        hamilt::Hamilt<std::complex<double>>* p_ham,
                                                        const std::string& SR_filename,
                                                        const bool& binary,
                                                        const double& sparse_thr);
