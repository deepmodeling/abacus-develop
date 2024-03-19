#ifndef EXX_LRI_INTERFACE_HPP
#define EXX_LRI_INTERFACE_HPP

#include "Exx_LRI_interface.h"
#include "module_ri/exx_abfs-jle.h"
#include "module_ri/exx_opt_orb.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"

#include <sys/time.h>

template<typename T, typename Tdata>
void Exx_LRI_Interface<T, Tdata>::write_Hexxs(const std::string& file_name) const
{
	ModuleBase::TITLE("Exx_LRI","write_Hexxs");
	ModuleBase::timer::tick("Exx_LRI", "write_Hexxs");
	std::ofstream ofs(file_name, std::ofstream::binary);
	cereal::BinaryOutputArchive oar(ofs);
	oar(this->exx_ptr->Hexxs);
	ModuleBase::timer::tick("Exx_LRI", "write_Hexxs");
}

template<typename T, typename Tdata>
void Exx_LRI_Interface<T, Tdata>::read_Hexxs(const std::string& file_name)
{
	ModuleBase::TITLE("Exx_LRI","read_Hexxs");
	ModuleBase::timer::tick("Exx_LRI", "read_Hexxs");
	std::ifstream ifs(file_name, std::ofstream::binary);
	cereal::BinaryInputArchive iar(ifs);
	iar(this->exx_ptr->Hexxs);
	ModuleBase::timer::tick("Exx_LRI", "read_Hexxs");
}
template<typename T, typename Tdata>
void Exx_LRI_Interface<T, Tdata>::exx_beforescf(const K_Vectors& kv, const Charge_Mixing& chgmix)
{
#ifdef __MPI
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        if (GlobalC::restart.info_load.load_H_finish && !GlobalC::restart.info_load.restart_exx) XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);
        else
        {
            if (GlobalC::ucell.atoms[0].ncpp.xc_func == "HF" || GlobalC::ucell.atoms[0].ncpp.xc_func == "PBE0" || GlobalC::ucell.atoms[0].ncpp.xc_func == "HSE")
            {
                XC_Functional::set_xc_type("pbe");
            }
            else if (GlobalC::ucell.atoms[0].ncpp.xc_func == "SCAN0")
            {
                XC_Functional::set_xc_type("scan");
            }
        }
        this->exx_ptr->cal_exx_ions();
    }

		if (Exx_Abfs::Jle::generate_matrix)
		{
			//program should be stopped after this judgement
			Exx_Opt_Orb exx_opt_orb;
			exx_opt_orb.generate_matrix(kv);
			ModuleBase::timer::tick("ESolver_KS_LCAO", "beforescf");
			return;
		}
		
		// set initial parameter for mix_DMk_2D
		if(GlobalC::exx_info.info_global.cal_exx)
		{
			this->mix_DMk_2D.set_nks(kv.nks, GlobalV::GAMMA_ONLY_LOCAL);
			if(GlobalC::exx_info.info_global.separate_loop)
                this->mix_DMk_2D.set_mixing(nullptr);
			else
				this->mix_DMk_2D.set_mixing(chgmix.get_mixing());
        }
        // for exx two_level scf
        this->two_level_step = 0;
#endif // __MPI
}

template<typename T, typename Tdata>
void Exx_LRI_Interface<T, Tdata>::exx_eachiterinit(const elecstate::DensityMatrix<T, double>& dm, const int& iter)
{
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        if (!GlobalC::exx_info.info_global.separate_loop && this->two_level_step)
        {
			const bool flag_restart = (iter==1) ? true : false;
            this->mix_DMk_2D.mix(dm.get_DMK_vector(), flag_restart);
			const std::vector<std::map<int,std::map<std::pair<int, std::array<int, 3>>,RI::Tensor<Tdata>>>>
				Ds = GlobalV::GAMMA_ONLY_LOCAL
					? RI_2D_Comm::split_m2D_ktoR<Tdata>(*this->exx_ptr->p_kv, this->mix_DMk_2D.get_DMk_gamma_out(), *dm.get_paraV_pointer())
					: RI_2D_Comm::split_m2D_ktoR<Tdata>(*this->exx_ptr->p_kv, this->mix_DMk_2D.get_DMk_k_out(), *dm.get_paraV_pointer());
            this->exx_ptr->cal_exx_elec(Ds, *dm.get_paraV_pointer());
        }
    }
}

template<typename T, typename Tdata>
void Exx_LRI_Interface<T, Tdata>::exx_hamilt2density(elecstate::ElecState& elec, const Parallel_Orbitals& pv, const  int iter)
{
    // Peize Lin add 2020.04.04
    if (XC_Functional::get_func_type() == 4 || XC_Functional::get_func_type() == 5)
    {
        // add exx
        // Peize Lin add 2016-12-03
        if (GlobalC::restart.info_load.load_H_finish && !GlobalC::restart.info_load.restart_exx
            && this->two_level_step == 0 && iter == 1)
        {
            if (GlobalV::MY_RANK == 0)GlobalC::restart.load_disk("Eexx", 0, 1, &this->exx_ptr->Eexx);
            Parallel_Common::bcast_double(this->exx_ptr->Eexx);
            this->exx_ptr->Eexx /= GlobalC::exx_info.info_global.hybrid_alpha;
        }
        elec.set_exx(this->get_Eexx());
    }
    else
    {
        elec.f_en.exx = 0.;
    }
}

template<typename T, typename Tdata>
bool Exx_LRI_Interface<T, Tdata>::exx_after_converge(
    hamilt::Hamilt<T>& hamilt,
    LCAO_Matrix& lm,
    const elecstate::DensityMatrix<T, double>& dm,
    const K_Vectors& kv,
    int& iter)
{
    auto restart_reset = [this]()
        { // avoid calling restart related procedure in the subsequent ion steps
            GlobalC::restart.info_load.restart_exx = true;
            this->exx_ptr->Eexx = 0;
        };
    if (GlobalC::exx_info.info_global.cal_exx)
    {
        // no separate_loop case
        if (!GlobalC::exx_info.info_global.separate_loop)
        {
            GlobalC::exx_info.info_global.hybrid_step = 1;

            // in no_separate_loop case, scf loop only did twice
            // in first scf loop, exx updated once in beginning,
            // in second scf loop, exx updated every iter

            if (this->two_level_step)
            {
                restart_reset();
                return true;
            }
            else
            {
                // update exx and redo scf
                XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);
                iter = 0;
                std::cout << " Entering 2nd SCF, where EXX is updated" << std::endl;
                this->two_level_step++;
                return false;
            }
        }
        // has separate_loop case
        // exx converged or get max exx steps
        else if (this->two_level_step == GlobalC::exx_info.info_global.hybrid_step
            || (iter == 1 && this->two_level_step != 0))
        {
            restart_reset();
            return true;
        }
        else
        {
            // update exx and redo scf
            if (this->two_level_step == 0)
            {
                XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].ncpp.xc_func);
            }

            std::cout << " Updating EXX " << std::flush;
            timeval t_start;       gettimeofday(&t_start, NULL);

            const bool flag_restart = (this->two_level_step == 0) ? true : false;
            this->mix_DMk_2D.mix(dm.get_DMK_vector(), flag_restart);

            // GlobalC::exx_lcao.cal_exx_elec(p_esolver->LOC, p_esolver->LOWF.wfc_k_grid);
			const std::vector<std::map<int,std::map<std::pair<int, std::array<int, 3>>,RI::Tensor<Tdata>>>>
				Ds = GlobalV::GAMMA_ONLY_LOCAL
					? RI_2D_Comm::split_m2D_ktoR<Tdata>(*this->exx_ptr->p_kv, this->mix_DMk_2D.get_DMk_gamma_out(), *dm.get_paraV_pointer())
					: RI_2D_Comm::split_m2D_ktoR<Tdata>(*this->exx_ptr->p_kv, this->mix_DMk_2D.get_DMk_k_out(), *dm.get_paraV_pointer());
            this->exx_ptr->cal_exx_elec(Ds, *dm.get_paraV_pointer());
            iter = 0;
            this->two_level_step++;
            
            timeval t_end;       gettimeofday(&t_end, NULL);
            std::cout << "and rerun SCF\t"
                << std::setprecision(3) << std::setiosflags(std::ios::scientific)
                << (double)(t_end.tv_sec-t_start.tv_sec) + (double)(t_end.tv_usec-t_start.tv_usec)/1000000.0 
                << std::defaultfloat << " (s)" << std::endl;
            return false;
        }
    }
    restart_reset();
    return true;
}
#endif