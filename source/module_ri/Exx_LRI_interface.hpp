#ifndef EXX_LRI_INTERFACE_HPP
#define EXX_LRI_INTERFACE_HPP

#include "Exx_LRI_interface.h"
#include "module_ri/exx_abfs-jle.h"
#include "module_ri/exx_opt_orb.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/op_exx_lcao.h"

#include <sys/time.h>
#include "module_io/csr_reader.h"
#include "module_io/write_HS_sparse.h"

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

template<typename Tdata>
inline void print_tensor(const RI::Tensor<Tdata>& t, const std::string& name, const double& threshold = 0.0)
{
    std::cout << name << ":\n";
    for (int i = 0;i < t.shape[0];++i)
    {
        for (int j = 0;j < t.shape[1];++j)
            std::cout << ((std::abs(t(i, j)) > threshold) ? t(i, j) : static_cast<Tdata>(0)) << " ";
        std::cout << std::endl;
    }
}
template<typename Tdata>
inline void print_HR(const std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>>& HR, const std::string name, const double& threshold)
{
    for (auto& HR_ia1 : HR)
    {
        int iat1 = HR_ia1.first;
        for (auto& HR_ia12R : HR_ia1.second)
        {
            int iat2 = HR_ia12R.first.first;
            std::array<int, 3> R = HR_ia12R.first.second;
            const RI::Tensor<Tdata>& HR_tensor = HR_ia12R.second;
            std::cout << "atom pair (" << iat1 << ", " << iat2 << "), R=(" << R[0] << "," << R[1] << "," << R[2] << "), ";
            print_tensor(HR_tensor, name, threshold);
        }
    }
}

template<typename T, typename Tdata>
void Exx_LRI_Interface<T, Tdata>::write_Hexxs(const std::string& file_name, const Parallel_Orbitals& pv) const
{
    ModuleBase::TITLE("Exx_LRI", "write_Hexxs");
    ModuleBase::timer::tick("Exx_LRI", "write_Hexxs");
    std::set<Abfs::Vector3_Order<int>> all_R_coor;
    double sparse_threshold = 1e-10;
    for (int is = 0;is < this->exx_ptr->Hexxs.size();++is)
    {
        for (const auto& HexxA : this->exx_ptr->Hexxs[is])
        {
            const int iat0 = HexxA.first;
            for (const auto& HexxB : HexxA.second)
            {
                const int iat1 = HexxB.first.first;
                const Abfs::Vector3_Order<int> R = RI_Util::array3_to_Vector3(HexxB.first.second);
                all_R_coor.insert(R);
            }
        }
        ModuleIO::save_sparse(
            this->calculate_RI_Tensor_sparse(sparse_threshold, this->exx_ptr->Hexxs[is], pv),
            all_R_coor,
            sparse_threshold,
            false, //binary
            file_name + "_" + std::to_string(is) + ".csr",
            pv,
            "Hexxs_" + std::to_string(is)
        );
    }

    // test: output Hexxs
    print_HR(this->exx_ptr->Hexxs[0], "Hexxs[0]_write", sparse_threshold);

    ModuleBase::timer::tick("Exx_LRI", "write_Hexxs");
}


template<typename T, typename Tdata>
std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, Tdata>>>
Exx_LRI_Interface<T, Tdata>::calculate_RI_Tensor_sparse(const double& sparse_threshold,
    const std::map<int, std::map<TAC, RI::Tensor<Tdata>>>& hR,
    const Parallel_Orbitals& paraV)const
{
    ModuleBase::TITLE("Exx_LRI_Interface", "calculate_HContainer_sparse_d");
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, Tdata>>> target;
    auto row_indexes = paraV.get_indexes_row();
    auto col_indexes = paraV.get_indexes_col();
    // for (int iap = 0;iap < hR.size_atom_pairs();++iap)
    for (auto& a1_a2R_data : hR)
    {
        int atom_i = a1_a2R_data.first;
        for (auto& a2R_data : a1_a2R_data.second)
        {
            int atom_j = a2R_data.first.first;
            int start_i = paraV.atom_begin_row[atom_i];
            int start_j = paraV.atom_begin_col[atom_j];
            int row_size = paraV.get_row_size(atom_i);
            int col_size = paraV.get_col_size(atom_j);
            const TC& R = a2R_data.first.second;
            auto& matrix = a2R_data.second;
            Abfs::Vector3_Order<int> dR(R[0], R[1], R[2]);
            for (int i = 0;i < row_size;++i)
            {
                int mu = row_indexes[start_i + i];
                for (int j = 0;j < col_size;++j)
                {
                    int nu = col_indexes[start_j + j];
                    const auto& value_tmp = matrix(i, j);
                    if (std::abs(value_tmp) > sparse_threshold)
                    {
                        target[dR][mu][nu] = value_tmp;
                    }
                }
            }
        }
    }
    return target;
}
template<typename T, typename Tdata>
void Exx_LRI_Interface<T, Tdata>::read_Hexxs(const std::string& file_name, const Parallel_Orbitals& pv, const UnitCell& ucell)
{
    ModuleBase::TITLE("Exx_LRI", "read_Hexxs");
    ModuleBase::timer::tick("Exx_LRI", "read_Hexxs");
    this->exx_ptr->Hexxs.resize(GlobalV::NSPIN);
    for (int is = 0;is < GlobalV::NSPIN;++is)
    {
        ModuleIO::csrFileReader<Tdata> csr(file_name + "_" + std::to_string(is) + ".csr");
        int nR = csr.getNumberOfR();
        int nbasis = csr.getMatrixDimension();
        assert(nbasis == GlobalV::NLOCAL);
        // allocate Hexxs[is]
        for (int iat1 = 0; iat1 < ucell.nat; ++iat1)
            for (int iat2 = 0;iat2 < ucell.nat;++iat2)
                for (int iR = 0;iR < nR;++iR)
                {
                    const std::vector<int>& R = csr.getRCoordinate(iR);
                    TC dR({ R[0], R[1], R[2] });
                    this->exx_ptr->Hexxs[is][iat1][{iat2, dR}] = RI::Tensor<Tdata>({ static_cast<size_t>(ucell.atoms[ucell.iat2it[iat1]].nw), static_cast<size_t>(ucell.atoms[ucell.iat2it[iat2]].nw) });
                }
        // read Hexxs[is]
        for (int i = 0;i < GlobalV::NLOCAL;++i)
            for (int j = 0;j < GlobalV::NLOCAL;++j)
                for (int iR = 0;iR < nR;++iR)
                {
                    int iat1 = ucell.iwt2iat[i];
                    int iat2 = ucell.iwt2iat[j];
                    const std::vector<int>& R = csr.getRCoordinate(iR);
                    const auto& matrix = csr.getMatrix(iR);
                    TC dR({ R[0], R[1], R[2] });
                    this->exx_ptr->Hexxs.at(is).at(iat1).at({ iat2, dR })(ucell.iwt2iw[i], ucell.iwt2iw[j]) = matrix(i, j);
                }
    }
    print_HR(this->exx_ptr->Hexxs[0], "Hexxs[0]_read", 1e-10);
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