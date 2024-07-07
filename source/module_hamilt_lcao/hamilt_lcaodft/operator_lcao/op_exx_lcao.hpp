#ifndef OPEXXLCAO_HPP
#define OPEXXLCAO_HPP

#ifdef __EXX

#include "op_exx_lcao.h"
#include "module_ri/RI_2D_Comm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

namespace hamilt
{

template <typename TK, typename TR>
OperatorEXX<OperatorLCAO<TK, TR>>::OperatorEXX(HS_Matrix_K<TK>* hsk_in,
	LCAO_Matrix* LM_in,
	hamilt::HContainer<TR>* hR_in,
	const K_Vectors& kv_in,
	std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd_in,
	std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc_in,
    Add_Hexx_Type add_hexx_type_in,
	int* two_level_step_in,
	const bool restart_in)
	: OperatorLCAO<TK, TR>(hsk_in, kv_in.kvec_d, hR_in),
    kv(kv_in),
    Hexxd(Hexxd_in),
    Hexxc(Hexxc_in),
    add_hexx_type(add_hexx_type_in),
    two_level_step(two_level_step_in),
    restart(restart_in)
{
    this->cal_type = calculation_type::lcao_exx;
    // if k points has no shift, use cell_nearest to reduce the memory cost
    this->use_cell_nearest = (ModuleBase::Vector3<double>(std::fmod(this->kv.get_koffset(0), 1.0),
        std::fmod(this->kv.get_koffset(1), 1.0), std::fmod(this->kv.get_koffset(2), 1.0)).norm() < 1e-10);

    if (this->add_hexx_type == Add_Hexx_Type::R)
    {
        // set cell_nearest
        std::map<int, std::array<double, 3>> atoms_pos;
        for (int iat = 0; iat < GlobalC::ucell.nat; ++iat) {
            atoms_pos[iat] = RI_Util::Vector3_to_array3(
                GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat]]
                .tau[GlobalC::ucell.iat2ia[iat]]);
        }
        const std::array<std::array<double, 3>, 3> latvec
            = { RI_Util::Vector3_to_array3(GlobalC::ucell.a1),
               RI_Util::Vector3_to_array3(GlobalC::ucell.a2),
               RI_Util::Vector3_to_array3(GlobalC::ucell.a3) };
        const std::array<int, 3> Rs_period = { this->kv.nmp[0], this->kv.nmp[1], this->kv.nmp[2] };
        this->cell_nearest.init(atoms_pos, latvec, Rs_period);

        // reallocate hR if needed (Hexxd temp)
        auto Rs = RI_Util::get_Born_von_Karmen_cells(Rs_period);
        bool need_allocate = false;
        for (int iat0 = 0;iat0 < GlobalC::ucell.nat;++iat0)
        {
            for (int iat1 = 0;iat1 < GlobalC::ucell.nat;++iat1)
            {
                for (auto& cell : Rs)
                {
                    const Abfs::Vector3_Order<int>& R = RI_Util::array3_to_Vector3(
                        (this->use_cell_nearest ?
                            cell_nearest.get_cell_nearest_discrete(iat0, iat1, cell)
                            : cell));
                    hamilt::BaseMatrix<TR>* HlocR = this->hR->find_matrix(iat0, iat1, R.x, R.y, R.z);
                    if (HlocR == nullptr)
                    { // add R to HContainer
                        need_allocate = true;
                        hamilt::AtomPair<TR> tmp(iat0, iat1, R.x, R.y, R.z, this->LM->ParaV);
                        this->hR->insert_pair(tmp);
                    }
                }
            }
        }
        if (need_allocate) this->hR->allocate(nullptr, true);
    }

    if (this->restart)
	{///  Now only Hexx depends on DM, so we can directly read Hexx to reduce the computational cost.
	/// If other operators depends on DM, we can also read DM and then calculate the operators to save the memory to store operator terms.
        assert(this->two_level_step != nullptr);

        if (this->add_hexx_type == Add_Hexx_Type::k)
        {
            /// read in Hexx(k)
            if (std::is_same<TK, double>::value)
            {
                this->LM->Hexxd_k_load.resize(this->kv.get_nks());
                for (int ik = 0; ik < this->kv.get_nks(); ik++)
                {
                    this->LM->Hexxd_k_load[ik].resize(this->LM->ParaV->get_local_size(), 0.0);
                    this->restart = GlobalC::restart.load_disk(
                        "Hexx", ik,
                        this->LM->ParaV->get_local_size(), this->LM->Hexxd_k_load[ik].data(), false);
                    if (!this->restart) break;
                }
            }
            else
            {
                this->LM->Hexxc_k_load.resize(this->kv.get_nks());
                for (int ik = 0; ik < this->kv.get_nks(); ik++)
                {
                    this->LM->Hexxc_k_load[ik].resize(this->LM->ParaV->get_local_size(), 0.0);
                    this->restart = GlobalC::restart.load_disk(
                        "Hexx", ik,
                        this->LM->ParaV->get_local_size(), this->LM->Hexxc_k_load[ik].data(), false);
                    if (!this->restart) break;
                }
            }
        }
        else if (this->add_hexx_type == Add_Hexx_Type::R)
        {
            // refactor IO-csr functions
        }

		if (!this->restart)
			std::cout << "WARNING: Hexx not found, restart from the non-exx loop." << std::endl
			          << "If the loaded charge density is EXX-solved, this may lead to poor convergence." << std::endl;
		GlobalC::restart.info_load.load_H_finish = this->restart;
	}
}

template<typename TK, typename TR>
void OperatorEXX<OperatorLCAO<TK, TR>>::contributeHR()
{
    // Peize Lin add 2016-12-03
    if (GlobalV::CALCULATION != "nscf" && this->two_level_step != nullptr && *this->two_level_step == 0 && !this->restart) return;  //in the non-exx loop, do nothing 
    if (XC_Functional::get_func_type() == 4 || XC_Functional::get_func_type() == 5)
    {
        if (this->restart && this->two_level_step != nullptr)
        {
            if (*this->two_level_step == 0)
            {
                this->add_loaded_HexxR();
                return;
            }
            else // clear loaded Hexx and release memory
            {
                this->clear_loaded_HexxR();
            }
        }
        // cal H(k) from H(R) normally
        if (GlobalC::exx_info.info_ri.real_number)
            RI_2D_Comm::add_HexxR(
                this->current_spin,
                GlobalC::exx_info.info_global.hybrid_alpha,
                this->Hexxd == nullptr ? *this->LM->Hexxd : *this->Hexxd,
                *this->LM->ParaV,
                GlobalV::NPOL,
                *this->hR,
                this->use_cell_nearest ? &this->cell_nearest : nullptr);
        else
            RI_2D_Comm::add_HexxR(
                this->current_spin,
                GlobalC::exx_info.info_global.hybrid_alpha,
                this->Hexxc == nullptr ? *this->LM->Hexxc : *this->Hexxc,
                *this->LM->ParaV,
                GlobalV::NPOL,
                *this->hR,
                this->use_cell_nearest ? &this->cell_nearest : nullptr);
    }
    if (GlobalV::NSPIN == 2) this->current_spin = 1 - this->current_spin;
}

template<typename TK, typename TR>
void OperatorEXX<OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    // Peize Lin add 2016-12-03
    if (GlobalV::CALCULATION != "nscf" && this->two_level_step != nullptr && *this->two_level_step == 0 && !this->restart) return;  //in the non-exx loop, do nothing 
    if (XC_Functional::get_func_type() == 4 || XC_Functional::get_func_type() == 5)
    {
        if (this->restart && this->two_level_step != nullptr)
        {
            if (*this->two_level_step == 0)
            {
                this->add_loaded_Hexx(ik);
                return;
            }
            else // clear loaded Hexx and release memory
            {
                if (this->LM->Hexxd_k_load.size() > 0)
                {
                    this->LM->Hexxd_k_load.clear();
                    this->LM->Hexxd_k_load.shrink_to_fit();
                }
                else if (this->LM->Hexxc_k_load.size() > 0)
                {
                    this->LM->Hexxc_k_load.clear();
                    this->LM->Hexxc_k_load.shrink_to_fit();
                }
            }
        }
        // cal H(k) from H(R) normally

        if (GlobalC::exx_info.info_ri.real_number)
            RI_2D_Comm::add_Hexx(
                this->kv,
                ik,
                GlobalC::exx_info.info_global.hybrid_alpha,
                this->Hexxd == nullptr ? *this->LM->Hexxd : *this->Hexxd,
                *this->LM->ParaV,
                this->hsk->get_hk());
        else
            RI_2D_Comm::add_Hexx(
                this->kv,
                ik,
                GlobalC::exx_info.info_global.hybrid_alpha,
                this->Hexxc == nullptr ? *this->LM->Hexxc : *this->Hexxc,
                *this->LM->ParaV,
                this->hsk->get_hk());
    }
}

} // namespace hamilt
#endif // __EXX
#endif // OPEXXLCAO_HPP