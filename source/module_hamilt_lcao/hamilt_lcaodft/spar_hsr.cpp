#include "spar_hsr.h"
#include "spar_dh.h"
#include "spar_u.h"
#include "spar_exx.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

void sparse_format::cal_HSR(
        const Parallel_Orbitals &pv,
		LCAO_Matrix &lm,
		Grid_Driver &grid,
		const int &current_spin, 
		const double &sparse_thr, 
		const int (&nmp)[3], 
		hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("sparse_format","cal_HSR");

    sparse_format::set_R_range(lm.all_R_coor, grid);

    const int nspin = GlobalV::NSPIN;

    //cal_STN_R_sparse(current_spin, sparse_thr);
    if(nspin==1 || nspin==2)
    {
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao = 
        dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham);

		sparse_format::cal_HContainer_d(
                pv,
                current_spin, 
				sparse_thr, 
				*(p_ham_lcao->getHR()), 
				lm.HR_sparse[current_spin]);

		sparse_format::cal_HContainer_d(
				pv,
				current_spin, 
				sparse_thr, 
				*(p_ham_lcao->getSR()), 
				lm.SR_sparse);
    }
	else if(nspin==4)
	{
		hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_ham_lcao = 
			dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham);

		sparse_format::cal_HContainer_cd(
				pv,
				current_spin, 
				sparse_thr, 
				*(p_ham_lcao->getHR()), 
				lm.HR_soc_sparse);

		sparse_format::cal_HContainer_cd(
				pv,
				current_spin, 
				sparse_thr, 
				*(p_ham_lcao->getSR()), 
				lm.SR_soc_sparse);
	}
	else
	{
        ModuleBase::WARNING_QUIT("cal_HSR","check the value of nspin.");       
    }

    // only old DFT+U method need to cal extra contribution to HR
    if (GlobalV::dft_plus_u == 2)
    {
        if(nspin==1 || nspin==2)
        {
			cal_HR_dftu(
					pv,
					lm.all_R_coor,
					lm.SR_sparse,
					lm.HR_sparse,
					current_spin, 
					sparse_thr);
		}
		else if(nspin==4)
		{
			cal_HR_dftu_soc(
					pv,
					lm.all_R_coor,
					lm.SR_soc_sparse,
					lm.HR_soc_sparse,
					current_spin, 
					sparse_thr);
        }
		else
		{
			ModuleBase::WARNING_QUIT("cal_HSR","check the value of nspin.");       
		}
	}

#ifdef __EXX
#ifdef __MPI
    // if EXX is considered
    if( GlobalC::exx_info.info_global.cal_exx )
    {
		if(GlobalC::exx_info.info_ri.real_number)
		{
			sparse_format::cal_HR_exx(lm, current_spin, sparse_thr, nmp, *lm.Hexxd);
		}
		else
		{
			sparse_format::cal_HR_exx(lm, current_spin, sparse_thr, nmp, *lm.Hexxc);
		}
	}
#endif // __MPI
#endif // __EXX

    sparse_format::clear_zero_elements(lm, current_spin, sparse_thr);

    return;
}


void sparse_format::cal_HContainer_d(
        const Parallel_Orbitals &pv,
		const int &current_spin, 
		const double &sparse_thr, 
		const hamilt::HContainer<double>& hR, 
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& target)
{
    ModuleBase::TITLE("sparse_format","cal_HContainer_d");

    auto row_indexes = pv.get_indexes_row();
    auto col_indexes = pv.get_indexes_col();
    for(int iap=0;iap<hR.size_atom_pairs();++iap)
    {
        int atom_i = hR.get_atom_pair(iap).get_atom_i();
        int atom_j = hR.get_atom_pair(iap).get_atom_j();
        int start_i = pv.atom_begin_row[atom_i];
        int start_j = pv.atom_begin_col[atom_j];
        int row_size = pv.get_row_size(atom_i);
        int col_size = pv.get_col_size(atom_j);
        for(int iR=0;iR<hR.get_atom_pair(iap).get_R_size();++iR)
        {
            auto& matrix = hR.get_atom_pair(iap).get_HR_values(iR);
            int* r_index = hR.get_atom_pair(iap).get_R_index(iR);
            Abfs::Vector3_Order<int> dR(r_index[0], r_index[1], r_index[2]);
            for(int i=0;i<row_size;++i)
            {
                int mu = row_indexes[start_i+i];
                for(int j=0;j<col_size;++j)
                {
                    int nu = col_indexes[start_j+j];
                    const auto& value_tmp = matrix.get_value(i,j);
                    if(std::abs(value_tmp)>sparse_thr)
                    {
                        target[dR][mu][nu] = value_tmp;
                    }
                }
            }
        }
    }

    return;
}

void sparse_format::cal_HContainer_cd(
        const Parallel_Orbitals &pv,
		const int &current_spin, 
		const double &sparse_thr, 
		const hamilt::HContainer<std::complex<double>>& hR, 
		std::map<Abfs::Vector3_Order<int>, 
		std::map<size_t, std::map<size_t, std::complex<double>>>>& target)
{
    ModuleBase::TITLE("sparse_format","cal_HContainer_cd");

    auto row_indexes = pv.get_indexes_row();
    auto col_indexes = pv.get_indexes_col();
    for(int iap=0;iap<hR.size_atom_pairs();++iap)
    {
        int atom_i = hR.get_atom_pair(iap).get_atom_i();
        int atom_j = hR.get_atom_pair(iap).get_atom_j();
        int start_i = pv.atom_begin_row[atom_i];
        int start_j = pv.atom_begin_col[atom_j];
        int row_size = pv.get_row_size(atom_i);
        int col_size = pv.get_col_size(atom_j);
        for(int iR=0;iR<hR.get_atom_pair(iap).get_R_size();++iR)
        {
            auto& matrix = hR.get_atom_pair(iap).get_HR_values(iR);
            int* r_index = hR.get_atom_pair(iap).get_R_index(iR);
            Abfs::Vector3_Order<int> dR(r_index[0], r_index[1], r_index[2]);
            for(int i=0;i<row_size;++i)
            {
                int mu = row_indexes[start_i+i];
                for(int j=0;j<col_size;++j)
                {
                    int nu = col_indexes[start_j+j];
                    const auto& value_tmp = matrix.get_value(i,j);
                    if(std::abs(value_tmp)>sparse_thr)
                    {
                        target[dR][mu][nu] = value_tmp;
                    }
                }
            }
        }
    }

    return;
}


// in case there are elements smaller than the threshold
void sparse_format::clear_zero_elements(
        LCAO_Matrix &lm,
		const int &current_spin, 
		const double &sparse_thr)
{
    ModuleBase::TITLE("sparse_format","clear_zero_elements");

    if(GlobalV::NSPIN != 4)
    {
        for (auto &R_loop : lm.HR_sparse[current_spin])
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_thr)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }

        for (auto &R_loop : lm.SR_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_thr)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }
            }
        }

    }
    else
    {
        for (auto &R_loop : lm.HR_soc_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_thr)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }// end while iter
            }// end row loop
        }// end R loop

        for (auto &R_loop : lm.SR_soc_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second;
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_thr)
                    {
                        col_map.erase(iter++);
                    }
                    else
                    {
                        iter++;
                    }
                }// end while iter
            }// end row_loop
        }// end R_loop
    }

    return;
}
