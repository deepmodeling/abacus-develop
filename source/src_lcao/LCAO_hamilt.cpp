#include "../src_pw/global.h"
#include "LCAO_hamilt.h"
#include "build_st_pw.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "global_fp.h" // mohan add 2021-01-30
#include "dftu.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_xc/xc_functional.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"	//caoyu add 2021-07-26
#endif
#include "../module_base/timer.h"

LCAO_Hamilt::LCAO_Hamilt()
{
}

LCAO_Hamilt::~LCAO_Hamilt()
{
    if(GlobalV::test_deconstructor)
    {
        std::cout << " ~LCAO_Hamilt()" << std::endl;	
    }
}

//--------------------------------------------
// prepare grid network for Gint(grid integral)
//--------------------------------------------
void LCAO_Hamilt::grid_prepare(void)
{
    ModuleBase::TITLE("LCAO_Hamilt","grid_prepare");
    ModuleBase::timer::tick("LCAO_Hamilt","grid_prepare");

    if(GlobalV::GAMMA_ONLY_LOCAL)
    {   
        this->GG.prep_grid(GlobalC::bigpw->nbx, GlobalC::bigpw->nby, GlobalC::bigpw->nbzp, GlobalC::bigpw->nbzp_start, GlobalC::rhopw->nxyz);	

    }
    else // multiple k-points
    {
        // calculate the grid integration of 'Vl' matrix for l-points algorithms.
        this->GK.prep_grid(GlobalC::bigpw->nbx, GlobalC::bigpw->nby, GlobalC::bigpw->nbzp, GlobalC::bigpw->nbzp_start, GlobalC::rhopw->nxyz);
    }

    ModuleBase::timer::tick("LCAO_Hamilt","grid_prepare");
    return;
}

void LCAO_Hamilt::set_R_range_sparse()
{
    int R_minX = int(GlobalC::GridD.getD_minX());
    int R_minY = int(GlobalC::GridD.getD_minY());
    int R_minZ = int(GlobalC::GridD.getD_minZ());

    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    for(int ix = 0; ix < R_x; ix++)
    {
        for(int iy = 0; iy < R_y; iy++)
        {
            for(int iz = 0; iz < R_z; iz++)
            {
                Abfs::Vector3_Order<int> temp_R(ix+R_minX, iy+R_minY, iz+R_minZ);
                this->LM->all_R_coor.insert(temp_R);
            }
        }
    }

    return;
}

void LCAO_Hamilt::calculate_STN_R_sparse(const int &current_spin, const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_STN_R_sparse");

    int index = 0;
    ModuleBase::Vector3<double> dtau, tau1, tau2;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;

    double temp_value_double;
    std::complex<double> temp_value_complex;

    for(int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
            Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int start = GlobalC::ucell.itiaiw2iwt(T1,I1,0);

            for(int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
                const int I2 = GlobalC::GridD.getNatom(ad);
                Atom* atom2 = &GlobalC::ucell.atoms[T2];

                tau2 = GlobalC::GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                bool adj = false;

                if(distance < rcut) adj = true;
                else if(distance >= rcut)
                {
                    for(int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
                    {
                        const int T0 = GlobalC::GridD.getType(ad0);

                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                        double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                        double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        if( distance1 < rcut1 && distance2 < rcut2 )
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if(adj)
                {
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);

                    Abfs::Vector3_Order<int> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z);

                    for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = this->LM->ParaV->trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = this->LM->ParaV->trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            if(GlobalV::NSPIN!=4)
                            {
                                if (current_spin == 0)
                                {
                                    temp_value_double = this->LM->SlocR[index];
                                    if (std::abs(temp_value_double) > sparse_threshold)
                                    {
                                        this->LM->SR_sparse[dR][iw1_all][iw2_all] = temp_value_double;
                                    }
                                }

                                temp_value_double = this->LM->Hloc_fixedR[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    this->LM->HR_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                            }
                            else
                            {
                                temp_value_complex = this->LM->SlocR_soc[index];
                                if(std::abs(temp_value_complex) > sparse_threshold)
                                {
                                    this->LM->SR_soc_sparse[dR][iw1_all][iw2_all] = temp_value_complex;
                                }

                                temp_value_complex = this->LM->Hloc_fixedR_soc[index];
                                if(std::abs(temp_value_complex) > sparse_threshold)
                                {
                                    this->LM->HR_soc_sparse[dR][iw1_all][iw2_all] = temp_value_complex;
                                }
                            }

                            ++index;
                        }
                    }
                }
            }
        }
    }

    return;
}


void LCAO_Hamilt::calculate_STN_R_sparse_for_S(const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_STN_R_sparse_for_S");

    int index = 0;
    ModuleBase::Vector3<double> dtau, tau1, tau2;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;

    double temp_value_double;
    std::complex<double> temp_value_complex;

    for(int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
            Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int start = GlobalC::ucell.itiaiw2iwt(T1,I1,0);

            for(int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
                const int I2 = GlobalC::GridD.getNatom(ad);
                Atom* atom2 = &GlobalC::ucell.atoms[T2];

                tau2 = GlobalC::GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                bool adj = false;

                if(distance < rcut) adj = true;
                else if(distance >= rcut)
                {
                    for(int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
                    {
                        const int T0 = GlobalC::GridD.getType(ad0);

                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;

                        double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                        double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                        double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        if( distance1 < rcut1 && distance2 < rcut2 )
                        {
                            adj = true;
                            break;
                        }
                    }
                }

                if(adj)
                {
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);

                    Abfs::Vector3_Order<int> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z);

                    for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
                    {
                        const int iw1_all = start + ii;
                        const int mu = this->LM->ParaV->trace_loc_row[iw1_all];

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = this->LM->ParaV->trace_loc_col[iw2_all];

                            if(nu<0)continue;

                            if(GlobalV::NSPIN!=4)
                            {
                                temp_value_double = this->LM->SlocR[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    this->LM->SR_sparse[dR][iw1_all][iw2_all] = temp_value_double;
                                }
                            }
                            else
                            {
                                temp_value_complex = this->LM->SlocR_soc[index];
                                if(std::abs(temp_value_complex) > sparse_threshold)
                                {
                                    this->LM->SR_soc_sparse[dR][iw1_all][iw2_all] = temp_value_complex;
                                }
                            }

                            ++index;
                        }
                    }
                }
            }
        }
    }

    return;
}

void LCAO_Hamilt::calculate_HSR_sparse(const int &current_spin, const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_HSR_sparse");

    set_R_range_sparse();

    calculate_STN_R_sparse(current_spin, sparse_threshold);

    GK.cal_vlocal_R_sparseMatrix(current_spin, sparse_threshold, this->LM);

    if (GlobalV::dft_plus_u)
    {
        if (GlobalV::NSPIN != 4)
        {
            calculat_HR_dftu_sparse(current_spin, sparse_threshold);
        }
        else
        {
            calculat_HR_dftu_soc_sparse(current_spin, sparse_threshold);
        }
    }

#ifdef __EXX
#ifdef __MPI
    if( GlobalC::exx_info.info_global.cal_exx )
    {
        //calculate_HR_exx_sparse(current_spin, sparse_threshold);
        if(GlobalV::GAMMA_ONLY_LOCAL)
            this->calculate_HR_exx_sparse(current_spin, sparse_threshold, GlobalC::exx_lri_double.Hexxs);
        else
            this->calculate_HR_exx_sparse(current_spin, sparse_threshold, GlobalC::exx_lri_complex.Hexxs);
    }
#endif // __MPI
#endif // __EXX

    clear_zero_elements(current_spin, sparse_threshold);
}

void LCAO_Hamilt::calculate_SR_sparse(const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_SR_sparse");
    set_R_range_sparse();
    calculate_STN_R_sparse_for_S(sparse_threshold);
}

void LCAO_Hamilt::calculat_HR_dftu_sparse(const int &current_spin, const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculat_HR_dftu_sparse");
    ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_sparse");

    int total_R_num = this->LM->all_R_coor.size();
    int *nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    int count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        auto iter = this->LM->SR_sparse.find(R_coor);
        if (iter != this->LM->SR_sparse.end())
        {
            for (auto &row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }
        count++;
    }

    Parallel_Reduce::reduce_int_all(nonzero_num, total_R_num);

    double *HR_tmp = new double[this->LM->ParaV->nloc];
    double *SR_tmp = new double[this->LM->ParaV->nloc];

    int ir;
    int ic;
    int iic;
    auto &temp_HR_sparse = this->LM->HR_sparse[current_spin];

    count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        if (nonzero_num[count] != 0)
        {
            ModuleBase::GlobalFunc::ZEROS(HR_tmp, this->LM->ParaV->nloc);
            ModuleBase::GlobalFunc::ZEROS(SR_tmp, this->LM->ParaV->nloc);

            auto iter = this->LM->SR_sparse.find(R_coor);
            if (iter != this->LM->SR_sparse.end())
            {
                for (auto &row_loop : iter->second)
                {
                    ir = this->LM->ParaV->trace_loc_row[row_loop.first];
                    for (auto &col_loop : row_loop.second)
                    {
                        ic = this->LM->ParaV->trace_loc_col[col_loop.first];
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * this->LM->ParaV->nrow;
                        }
                        else
                        {
                            iic = ir * this->LM->ParaV->ncol + ic;
                        }
                        SR_tmp[iic] = col_loop.second;
                    }
                }
            }

            GlobalC::dftu.cal_eff_pot_mat_R_double(current_spin, SR_tmp, HR_tmp);

            for (int i = 0; i < GlobalV::NLOCAL; ++i)
            {
                ir = this->LM->ParaV->trace_loc_row[i];
                if (ir >= 0)
                {
                    for (int j = 0; j < GlobalV::NLOCAL; ++j)
                    {
                        ic = this->LM->ParaV->trace_loc_col[j];
                        if (ic >= 0)
                        {
                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                            {
                                iic = ir + ic * this->LM->ParaV->nrow;
                            }
                            else
                            {
                                iic = ir * this->LM->ParaV->ncol + ic;
                            }

                            if (std::abs(HR_tmp[iic]) > sparse_threshold)
                            {
                                double &value = temp_HR_sparse[R_coor][i][j];
                                value += HR_tmp[iic];
                                if (std::abs(value) <= sparse_threshold)
                                {
                                    temp_HR_sparse[R_coor][i].erase(j);
                                }
                            }
                        }
                    }
                }
            }

        }
        
        count++;
    }

    delete[] nonzero_num;
    delete[] HR_tmp;
    delete[] SR_tmp;
    nonzero_num = nullptr;
    HR_tmp = nullptr;
    SR_tmp = nullptr;

    ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_sparse");

}

void LCAO_Hamilt::calculat_HR_dftu_soc_sparse(const int &current_spin, const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");
    ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");

    int total_R_num = this->LM->all_R_coor.size();
    int *nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    int count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        auto iter = this->LM->SR_soc_sparse.find(R_coor);
        if (iter != this->LM->SR_soc_sparse.end())
        {
            for (auto &row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }
        count++;
    }

    Parallel_Reduce::reduce_int_all(nonzero_num, total_R_num);

    std::complex<double> *HR_soc_tmp = new std::complex<double>[this->LM->ParaV->nloc];
    std::complex<double> *SR_soc_tmp = new std::complex<double>[this->LM->ParaV->nloc];

    int ir;
    int ic;
    int iic;

    count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        if (nonzero_num[count] != 0)
        {
            ModuleBase::GlobalFunc::ZEROS(HR_soc_tmp, this->LM->ParaV->nloc);
            ModuleBase::GlobalFunc::ZEROS(SR_soc_tmp, this->LM->ParaV->nloc);

            auto iter = this->LM->SR_soc_sparse.find(R_coor);
            if (iter != this->LM->SR_soc_sparse.end())
            {
                for (auto &row_loop : iter->second)
                {
                    ir = this->LM->ParaV->trace_loc_row[row_loop.first];
                    for (auto &col_loop : row_loop.second)
                    {
                        ic = this->LM->ParaV->trace_loc_col[col_loop.first];
                        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            iic = ir + ic * this->LM->ParaV->nrow;
                        }
                        else
                        {
                            iic = ir * this->LM->ParaV->ncol + ic;
                        }
                        SR_soc_tmp[iic] = col_loop.second;
                    }
                }
            }

            GlobalC::dftu.cal_eff_pot_mat_R_complex_double(current_spin, SR_soc_tmp, HR_soc_tmp);

            for (int i = 0; i < GlobalV::NLOCAL; ++i)
            {
                ir = this->LM->ParaV->trace_loc_row[i];
                if (ir >= 0)
                {
                    for (int j = 0; j < GlobalV::NLOCAL; ++j)
                    {
                        ic = this->LM->ParaV->trace_loc_col[j];
                        if (ic >= 0)
                        {
                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                            {
                                iic = ir + ic * this->LM->ParaV->nrow;
                            }
                            else
                            {
                                iic = ir * this->LM->ParaV->ncol + ic;
                            }

                            if (std::abs(HR_soc_tmp[iic]) > sparse_threshold)
                            {
                                std::complex<double> &value = this->LM->HR_soc_sparse[R_coor][i][j];
                                value += HR_soc_tmp[iic];
                                if (std::abs(value) <= sparse_threshold)
                                {
                                    this->LM->HR_soc_sparse[R_coor][i].erase(j);
                                }
                            }
                        }
                    }
                }
            }

        }
        
        count++;
    }

    delete[] nonzero_num;
    delete[] HR_soc_tmp;
    delete[] SR_soc_tmp;
    nonzero_num = nullptr;
    HR_soc_tmp = nullptr;
    SR_soc_tmp = nullptr;

    ModuleBase::timer::tick("LCAO_Hamilt","calculat_HR_dftu_soc_sparse");

}

#ifdef __EXX
// Peize Lin add 2021.11.16
void LCAO_Hamilt::calculate_HR_exx_sparse(const int &current_spin, const double &sparse_threshold)
{
	ModuleBase::TITLE("LCAO_Hamilt","calculate_HR_exx_sparse");
	ModuleBase::timer::tick("LCAO_Hamilt","calculate_HR_exx_sparse");	

	const Abfs::Vector3_Order<int> Rs_period(GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]);
	if(Rs_period.x<=0 || Rs_period.y<=0 || Rs_period.z<=0)
		throw std::invalid_argument("Rs_period = ("+ModuleBase::GlobalFunc::TO_STRING(Rs_period.x)+","+ModuleBase::GlobalFunc::TO_STRING(Rs_period.y)+","+ModuleBase::GlobalFunc::TO_STRING(Rs_period.z)+").\n"
			+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	const std::vector<Abfs::Vector3_Order<int>> Rs = Abfs::get_Born_von_Karmen_boxes( Rs_period );

	const int ik_begin = (GlobalV::NSPIN==2) ? (current_spin*GlobalC::kv.nks/2) : 0;
	const int ik_end = (GlobalV::NSPIN==2) ? ((current_spin+1)*GlobalC::kv.nks/2) : GlobalC::kv.nks;
	const double frac = (GlobalV::NSPIN==1) ? 0.5 : 1.0;                        // Peize Lin add 2022.07.09

	for(const Abfs::Vector3_Order<int> &R : Rs)
	{
		ModuleBase::matrix HexxR;
		for(int ik=ik_begin; ik<ik_end; ++ik)
		{
			ModuleBase::matrix HexxR_tmp;
			if(GlobalV::GAMMA_ONLY_LOCAL)
				HexxR_tmp = GlobalC::exx_info.info_global.hybrid_alpha
					* GlobalC::exx_lcao.Hexx_para.HK_Gamma_m2D[ik]
					* (GlobalC::kv.wk[ik] * frac);
			else
				HexxR_tmp = GlobalC::exx_info.info_global.hybrid_alpha
					* (GlobalC::exx_lcao.Hexx_para.HK_K_m2D[ik]
					* std::exp( ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (GlobalC::kv.kvec_c[ik] * (R*GlobalC::ucell.latvec)) )).real()
					* (GlobalC::kv.wk[ik] * frac);

			if(HexxR.c)
				HexxR += HexxR_tmp;
			else
				HexxR = std::move(HexxR_tmp);
		}

		for(int iwt1_local=0; iwt1_local<HexxR.nr; ++iwt1_local)
		{
			const int iwt1_global = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER()
				? this->LM->ParaV->MatrixInfo.col_set[iwt1_local]
				: this->LM->ParaV->MatrixInfo.row_set[iwt1_local];
			for(int iwt2_local=0; iwt2_local<HexxR.nc; ++iwt2_local)
			{
			    const int iwt2_global = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER()
					? this->LM->ParaV->MatrixInfo.row_set[iwt2_local]
					: this->LM->ParaV->MatrixInfo.col_set[iwt2_local];
				if(std::abs(HexxR(iwt1_local,iwt2_local)) > sparse_threshold)
				{
					if(GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
					{
						auto &HR_sparse_ptr = this->LM->HR_sparse[current_spin][R][iwt1_global];
						auto &HR_sparse = HR_sparse_ptr[iwt2_global];
						HR_sparse += HexxR(iwt1_local,iwt2_local);
						if(std::abs(HR_sparse) < sparse_threshold)
							HR_sparse_ptr.erase(iwt2_global);
					}
					else
					{
						auto &HR_sparse_ptr = this->LM->HR_soc_sparse[R][iwt1_global];
						auto &HR_sparse = HR_sparse_ptr[iwt2_global];
						HR_sparse += HexxR(iwt1_local,iwt2_local);
						if(std::abs(HR_sparse) < sparse_threshold)
							HR_sparse_ptr.erase(iwt2_global);
					}
				}
			}
		}
	}

    // In the future it should be changed to mpi communication, since some Hexx(R) of R in Rs may be zeros
    this->LM->all_R_coor.insert(Rs.begin(),Rs.end());
    
	ModuleBase::timer::tick("LCAO_Hamilt","calculate_HR_exx_sparse");	
}
#endif // __EXX

// in case there are elements smaller than the threshold
void LCAO_Hamilt::clear_zero_elements(const int &current_spin, const double &sparse_threshold)
{
    if(GlobalV::NSPIN != 4)
    {
        for (auto &R_loop : this->LM->HR_sparse[current_spin])
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second; 
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
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

        for (auto &R_loop : this->LM->SR_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second; 
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
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
        for (auto &R_loop : this->LM->HR_soc_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second; 
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
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

        for (auto &R_loop : this->LM->SR_soc_sparse)
        {
            for (auto &row_loop : R_loop.second)
            {
                auto &col_map = row_loop.second; 
                auto iter = col_map.begin();
                while (iter != col_map.end())
                {
                    if (std::abs(iter->second) <= sparse_threshold)
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
}

void LCAO_Hamilt::destroy_all_HSR_sparse(void)
{
	this->LM->destroy_HS_R_sparse();
}
