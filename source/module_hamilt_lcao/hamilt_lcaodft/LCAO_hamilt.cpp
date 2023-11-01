#include "LCAO_hamilt.h"

#include "module_base/parallel_reduce.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"	//caoyu add 2021-07-26
#endif
#include "module_base/timer.h"

#ifdef __EXX
#include "LCAO_hamilt.hpp"
#endif

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
void LCAO_Hamilt::grid_prepare(const Grid_Technique& gt, const ModulePW::PW_Basis& rhopw, const ModulePW::PW_Basis_Big& bigpw)
{
    ModuleBase::TITLE("LCAO_Hamilt","grid_prepare");
    ModuleBase::timer::tick("LCAO_Hamilt","grid_prepare");

    if(GlobalV::GAMMA_ONLY_LOCAL)
    {
        this->GG.prep_grid(gt, bigpw.nbx, bigpw.nby, bigpw.nbzp, bigpw.nbzp_start,
            rhopw.nxyz, bigpw.bx, bigpw.by, bigpw.bz, bigpw.bxyz, bigpw.nbxx,
            rhopw.ny, rhopw.nplane, rhopw.startz_current);

    }
    else // multiple k-points
    {
        // calculate the grid integration of 'Vl' matrix for l-points algorithms.
        this->GK.prep_grid(gt, bigpw.nbx, bigpw.nby, bigpw.nbzp, bigpw.nbzp_start,
            rhopw.nxyz, bigpw.bx, bigpw.by, bigpw.bz, bigpw.bxyz, bigpw.nbxx,
            rhopw.ny, rhopw.nplane, rhopw.startz_current);
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
                        const int mu = this->LM->ParaV->global2local_row(iw1_all);

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = this->LM->ParaV->global2local_col(iw2_all);

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

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
void LCAO_Hamilt::calculate_HContainer_sparse_d(const int &current_spin, const double &sparse_threshold, const hamilt::HContainer<double>& hR, std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& target)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_HContainer_sparse_d");

    const Parallel_Orbitals* paraV = this->LM->ParaV;
    auto row_indexes = paraV->get_indexes_row();
    auto col_indexes = paraV->get_indexes_col();
    for(int iap=0;iap<hR.size_atom_pairs();++iap)
    {
        int atom_i = hR.get_atom_pair(iap).get_atom_i();
        int atom_j = hR.get_atom_pair(iap).get_atom_j();
        int start_i = paraV->atom_begin_row[atom_i];
        int start_j = paraV->atom_begin_col[atom_j];
        int row_size = paraV->get_row_size(atom_i);
        int col_size = paraV->get_col_size(atom_j);
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
                    if(std::abs(value_tmp)>sparse_threshold)
                    {
                        target[dR][mu][nu] = value_tmp;
                    }
                }
            }
        }
    }

    return;
}

void LCAO_Hamilt::calculate_HContainer_sparse_cd(const int &current_spin, const double &sparse_threshold, const hamilt::HContainer<std::complex<double>>& hR, std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>>& target)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_HContainer_sparse_cd");

    const Parallel_Orbitals* paraV = this->LM->ParaV;
    auto row_indexes = paraV->get_indexes_row();
    auto col_indexes = paraV->get_indexes_col();
    for(int iap=0;iap<hR.size_atom_pairs();++iap)
    {
        int atom_i = hR.get_atom_pair(iap).get_atom_i();
        int atom_j = hR.get_atom_pair(iap).get_atom_j();
        int start_i = paraV->atom_begin_row[atom_i];
        int start_j = paraV->atom_begin_col[atom_j];
        int row_size = paraV->get_row_size(atom_i);
        int col_size = paraV->get_col_size(atom_j);
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
                    if(std::abs(value_tmp)>sparse_threshold)
                    {
                        target[dR][mu][nu] = value_tmp;
                    }
                }
            }
        }
    }

    return;
}

void LCAO_Hamilt::calculate_HSR_sparse(const int &current_spin, const double &sparse_threshold, const int (&nmp)[3], hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_HSR_sparse");

    set_R_range_sparse();

    //calculate_STN_R_sparse(current_spin, sparse_threshold);
    if(GlobalV::NSPIN!=4)
    {
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham);
        this->calculate_HContainer_sparse_d(current_spin, sparse_threshold, *(p_ham_lcao->getHR()), this->LM->HR_sparse[current_spin]);
        this->calculate_HContainer_sparse_d(current_spin, sparse_threshold, *(p_ham_lcao->getSR()), this->LM->SR_sparse);
    }
    else
    {
        hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_ham_lcao = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham);
        this->calculate_HContainer_sparse_cd(current_spin, sparse_threshold, *(p_ham_lcao->getHR()), this->LM->HR_soc_sparse);
        this->calculate_HContainer_sparse_cd(current_spin, sparse_threshold, *(p_ham_lcao->getSR()), this->LM->SR_soc_sparse);
    }

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
        if(GlobalC::exx_info.info_ri.real_number)
            this->calculate_HR_exx_sparse(current_spin, sparse_threshold, nmp, *this->LM->Hexxd);
        else
            this->calculate_HR_exx_sparse(current_spin, sparse_threshold, nmp, *this->LM->Hexxc);
    }
#endif // __MPI
#endif // __EXX

    clear_zero_elements(current_spin, sparse_threshold);
}

void LCAO_Hamilt::calculate_dH_sparse(const int &current_spin, const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_dH_sparse");

    set_R_range_sparse();

    const int nnr = this->LM->ParaV->nnr;
    this->LM->DHloc_fixedR_x = new double[nnr];
    this->LM->DHloc_fixedR_y = new double[nnr];
    this->LM->DHloc_fixedR_z = new double[nnr];

    ModuleBase::GlobalFunc::ZEROS(this->LM->DHloc_fixedR_x, this->LM->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(this->LM->DHloc_fixedR_y, this->LM->ParaV->nloc);
    ModuleBase::GlobalFunc::ZEROS(this->LM->DHloc_fixedR_z, this->LM->ParaV->nloc);
    // calculate dT=<phi|kin|dphi> in LCAO
    // calculate T + VNL(P1) in LCAO basis
    if(GlobalV::CAL_STRESS)
	{
        GlobalV::CAL_STRESS = false;
        this->genH.build_ST_new('T', true, GlobalC::ucell, this->LM->Hloc_fixedR.data());
        GlobalV::CAL_STRESS = true;
    }
    else
    {
        this->genH.build_ST_new('T', true, GlobalC::ucell, this->LM->Hloc_fixedR.data());
    }
    this->genH.build_Nonlocal_mu_new (this->LM->Hloc_fixed.data(), true);
    
    calculate_dSTN_R_sparse(current_spin, sparse_threshold);

    delete[] this->LM->DHloc_fixedR_x;
    delete[] this->LM->DHloc_fixedR_y;
    delete[] this->LM->DHloc_fixedR_z;

    GK.cal_dvlocal_R_sparseMatrix(current_spin, sparse_threshold, this->LM);
}

void LCAO_Hamilt::calculate_STN_R_sparse_for_T(const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_STN_R_sparse_for_T");

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
                        const int mu = this->LM->ParaV->global2local_row(iw1_all);

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = this->LM->ParaV->global2local_col(iw2_all);

                            if(nu<0)continue;

                            if(GlobalV::NSPIN!=4)
                            {
                                temp_value_double = this->LM->Hloc_fixedR[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    this->LM->TR_sparse[dR][iw1_all][iw2_all] = temp_value_double;
                                }
                            }
                            else
                            {
                                temp_value_complex = this->LM->Hloc_fixedR_soc[index];
                                if(std::abs(temp_value_complex) > sparse_threshold)
                                {
                                    this->LM->TR_soc_sparse[dR][iw1_all][iw2_all] = temp_value_complex;
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

void LCAO_Hamilt::calculate_SR_sparse(const double &sparse_threshold, hamilt::Hamilt<std::complex<double>>* p_ham)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_SR_sparse");
    set_R_range_sparse();
    //calculate_STN_R_sparse(current_spin, sparse_threshold);
    if(GlobalV::NSPIN!=4)
    {
        hamilt::HamiltLCAO<std::complex<double>, double>* p_ham_lcao = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(p_ham);
        this->calculate_HContainer_sparse_d(0, sparse_threshold, *(p_ham_lcao->getSR()), this->LM->SR_sparse);
    }
    else
    {
        hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>* p_ham_lcao = dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(p_ham);
        this->calculate_HContainer_sparse_cd(0, sparse_threshold, *(p_ham_lcao->getSR()), this->LM->SR_soc_sparse);
    }
}

void LCAO_Hamilt::calculate_TR_sparse(const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_TR_sparse");
    
    //need to rebuild T(R)
    this->LM->Hloc_fixedR.resize(this->LM->ParaV->nnr);
    this->LM->zeros_HSR('T');
    this->genH.build_ST_new('T', 0, GlobalC::ucell, this->LM->Hloc_fixedR.data());
    set_R_range_sparse();
    calculate_STN_R_sparse_for_T(sparse_threshold);
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

    Parallel_Reduce::reduce_all(nonzero_num, total_R_num);

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
                    ir = this->LM->ParaV->global2local_row(row_loop.first);
                    for (auto &col_loop : row_loop.second)
                    {
                        ic = this->LM->ParaV->global2local_col(col_loop.first);
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
                ir = this->LM->ParaV->global2local_row(i);
                if (ir >= 0)
                {
                    for (int j = 0; j < GlobalV::NLOCAL; ++j)
                    {
                        ic = this->LM->ParaV->global2local_col(j);
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

    Parallel_Reduce::reduce_all(nonzero_num, total_R_num);

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
                    ir = this->LM->ParaV->global2local_row(row_loop.first);
                    for (auto &col_loop : row_loop.second)
                    {
                        ic = this->LM->ParaV->global2local_col(col_loop.first);
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
                ir = this->LM->ParaV->global2local_row(i);
                if (ir >= 0)
                {
                    for (int j = 0; j < GlobalV::NLOCAL; ++j)
                    {
                        ic = this->LM->ParaV->global2local_col(j);
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

void LCAO_Hamilt::calculate_dSTN_R_sparse(const int &current_spin, const double &sparse_threshold)
{
    ModuleBase::TITLE("LCAO_Hamilt","calculate_dSTN_R_sparse");

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
                        const int mu = this->LM->ParaV->global2local_row(iw1_all);

                        if(mu<0)continue;

                        for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
                        {
                            int iw2_all = start2 + jj;
                            const int nu = this->LM->ParaV->global2local_col(iw2_all);

                            if(nu<0)continue;

                            if(GlobalV::NSPIN!=4)
                            {
                                temp_value_double = this->LM->DHloc_fixedR_x[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    this->LM->dHRx_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                                temp_value_double = this->LM->DHloc_fixedR_y[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    this->LM->dHRy_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                                temp_value_double = this->LM->DHloc_fixedR_z[index];
                                if (std::abs(temp_value_double) > sparse_threshold)
                                {
                                    this->LM->dHRz_sparse[current_spin][dR][iw1_all][iw2_all] = temp_value_double;
                                }
                            }
                            else
                            {
                                ModuleBase::WARNING_QUIT("calculate_dSTN_R_sparse","soc not supported!");
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

void LCAO_Hamilt::destroy_all_HSR_sparse(void)
{
	this->LM->destroy_HS_R_sparse();
}

void LCAO_Hamilt::destroy_TR_sparse(void)
{
	this->LM->destroy_T_R_sparse();
}

void LCAO_Hamilt::destroy_dH_R_sparse(void)
{
	this->LM->destroy_dH_R_sparse();
}
