//wenfei 2022-1-11
//This file contains subroutines for calculating pdm,
//which is defind as sum_mu,nu rho_mu,nu (<chi_mu|alpha><alpha|chi_nu>);
//as well as gdmx, which is the gradient of pdm, defined as
//sum_mu,nu rho_mu,nu d/dX(<chi_mu|alpha><alpha|chi_nu>)

//It also contains subroutines for printing pdm and gdmx
//for checking purpose

//There are 6 subroutines in this file:
//1. cal_projected_DM, which is used for calculating pdm for gamma point calculation
//2. cal_projected_DM_k, counterpart of 1, for multi-k
//3. check_projected_dm, which prints pdm to descriptor.dat

//4. cal_gdmx, calculating gdmx (and optionally gdm_epsl for stress) for gamma point
//5. cal_gdmx_k, counterpart of 3, for multi-k
//6. check_gdmx, which prints gdmx to a series of .dat files

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "module_base/vector3.h"
#include "module_base/timer.h"
#include "module_base/constants.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_base/libm/libm.h"

//this subroutine performs the calculation of projected density matrices
//pdm_m,m'=\sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
void LCAO_Deepks::cal_projected_DM(const elecstate::DensityMatrix<double, double>* dm, 
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_projected_DM");

    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    if (GlobalV::init_chg == "file" && !this->init_pdm) //for DeePKS NSCF calculation 
    {
        std::ifstream ifs("pdm.dat");
        if (!ifs)
        {
            ModuleBase::WARNING_QUIT("LCAO_Deepks::cal_projected_DM", "Can not find the file pdm.dat . Please do DeePKS SCF calculation first.");
        }
        for(int inl=0;inl<this->inlmax;inl++)
        {
            for(int ind=0;ind<pdm_size;ind++)
            {
                double c;
			    ifs >> c;
                pdm[inl][ind] = c;
            }
        }
        this->init_pdm = true;
        return;
    }

    //if(dm.size() == 0 || dm[0].size() == 0)
    //{
    //    return;
    //}
    ModuleBase::timer::tick("LCAO_Deepks","cal_projected_DM");

    for(int inl=0;inl<inlmax;inl++)
    {
        ModuleBase::GlobalFunc::ZEROS(pdm[inl],pdm_size);
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    int nrow = this->pv->nrow;
    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            //trace alpha orbital
            std::vector<int> trace_alpha_row;
            std::vector<int> trace_alpha_col;
            int ib=0;
            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
            {
                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                {
                    const int inl = this->inl_index[T0](I0, L0, N0);
                    const int nm = 2*L0+1;
            
                    for (int m1=0; m1<nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                    {
                        for (int m2=0; m2<nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                        {
                            trace_alpha_row.push_back(ib+m1);
                            trace_alpha_col.push_back(ib+m2);
                        }
                    }
                    ib+=nm;
                }
            }
            const int trace_alpha_size = trace_alpha_row.size();

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1,I1);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut(); 

                auto row_indexes = pv->get_indexes_row(ibt1);
                const int row_size = row_indexes.size();
                if(row_size == 0) continue;

                // no possible to unexist key
                std::vector<double> s_1t(trace_alpha_size * row_size);
                std::vector<double> g_1dmt(trace_alpha_size * row_size, 0.0);
                for(int irow=0;irow<row_size;irow++)
                {
                    const double* row_ptr = this->nlm_save[iat][ad1][row_indexes[irow]][0].data();
                    
                    for(int i=0;i<trace_alpha_size;i++)
                    {
                        s_1t[i * row_size + irow] = row_ptr[trace_alpha_row[i]];
                    }
                }

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
                    const int ibt2 = ucell.itia2iat(T2,I2);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

                    auto col_indexes = pv->get_indexes_col(ibt2);
                    const int col_size = col_indexes.size();
                    if(col_size == 0) continue;

                    std::vector<double> s_2t(trace_alpha_size * col_size);
                    // no possible to unexist key
                    for(int icol=0;icol<col_size;icol++)
                    {
                        const double* col_ptr = this->nlm_save[iat][ad2][col_indexes[icol]][0].data();
                        for(int i=0;i<trace_alpha_size;i++)
                        {
                            s_2t[i * col_size + icol] = col_ptr[trace_alpha_col[i]];
                        }
                    }
                    // prepare DM_gamma from DMR
                    std::vector<double> dm_array(row_size*col_size, 0.0);
                    const double* dm_current;
                    for(int is=0;is<dm->get_DMR_vector().size();is++)
                    {
                        dm_current = dm->get_DMR_vector()[is]->find_matrix(ibt1, ibt2, 0, 0, 0)->get_pointer();
                        for(int idm=0;idm<row_size*col_size;idm++)
                        {
                            dm_array[idm] += dm_current[idm];
                        }
                    }
                    dm_current = dm_array.data();
                    //dgemm for s_2t and dm_current to get g_1dmt
                    constexpr char transa='T', transb='N';
                    const double gemm_alpha = 1.0, gemm_beta = 1.0;
                    dgemm_(
                        &transa, &transb, 
                        &row_size, 
                        &trace_alpha_size, 
                        &col_size, 
                        &gemm_alpha, 
                        dm_current, 
                        &col_size,
                        s_2t.data(),  
                        &col_size,     
                        &gemm_beta,      
                        g_1dmt.data(),    
                        &row_size);
				}//ad2
                // do dot of g_1dmt and s_1t to get orbital_pdm_shell
                int ib=0, index=0, inc=1;
                for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                {
                    for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                    {
                        const int inl = this->inl_index[T0](I0, L0, N0);
                        const int nm = 2*L0+1;
                
                        for (int m1=0; m1<nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                        {
                            for (int m2=0; m2<nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                            {
                                int ind = m1*nm + m2;
                                pdm[inl][ind] += 
                                    ddot_(&row_size, g_1dmt.data()+index*row_size, &inc, s_1t.data()+index*row_size, &inc);
                                index++;
                            }
                        }
                        ib+=nm;
                    }
                }
			}//ad1
        }//I0
    }//T0

#ifdef __MPI
    allsum_deepks(this->inlmax,pdm_size,this->pdm);
#endif
    ModuleBase::timer::tick("LCAO_Deepks","cal_projected_DM"); 
    return;
}

void LCAO_Deepks::cal_projected_DM_k(const elecstate::DensityMatrix<std::complex<double>, double>* dm,
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD,
    const int nks,
    const std::vector<ModuleBase::Vector3<double>> &kvec_d)
{
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);

    if (GlobalV::init_chg == "file" && !this->init_pdm) //for DeePKS NSCF calculation 
    {
        std::ifstream ifs("pdm.dat");
        if (!ifs)
        {
            ModuleBase::WARNING_QUIT("LCAO_Deepks::cal_projected_DM_k","Can not find the file pdm.dat . Please do DeePKS SCF calculation first.");
        }
        for(int inl=0;inl<this->inlmax;inl++)
        {
            for(int ind=0;ind<pdm_size;ind++)
            {
                double c;
			    ifs >> c;
                pdm[inl][ind] = c;
            }
        }
        this->init_pdm = true;
        return;
    }

    //check for skipping
    //if(dm.size() == 0 || dm[0].size() == 0)
    //{
    //    return;
    //}
    ModuleBase::timer::tick("LCAO_Deepks","cal_projected_DM_k");

    for(int inl=0;inl<inlmax;inl++)
    {
        ModuleBase::GlobalFunc::ZEROS(pdm[inl],pdm_size);
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            //trace alpha orbital
            std::vector<int> trace_alpha_row;
            std::vector<int> trace_alpha_col;
            int ib=0;
            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
            {
                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                {
                    const int inl = this->inl_index[T0](I0, L0, N0);
                    const int nm = 2*L0+1;
            
                    for (int m1=0; m1<nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                    {
                        for (int m2=0; m2<nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                        {
                            trace_alpha_row.push_back(ib+m1);
                            trace_alpha_col.push_back(ib+m2);
                        }
                    }
                    ib+=nm;
                }
            }
            const int trace_alpha_size = trace_alpha_row.size();

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1,I1);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut();

                ModuleBase::Vector3<double> dR1(GridD.getBox(ad1).x, GridD.getBox(ad1).y, GridD.getBox(ad1).z); 

                auto row_indexes = pv->get_indexes_row(ibt1);
                const int row_size = row_indexes.size();
                if(row_size == 0) continue;

                key_tuple key_1(ibt1,dR1.x,dR1.y,dR1.z);
                if(this->nlm_save_k[iat].find(key_1) == this->nlm_save_k[iat].end()) continue;
                std::vector<double> s_1t(trace_alpha_size * row_size);
                std::vector<double> g_1dmt(trace_alpha_size * row_size, 0.0);
                for(int irow=0;irow<row_size;irow++)
                {
                    const double* row_ptr = this->nlm_save_k[iat][key_1][row_indexes[irow]][0].data();
                    for(int i=0;i<trace_alpha_size;i++)
                    {
                        s_1t[i * row_size + irow] = row_ptr[trace_alpha_row[i]];
                    }
                }

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
                    const int ibt2 = ucell.itia2iat(T2,I2);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    ModuleBase::Vector3<double> dR2(GridD.getBox(ad2).x, GridD.getBox(ad2).y, GridD.getBox(ad2).z);
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

                    auto col_indexes = pv->get_indexes_col(ibt2);
                    const int col_size = col_indexes.size();
                    if(col_size == 0) continue;

                    std::vector<double> s_2t(trace_alpha_size * col_size);
                    key_tuple key_2(ibt2,dR2.x,dR2.y,dR2.z);
                    if(this->nlm_save_k[iat].find(key_2) == this->nlm_save_k[iat].end()) continue;
                    for(int icol=0;icol<col_size;icol++)
                    {
                        const double* col_ptr = this->nlm_save_k[iat][key_2][col_indexes[icol]][0].data();
                        for(int i=0;i<trace_alpha_size;i++)
                        {
                            s_2t[i * col_size + icol] = col_ptr[trace_alpha_col[i]];
                        }
                    }
                    // prepare DM_gamma from DMR
                    std::vector<double> dm_array(row_size*col_size, 0.0);
                    const double* dm_current;
                    for(int is=0;is<dm->get_DMR_vector().size();is++)
                    {
                        auto tmp_matrix = dm->get_DMR_vector()[is]->find_matrix(ibt1, ibt2, (dR2-dR1).x, (dR2-dR1).y, (dR2-dR1).z);
                        if(tmp_matrix == nullptr)
                        {
                            dm_current = nullptr;
                            break;
                        }
                        dm_current = tmp_matrix->get_pointer();
                        for(int idm=0;idm<row_size*col_size;idm++)
                        {
                            dm_array[idm] += dm_current[idm];
                        }
                    }
                    if(dm_current == nullptr) continue;
                    dm_current = dm_array.data();
                    //dgemm for s_2t and dm_current to get g_1dmt
                    constexpr char transa='T', transb='N';
                    const double gemm_alpha = 1.0, gemm_beta = 1.0;
                    dgemm_(
                        &transa, &transb, 
                        &row_size, 
                        &trace_alpha_size, 
                        &col_size, 
                        &gemm_alpha, 
                        dm_current, 
                        &col_size,
                        s_2t.data(),  
                        &col_size,     
                        &gemm_beta,      
                        g_1dmt.data(),    
                        &row_size);
				}//ad2
                // do dot of g_1dmt and s_1t to get orbital_pdm_shell
                int ib=0, index=0, inc=1;
                for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                {
                    for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                    {
                        const int inl = this->inl_index[T0](I0, L0, N0);
                        const int nm = 2*L0+1;
                
                        for (int m1=0; m1<nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                        {
                            for (int m2=0; m2<nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                            {
                                int ind = m1*nm + m2;
                                pdm[inl][ind] += 
                                    ddot_(&row_size, g_1dmt.data()+index*row_size, &inc, s_1t.data()+index*row_size, &inc);
                                index++;
                            }
                        }
                        ib+=nm;
                    }
                }
			}//ad1
        }//I0
    }//T0

#ifdef __MPI
    allsum_deepks(this->inlmax,pdm_size,this->pdm);
#endif
    ModuleBase::timer::tick("LCAO_Deepks","cal_projected_DM_k");
    return;
    
}

void LCAO_Deepks::check_projected_dm(void)
{
    std::ofstream ofs("pdm.dat");
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    ofs<<std::setprecision(10);
    for(int inl=0;inl<inlmax;inl++)
    {
        for(int ind=0;ind<pdm_size;ind++)
        {
            ofs << pdm[inl][ind] << " ";
        }
        ofs << std::endl;
    }
}

//this subroutine calculates the gradient of projected density matrices
//gdmx_m,m = d/dX sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
//if stress label is enabled, the gradient of PDM wrt strain tensor will 
//be calculated: 
//gdm_epsl = d/d\epsilon_{ab} * 
//           sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
void LCAO_Deepks::cal_gdmx(const std::vector<double>& dm,
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD,
    const bool isstress)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_gdmx");
    ModuleBase::timer::tick("LCAO_Deepks","cal_gdmx");
    //get DS_alpha_mu and S_nu_beta

    int size = (2 * lmaxd + 1) * (2 * lmaxd + 1);
    int nrow = this->pv->nrow;
    for (int iat = 0;iat < ucell.nat;iat++)
    {
        for (int inl = 0;inl < inlmax;inl++)
        {
            ModuleBase::GlobalFunc::ZEROS(gdmx[iat][inl], size);
            ModuleBase::GlobalFunc::ZEROS(gdmy[iat][inl], size);
            ModuleBase::GlobalFunc::ZEROS(gdmz[iat][inl], size);
        }
    }

    if (isstress)
    {
        for (int ipol = 0;ipol < 6;ipol++)
        {
            for (int inl = 0;inl < inlmax;inl++)
            {
                ModuleBase::GlobalFunc::ZEROS(gdm_epsl[ipol][inl], size);
            }
        }
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();

    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);//on which alpha is located
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1,I1); //on which chi_mu is located
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut(); 

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    const int ibt2 = ucell.itia2iat(T2,I2);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

                    double r0[3];
                    double r1[3];
                    if(isstress)
                    {
                        r1[0] = ( tau1.x - tau0.x) ;
                        r1[1] = ( tau1.y - tau0.y) ;
                        r1[2] = ( tau1.z - tau0.z) ;
                        r0[0] = ( tau2.x - tau0.x) ;
                        r0[1] = ( tau2.y - tau0.y) ;
                        r0[2] = ( tau2.z - tau0.z) ;
                    }
                    auto row_indexes = pv->get_indexes_row(ibt1);
                    auto col_indexes = pv->get_indexes_col(ibt2);
                    if(row_indexes.size() * col_indexes.size() == 0) continue;

                    hamilt::AtomPair<double> dm_pair(ibt1, ibt2, 0, 0, 0, pv);
                    dm_pair.allocate(nullptr, 1);
                    if(ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                    {
                        dm_pair.add_from_matrix(dm.data(), pv->get_row_size(), 1.0, 1);
                    }
                    else
                    {
                        dm_pair.add_from_matrix(dm.data(), pv->get_col_size(), 1.0, 0);
                    }
                    const double* dm_current = dm_pair.get_pointer();

					for (int iw1=0; iw1<row_indexes.size(); ++iw1)
                    {
                        for (int iw2=0; iw2<col_indexes.size(); ++iw2)
                        {
                            std::vector<double> nlm1 = this->nlm_save[iat][ad1][row_indexes[iw1]][0];
                            std::vector<std::vector<double>> nlm2 = this->nlm_save[iat][ad2][col_indexes[iw2]];

                            assert(nlm1.size()==nlm2[0].size());

                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < nm;++m1)
                                    {
                                        for (int m2 = 0; m2 <nm; ++m2)
                                        {
                                            //(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'>
                                            gdmx[iat][inl][m1*nm+m2] += nlm2[1][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);
                                            gdmy[iat][inl][m1*nm+m2] += nlm2[2][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);
                                            gdmz[iat][inl][m1*nm+m2] += nlm2[3][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);

                                            //(<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                            gdmx[iat][inl][m2*nm+m1] += nlm2[1][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);
                                            gdmy[iat][inl][m2*nm+m1] += nlm2[2][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);
                                            gdmz[iat][inl][m2*nm+m1] += nlm2[3][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);                                            

                                            //(<chi_mu|d/dX alpha_m>)<chi_nu|alpha_m'> = -(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'>
                                            gdmx[ibt2][inl][m1*nm+m2] -= nlm2[1][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);                                               
                                            gdmy[ibt2][inl][m1*nm+m2] -= nlm2[2][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);                                               
                                            gdmz[ibt2][inl][m1*nm+m2] -= nlm2[3][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);

                                            //(<chi_nu|d/dX alpha_m'>)<chi_mu|alpha_m> = -(<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                            gdmx[ibt2][inl][m2*nm+m1] -= nlm2[1][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);                                               
                                            gdmy[ibt2][inl][m2*nm+m1] -= nlm2[2][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);                                               
                                            gdmz[ibt2][inl][m2*nm+m1] -= nlm2[3][ib+m2] * nlm1[ib+m1] * *dm_current; //dm(iw1_local,iw2_local);                                            

                                            if (isstress)
                                            {
                                                int mm = 0;
                                                for(int ipol=0;ipol<3;ipol++)
                                                {
                                                    for(int jpol=ipol;jpol<3;jpol++)
                                                    {
                                                        //gdm_epsl[mm][inl][m2*nm+m1] += ucell.lat0 * dm(iw1_local, iw2_local) * (nlm2[jpol+1][ib+m2] * nlm1[ib+m1] * r0[ipol]);
                                                        gdm_epsl[mm][inl][m2*nm+m1] += ucell.lat0 * *dm_current * (nlm2[jpol+1][ib+m2] * nlm1[ib+m1] * r0[ipol]);
                                                        mm++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());
                            if  (isstress)
                            {
                                nlm1 = this->nlm_save[iat][ad2][col_indexes[iw2]][0];
                                nlm2 = this->nlm_save[iat][ad1][row_indexes[iw1]];

                                assert(nlm1.size()==nlm2[0].size());  
                                int ib=0;
                                for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                                {
                                    for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                    {
                                        const int inl = this->inl_index[T0](I0, L0, N0);
                                        const int nm = 2*L0+1;
                                        for (int m1 = 0;m1 < nm; ++m1)
                                        {
                                            for (int m2 = 0; m2 < nm; ++m2)
                                            {
                                                int mm = 0;
                                                for(int ipol=0;ipol<3;ipol++)
                                                {
                                                    for(int jpol=ipol;jpol<3;jpol++)
                                                    {
                                                        gdm_epsl[mm][inl][m2*nm+m1]  += ucell.lat0 * *dm_current * (nlm1[ib+m1] * nlm2[jpol+1][ib+m2] * r1[ipol]);
                                                        mm++;
                                                    }
                                                }
                                            }
                                        }
                                        ib+=nm;
                                    }
                                }
                            }
                            dm_current++;
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }//I0
    }//T0

#ifdef __MPI
    for(int iat=0;iat<ucell.nat;iat++)
    {
        allsum_deepks(this->inlmax,size,this->gdmx[iat]);
        allsum_deepks(this->inlmax,size,this->gdmy[iat]);
        allsum_deepks(this->inlmax,size,this->gdmz[iat]);
    }
    if (isstress)
    {
        for(int ipol=0;ipol<6;ipol++)
        {
            allsum_deepks(this->inlmax,size,this->gdm_epsl[ipol]);
        }
    }
#endif
    ModuleBase::timer::tick("LCAO_Deepks","cal_gdmx");
    return;
}

void LCAO_Deepks::cal_gdmx_k(const std::vector<std::vector<std::complex<double>>>& dm,
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver& GridD,
    const int nks,
    const std::vector<ModuleBase::Vector3<double>> &kvec_d,
    const bool isstress)
{
    ModuleBase::TITLE("LCAO_Deepks", "cal_gdmx_k");
    ModuleBase::timer::tick("LCAO_Deepks","cal_gdmx_k");

    int size = (2 * lmaxd + 1) * (2 * lmaxd + 1);

    for (int iat = 0;iat < ucell.nat;iat++)
    {
        for (int inl = 0;inl < inlmax;inl++)
        {
            ModuleBase::GlobalFunc::ZEROS(gdmx[iat][inl], size);
            ModuleBase::GlobalFunc::ZEROS(gdmy[iat][inl], size);
            ModuleBase::GlobalFunc::ZEROS(gdmz[iat][inl], size);
        }
    }

    if (isstress)
    {
        for (int ipol = 0;ipol < 6;ipol++)
        {
            for (int inl = 0;inl < inlmax;inl++)
            {
                ModuleBase::GlobalFunc::ZEROS(gdm_epsl[ipol][inl], size);
            }
        }
    }

    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    int nrow = this->pv->nrow;
    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);//on which alpha is located
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int ibt1 = ucell.itia2iat(T1,I1); //on which chi_mu is located
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = orb.Phi[T1].getRcut();

                ModuleBase::Vector3<double> dR1(GridD.getBox(ad1).x, GridD.getBox(ad1).y, GridD.getBox(ad1).z); 

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                    const int ibt2 = ucell.itia2iat(T2,I2);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    ModuleBase::Vector3<double> dR2(GridD.getBox(ad2).x, GridD.getBox(ad2).y, GridD.getBox(ad2).z);
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

                    double r0[3];
                    double r1[3];
                    if(isstress)
                    {
                        r1[0] = ( tau1.x - tau0.x) ;
                        r1[1] = ( tau1.y - tau0.y) ;
                        r1[2] = ( tau1.z - tau0.z) ;
                        r0[0] = ( tau2.x - tau0.x) ;
                        r0[1] = ( tau2.y - tau0.y) ;
                        r0[2] = ( tau2.z - tau0.z) ;
                    }

                    auto row_indexes = pv->get_indexes_row(ibt1);
                    auto col_indexes = pv->get_indexes_col(ibt2);
                    if(row_indexes.size() * col_indexes.size() == 0) continue;

                    hamilt::AtomPair<double> dm_pair(ibt1, ibt2, (dR2-dR1).x, (dR2-dR1).y, (dR2-dR1).z, pv);
                    dm_pair.allocate(nullptr, 1);
                    for(int ik=0;ik<nks;ik++)
                    {
                        const double arg = - (kvec_d[ik] * (dR2-dR1) ) * ModuleBase::TWO_PI;
                        double sinp, cosp;
                        ModuleBase::libm::sincos(arg, &sinp, &cosp);
                        const std::complex<double> kphase = std::complex<double>(cosp, sinp);
                        if(ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
                        {
                            dm_pair.add_from_matrix(dm[ik].data(), pv->get_row_size(), kphase, 1);
                        }
                        else
                        {
                            dm_pair.add_from_matrix(dm[ik].data(), pv->get_col_size(), kphase, 0);
                        }
                    }
                    const double* dm_current = dm_pair.get_pointer();

                    key_tuple key_1(ibt1,dR1.x,dR1.y,dR1.z);
                    key_tuple key_2(ibt2,dR2.x,dR2.y,dR2.z);
					for (int iw1l = 0; iw1l < row_indexes.size(); ++iw1l)
                    {
                        std::vector<double> nlm1 = this->nlm_save_k[iat][key_1][row_indexes[iw1l]][0];
                        for (int iw2l = 0; iw2l < col_indexes.size(); ++iw2l)
                        {
                            std::vector<std::vector<double>> nlm2 = this->nlm_save_k[iat][key_2][col_indexes[iw2l]];

                            assert(nlm1.size()==nlm2[0].size());

                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                    {
                                        for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                        {
                                            //(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'>
                                            gdmx[iat][inl][m1*nm+m2] += nlm2[1][ib+m2] * nlm1[ib+m1] * *dm_current;
                                            gdmy[iat][inl][m1*nm+m2] += nlm2[2][ib+m2] * nlm1[ib+m1] * *dm_current;
                                            gdmz[iat][inl][m1*nm+m2] += nlm2[3][ib+m2] * nlm1[ib+m1] * *dm_current;

                                            //(<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                            gdmx[iat][inl][m2*nm+m1] += nlm2[1][ib+m2] * nlm1[ib+m1] * *dm_current;
                                            gdmy[iat][inl][m2*nm+m1] += nlm2[2][ib+m2] * nlm1[ib+m1] * *dm_current;
                                            gdmz[iat][inl][m2*nm+m1] += nlm2[3][ib+m2] * nlm1[ib+m1] * *dm_current;                                            

                                            //(<chi_mu|d/dX alpha_m>)<chi_nu|alpha_m'> = -(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'>
                                            gdmx[ibt2][inl][m1*nm+m2] -= nlm2[1][ib+m2] * nlm1[ib+m1] * *dm_current;                                               
                                            gdmy[ibt2][inl][m1*nm+m2] -= nlm2[2][ib+m2] * nlm1[ib+m1] * *dm_current;                                               
                                            gdmz[ibt2][inl][m1*nm+m2] -= nlm2[3][ib+m2] * nlm1[ib+m1] * *dm_current;

                                            //(<chi_nu|d/dX alpha_m'>)<chi_mu|alpha_m> = -(<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                            gdmx[ibt2][inl][m2*nm+m1] -= nlm2[1][ib+m2] * nlm1[ib+m1] * *dm_current;                                               
                                            gdmy[ibt2][inl][m2*nm+m1] -= nlm2[2][ib+m2] * nlm1[ib+m1] * *dm_current;                                               
                                            gdmz[ibt2][inl][m2*nm+m1] -= nlm2[3][ib+m2] * nlm1[ib+m1] * *dm_current;     


                                            if (isstress)
                                            {
                                                int mm = 0;
                                                for(int ipol=0;ipol<3;ipol++)
                                                {
                                                    for(int jpol=ipol;jpol<3;jpol++)
                                                    {
                                                        gdm_epsl[mm][inl][m2*nm+m1] += ucell.lat0 * *dm_current * (nlm2[jpol+1][ib+m2] * nlm1[ib+m1] * r0[ipol]);
                                                        mm++;
                                                    }
                                                }
                                            }

                                        }
                                    }
                                    ib+=nm;
                                }
                            }
                            assert(ib==nlm1.size());

                            if  (isstress)
                            {
                                nlm1 = this->nlm_save_k[iat][key_2][col_indexes[iw2l]][0];
                                nlm2 = this->nlm_save_k[iat][key_1][row_indexes[iw1l]];

                                assert(nlm1.size()==nlm2[0].size());  
                                int ib=0;
                                for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                                {
                                    for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                    {
                                        const int inl = this->inl_index[T0](I0, L0, N0);
                                        const int nm = 2*L0+1;
                                        for (int m1 = 0;m1 < nm; ++m1)
                                        {
                                            for (int m2 = 0; m2 < nm; ++m2)
                                            {
                                                int mm = 0;
                                                for(int ipol=0;ipol<3;ipol++)
                                                {
                                                    for(int jpol=ipol;jpol<3;jpol++)
                                                    {
                                                        gdm_epsl[mm][inl][m2*nm+m1]  += ucell.lat0 * *dm_current * (nlm1[ib+m1] * nlm2[jpol+1][ib+m2] * r1[ipol]);
                                                        mm++;
                                                    }
                                                }
                                            }
                                        }
                                        ib+=nm;
                                    }
                                }
                            }
                            dm_current++;
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }//I0
    }//T0

#ifdef __MPI
    for(int iat=0;iat<ucell.nat;iat++)
    {
        allsum_deepks(this->inlmax,size,this->gdmx[iat]);
        allsum_deepks(this->inlmax,size,this->gdmy[iat]);
        allsum_deepks(this->inlmax,size,this->gdmz[iat]);
    }
    if (isstress)
    {
        for(int ipol=0;ipol<6;ipol++)
        {
            allsum_deepks(this->inlmax,size,this->gdm_epsl[ipol]);
        }
    }
#endif
    ModuleBase::timer::tick("LCAO_Deepks","cal_gdmx_k");
    return;
}

void LCAO_Deepks::check_gdmx(const int nat)
{
    std::stringstream ss;
    std::ofstream ofs_x;
    std::ofstream ofs_y;
    std::ofstream ofs_z;

    ofs_x<<std::setprecision(10);
    ofs_y<<std::setprecision(10);
    ofs_z<<std::setprecision(10);

    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    for(int ia=0;ia<nat;ia++)
    {
        ss.str("");
        ss<<"gdmx_"<<ia<<".dat";
        ofs_x.open(ss.str().c_str());
        ss.str("");
        ss<<"gdmy_"<<ia<<".dat";
        ofs_y.open(ss.str().c_str());
        ss.str("");
        ss<<"gdmz_"<<ia<<".dat";
        ofs_z.open(ss.str().c_str());

        for(int inl=0;inl<inlmax;inl++)
        {
            for(int ind=0;ind<pdm_size;ind++)
            {
                ofs_x << gdmx[ia][inl][ind] << " ";
                ofs_y << gdmy[ia][inl][ind] << " ";
                ofs_z << gdmz[ia][inl][ind] << " ";
            }
            ofs_x << std::endl;
            ofs_y << std::endl;
            ofs_z << std::endl;
        }
        ofs_x.close();
        ofs_y.close();
        ofs_z.close();
    }
}

#endif