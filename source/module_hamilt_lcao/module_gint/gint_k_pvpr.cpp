#include "gint_k.h"
#include "grid_technique.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/libm/libm.h"
#include "module_base/memory.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void Gint_k::allocate_pvpR(void)
{
    ModuleBase::TITLE("Gint_k","allocate_pvpR");

    if(this->pvpR_alloc_flag)
    {
        return; //Liuxh add, 20181012
        ModuleBase::WARNING_QUIT("Gint_k::allocate_pvpR","pvpR has been allocated!");
    }

    //xiaohui modify 2015-05-30
    // the number of matrix element <phi_0 | V | phi_R> is nnrg.
    this->pvpR_reduced = new double*[GlobalV::NSPIN];
    for(int is =0;is<GlobalV::NSPIN;is++)
    {
        this->pvpR_reduced[is] = new double[this->gridt->nnrg];	
        ModuleBase::GlobalFunc::ZEROS( pvpR_reduced[is], this->gridt->nnrg);
    }

    ModuleBase::Memory::record("pvpR_reduced", sizeof(double) * this->gridt->nnrg * GlobalV::NSPIN);

    this->pvpR_alloc_flag = true;
    return;
}

void Gint_k::destroy_pvpR(void)
{
    ModuleBase::TITLE("Gint_k","destroy_pvpR");
    
    if(!pvpR_alloc_flag)
    {
        return;
    }
    
    for(int is =0;is<GlobalV::NSPIN;is++) delete[] pvpR_reduced[is];
    delete[] pvpR_reduced;

    this->pvpR_alloc_flag = false;
    return;
}

// folding the matrix for 'ik' k-point.
// H(k)=\sum{R} H(R)exp(ikR) 
void Gint_k::folding_vl_k(const int &ik, 
                        LCAO_Matrix *LM, 
                        const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("Gint_k","folding_vl_k");
    ModuleBase::timer::tick("Gint_k","folding_vl_k");

    if(!pvpR_alloc_flag)
    {
        ModuleBase::WARNING_QUIT("Gint_k::destroy_pvpR","pvpR hasnot been allocated yet!");
    }

    //####################### EXPLAIN #################################
    // 1. what is gridt->lgd ?
    // gridt->lgd is the number of orbitals in each processor according
    // to the division of real space FFT grid.
    // 
    // 2. why the folding of vlocal is different from folding of 
    // < phi_0i | T+Vnl | phi_Rj > ?
    // Because the (i,j) is different for T+Vnl and Vlocal
    // The first part is due to 2D division of H and S matrix,
    // The second part is due to real space division. 
    // 
    // here we construct a temporary matrix to store the
    // matrix element < phi_0 | Vlocal | phi_R >
    //#################################################################

    int lgd = this->gridt->lgd;
    std::complex<double>** pvp = new std::complex<double>*[lgd];
    std::complex<double>* pvp_base = new std::complex<double>[lgd * lgd];
    for(int i=0; i<lgd; i++)
    {
        pvp[i] = pvp_base + i * lgd;
    }

    std::complex<double>*** pvp_nc;
    std::complex<double>* pvp_nc_base;
    if(GlobalV::NSPIN==4)
    {
        pvp_nc_base = new std::complex<double>[4 * lgd * lgd];
        pvp_nc=new std::complex<double>**[4];
        for(int spin=0;spin<4;spin++)
        {
            pvp_nc[spin] = new std::complex<double>*[lgd];
            for(int i=0; i<lgd; i++)
            {
                pvp_nc[spin][i] = pvp_nc_base + spin * lgd * lgd + i * lgd;
            }
        }
    }

    auto init_pvp = [&](int num_threads, int thread_id)
    {
        int beg, len;
        ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, lgd * lgd, 256, beg, len);
        ModuleBase::GlobalFunc::ZEROS(pvp_base + beg, len);
        if(GlobalV::NSPIN==4)
        {
            ModuleBase::GlobalFunc::ZEROS(pvp_nc_base + 4 * beg, 4 * len);
        }
    };
    ModuleBase::OMP_PARALLEL(init_pvp);

#ifdef _OPENMP
#pragma omp parallel
{
#endif
    ModuleBase::Vector3<double> tau1, dtau, dR;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
    {
        const int T1 = GlobalC::ucell.iat2it[iat];
        const int I1 = GlobalC::ucell.iat2ia[iat];
        {
            // atom in this grid piece.
            if(this->gridt->in_this_processor[iat])
            {
                Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

                // get the start positions of elements.
                const int DM_start = this->gridt->nlocstartg[iat];

                // get the coordinates of adjacent atoms.
                tau1 = GlobalC::ucell.atoms[T1].tau[I1];
                //GlobalC::GridD.Find_atom(tau1);	
                AdjacentAtomInfo adjs;
                GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1, &adjs);	
                // search for the adjacent atoms.
                int nad = 0;

                for (int ad = 0; ad < adjs.adj_num+1; ad++)
                {
                    // get iat2
                    const int T2 = adjs.ntype[ad];
                    const int I2 = adjs.natom[ad];
                    const int iat2 = GlobalC::ucell.itia2iat(T2, I2);


                    // adjacent atom is also on the grid.
                    if(this->gridt->in_this_processor[iat2])
                    {
                        Atom* atom2 = &GlobalC::ucell.atoms[T2];
                        dtau = adjs.adjacent_tau[ad] - tau1;
                        double distance = dtau.norm() * GlobalC::ucell.lat0;
                        double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                        // for the local part, only need to calculate <phi_i | phi_j> within range
                        // mohan note 2012-07-06
                        if(distance < rcut)
                        {
                            const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0); 

                            // calculate the distance between iat1 and iat2.
                            // ModuleBase::Vector3<double> dR = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                            dR.x = adjs.box[ad].x;
                            dR.y = adjs.box[ad].y;
                            dR.z = adjs.box[ad].z;

                            // calculate the phase factor exp(ikR).
                            const double arg = (kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                            double sinp, cosp;
                            ModuleBase::libm::sincos(arg, &sinp, &cosp);
                            const std::complex<double> phase = std::complex<double>(cosp, sinp);
                            int ixxx = DM_start + this->gridt->find_R2st[iat][nad];
                            
                            if(GlobalV::NSPIN!=4)
                            {
                                for(int iw=0; iw<atom1->nw; iw++)
                                {
                                    std::complex<double> *vij = pvp[this->gridt->trace_lo[start1+iw]];
                                    const int* iw2_lo = &this->gridt->trace_lo[start2];
                                    // get the <phi | V | phi>(R) Hamiltonian.
                                    const double *vijR = &pvpR_reduced[0][ixxx];
                                    for(int iw2 = 0; iw2<atom2->nw; ++iw2)
                                    {
                                        vij[iw2_lo[iw2]] += vijR[iw2] * phase; 
                                        //if(((start1+iw == 238 ) && ( start2+iw2 == 1089 )))
                                        //    GlobalV::ofs_running<<__FILE__<<__LINE__<<" "<<iat<<" "<<iat2<<" "<<ixxx<<" "<<start1+iw<<" "<<start2+iw2<<" "<<vijR[iw2]<<" "<<iw2_lo[iw2]<<" "<<vij[iw2_lo[iw2]]<<std::endl;
                                    }
                                    ixxx += atom2->nw;
                                }
                            }
                            else
                            {
                                for(int iw=0; iw<atom1->nw; iw++)
                                {
                                    int iw2_lo = this->gridt->trace_lo[start2]/GlobalV::NPOL;
                                    for(int spin = 0;spin<4;spin++) 
                                    {
                                        auto vij = pvp_nc[spin][this->gridt->trace_lo[start1]/GlobalV::NPOL + iw];
                                        auto vijR = &pvpR_reduced[spin][ixxx];
                                        auto vijs = &vij[iw2_lo];
                                        for(int iw2 = 0; iw2<atom2->nw; ++iw2)
                                        {
                                            vijs[iw2] += vijR[iw2] * phase; 
                                        }
                                    }
                                    ixxx += atom2->nw;
                                }
                            }
                            ++nad;
                        }// end distane<rcut
                    }
                }// end ad
            }
        }// end ia
    }// end it
#ifdef _OPENMP
}
#endif

    // Distribution of data.
    ModuleBase::timer::tick("Gint_k","Distri");
    const int nlocal = GlobalV::NLOCAL;
    std::vector<std::complex<double>> tmp(nlocal);
    const double sign_table[2] = {1.0, -1.0};
#ifdef _OPENMP
#pragma omp parallel
{
#endif
    //loop each row with index i, than loop each col with index j 
    for (int i=0; i<nlocal; i++)
    {
#ifdef _OPENMP
#pragma omp for
#endif
        for (int j=0; j<nlocal; j++)
        {
            tmp[j] = std::complex<double>(0.0, 0.0);
        }
        int i_flag = i & 1; // i % 2 == 0
        const int mug = this->gridt->trace_lo[i];
        // if the row element is on this processor.
        if (mug >= 0)
        {
            if(GlobalV::NSPIN!=4)
            {
#ifdef _OPENMP
#pragma omp for
#endif
                for (int j=0; j<nlocal; j++)
                {
                    const int nug = this->gridt->trace_lo[j];
                    // if the col element is on this processor.
                    if (nug >=0)
                    {
                        if (mug <= nug)
                        {
                            // pvp is symmetric, only half is calculated.
                            tmp[j] = pvp[mug][nug];
                        }
                        else
                        {
                            // need to get elements from the other half.
                            // I have question on this! 2011-02-22
                            tmp[j] = conj(pvp[nug][mug]);
                        }
                    }
                }
            }
            else
            {
                const int mug0 = mug/GlobalV::NPOL;
                if (GlobalV::DOMAG)
                {
#ifdef _OPENMP
#pragma omp for
#endif
                    for (int j=0; j<nlocal; j++)
                    {
                        int j_flag = j & 1; // j % 2 == 0
                        int ij_same = i_flag ^ j_flag ? 0 : 1;
                        const int nug = this->gridt->trace_lo[j];
                        const int nug0 = nug/GlobalV::NPOL;
                        double sign = sign_table[j_flag];
                        // if the col element is on this processor.
                        if (nug >=0)
                        {
                            if (mug <= nug)
                            {
                                if (ij_same)
                                {
                                    //spin = 0;
                                    //spin = 3;
                                    tmp[j] = pvp_nc[0][mug0][nug0]+sign*pvp_nc[3][mug0][nug0];
                                }
                                else
                                {
                                    // spin = 1;
                                    // spin = 2;
                                    tmp[j] = pvp_nc[1][mug0][nug0] + sign*std::complex<double>(0.0,1.0) * pvp_nc[2][mug0][nug0];
                                }
                            }
                            else
                            {
                                if (ij_same)
                                {
                                    //spin = 0;
                                    //spin = 3;
                                    tmp[j] = conj(pvp_nc[0][nug0][mug0]+sign*pvp_nc[3][nug0][mug0]);
                                }
                                else
                                {
                                    // spin = 1;
                                    //spin = 2;
                                    tmp[j] = conj(pvp_nc[1][nug0][mug0] + sign*std::complex<double>(0.0,1.0) * pvp_nc[2][nug0][mug0]);
                                }
                            }
                        }
                    }
                }
                else
                {
#ifdef _OPENMP
#pragma omp for
#endif
                    for (int j=0; j<nlocal; j++)
                    {
                        int j_flag = j & 1; // j % 2 == 0
                        int ij_same = i_flag ^ j_flag ? 0 : 1;

                        if (!ij_same)
                            continue;

                        const int nug = this->gridt->trace_lo[j];
                        const int nug0 = nug/GlobalV::NPOL;
                        double sign = sign_table[j_flag];
                        // if the col element is on this processor.
                        if (nug >=0)
                        {
                            if (mug <= nug)
                            {
                                //spin = 0;
                                //spin = 3;
                                tmp[j] = pvp_nc[0][mug0][nug0]+sign*pvp_nc[3][mug0][nug0];
                            }
                            else
                            {
                                //spin = 0;
                                //spin = 3;
                                tmp[j] = conj(pvp_nc[0][nug0][mug0]+sign*pvp_nc[3][nug0][mug0]);
                            }
                        }
                    }
                }
            }
        }
#ifdef _OPENMP
#pragma omp single
{
#endif
        // collect the matrix after folding.
            Parallel_Reduce::reduce_pool(tmp.data(), tmp.size());
#ifdef _OPENMP
}
#endif

        //-----------------------------------------------------
        // NOW! Redistribute the Hamiltonian matrix elements
        // according to the HPSEPS's 2D distribution methods.
        //-----------------------------------------------------
#ifdef _OPENMP
#pragma omp for
#endif
        for (int j=0; j<nlocal; j++)
        {
            if (!LM->ParaV->in_this_processor(i,j))
            {
                continue;
            }
            // set the matrix value.
            LM->set_HSk(i,j,tmp[j],'L');
        }
    }
#ifdef _OPENMP
}
#endif

    // delete the tmp matrix.
    delete[] pvp;
    delete[] pvp_base;
    if(GlobalV::NSPIN==4)
    {
        for(int spin =0;spin<4;spin++)
        {
            delete[] pvp_nc[spin];
        }
        delete[] pvp_nc;
        delete[] pvp_nc_base;
    }
    ModuleBase::timer::tick("Gint_k","Distri");

    ModuleBase::timer::tick("Gint_k","folding_vl_k");
    return;
}

#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"

//transfer_pvpR, NSPIN = 1 or 2
void Gint_k::transfer_pvpR(hamilt::HContainer<double> *hR)
{
    ModuleBase::TITLE("Gint_k","transfer_pvpR");
    ModuleBase::timer::tick("Gint_k","transfer_pvpR");

    if(!pvpR_alloc_flag || this->hRGint == nullptr)
    {
        ModuleBase::WARNING_QUIT("Gint_k::destroy_pvpR","pvpR hasnot been allocated yet!");
    }
    this->hRGint->set_zero();

    const int npol = GlobalV::NPOL;

    for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
    {
        const int T1 = GlobalC::ucell.iat2it[iat];
        const int I1 = GlobalC::ucell.iat2ia[iat];
        {
            // atom in this grid piece.
            if(this->gridt->in_this_processor[iat])
            {
                Atom* atom1 = &GlobalC::ucell.atoms[T1];

                // get the start positions of elements.
                const int DM_start = this->gridt->nlocstartg[iat];

                // get the coordinates of adjacent atoms.
                auto& tau1 = GlobalC::ucell.atoms[T1].tau[I1];
                //GlobalC::GridD.Find_atom(tau1);	
                AdjacentAtomInfo adjs;
                GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1, &adjs);	
                // search for the adjacent atoms.
                int nad = 0;

                for (int ad = 0; ad < adjs.adj_num+1; ad++)
                {
                    // get iat2
                    const int T2 = adjs.ntype[ad];
                    const int I2 = adjs.natom[ad];
                    const int iat2 = GlobalC::ucell.itia2iat(T2, I2);


                    // adjacent atom is also on the grid.
                    if(this->gridt->in_this_processor[iat2])
                    {
                        Atom* atom2 = &GlobalC::ucell.atoms[T2];
                        auto dtau = adjs.adjacent_tau[ad] - tau1;
                        double distance = dtau.norm() * GlobalC::ucell.lat0;
                        double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                        if(distance < rcut)
                        {
                            if( iat > iat2)
                            {// skip the lower triangle.
                                nad++;
                                continue;
                            }
                            // calculate the distance between iat1 and iat2.
                            // ModuleBase::Vector3<double> dR = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                            auto& dR = adjs.box[ad];
                            //dR.x = adjs.box[ad].x;
                            //dR.y = adjs.box[ad].y;
                            //dR.z = adjs.box[ad].z;

                            int ixxx = DM_start + this->gridt->find_R2st[iat][nad];
                            
                            hamilt::BaseMatrix<double>* tmp_matrix = this->hRGint->find_matrix(iat, iat2, dR.x, dR.y, dR.z);
#ifdef __DEBUG
                            assert(tmp_matrix != nullptr);
#endif
                            double* tmp_pointer = tmp_matrix->get_pointer();
                            const double *vijR = &pvpR_reduced[0][ixxx];
                            for(int iw=0; iw<atom1->nw; iw++)
                            {
                                for(int iw2 = 0; iw2<atom2->nw; ++iw2)
                                {
                                    *tmp_pointer++ = *vijR++;
                                }
                            }
                            // save the lower triangle.
                            if(iat < iat2)//skip iat == iat2
                            {
                                hamilt::BaseMatrix<double>* conj_matrix = this->hRGint->find_matrix(iat2, iat, -dR.x, -dR.y, -dR.z);
#ifdef __DEBUG
                                assert(conj_matrix != nullptr);
#endif
                                tmp_pointer = tmp_matrix->get_pointer();
                                for(int iw=0; iw<atom1->nw; iw++)
                                {
                                    for(int iw2 = 0; iw2<atom2->nw; ++iw2)
                                    {
                                        conj_matrix->get_value(iw2, iw) = *tmp_pointer++; 
                                    }
                                }
                            }
                            ++nad;
                        }// end distane<rcut
                    }
                }// end ad
            }
        }// end ia
    }// end it

    // ===================================
    // transfer HR from Gint to Veff<OperatorLCAO<std::complex<double>, double>>
    // ===================================
#ifdef __MPI
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size == 1)
    {
        hR->add(*this->hRGint);
    }
    else
    {
        hamilt::transferSerials2Parallels(*this->hRGint, hR);
    }
#else
    hR->add(*this->hRGint);
#endif
    ModuleBase::timer::tick("Gint_k","transfer_pvpR");
    
    return;
}

//transfer_pvpR, NSPIN = 4
void Gint_k::transfer_pvpR(hamilt::HContainer<std::complex<double>> *hR)
{
    ModuleBase::TITLE("Gint_k","transfer_pvpR");
    ModuleBase::timer::tick("Gint_k","transfer_pvpR");

    if(!pvpR_alloc_flag || this->hRGintCd == nullptr)
    {
        ModuleBase::WARNING_QUIT("Gint_k::destroy_pvpR","pvpR hasnot been allocated yet!");
    }
    this->hRGintCd->set_zero();

    const int npol = GlobalV::NPOL;

    for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
    {
        const int T1 = GlobalC::ucell.iat2it[iat];
        const int I1 = GlobalC::ucell.iat2ia[iat];
        {
            // atom in this grid piece.
            if(this->gridt->in_this_processor[iat])
            {
                Atom* atom1 = &GlobalC::ucell.atoms[T1];

                // get the start positions of elements.
                const int DM_start = this->gridt->nlocstartg[iat];

                // get the coordinates of adjacent atoms.
                auto& tau1 = GlobalC::ucell.atoms[T1].tau[I1];
                //GlobalC::GridD.Find_atom(tau1);	
                AdjacentAtomInfo adjs;
                GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1, &adjs);	
                // search for the adjacent atoms.
                int nad = 0;

                for (int ad = 0; ad < adjs.adj_num+1; ad++)
                {
                    // get iat2
                    const int T2 = adjs.ntype[ad];
                    const int I2 = adjs.natom[ad];
                    const int iat2 = GlobalC::ucell.itia2iat(T2, I2);

                    // adjacent atom is also on the grid.
                    if(this->gridt->in_this_processor[iat2])
                    {
                        Atom* atom2 = &GlobalC::ucell.atoms[T2];
                        auto dtau = adjs.adjacent_tau[ad] - tau1;
                        double distance = dtau.norm() * GlobalC::ucell.lat0;
                        double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                        if(distance < rcut)
                        {
                            if( iat > iat2)
                            {// skip the lower triangle.
                                nad++;
                                continue;
                            }
                            // calculate the distance between iat1 and iat2.
                            // ModuleBase::Vector3<double> dR = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                            auto& dR = adjs.box[ad];
                            //dR.x = adjs.box[ad].x;
                            //dR.y = adjs.box[ad].y;
                            //dR.z = adjs.box[ad].z;

                            int ixxx = DM_start + this->gridt->find_R2st[iat][nad];
                            
                            hamilt::BaseMatrix<std::complex<double>>* tmp_matrix = this->hRGintCd->find_matrix(iat, iat2, dR.x, dR.y, dR.z);
#ifdef __DEBUG
                            assert(tmp_matrix != nullptr);
#endif
                            std::complex<double>* tmp_pointer = tmp_matrix->get_pointer();
                            std::vector<int> step_trace(4, 0);
                            for (int is = 0; is < 2; is++)
                            {
                                for (int is2 = 0; is2 < 2; is2++)
                                {
                                    step_trace[is * 2 + is2] = atom2->nw * 2 * is + is2;
                                }
                            }
                            const double* vijR[4];
                            for(int spin = 0; spin < 4; spin++) 
                            {
                                vijR[spin] = &pvpR_reduced[spin][ixxx];
                            }
                            for(int iw=0; iw<atom1->nw; iw++)
                            {
                                for(int iw2 = 0; iw2<atom2->nw; ++iw2)
                                {
                                    tmp_pointer[step_trace[0]] = *vijR[0] + *vijR[3]; 
                                    tmp_pointer[step_trace[3]] = *vijR[0] - *vijR[3];
                                    tmp_pointer += 2;
                                    vijR[0]++;
                                    vijR[3]++;
                                }
                                tmp_pointer += 2 * atom2->nw;
                            }
                            if(GlobalV::DOMAG)
                            {
                                tmp_pointer = tmp_matrix->get_pointer();
                                for(int iw=0; iw<atom1->nw; iw++)
                                {
                                    for(int iw2 = 0; iw2<atom2->nw; ++iw2)
                                    {
                                        tmp_pointer[step_trace[1]] = *vijR[1] + std::complex<double>(0.0,1.0) * *vijR[2]; 
                                        tmp_pointer[step_trace[2]] = *vijR[1] - std::complex<double>(0.0,1.0) * *vijR[2];
                                        tmp_pointer += 2;
                                        vijR[1]++;
                                        vijR[2]++;
                                    }
                                    tmp_pointer += 2 * atom2->nw;
                                }
                            }
                            // save the lower triangle.
                            if(iat < iat2)
                            {
                                hamilt::BaseMatrix<std::complex<double>>* conj_matrix = this->hRGintCd->find_matrix(iat2, iat, -dR.x, -dR.y, -dR.z);
#ifdef __DEBUG
                                assert(conj_matrix != nullptr);
#endif
                                tmp_pointer = tmp_matrix->get_pointer();
                                for(int iw1=0; iw1<atom1->nw * 2; ++iw1)
                                {
                                    for(int iw2 = 0; iw2<atom2->nw * 2; ++iw2)
                                    {
                                        conj_matrix->get_value(iw2, iw1) = conj(*tmp_pointer);
                                        tmp_pointer++;
                                    }
                                }
                            }
                            ++nad;
                        }// end distane<rcut
                    }
                }// end ad
            }
        }// end ia
    }// end it

    // ===================================
    // transfer HR from Gint to Veff<OperatorLCAO<std::complex<double>, std::complex<double>>>
    // ===================================
#ifdef __MPI
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size == 1)
    {
        hR->add(*this->hRGintCd);
    }
    else
    {
        hamilt::transferSerials2Parallels<std::complex<double>>(*this->hRGintCd, hR);
    }
#else
    hR->add(*this->hRGintCd);
#endif
    
    ModuleBase::timer::tick("Gint_k","transfer_pvpR");
    return;
}