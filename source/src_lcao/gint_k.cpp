#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_reduce.h"
#include "gint_k.h"
#include "../module_orbital/ORB_read.h"
#include "grid_technique.h"
#include "../module_base/ylm.h"
#include "../src_pw/global.h"
#include "global_fp.h" // mohan add 2021-01-30
#include "../module_base/memory.h"
#include "../module_base/timer.h"

Gint_k::Gint_k()
{
    ik_now = 0;	
    pvpR_alloc_flag = false;
    spin_now = -1; // for a start value, must not equal 1,2 or 4.
    reduced = true;// use reduced memory for H storage.
}

Gint_k::~Gint_k()
{
    
}

void Gint_k::reset_spin(const int &spin_now_in)
{
    this->spin_now = spin_now_in;
    return;
}




void Gint_k::allocate_pvpR(void)
{
    ModuleBase::TITLE("Gint_k","allocate_pvpR");

    if(this->pvpR_alloc_flag)
    {
        return; //Liuxh add, 20181012
        ModuleBase::WARNING_QUIT("Gint_k::allocate_pvpR","pvpR has been allocated!");
    }

    //	reduced = GlobalV::NURSE; 
    //xiaohui modify 2015-05-30
    //std::cout << " reduced algorithm for grid integration = " << reduced << std::endl;

    if(this->reduced)
    {
        // the number of matrix element <phi_0 | V | phi_R> is GlobalC::GridT.nnrg.
        this->pvpR_reduced = new double*[GlobalV::NSPIN];
        for(int is =0;is<GlobalV::NSPIN;is++)
        {
            this->pvpR_reduced[is] = new double[GlobalC::GridT.nnrg];	
            ModuleBase::GlobalFunc::ZEROS( pvpR_reduced[is], GlobalC::GridT.nnrg);
        }

        double mem = ModuleBase::Memory::record("allocate_pvpR", "pvpR_reduced", GlobalC::GridT.nnrg * GlobalV::NSPIN , "double");

        if(GlobalV::OUT_LEVEL != "m") 
        {
            GlobalV::ofs_running << " Memory of pvpR : " << mem << " MB" << std::endl;
        }

        if( mem > 800 )
        {
            GlobalV::ofs_warning << " memory for pvpR = " << mem << std::endl;
            GlobalV::ofs_warning << " which is larger than 800 MB ! " << std::endl;
            ModuleBase::WARNING_QUIT("Gint_k","allocate_pvpR");
        }

    }
    else
    {
        double mem = ModuleBase::Memory::record("allocate_pvpR", "pvpR", GlobalC::GridT.lgd * GlobalC::GridT.nutot
                * GlobalC::GridT.lgd * GlobalC::GridT.nutot , "double");

        if(GlobalV::OUT_LEVEL != "m") 
        {
            GlobalV::ofs_running << " Memory of pvpR : " << mem << " MB" << std::endl;
        }

        if( mem > 800 )
        {
            GlobalV::ofs_warning << " memory for pvpR = " << mem << std::endl;
            GlobalV::ofs_warning << " which is larger than 800 MB ! " << std::endl;
            ModuleBase::WARNING_QUIT("Gint_k","allocate_pvpR");
        }

        //----------------------------------------------
        // allocate the std::complex matrix !!
        // nutot : total number of unitcells involved.
        // this may be very large, at least 
        // 3*3*3 = 27.
        //----------------------------------------------
        const int LDIM=GlobalC::GridT.lgd*GlobalC::GridT.nutot;

        this->pvpR_pool = new double[LDIM*LDIM];
        ModuleBase::GlobalFunc::ZEROS(pvpR_pool, LDIM*LDIM);

        this->pvpR = new double*[LDIM];
        for(int i=0; i<LDIM; i++)
        {
            pvpR[i] = &pvpR_pool[i*LDIM];
        }
    }

    this->pvpR_alloc_flag = true;
    return;
}



void Gint_k::destroy_pvpR(void)
{
    ModuleBase::TITLE("Gint_k","destroy_pvpR");
    
    if(!pvpR_alloc_flag)
    {
        ModuleBase::WARNING_QUIT("Gint_k::destroy_pvpR","<phi_0i | V | phi_Rj> matrix has not been allocated yet!");
    }
    
    if(this->reduced)
    {
        for(int is =0;is<GlobalV::NSPIN;is++) delete[] pvpR_reduced[is];
        delete[] pvpR_reduced;
    }	
    else
    {
        delete[] pvpR;
        delete[] pvpR_pool;
    }

    this->pvpR_alloc_flag = false;
    return;
}





// fold the <phi | vl |dphi(R)> * DM(R) to 
// calculate the force.
void Gint_k::folding_force(
    ModuleBase::matrix& fvl_dphi,
    double* pvdpx, 
    double* pvdpy, 
    double* pvdpz)
{
    ModuleBase::TITLE("Gint_k","folding_force");
    ModuleBase::timer::tick("Gint_k","folding_force");

    //xiaohui modify 2013-12-17, test
//	assert(GlobalC::GridT.lgd > 0); //mohan add 2012-06-10

    // mohan add 2014-01-20
    const int lgd = GlobalC::GridT.lgd;

    double** ppx;
    double** ppy;
    double** ppz;

    if(GlobalC::GridT.lgd>0)
    {
        ppx = new double*[lgd];
          ppy = new double*[lgd];
          ppz = new double*[lgd];
          for(int i=0; i<lgd; i++)
          {
            ppx[i] = new double[lgd];
            ppy[i] = new double[lgd];
            ppz[i] = new double[lgd];
            ModuleBase::GlobalFunc::ZEROS( ppx[i], lgd);
            ModuleBase::GlobalFunc::ZEROS( ppy[i], lgd);
            ModuleBase::GlobalFunc::ZEROS( ppz[i], lgd);
        }
    }
    
    ModuleBase::Vector3<double> tau1, dtau;
    for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1=0; I1< atom1->na; ++I1)
        {
            const int iat = GlobalC::ucell.itia2iat(T1,I1);
            if(GlobalC::GridT.in_this_processor[iat])
            {
                assert( lgd > 0 );

                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                // get the start positions of elements.
                const int DM_start = GlobalC::GridT.nlocstartg[iat];
                // get the coordinates of adjacent atoms.
                tau1 = atom1->tau[I1];
                //GlobalC::GridD.Find_atom(tau1);
                GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
                // search for the adjacent atoms.
                int nad = 0;
                for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
                {
                    // get iat2
                    const int T2 = GlobalC::GridD.getType(ad);
                    const int I2 = GlobalC::GridD.getNatom(ad);
                    const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                    if(GlobalC::GridT.in_this_processor[iat2])
                    {
                        Atom* atom2 = &GlobalC::ucell.atoms[T2];
                        dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                        double distance = dtau.norm() * GlobalC::ucell.lat0;
                        double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
                        if(distance < rcut)
                        {
                            const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                            int ixxx = DM_start + GlobalC::GridT.find_R2st[iat][nad];
                            for(int iw=0; iw<atom1->nw; iw++)
                            {
                                const int iw_all = start1+iw;
                                const int iw_local = GlobalC::GridT.trace_lo[iw_all];
                                // iw1_lo
                                double *vijx = ppx[iw_local];
                                double *vijy = ppy[iw_local];
                                double *vijz = ppz[iw_local];

                                double *vRx = &pvdpx[ixxx]; //just fold R to normal matrix.
                                double *vRy = &pvdpy[ixxx];
                                double *vRz = &pvdpz[ixxx];

                                int* iw2_lo = &GlobalC::GridT.trace_lo[start2];
                                int* iw2_end = iw2_lo + atom2->nw;

                                for(; iw2_lo<iw2_end; ++iw2_lo, ++vRx, ++vRy, ++vRz)
                                {
                                    vijx[iw2_lo[0]] += vRx[0] ;
                                    vijy[iw2_lo[0]] += vRy[0] ;
                                    vijz[iw2_lo[0]] += vRz[0] ;
                                }
                                ixxx += atom2->nw;
                            }
                            ++nad;
                        }//end distance<rcut
                    }
                }//end ad
            }
        }//end ia
    }//end it



    double* tmp = new double[GlobalV::NLOCAL*3];
    for(int i=0; i<GlobalV::NLOCAL; ++i)
    {
        ModuleBase::GlobalFunc::ZEROS(tmp, 3*GlobalV::NLOCAL);
        const int mug = GlobalC::GridT.trace_lo[i];
        // if the row element is on this processor
        if(mug>=0)
        {
            //GlobalV::ofs_running << " i=" << i << " mug=" << mug << std::endl;
            for(int j=0; j<GlobalV::NLOCAL; ++j)
            {
                const int nug = GlobalC::GridT.trace_lo[j];
                // if the col element is on this processor
                if(nug>=0)
                {
    //				if(mug<nug)
    //				{
                        const int index = 3*j;
                        tmp[index] = ppx[mug][nug];
                        tmp[index+1] = ppy[mug][nug];
                        tmp[index+2] = ppz[mug][nug];
    //				}
    //				else
    //				{
                    //	tmpx[j] = 0.0;
                    //	tmpy[j] = 0.0;
                    //	tmpz[j] = 0.0;
    //				}
                }
            }
        }
        // collect the matrix after folding.
        Parallel_Reduce::reduce_double_pool( tmp, GlobalV::NLOCAL*3 );
        for (int j=0; j<GlobalV::NLOCAL; j++)
        {
            if (!this->LM->ParaV->in_this_processor(i,j))
            {
                continue;
            }
            const int iat = GlobalC::ucell.iwt2iat[i];
            const int index = 3*j;
            fvl_dphi(iat,0) += 2.0*tmp[index];	
            fvl_dphi(iat,1) += 2.0*tmp[index+1];	
            fvl_dphi(iat,2) += 2.0*tmp[index+2];	
        }
    }
    delete[] tmp;

    // mohan add 2014-01-20
    if(GlobalC::GridT.lgd > 0)
    {
        //-------------------------
        // delete the tmp matrix.
        //-------------------------
        for(int i=0; i<GlobalC::GridT.lgd; i++)
        {
            delete[] ppx[i];
            delete[] ppy[i];
            delete[] ppz[i];
        }
        delete[] ppx;
        delete[] ppy;
        delete[] ppz;
    }
    ModuleBase::timer::tick("Gint_k","folding_force");
    return;
}

// fold the <phi | vl * R_beta|dphi(R_alpha)> * DM(R) to 
// calculate the stress.
void Gint_k::folding_stress(
    const bool isforce,
    const bool isstress,
    ModuleBase::matrix& fvl_dphi, 
    ModuleBase::matrix& svl_dphi,
    double* pvdpx, 
    double* pvdpy, 
    double* pvdpz,
    double* pvdp11, 
    double* pvdp22, 
    double* pvdp33,
    double* pvdp12, 
    double* pvdp13, 
    double* pvdp23)
{
    ModuleBase::TITLE("Gint_k","folding_stress");
    if(!isforce&&!isstress) return;

    const int lgd = GlobalC::GridT.lgd;

    if(isforce)
    {
        ModuleBase::timer::tick("Gint_k","folding_force");

        double** ppx;
        double** ppy;
        double** ppz;
        if(lgd>0)
        {
            ppx = new double*[lgd];
            ppy = new double*[lgd];
            ppz = new double*[lgd];
            for(int i=0; i<lgd; i++)
            {
                ppx[i] = new double[lgd];
                ppy[i] = new double[lgd];
                ppz[i] = new double[lgd];
                ModuleBase::GlobalFunc::ZEROS( ppx[i], lgd);
                ModuleBase::GlobalFunc::ZEROS( ppy[i], lgd);
                ModuleBase::GlobalFunc::ZEROS( ppz[i], lgd);
            }
        }
        ModuleBase::Vector3<double> tau1, dtau;
        for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
        {
            const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                    
            for (int I1 =0; I1< atom1->na; ++I1)
            {
                const int iat = GlobalC::ucell.itia2iat(T1,I1);
                if(GlobalC::GridT.in_this_processor[iat])
                {
                    assert( lgd > 0 );

                    const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                    // get the start positions of elements.
                    const int DM_start = GlobalC::GridT.nlocstartg[iat];
                    // get the coordinates of adjacent atoms.
                    tau1 = atom1->tau[I1];
                    //GlobalC::GridD.Find_atom(tau1);
                    GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
                    // search for the adjacent atoms.
                    int nad = 0;
                    for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
                    {
                        // get iat2
                        const int T2 = GlobalC::GridD.getType(ad);
                        const int I2 = GlobalC::GridD.getNatom(ad);

                        const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad);
                        const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                        if(GlobalC::GridT.in_this_processor[iat2])
                        {
                            Atom* atom2 = &GlobalC::ucell.atoms[T2];
                            dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                            double distance = dtau.norm() * GlobalC::ucell.lat0;
                            double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
                            if(distance < rcut)
                            {
                                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

                                //index for iteration below
                                const int* iw2_lo = &GlobalC::GridT.trace_lo[start2];
                                const int* iw2_end = iw2_lo + atom2->nw;

                                int ixxx = DM_start + GlobalC::GridT.find_R2st[iat][nad];
                                for(int iw=0; iw<atom1->nw; iw++)
                                {
                                    const int iw_all = start1+iw;
                                    const int iw_local = GlobalC::GridT.trace_lo[iw_all];
                                    // iw1_lo
                                    double *vijx = ppx[iw_local];
                                    double *vijy = ppy[iw_local];
                                    double *vijz = ppz[iw_local];

                                    //just fold R to normal matrix.
                                    double *vRx = &pvdpx[ixxx];
                                    double *vRy = &pvdpy[ixxx];
                                    double *vRz = &pvdpz[ixxx];

                                    for(auto iteration=iw2_lo; iteration<iw2_end; ++iteration, ++vRx, ++vRy, ++vRz)
                                    {
                                        vijx[iteration[0]] += vRx[0];
                                        vijy[iteration[0]] += vRy[0];
                                        vijz[iteration[0]] += vRz[0];
                                    }

                                    ixxx += atom2->nw;
                                }
                                ++nad;
                            }//end distance<rcut
                        }
                    }//end ad
                }
            }//end ia
        }//end it

        double* tmp = new double[GlobalV::NLOCAL*3];
        for(int i=0; i<GlobalV::NLOCAL; ++i)
        {
            ModuleBase::GlobalFunc::ZEROS(tmp, 3*GlobalV::NLOCAL);
            const int mug = GlobalC::GridT.trace_lo[i];
            // if the row element is on this processor
            if(mug>=0)
            {
                //GlobalV::ofs_running << " i=" << i << " mug=" << mug << std::endl;
                for(int j=0; j<GlobalV::NLOCAL; ++j)
                {
                    const int nug = GlobalC::GridT.trace_lo[j];
                    // if the col element is on this process or
                    if(nug>=0)
                    {
                        const int index = 3*j;
                        tmp[index] = ppx[mug][nug];
                        tmp[index+1] = ppy[mug][nug];
                        tmp[index+2] = ppz[mug][nug];
                    }
                }
            }
            // collect the matrix after folding.
            Parallel_Reduce::reduce_double_pool( tmp, GlobalV::NLOCAL*3 );
            for (int j=0; j<GlobalV::NLOCAL; j++)
            {
                if (!this->LM->ParaV->in_this_processor(i,j))
                {
                    continue;
                }
                const int iat = GlobalC::ucell.iwt2iat[i];
                const int index = 3*j;
                fvl_dphi(iat,0) += 2.0*tmp[index];
                fvl_dphi(iat,1) += 2.0*tmp[index+1];
                fvl_dphi(iat,2) += 2.0*tmp[index+2];
            }
        }
        delete[] tmp;

        if(lgd > 0)
        {
            //-------------------------
            // delete the tmp matrix.
            //-------------------------
            for(int i=0; i<lgd; i++)
            {
                delete[] ppx[i];
                delete[] ppy[i];
                delete[] ppz[i];
            }
            delete[] ppx;
            delete[] ppy;
            delete[] ppz;
        }
        ModuleBase::timer::tick("Gint_k","folding_force");
    }
    if(isstress)
    {
        ModuleBase::timer::tick("Gint_k","folding_stress");

        double** pp11;
        double** pp22;
        double** pp33;
        double** pp12;
        double** pp13;
        double** pp23;

        if(lgd>0)
        {
            pp11 = new double*[lgd];
            pp22 = new double*[lgd];
            pp33 = new double*[lgd];
            pp12 = new double*[lgd];
            pp13 = new double*[lgd];
            pp23 = new double*[lgd];
            for(int i=0; i<lgd; i++)
            {
                pp11[i] = new double[lgd];
                pp22[i] = new double[lgd];
                pp33[i] = new double[lgd];
                ModuleBase::GlobalFunc::ZEROS( pp11[i], lgd);
                ModuleBase::GlobalFunc::ZEROS( pp22[i], lgd);
                ModuleBase::GlobalFunc::ZEROS( pp33[i], lgd);
                pp12[i] = new double[lgd];
                pp13[i] = new double[lgd];
                pp23[i] = new double[lgd];
                ModuleBase::GlobalFunc::ZEROS( pp12[i], lgd);
                ModuleBase::GlobalFunc::ZEROS( pp13[i], lgd);
                ModuleBase::GlobalFunc::ZEROS( pp23[i], lgd);
            }
        }
        ModuleBase::Vector3<double> tau1, dtau;
        for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
        {
            const Atom* atom1 = &GlobalC::ucell.atoms[T1];
                    
            for (int I1 =0; I1< atom1->na; ++I1)
            {
                const int iat = GlobalC::ucell.itia2iat(T1,I1);
                if(GlobalC::GridT.in_this_processor[iat])
                {
                    assert( lgd > 0 );

                    const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                    // get the start positions of elements.
                    const int DM_start = GlobalC::GridT.nlocstartg[iat];
                    // get the coordinates of adjacent atoms.
                    tau1 = atom1->tau[I1];
                    //GlobalC::GridD.Find_atom(tau1);
                    GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
                    // search for the adjacent atoms.
                    int nad = 0;
                    for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
                    {
                        // get iat2
                        const int T2 = GlobalC::GridD.getType(ad);
                        const int I2 = GlobalC::GridD.getNatom(ad);

                        const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad);
                        const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                        if(GlobalC::GridT.in_this_processor[iat2])
                        {
                            Atom* atom2 = &GlobalC::ucell.atoms[T2];
                            dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                            double distance2 = dtau.norm2() * pow(GlobalC::ucell.lat0, 2.0);
                            double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();
                            if(distance2 < pow(rcut,2.0))
                            {
                                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

                                //index for iteration below
                                const int* iw2_lo = &GlobalC::GridT.trace_lo[start2];
                                const int* iw2_end = iw2_lo + atom2->nw;

                                int ixxx = DM_start + GlobalC::GridT.find_R2st[iat][nad];
                                for(int iw=0; iw<atom1->nw; iw++)
                                {
                                    const int iw_all = start1+iw;
                                    const int iw_local = GlobalC::GridT.trace_lo[iw_all];

                                    double *vij11 = pp11[iw_local];
                                    double *vij22 = pp22[iw_local];
                                    double *vij33 = pp33[iw_local];
                                    double *vij12 = pp12[iw_local];
                                    double *vij13 = pp13[iw_local];
                                    double *vij23 = pp23[iw_local];

                                    double *vR11 = &pvdp11[ixxx]; 
                                    double *vR22 = &pvdp22[ixxx];
                                    double *vR33 = &pvdp33[ixxx];
                                    double *vR12 = &pvdp12[ixxx]; 
                                    double *vR13 = &pvdp13[ixxx];
                                    double *vR23 = &pvdp23[ixxx];

                                    for(auto iteration=iw2_lo; iteration<iw2_end; ++iteration, ++vR11, 
                                        ++vR22,++vR33,++vR12,++vR13,++vR23)
                                    {
                                        vij11[iteration[0]] += vR11[0];
                                        vij22[iteration[0]] += vR22[0];
                                        vij33[iteration[0]] += vR33[0];
                                        vij12[iteration[0]] += vR12[0];
                                        vij13[iteration[0]] += vR13[0];
                                        vij23[iteration[0]] += vR23[0];
                                    }

                                    ixxx += atom2->nw;
                                }
                                ++nad;
                            }//end distance<rcut
                        }
                    }//end ad
                }
            }//end ia
        }//end it

        double* tmp1 = new double[GlobalV::NLOCAL*6];
        for(int i=0; i<GlobalV::NLOCAL; ++i)
        {
            ModuleBase::GlobalFunc::ZEROS(tmp1, 6*GlobalV::NLOCAL);
            const int mug = GlobalC::GridT.trace_lo[i];
            // if the row element is on this processor
            if(mug>=0)
            {
                //GlobalV::ofs_running << " i=" << i << " mug=" << mug << std::endl;
                for(int j=0; j<GlobalV::NLOCAL; ++j)
                {
                    const int nug = GlobalC::GridT.trace_lo[j];
                    // if the col element is on this process or
                    if(nug>=0)
                    {
                            const int index1 = 6*j;
                            tmp1[index1] = pp11[mug][nug];
                            tmp1[index1+1] = pp22[mug][nug];
                            tmp1[index1+2] = pp33[mug][nug];
                            tmp1[index1+3] = pp12[mug][nug];
                            tmp1[index1+4] = pp13[mug][nug];
                            tmp1[index1+5] = pp23[mug][nug];

                    }
                }
            }
            // collect the matrix after folding.
            Parallel_Reduce::reduce_double_pool( tmp1, GlobalV::NLOCAL*6 );
            for (int j=0; j<GlobalV::NLOCAL; j++)
            {
                if (!this->LM->ParaV->in_this_processor(i,j))
                {
                    continue;
                }
                const int index1 = 6*j;
                svl_dphi(0,0) -= 2.0*tmp1[index1];
                svl_dphi(1,1) -= 2.0*tmp1[index1+1];
                svl_dphi(2,2) -= 2.0*tmp1[index1+2];
                svl_dphi(0,1) -= 2.0*tmp1[index1+3];
                svl_dphi(0,2) -= 2.0*tmp1[index1+4];
                svl_dphi(1,2) -= 2.0*tmp1[index1+5];
            }
        }

        delete[] tmp1;

        if(lgd > 0)
        {
            //-------------------------
            // delete the tmp matrix.
            //-------------------------
            for(int i=0; i<lgd; i++)
            {
                delete[] pp11[i];
                delete[] pp22[i];
                delete[] pp33[i];
                delete[] pp12[i];
                delete[] pp13[i];
                delete[] pp23[i];
            }
            delete[] pp11;
            delete[] pp22;
            delete[] pp33;
            delete[] pp12;
            delete[] pp13;
            delete[] pp23;
        }
        ModuleBase::timer::tick("Gint_k","folding_stress");
    }
    return;
}

// folding the matrix for 'ik' k-point.
// H(k)=\sum{R} H(R)exp(ikR) 
void Gint_k::folding_vl_k(const int &ik)
{
    ModuleBase::TITLE("Gint_k","folding_vl_k");
    ModuleBase::timer::tick("Gint_k","folding_vl_k");

    if(!pvpR_alloc_flag)
    {
        ModuleBase::WARNING_QUIT("Gint_k::destroy_pvpR","pvpR hasnot been allocated yet!");
    }

    //####################### EXPLAIN #################################
    // 1. what is GlobalC::GridT.lgd ?
    // GlobalC::GridT.lgd is the number of orbitals in each processor according
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
    this->ik_now = ik;
    this->pvp = new std::complex<double>*[GlobalC::GridT.lgd];
    for(int i=0; i<GlobalC::GridT.lgd; i++)
    {
        this->pvp[i] = new std::complex<double>[GlobalC::GridT.lgd];
        ModuleBase::GlobalFunc::ZEROS( this->pvp[i], GlobalC::GridT.lgd);
    }

    if(!reduced)
    {	
        ModuleBase::Vector3<double> dR;
        double arg;
        std::complex<double> phase;
        std::complex<double> *pp1;
        double *pp2;
        int count;
        for(int k=0; k<GlobalC::GridT.nutot; k++)
        {
            const int R1x = GlobalC::GridT.ucell_index2x[k];
            const int R1y = GlobalC::GridT.ucell_index2y[k];
            const int R1z = GlobalC::GridT.ucell_index2z[k];

            const int dimk = GlobalC::GridT.lgd*k;
            for(int m=0; m<GlobalC::GridT.nutot; m++)
            {
                //------------------------------------------------
                // exp(k dot dR)
                // dR is the index of box in Crystal coordinates
                //------------------------------------------------
                dR.x = GlobalC::GridT.ucell_index2x[m] - R1x;
                dR.y = GlobalC::GridT.ucell_index2y[m] - R1y;
                dR.z = GlobalC::GridT.ucell_index2z[m] - R1z;

                arg = (GlobalC::kv.kvec_d[ this->ik_now ] * dR) * ModuleBase::TWO_PI;
                phase = std::complex<double>(cos(arg), sin(arg));
                for(int i=0; i<GlobalC::GridT.lgd; i++)
                {
                    pp1 = this->pvp[i];
                    pp2 = this->pvpR[i+GlobalC::GridT.lgd*m];
                    count = dimk;
                    for(int j=0; j<GlobalC::GridT.lgd; j++)
                    {
                        // folding matrix
                        pp1[j] += pp2[count] * phase;
                        ++count;
                    }
                }
            }
        }
    }
    else
    {
        int lgd = 0;
        ModuleBase::Vector3<double> tau1, dtau, dR;
        for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
        {
            for(int I1=0; I1<GlobalC::ucell.atoms[T1].na; ++I1)
            {
                // get iat
                const int iat = GlobalC::ucell.itia2iat(T1,I1);
                // atom in this grid piece.
                if(GlobalC::GridT.in_this_processor[iat])
                {
                    Atom* atom1 = &GlobalC::ucell.atoms[T1];
                    const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

                    // get the start positions of elements.
                    const int DM_start = GlobalC::GridT.nlocstartg[iat];

                    // get the coordinates of adjacent atoms.
                    tau1 = GlobalC::ucell.atoms[T1].tau[I1];
                    //GlobalC::GridD.Find_atom(tau1);	
                    GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);	
                    // search for the adjacent atoms.
                    int nad = 0;


                    for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ad++)
                    {
                        // get iat2
                        const int T2 = GlobalC::GridD.getType(ad);
                        const int I2 = GlobalC::GridD.getNatom(ad);
                        const int iat2 = GlobalC::ucell.itia2iat(T2, I2);


                        // adjacent atom is also on the grid.
                        if(GlobalC::GridT.in_this_processor[iat2])
                        {
                            Atom* atom2 = &GlobalC::ucell.atoms[T2];
                            dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                            double distance = dtau.norm() * GlobalC::ucell.lat0;
                            double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                            // for the local part, only need to calculate <phi_i | phi_j> within range
                            // mohan note 2012-07-06
                            if(distance < rcut)
                            {
                                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0); 

                                // calculate the distance between iat1 and iat2.
                                // ModuleBase::Vector3<double> dR = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                                dR.x = GlobalC::GridD.getBox(ad).x;
                                dR.y = GlobalC::GridD.getBox(ad).y;
                                dR.z = GlobalC::GridD.getBox(ad).z;

                                // calculate the phase factor exp(ikR).
                                const double arg = (GlobalC::kv.kvec_d[ this->ik_now ] * dR) * ModuleBase::TWO_PI;
                                std::complex<double> phase = std::complex<double>(cos(arg), sin(arg));
                                int ixxx = DM_start + GlobalC::GridT.find_R2st[iat][nad];
                                for(int iw=0; iw<atom1->nw; iw++)
                                {
                                    // iw1_lo
                                    std::complex<double> *vij = this->pvp[GlobalC::GridT.trace_lo[start1+iw]];


                                    int* iw2_lo = &GlobalC::GridT.trace_lo[start2];
                                    int* iw2_end = iw2_lo + atom2->nw;

                                    // get the <phi | V | phi>(R) Hamiltonian.
                                    double *vijR = &pvpR_reduced[0][ixxx];
                                    // std::complex<double> *vijR_soc = &pvpR_reduced_soc[ixxx];
                                    for(; iw2_lo<iw2_end; ++iw2_lo, ++vijR)
                                    {
                                        vij[iw2_lo[0]] += vijR[0] * phase; 
                                    }
                                    ixxx += atom2->nw;
                                    ++lgd;
                                }

                                ++nad;
                            }// end distane<rcut
                        }
                    }// end ad
                }
            }// end ia
        }// end it

        //------------------
        // To test the pvpR
        //------------------
/*
        for(int i=0; i<GlobalC::GridT.nlocdimg[0]; i++)
        {
            const int DM_start = GlobalC::GridT.nlocstartg[0];
            const int j = i + DM_start;
            if( abs(pvpR_reduced[j]) > 1.0e-5  )
            {
//				std::cout << " pvpR_reduced[" << i <<"] = " << pvpR_reduced[j] << std::endl;
            }
        }
*/

    }

    //----------------------
    // Print the pvp matrix
    //----------------------
/*
    std::cout << " pvp matrix:" << std::endl;
    for(int i=0; i<GlobalC::GridT.lgd; i++)
    {
        for(int j=0; j<GlobalC::GridT.lgd; j++)
        {
            std::cout << std::setw(15) << pvp[i][j].real();
        }
        std::cout << std::endl;
    }
    */

    // Distribution of data.
    ModuleBase::timer::tick("Gint_k","Distri");
    std::complex<double>* tmp;
    for (int i=0; i<GlobalV::NLOCAL; i++)
    {
        tmp = new std::complex<double>[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(tmp, GlobalV::NLOCAL);
        const int mug = GlobalC::GridT.trace_lo[i];
        // if the row element is on this processor.
        if (mug >= 0)
        {
            for (int j=0; j<GlobalV::NLOCAL; j++)
            {
                const int nug = GlobalC::GridT.trace_lo[j];
                // if the col element is on this processor.
                if (nug >=0)
                {
                    if (mug <= nug)
                    {
                        // pvp is symmetric, only half is calculated.
                        tmp[j] = this->pvp[mug][nug];
                    }
                    else
                    {
                        // need to get elements from the other half.
                        // I have question on this! 2011-02-22
                        tmp[j] = conj(this->pvp[nug][mug]);
                    }
                }
            }
        }
        // collect the matrix after folding.
        Parallel_Reduce::reduce_complex_double_pool( tmp, GlobalV::NLOCAL );

        //-----------------------------------------------------
        // NOW! Redistribute the Hamiltonian matrix elements
        // according to the HPSEPS's 2D distribution methods.
        //-----------------------------------------------------
        for (int j=0; j<GlobalV::NLOCAL; j++)
        {
            if (!this->LM->ParaV->in_this_processor(i,j))
            {
                continue;
            }
            // set the matrix value.
            this->LM->set_HSk(i,j,tmp[j],'L');
        }
        delete[] tmp;
    }

    // delete the tmp matrix.
    for(int i=0; i<GlobalC::GridT.lgd; i++)
    {
        delete[] pvp[i];
    }
    delete[] pvp;
    ModuleBase::timer::tick("Gint_k","Distri");

    ModuleBase::timer::tick("Gint_k","folding_vl_k");
    return;
}

// folding the matrix for 'ik' k-point.
// H(k)=\sum{R} H(R)exp(ikR) 
// for non-collinear case  
void Gint_k::folding_vl_k_nc(const int &ik)
{
    ModuleBase::TITLE("Gint_k","folding_vl_k_nc");
    ModuleBase::timer::tick("Gint_k","folding_vl_k_nc");

    if(!pvpR_alloc_flag)
    {
        ModuleBase::WARNING_QUIT("Gint_k::destroy_pvpR","pvpR hasnot been allocated yet!");
    }

    this->ik_now = ik;
//	std::complex<double>** pvp_nc[4];
    for(int spin=0;spin<4;spin++)
    {
        pvp_nc[spin] = new std::complex<double>*[GlobalC::GridT.lgd];
        for(int i=0; i<GlobalC::GridT.lgd; i++)
        {
            pvp_nc[spin][i] = new std::complex<double>[GlobalC::GridT.lgd];
            ModuleBase::GlobalFunc::ZEROS( this->pvp_nc[spin][i], GlobalC::GridT.lgd);
        }
    }

    if(!reduced)
    {	
        ModuleBase::Vector3<double> dR;
        double arg;
        std::complex<double> phase;
        std::complex<double> *pp1;
        double *pp2;
        int count;
        for(int k=0; k<GlobalC::GridT.nutot; k++)
        {
            const int R1x = GlobalC::GridT.ucell_index2x[k];
            const int R1y = GlobalC::GridT.ucell_index2y[k];
            const int R1z = GlobalC::GridT.ucell_index2z[k];

            const int dimk = GlobalC::GridT.lgd*k;
            for(int m=0; m<GlobalC::GridT.nutot; m++)
            {
                //------------------------------------------------
                // exp(k dot dR)
                // dR is the index of box in Crystal coordinates
                //------------------------------------------------
                dR.x = GlobalC::GridT.ucell_index2x[m] - R1x;
                dR.y = GlobalC::GridT.ucell_index2y[m] - R1y;
                dR.z = GlobalC::GridT.ucell_index2z[m] - R1z;

                arg = (GlobalC::kv.kvec_d[ this->ik_now ] * dR) * ModuleBase::TWO_PI;
                phase = std::complex<double>(cos(arg), sin(arg));
                for(int i=0; i<GlobalC::GridT.lgd; i++)
                {
                    pp1 = this->pvp_nc[0][i];
                    pp2 = this->pvpR[i+GlobalC::GridT.lgd*m];
                    count = dimk;
                    for(int j=0; j<GlobalC::GridT.lgd; j++)
                    {
                        // folding matrix
                        pp1[j] += pp2[count] * phase;
                        ++count;
                    }
                }
            }
        }
    }
    else
    {
        int lgd = 0;
        ModuleBase::Vector3<double> tau1, dtau, dR;
        for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
        {
            for(int I1=0; I1<GlobalC::ucell.atoms[T1].na; ++I1)
            {
                // get iat
                const int iat = GlobalC::ucell.itia2iat(T1,I1);
                // atom in this grid piece.
                if(GlobalC::GridT.in_this_processor[iat])
                {
                    Atom* atom1 = &GlobalC::ucell.atoms[T1];
                    const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

                    // get the start positions of elements.
                    const int DM_start = GlobalC::GridT.nlocstartg[iat];

                    // get the coordinates of adjacent atoms.
                    tau1 = GlobalC::ucell.atoms[T1].tau[I1];
                    //GlobalC::GridD.Find_atom(tau1);	
                    GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);	
                    // search for the adjacent atoms.
                    int nad = 0;


                    for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ad++)
                    {
                        // get iat2
                        const int T2 = GlobalC::GridD.getType(ad);
                        const int I2 = GlobalC::GridD.getNatom(ad);
                        const int iat2 = GlobalC::ucell.itia2iat(T2, I2);


                        // adjacent atom is also on the grid.
                        if(GlobalC::GridT.in_this_processor[iat2])
                        {
                            Atom* atom2 = &GlobalC::ucell.atoms[T2];
                            dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                            double distance = dtau.norm() * GlobalC::ucell.lat0;
                            double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                            // for the local part, only need to calculate <phi_i | phi_j> within range
                            // mohan note 2012-07-06
                            if(distance < rcut)
                            {
                                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

                                // calculate the distance between iat1 and iat2.
                                // ModuleBase::Vector3<double> dR = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                                dR.x = GlobalC::GridD.getBox(ad).x;
                                dR.y = GlobalC::GridD.getBox(ad).y;
                                dR.z = GlobalC::GridD.getBox(ad).z;

                                // calculate the phase factor exp(ikR).
                                const double arg = (GlobalC::kv.kvec_d[ this->ik_now ] * dR) * ModuleBase::TWO_PI;
                                std::complex<double> phase = std::complex<double>(cos(arg), sin(arg));
                                int ixxx = DM_start + GlobalC::GridT.find_R2st[iat][nad];
                                for(int iw=0; iw<atom1->nw; iw++)
                                {
                                    // iw1_lo
                                    std::complex<double> *vij[4];
                                    for(int spin=0;spin<4;spin++)
                                        vij[spin] = this->pvp_nc[spin][GlobalC::GridT.trace_lo[start1]/GlobalV::NPOL + iw];


                                    int iw2_lo = GlobalC::GridT.trace_lo[start2]/GlobalV::NPOL;
                                    int iw2_end = iw2_lo + atom2->nw;

                                    // get the <phi | V | phi>(R) Hamiltonian.
//									std::complex<double> *vijR_soc = &pvpR_reduced_soc[ixxx];
                                    double *vijR[4];
                                    for(int spin = 0;spin<4;spin++) 
                                    {
                                        vijR[spin] = &pvpR_reduced[spin][ixxx];
                                    }
                                    for(; iw2_lo<iw2_end; ++iw2_lo, ++vijR[0], ++vijR[1],++vijR[2],++vijR[3])
                                    {
                                        for(int spin =0;spin<4;spin++)
                                        {
                                            vij[spin][iw2_lo] += vijR[spin][0] * phase; 
                                        }
//										else vij[iw2_lo[0]] += vijR_soc[0] * phase;
                                    }
                                    ixxx += atom2->nw;
                                    ++lgd;
                                }

                                ++nad;
                            }// end distane<rcut
                        }
                    }// end ad
                }
            }// end ia
        }// end it

        //------------------
        // To test the pvpR
        //------------------
/*
        for(int i=0; i<GlobalC::GridT.nlocdimg[0]; i++)
        {
            const int DM_start = GlobalC::GridT.nlocstartg[0];
            const int j = i + DM_start;
            if( abs(pvpR_reduced[j]) > 1.0e-5  )
            {
//				std::cout << " pvpR_reduced[" << i <<"] = " << pvpR_reduced[j] << std::endl;
            }
        }
*/

    }

    //----------------------
    // Print the pvp matrix
    //----------------------
/*
    std::cout << " pvp matrix:" << std::endl;
    for(int i=0; i<GlobalC::GridT.lgd; i++)
    {
        for(int j=0; j<GlobalC::GridT.lgd; j++)
        {
            std::cout << std::setw(15) << pvp[i][j].real();
        }
        std::cout << std::endl;
    }
    */

    // Distribution of data.
    ModuleBase::timer::tick("Gint_k","Distri");
    std::complex<double>* tmp;
    for (int i=0; i<GlobalV::NLOCAL; i++)
    {
        tmp = new std::complex<double>[GlobalV::NLOCAL];
        ModuleBase::GlobalFunc::ZEROS(tmp, GlobalV::NLOCAL);
        const int mug = GlobalC::GridT.trace_lo[i];
        const int mug0 = mug/GlobalV::NPOL;
        // if the row element is on this processor.
        if (mug >= 0)
        {
            for (int j=0; j<GlobalV::NLOCAL; j++)
            {
                const int nug = GlobalC::GridT.trace_lo[j];
                const int nug0 = nug/GlobalV::NPOL;
                // if the col element is on this processor.
                if (nug >=0)
                {
                    if (mug <= nug)
                    {
                        // pvp is symmetric, only half is calculated.
                        //int spin=0;
                        if(i%2==0&&j%2==0)
                        {
                            //spin = 0;
                            tmp[j] = this->pvp_nc[0][mug0][nug0]+this->pvp_nc[3][mug0][nug0];
                        }	
                        else if(i%2==1&&j%2==1)
                        {
                            //spin = 3;
                            tmp[j] = this->pvp_nc[0][mug0][nug0]-this->pvp_nc[3][mug0][nug0];
                        }
                        else if(i%2==0&&j%2==1)
                        {
                            // spin = 1;
                            if(!GlobalV::DOMAG) tmp[j] = 0;
                            else tmp[j] = this->pvp_nc[1][mug0][nug0] - std::complex<double>(0.0,1.0) * this->pvp_nc[2][mug0][nug0];
                        }
                        else if(i%2==1&&j%2==0) 
                        {
                            //spin = 2;
                            if(!GlobalV::DOMAG) tmp[j] = 0;
                            else tmp[j] = this->pvp_nc[1][mug0][nug0] + std::complex<double>(0.0,1.0) * this->pvp_nc[2][mug0][nug0];
                        }
                        else
                        {
                            ModuleBase::WARNING_QUIT("Gint_k::folding_vl_k_nc","index is wrong!");
                        }
                        //tmp[j] = this->pvp[spin][mug][nug];
                    }
                    else
                    {
                        // need to get elements from the other half.
                        // I have question on this! 2011-02-22
                        if(i%2==0&&j%2==0)
                        {
                            //spin = 0;
                            tmp[j] = conj(this->pvp_nc[0][nug0][mug0]+this->pvp_nc[3][nug0][mug0]);
                        }	
                        else if(i%2==1&&j%2==1)
                        {
                            //spin = 3;
                            tmp[j] = conj(this->pvp_nc[0][nug0][mug0]-this->pvp_nc[3][nug0][mug0]);
                        }
                        else if(i%2==1&&j%2==0)
                        {
                            // spin = 1;
                            if(!GlobalV::DOMAG) tmp[j] = 0;
                            else tmp[j] = conj(this->pvp_nc[1][nug0][mug0] - std::complex<double>(0.0,1.0) * this->pvp_nc[2][nug0][mug0]);
                        }
                        else if(i%2==0&&j%2==1) 
                        {
                            //spin = 2;
                            if(!GlobalV::DOMAG) tmp[j] = 0;
                            else tmp[j] = conj(this->pvp_nc[1][nug0][mug0] + std::complex<double>(0.0,1.0) * this->pvp_nc[2][nug0][mug0]);
                        }
                        else
                        {
                            ModuleBase::WARNING_QUIT("Gint_k::folding_vl_k_nc","index is wrong!");
                        }
                        //tmp[j] = conj(this->pvp[spin][nug][mug]);
                        
                    }
                }
            }
        }
        // collect the matrix after folding.
        Parallel_Reduce::reduce_complex_double_pool( tmp, GlobalV::NLOCAL );

        //-----------------------------------------------------
        // NOW! Redistribute the Hamiltonian matrix elements
        // according to the HPSEPS's 2D distribution methods.
        //-----------------------------------------------------
        for (int j=0; j<GlobalV::NLOCAL; j++)
        {
            if (!this->LM->ParaV->in_this_processor(i, j))
            {
                continue;
            }
            // set the matrix value.
            this->LM->set_HSk(i,j,tmp[j],'L');
        }
        delete[] tmp;
    }

    // delete the tmp matrix.
    for(int spin =0;spin<4;spin++)
    {
        for(int i=0; i<GlobalC::GridT.lgd; i++)
        {
            delete[] this->pvp_nc[spin][i];
        }
        delete[] this->pvp_nc[spin];
    }
    ModuleBase::timer::tick("Gint_k","Distri");

    ModuleBase::timer::tick("Gint_k","folding_vl_k_nc");
    return;
}

void Gint_k::set_ijk_atom(
    const int &grid_index, 
    const int &size,
    double*** psir_ylm, 
    double*** dr, 
    bool** cal_flag, 
    double** distance,
    const double &delta_r)
{
    const Numerical_Orbital_Lm* pointer;
    double mt[3];
    for (int id=0; id<size; id++)
    {
        // (2.1) get the atom type and atom index.
        const int mcell_index = GlobalC::GridT.bcell_start[grid_index] + id;	
        const int imcell = GlobalC::GridT.which_bigcell[mcell_index];
        const int iat = GlobalC::GridT.which_atom[mcell_index];
        const int it = GlobalC::ucell.iat2it[ iat ];
        const int ia = GlobalC::ucell.iat2ia[ iat ];

        // (2.2) get the distance between the grid and the atom.
        mt[0] = GlobalC::GridT.meshball_positions[imcell][0] - GlobalC::GridT.tau_in_bigcell[iat][0];
        mt[1] = GlobalC::GridT.meshball_positions[imcell][1] - GlobalC::GridT.tau_in_bigcell[iat][1];
        mt[2] = GlobalC::GridT.meshball_positions[imcell][2] - GlobalC::GridT.tau_in_bigcell[iat][2];

        for(int ib=0; ib<GlobalC::pw.bxyz; ib++)
        {
            // meshcell_pos: z is the fastest
            dr[ib][id][0] = GlobalC::GridT.meshcell_pos[ib][0] + mt[0];
            dr[ib][id][1] = GlobalC::GridT.meshcell_pos[ib][1] + mt[1];
            dr[ib][id][2] = GlobalC::GridT.meshcell_pos[ib][2] + mt[2];

            distance[ib][id] = std::sqrt(dr[ib][id][0]*dr[ib][id][0] 
            + dr[ib][id][1]*dr[ib][id][1] + dr[ib][id][2]*dr[ib][id][2]);

            //if(distance[ib][id] <= GlobalC::ORB.Phi[it].getRcut())
            // mohan reset this 2012-06-27
            // if(distance[ib][id] < GlobalC::ORB.Phi[it].getRcut() )
            // mohan reset this 2013-07-02 in Princeton
            // we should make absolutely sure that the distance is smaller than GlobalC::ORB.Phi[it].getRcut
            // this should be consistant with LCAO_nnr::cal_nnrg function 
            // typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
            // distance = 7.0000000000000000
            // GlobalC::ORB.Phi[it].getRcut = 7.0000000000000008
            if(distance[ib][id] < (GlobalC::ORB.Phi[it].getRcut() - 1.0e-15) )
            {
                cal_flag[ib][id]=true;
            }
            else
            {
                cal_flag[ib][id]=false;
                continue;
            }

            if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;

            std::vector<double> ylma;
            ModuleBase::Ylm::sph_harm ( GlobalC::ucell.atoms[it].nwl,
                    dr[ib][id][0] / distance[ib][id],
                    dr[ib][id][1] / distance[ib][id],
                    dr[ib][id][2] / distance[ib][id],
                    ylma);

            // these parameters are about interpolation
            // because once we know the distance from atom to grid point,
            // we can get the parameters we need to do interpolation and
            // store them first!! these can save a lot of effort.
            const double position = distance[ib][id] / delta_r;

            const int ip = static_cast<int>(position);

            const double dx = position - ip;
            const double dx2 = dx * dx;
            const double dx3 = dx2 * dx;

            const double c3 = 3.0*dx2-2.0*dx3;
            const double c1 = 1.0-c3;
            const double c2 = (dx-2.0*dx2+dx3)*delta_r;
            const double c4 = (dx3-dx2)*delta_r;

            Atom* atom1 = &GlobalC::ucell.atoms[it];
            double tmp=0.0;//mohan fix bug 2011-05-04
            for (int iw=0; iw< atom1->nw; iw++)
            {
                if ( atom1->iw2_new[iw] )
                {
                    pointer = &GlobalC::ORB.Phi[it].PhiLN(
                            atom1->iw2l[iw],
                            atom1->iw2n[iw]);

                    // Efficient!! to get the orbital value at this point.
                    tmp = c1*pointer->psi_uniform[ip] 
                        + c2*pointer->dpsi_uniform[ip]
                        + c3*pointer->psi_uniform[ip+1] 
                        + c4*pointer->dpsi_uniform[ip+1];
                }
                psir_ylm[ib][id][iw] = tmp * ylma[atom1->iw2_ylm[iw]];
            }// end iw.
        }//end ib
    }// int id
    return;
}




void Gint_k::allocate_pvpR_tr(void)
{
    ModuleBase::TITLE("Gint_k","allocate_pvpR_tr");

    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

//std::cout<<"R_x: "<<R_x<<std::endl;
//std::cout<<"R_y: "<<R_y<<std::endl;
//std::cout<<"R_z: "<<R_z<<std::endl;
//std::cout<<"GlobalC::GridT.lgd: "<<GlobalC::GridT.lgd<<std::endl;
    if(GlobalV::NSPIN!=4)
    {
        pvpR_tr = new double****[R_x];
        for(int ix=0; ix<R_x; ix++)
        {
            pvpR_tr[ix] = new double***[R_y];
            for(int iy=0; iy<R_y; iy++)
            {
                pvpR_tr[ix][iy] = new double**[R_z];
                for(int iz=0; iz<R_z; iz++)
                {
                    pvpR_tr[ix][iy][iz] = new double*[GlobalC::GridT.lgd];
                    for(int iw=0; iw<GlobalC::GridT.lgd; iw++)
                    {
                        pvpR_tr[ix][iy][iz][iw] = new double[GlobalC::GridT.lgd];
//do    uble mem = Memory::record("allocate_pvpR_tr", "pvpR_tr[ix][iy][iz][iw]", GlobalC::GridT.lgd , "double");
//co    ut<<" Memory of pvpR_tr[ix][iy][iz][iw]: "<<mem<<" MB"<<std::endl;
                        ModuleBase::GlobalFunc::ZEROS(pvpR_tr[ix][iy][iz][iw], GlobalC::GridT.lgd);
                    }
                }
            }
        }
    }
    else
    {
        pvpR_tr_soc = new std::complex<double>****[R_x];
        for(int ix=0; ix<R_x; ix++)
        {
            pvpR_tr_soc[ix] = new std::complex<double>***[R_y];
            for(int iy=0; iy<R_y; iy++)
            {
                pvpR_tr_soc[ix][iy] = new std::complex<double>**[R_z];
                for(int iz=0; iz<R_z; iz++)
                {
                    pvpR_tr_soc[ix][iy][iz] = new std::complex<double>*[GlobalC::GridT.lgd];
                    for(int iw=0; iw<GlobalC::GridT.lgd; iw++)
                    {
                        pvpR_tr_soc[ix][iy][iz][iw] = new std::complex<double>[GlobalC::GridT.lgd];
//do    uble mem = Memory::record("allocate_pvpR_tr", "pvpR_tr[ix][iy][iz][iw]", GlobalC::GridT.lgd , "double");
//co    ut<<" Memory of pvpR_tr[ix][iy][iz][iw]: "<<mem<<" MB"<<std::endl;
                        ModuleBase::GlobalFunc::ZEROS(pvpR_tr_soc[ix][iy][iz][iw], GlobalC::GridT.lgd);
                    }
                }
            }
        }
    }

    return;
}

void Gint_k::destroy_pvpR_tr(void)
{
    ModuleBase::TITLE("Gint_k","destroy_pvpR_tr");

    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    if(GlobalV::NSPIN!=4)
    {
        for(int ix=0; ix<R_x; ix++)
        {
            for(int iy=0; iy<R_y; iy++)
            {
                for(int iz=0; iz<R_z; iz++)
                {
                    for(int iw=0; iw<GlobalC::GridT.lgd; GlobalC::GridT.lgd++)
                    {
                        delete[] pvpR_tr[ix][iy][iz][iw];
                    }
                    delete[] pvpR_tr[ix][iy][iz];
                }
                delete[] pvpR_tr[ix][iy];
            }
            delete[] pvpR_tr[ix];
        }
        delete[] pvpR_tr;
    }
    else
    {
        for(int ix=0; ix<R_x; ix++)
        {
            for(int iy=0; iy<R_y; iy++)
            {
                for(int iz=0; iz<R_z; iz++)
                {
                    for(int iw=0; iw<GlobalC::GridT.lgd; GlobalC::GridT.lgd++)
                    {
                        delete[] pvpR_tr_soc[ix][iy][iz][iw];
                    }
                    delete[] pvpR_tr_soc[ix][iy][iz];
                }
                delete[] pvpR_tr_soc[ix][iy];
            }
            delete[] pvpR_tr_soc[ix];
        }
        delete[] pvpR_tr_soc;
    }

    return;
}

void Gint_k::distribute_pvpR_tr(void)
{
    ModuleBase::TITLE("Gint_k","distribute_pvpR_tr");

    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    for(int ix=0; ix<R_x; ix++)
    {
        for(int iy=0; iy<R_y; iy++)
        {
            for(int iz=0; iz<R_z; iz++)
            {
                double* tmp;
                std::complex<double>* tmp_soc;
                for(int i=0; i<GlobalV::NLOCAL; i++)
                {
                    if(GlobalV::NSPIN!=4)
                    {
                        tmp = new double[GlobalV::NLOCAL];
                        ModuleBase::GlobalFunc::ZEROS(tmp, GlobalV::NLOCAL);
                    }
                    else
                    {
                        tmp_soc = new std::complex<double>[GlobalV::NLOCAL];
                        ModuleBase::GlobalFunc::ZEROS(tmp_soc, GlobalV::NLOCAL);
                    }

                    const int mug = GlobalC::GridT.trace_lo[i];
//std::cout<<"mug: "<<mug<<std::endl;
                    if(mug >= 0)
                    {
                        for(int j=0; j<GlobalV::NLOCAL; j++)
                        {
                            const int nug = GlobalC::GridT.trace_lo[j];
//std::cout<<"nug: "<<nug<<std::endl;
                            if(nug >= 0)
                            {
                                //if(mug <= nug)
                                //{
//std::cout<<"ix: "<<ix<<std::endl;
//std::cout<<"iy: "<<iy<<std::endl;
//std::cout<<"iz: "<<iz<<std::endl;
//std::cout<<"mug: "<<mug<<std::endl;
//std::cout<<"nug: "<<nug<<std::endl;
//std::cout<<"pvpR_tr: "<<pvpR_tr[ix][iy][iz][mug][nug]<<std::endl;
                                    if(GlobalV::NSPIN!=4) tmp[j] = pvpR_tr[ix][iy][iz][mug][nug];
                                    else tmp_soc[j] = pvpR_tr_soc[ix][iy][iz][mug][nug];
//std::cout<<"tmp["<<j<<"]: "<<tmp[j]<<std::endl;
                                //}
                                //else
                                //{
                                    //tmp[j] = pvpR_tr[ix][iy][iz][nug][mug];
//std::cout<<"tmp["<<j<<"]: "<<tmp[j]<<std::endl;
                                //}
                            }
                        }
                    }
                    // collect the matrix after folding.
                    if(GlobalV::NSPIN!=4) Parallel_Reduce::reduce_double_pool( tmp, GlobalV::NLOCAL );
                    else Parallel_Reduce::reduce_complex_double_pool( tmp_soc, GlobalV::NLOCAL );
                    for(int j=0; j<GlobalV::NLOCAL; j++)
                    {
                        if(!this->LM->ParaV->in_this_processor(i,j))
                        {
                            continue;
                        }
                        else
                        {
                            //this->LM->set_HSk(i,j,tmp[j],'L');
                            if(GlobalV::NSPIN!=4) this->LM->set_HR_tr(ix,iy,iz,i,j,tmp[j]);
                            else this->LM->set_HR_tr_soc(ix,iy,iz,i,j,tmp_soc[j]);
                        }
                    }
                    if(GlobalV::NSPIN!=4) delete[] tmp;
                    else delete[] tmp_soc;
                }
            }
        }
    }

    return;
}


void Gint_k::cal_vlocal_R(const int current_spin)
{
    ModuleBase::TITLE("Gint_k","cal_vlocal_R");

    allocate_pvpR_tr();

    int lgd = 0;

    double R_minX = GlobalC::GridD.getD_minX();
    double R_minY = GlobalC::GridD.getD_minY();
    double R_minZ = GlobalC::GridD.getD_minZ();

    int R_x;
    int R_y;
    int R_z;

    ModuleBase::Vector3<double> tau1, dtau, dR;
    for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
    {
        for(int I1=0; I1<GlobalC::ucell.atoms[T1].na; ++I1)
        {
            const int iat = GlobalC::ucell.itia2iat(T1,I1);
            if(GlobalC::GridT.in_this_processor[iat])
            {
                Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

                const int DM_start = GlobalC::GridT.nlocstartg[iat];
                tau1 = GlobalC::ucell.atoms[T1].tau[I1];
                //GlobalC::GridD.Find_atom(tau1);        
                GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
                int nad2 = 0;

                for(int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ad++)
                {
                    const int T2 = GlobalC::GridD.getType(ad);
                    const int I2 = GlobalC::GridD.getNatom(ad);
                    const int iat2 = GlobalC::ucell.itia2iat(T2, I2);

                    if(GlobalC::GridT.in_this_processor[iat2])
                    {
                        Atom* atom2 = &GlobalC::ucell.atoms[T2];
                        dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                        double distance = dtau.norm() * GlobalC::ucell.lat0;
                        double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                        if(distance < rcut)
                        {
                            const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

                            dR.x = GlobalC::GridD.getBox(ad).x;
                            dR.y = GlobalC::GridD.getBox(ad).y;
                            dR.z = GlobalC::GridD.getBox(ad).z;

                            R_x = (int) (dR.x -R_minX);
                            R_y = (int) (dR.y -R_minY);
                            R_z = (int) (dR.z -R_minZ);

                            int ixxx = DM_start + GlobalC::GridT.find_R2st[iat][nad2];
                            for(int iw=0; iw<atom1->nw * GlobalV::NPOL; iw++)
                            {
                                int* iw2_lo = &GlobalC::GridT.trace_lo[start2];
                                int* iw2_end = iw2_lo + atom2->nw;

                                //double *vijR = &pvpR_reduced[ixxx];
                                double *vijR = &pvpR_reduced[current_spin][ixxx];
                                for(int iw2=0;iw2<atom2->nw * GlobalV::NPOL; iw2++)
                                {
                                    double *HlocR;
                                    std::complex<double> *HlocR_soc;
                                    if(GlobalV::NSPIN!=4) HlocR = &pvpR_tr[R_x][R_y][R_z][GlobalC::GridT.trace_lo[start1+iw]][GlobalC::GridT.trace_lo[start2+iw2]];
                                    else    HlocR_soc = &pvpR_tr_soc[R_x][R_y][R_z][GlobalC::GridT.trace_lo[start1+iw]][GlobalC::GridT.trace_lo[start2+iw2]];
                                    const int nw = atom2->nw;
                                    const int mug0 = iw/GlobalV::NPOL;
                                    const int nug0 = iw2/GlobalV::NPOL;
                                    const int iw_nowg = ixxx + mug0*nw + nug0;
                                    const int iw_nowg1 = ixxx + nug0*nw + mug0;
                                    if(GlobalV::NSPIN==4)
                                    {
                                            
                                            // if the col element is on this processor.
                                            
                                                // pvp is symmetric, only half is calculated.
                                                //int spin=0;
                                                if(iw%2==0&&iw2%2==0)
                                                {
                                                    //spin = 0;
                                                    HlocR_soc[0] = std::complex<double>(1.0,0.0) * pvpR_reduced[0][iw_nowg] + std::complex<double>(1.0,0.0) * pvpR_reduced[3][iw_nowg];
                                                }	
                                                else if(iw%2==1&&iw2%2==1)
                                                {
                                                    //spin = 3;
                                                    HlocR_soc[0] = std::complex<double>(1.0,0.0) * pvpR_reduced[0][iw_nowg] - std::complex<double>(1.0,0.0) * pvpR_reduced[3][iw_nowg];
                                                }
                                                else if(iw%2==0&&iw2%2==1)
                                                {
                                                    // spin = 1;
                                                    if(!GlobalV::DOMAG) HlocR_soc[0] = std::complex<double>(0.0,0.0);
                                                    else HlocR_soc[0] = pvpR_reduced[1][iw_nowg] - std::complex<double>(0.0,1.0) * pvpR_reduced[2][iw_nowg];
                                                }	
                                                else if(iw%2==1&&iw2%2==0) 
                                                {
                                                    //spin = 2;
                                                    if(!GlobalV::DOMAG) HlocR_soc[0] = std::complex<double>(0.0,0.0);
                                                    else HlocR_soc[0] = pvpR_reduced[1][iw_nowg] + std::complex<double>(0.0,1.0) * pvpR_reduced[2][iw_nowg];
                                                }
                                                else
                                                {
                                                    ModuleBase::WARNING_QUIT("Gint_k::folding_vl_k_nc","index is wrong!");
                                                }
                                    }//endif NC
                                    else
                                    {
                                        HlocR[0] = pvpR_reduced[current_spin][iw_nowg];
                                            //pvpR_tr[R_x][R_y][R_z][GlobalC::GridT.trace_lo[start1+iw]][iw2_lo[0]] = vijR[0];
                                    }//endif normal
                                }
//                               for(; iw2_lo<iw2_end; ++iw2_lo, ++vijR)
//                               {
//                                   pvpR_tr[R_x][R_y][R_z][GlobalC::GridT.trace_lo[start1+iw]][iw2_lo[0]] = vijR[0];
//                               }
//                               ixxx += atom2->nw;
                                ++lgd;
                            }
                            ++nad2;
                        }
                    }
                }
            }
        }
    }

    return;
}

void Gint_k::distribute_pvpR_sparseMatrix(
    const int current_spin, 
    const double &sparse_threshold, 
    const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> &pvpR_sparseMatrix
)
{
    ModuleBase::TITLE("Gint_k","distribute_pvpR_sparseMatrix");

    int total_R_num = this->LM->all_R_coor.size();
    int *nonzero_num = new int[total_R_num];
    int *minus_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    ModuleBase::GlobalFunc::ZEROS(minus_nonzero_num, total_R_num);
    int count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        auto iter = pvpR_sparseMatrix.find(R_coor);
        if (iter != pvpR_sparseMatrix.end())
        {
            for (auto &row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }

        auto minus_R_coor = -1 * R_coor;

        iter = pvpR_sparseMatrix.find(minus_R_coor);
        if (iter != pvpR_sparseMatrix.end())
        {
            for (auto &row_loop : iter->second)
            {
                minus_nonzero_num[count] += row_loop.second.size();
            }
        }
        
        count++;
    }

    Parallel_Reduce::reduce_int_all(nonzero_num, total_R_num);
    Parallel_Reduce::reduce_int_all(minus_nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_int_pool(nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_int_pool(minus_nonzero_num, total_R_num);

    double* tmp = nullptr;
    tmp = new double[GlobalV::NLOCAL];

    count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        if (nonzero_num[count] != 0 || minus_nonzero_num[count] != 0)
        {
            auto minus_R_coor = -1 * R_coor;

            for(int row = 0; row < GlobalV::NLOCAL; ++row)
            {        
                ModuleBase::GlobalFunc::ZEROS(tmp, GlobalV::NLOCAL);
                
                auto iter = pvpR_sparseMatrix.find(R_coor);
                if (iter != pvpR_sparseMatrix.end())
                {
                    
                    if(GlobalC::GridT.trace_lo[row] >= 0)
                    {
                        auto row_iter = iter->second.find(row);
                        if (row_iter != iter->second.end())
                        {
                            for (auto &value : row_iter->second)
                            {
                                tmp[value.first] = value.second;
                            }
                        }
                    }
                }

                auto minus_R_iter = pvpR_sparseMatrix.find(minus_R_coor);
                if (minus_R_iter != pvpR_sparseMatrix.end())
                {
                    for (int col = 0; col < row; ++col)
                    {
                        if(GlobalC::GridT.trace_lo[col] >= 0)
                        {
                            auto row_iter = minus_R_iter->second.find(col);
                            if (row_iter != minus_R_iter->second.end())
                            {
                                auto col_iter = row_iter->second.find(row);
                                if (col_iter != row_iter->second.end())
                                {
                                    tmp[col] = col_iter->second;
                                }

                            }
                        }
                    }
                }
                
                Parallel_Reduce::reduce_double_pool(tmp, GlobalV::NLOCAL);

                if (this->LM->ParaV->trace_loc_row[row] >= 0)
                {
                    for(int col = 0; col < GlobalV::NLOCAL; ++col)
                    {
                        if(this->LM->ParaV->trace_loc_col[col] >= 0)
                        {
                            if (std::abs(tmp[col]) > sparse_threshold)
                            {
                                double &value = this->LM->HR_sparse[current_spin][R_coor][row][col];
                                value += tmp[col];
                                if (std::abs(value) <= sparse_threshold)
                                {
                                    this->LM->HR_sparse[current_spin][R_coor][row].erase(col);
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
    delete[] minus_nonzero_num;
    delete[] tmp;
    nonzero_num = nullptr;
    minus_nonzero_num = nullptr;
    tmp = nullptr;

    return;

}

void Gint_k::distribute_pvpR_soc_sparseMatrix(
    const double &sparse_threshold, 
    const std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> &pvpR_soc_sparseMatrix
)
{
    ModuleBase::TITLE("Gint_k","distribute_pvpR_soc_sparseMatrix");

    int total_R_num = this->LM->all_R_coor.size();
    int *nonzero_num = new int[total_R_num];
    int *minus_nonzero_num = new int[total_R_num];
    ModuleBase::GlobalFunc::ZEROS(nonzero_num, total_R_num);
    ModuleBase::GlobalFunc::ZEROS(minus_nonzero_num, total_R_num);
    int count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        auto iter = pvpR_soc_sparseMatrix.find(R_coor);
        if (iter != pvpR_soc_sparseMatrix.end())
        {
            for (auto &row_loop : iter->second)
            {
                nonzero_num[count] += row_loop.second.size();
            }
        }

        auto minus_R_coor = -1 * R_coor;

        iter = pvpR_soc_sparseMatrix.find(minus_R_coor);
        if (iter != pvpR_soc_sparseMatrix.end())
        {
            for (auto &row_loop : iter->second)
            {
                minus_nonzero_num[count] += row_loop.second.size();
            }
        }
        
        count++;
    }

    Parallel_Reduce::reduce_int_all(nonzero_num, total_R_num);
    Parallel_Reduce::reduce_int_all(minus_nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_int_pool(nonzero_num, total_R_num);
    // Parallel_Reduce::reduce_int_pool(minus_nonzero_num, total_R_num);

    std::complex<double>* tmp_soc = nullptr;
    tmp_soc = new std::complex<double>[GlobalV::NLOCAL];

    count = 0;
    for (auto &R_coor : this->LM->all_R_coor)
    {
        if (nonzero_num[count] != 0 || minus_nonzero_num[count] != 0)
        {
            auto minus_R_coor = -1 * R_coor;

            for(int row = 0; row < GlobalV::NLOCAL; ++row)
            {        
                ModuleBase::GlobalFunc::ZEROS(tmp_soc, GlobalV::NLOCAL);
                
                auto iter = pvpR_soc_sparseMatrix.find(R_coor);
                if (iter != pvpR_soc_sparseMatrix.end())
                {
                    if(GlobalC::GridT.trace_lo[row] >= 0)
                    {
                        auto row_iter = iter->second.find(row);
                        if (row_iter != iter->second.end())
                        {
                            for (auto &value : row_iter->second)
                            {
                                tmp_soc[value.first] = value.second;
                            }
                        }
                    }
                }

                auto minus_R_iter = pvpR_soc_sparseMatrix.find(minus_R_coor);
                if (minus_R_iter != pvpR_soc_sparseMatrix.end())
                {
                    for (int col = 0; col < row; ++col)
                    {
                        if(GlobalC::GridT.trace_lo[col] >= 0)
                        {
                            auto row_iter = minus_R_iter->second.find(col);
                            if (row_iter != minus_R_iter->second.end())
                            {
                                auto col_iter = row_iter->second.find(row);
                                if (col_iter != row_iter->second.end())
                                {
                                    tmp_soc[col] = conj(col_iter->second);
                                }

                            }
                        }
                    }
                }
                
                Parallel_Reduce::reduce_complex_double_pool(tmp_soc, GlobalV::NLOCAL);

                if (this->LM->ParaV->trace_loc_row[row] >= 0)
                {
                    for(int col = 0; col < GlobalV::NLOCAL; ++col)
                    {
                        if(this->LM->ParaV->trace_loc_col[col] >= 0)
                        {
                            if (std::abs(tmp_soc[col]) > sparse_threshold)
                            {
                                std::complex<double> &value = this->LM->HR_soc_sparse[R_coor][row][col];
                                value += tmp_soc[col];
                                if (std::abs(value) <= sparse_threshold)
                                {
                                    this->LM->HR_soc_sparse[R_coor][row].erase(col);
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
    delete[] minus_nonzero_num;
    delete[] tmp_soc;
    nonzero_num = nullptr;
    minus_nonzero_num = nullptr;
    tmp_soc = nullptr;

    return;

}

void Gint_k::cal_vlocal_R_sparseMatrix(const int &current_spin, const double &sparse_threshold)
{
    ModuleBase::TITLE("Gint_k","cal_vlocal_R_sparseMatrix");

    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> pvpR_sparseMatrix;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> pvpR_soc_sparseMatrix;

    int lgd = 0;
    double temp_value_double;
    std::complex<double> temp_value_complex;

    ModuleBase::Vector3<double> tau1, dtau;
    for(int T1=0; T1<GlobalC::ucell.ntype; ++T1)
    {
        for(int I1=0; I1<GlobalC::ucell.atoms[T1].na; ++I1)
        {
            const int iat = GlobalC::ucell.itia2iat(T1,I1);
            if(GlobalC::GridT.in_this_processor[iat])
            {
                Atom* atom1 = &GlobalC::ucell.atoms[T1];
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

                const int DM_start = GlobalC::GridT.nlocstartg[iat];
                tau1 = GlobalC::ucell.atoms[T1].tau[I1];
                GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
                int nad2 = 0;

                for(int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ad++)
                {
                    const int T2 = GlobalC::GridD.getType(ad);
                    const int I2 = GlobalC::GridD.getNatom(ad);
                    const int iat2 = GlobalC::ucell.itia2iat(T2, I2);

                    if(GlobalC::GridT.in_this_processor[iat2])
                    {
                        Atom* atom2 = &GlobalC::ucell.atoms[T2];
                        dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                        double distance = dtau.norm() * GlobalC::ucell.lat0;
                        double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                        if(distance < rcut)
                        {
                            const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                            Abfs::Vector3_Order<int> dR(GlobalC::GridD.getBox(ad).x, GlobalC::GridD.getBox(ad).y, GlobalC::GridD.getBox(ad).z);
                            int ixxx = DM_start + GlobalC::GridT.find_R2st[iat][nad2];
                            for(int iw=0; iw<atom1->nw * GlobalV::NPOL; iw++)
                            {
                                for(int iw2=0;iw2<atom2->nw * GlobalV::NPOL; iw2++)
                                {
                                    const int nw = atom2->nw;
                                    const int mug0 = iw/GlobalV::NPOL;
                                    const int nug0 = iw2/GlobalV::NPOL;
                                    const int iw_nowg = ixxx + mug0*nw + nug0;

                                    if(GlobalV::NSPIN == 4)
                                    {		
                                        // pvp is symmetric, only half is calculated.

                                        if(iw%2==0&&iw2%2==0)
                                        {
                                            //spin = 0;
                                            temp_value_complex = std::complex<double>(1.0,0.0) * pvpR_reduced[0][iw_nowg] + std::complex<double>(1.0,0.0) * pvpR_reduced[3][iw_nowg];
                                            if(std::abs(temp_value_complex) > sparse_threshold)
                                            {
                                                pvpR_soc_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_complex;
                                            }
                                        }	
                                        else if(iw%2==1&&iw2%2==1)
                                        {
                                            //spin = 3;
                                            temp_value_complex = std::complex<double>(1.0,0.0) * pvpR_reduced[0][iw_nowg] - std::complex<double>(1.0,0.0) * pvpR_reduced[3][iw_nowg];
                                            if(std::abs(temp_value_complex) > sparse_threshold)
                                            {
                                                pvpR_soc_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_complex;
                                            }
                                        }
                                        else if(iw%2==0&&iw2%2==1)
                                        {
                                            // spin = 1;
                                            if(!GlobalV::DOMAG)
                                            {
                                                // do nothing
                                            }
                                            else
                                            {
                                                temp_value_complex = pvpR_reduced[1][iw_nowg] - std::complex<double>(0.0,1.0) * pvpR_reduced[2][iw_nowg];
                                                if(std::abs(temp_value_complex) > sparse_threshold)
                                                {
                                                    pvpR_soc_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_complex;
                                                }
                                            }
                                        }	
                                        else if(iw%2==1&&iw2%2==0) 
                                        {
                                            //spin = 2;
                                            if(!GlobalV::DOMAG)
                                            {
                                                // do nothing
                                            }
                                            else
                                            {
                                                temp_value_complex = pvpR_reduced[1][iw_nowg] + std::complex<double>(0.0,1.0) * pvpR_reduced[2][iw_nowg];
                                                if(std::abs(temp_value_complex) > sparse_threshold)
                                                {
                                                    pvpR_soc_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_complex;
                                                }
                                            }
                                        }
                                        else
                                        {
                                            ModuleBase::WARNING_QUIT("Gint_k::folding_vl_k_nc","index is wrong!");
                                        }
                                    } //endif NC
                                    else
                                    {
                                        temp_value_double = pvpR_reduced[current_spin][iw_nowg];
                                        if (std::abs(temp_value_double) > sparse_threshold)
                                        {
                                            pvpR_sparseMatrix[dR][start1 + iw][start2 + iw2] = temp_value_double;
                                        }

                                    } //endif normal

                                }

                                ++lgd;
                            }
                            ++nad2;
                        }
                    }
                }
            }
        }
    }

    if (GlobalV::NSPIN != 4)
    {
        distribute_pvpR_sparseMatrix(current_spin, sparse_threshold, pvpR_sparseMatrix);
    }
    else
    {
        distribute_pvpR_soc_sparseMatrix(sparse_threshold, pvpR_soc_sparseMatrix);
    }

    return;
}

