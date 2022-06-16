#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_reduce.h"
#include "gint_k.h"
#include "../module_orbital/ORB_read.h"
#include "grid_technique.h"
#include "../module_base/ylm.h"
#include "../src_pw/global.h"
#include "../src_lcao/global_fp.h" // mohan add 2021-01-30
#include "../module_base/memory.h"
#include "../module_base/timer.h"

void Gint_k::allocate_pvpR(void)
{
    ModuleBase::TITLE("Gint_k","allocate_pvpR");

    if(this->pvpR_alloc_flag)
    {
        return; //Liuxh add, 20181012
        ModuleBase::WARNING_QUIT("Gint_k::allocate_pvpR","pvpR has been allocated!");
    }

    //xiaohui modify 2015-05-30
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
void Gint_k::folding_vl_k(const int &ik, LCAO_Matrix *LM)
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

    std::complex<double>** pvp = new std::complex<double>*[GlobalC::GridT.lgd];
    for(int i=0; i<GlobalC::GridT.lgd; i++)
    {
        pvp[i] = new std::complex<double>[GlobalC::GridT.lgd];
        ModuleBase::GlobalFunc::ZEROS( pvp[i], GlobalC::GridT.lgd);
    }

    std::complex<double>*** pvp_nc;
    if(GlobalV::NSPIN==4)
    {
        pvp_nc=new std::complex<double>**[4];
        for(int spin=0;spin<4;spin++)
        {
            pvp_nc[spin] = new std::complex<double>*[GlobalC::GridT.lgd];
            for(int i=0; i<GlobalC::GridT.lgd; i++)
            {
                pvp_nc[spin][i] = new std::complex<double>[GlobalC::GridT.lgd];
                ModuleBase::GlobalFunc::ZEROS( pvp_nc[spin][i], GlobalC::GridT.lgd);
            }
        }

    }

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
                            const double arg = (GlobalC::kv.kvec_d[ ik ] * dR) * ModuleBase::TWO_PI;
                            std::complex<double> phase = std::complex<double>(cos(arg), sin(arg));
                            int ixxx = DM_start + GlobalC::GridT.find_R2st[iat][nad];
                            for(int iw=0; iw<atom1->nw; iw++)
                            {
                                // iw1_lo
                                if(GlobalV::NSPIN!=4)
                                {
                                    std::complex<double> *vij = pvp[GlobalC::GridT.trace_lo[start1+iw]];

                                    int* iw2_lo = &GlobalC::GridT.trace_lo[start2];
                                    int* iw2_end = iw2_lo + atom2->nw;

                                    // get the <phi | V | phi>(R) Hamiltonian.
                                    double *vijR = &pvpR_reduced[0][ixxx];
                                    for(; iw2_lo<iw2_end; ++iw2_lo, ++vijR)
                                    {
                                        vij[iw2_lo[0]] += vijR[0] * phase; 
                                    }
                                }
                                else
                                {
                                    std::complex<double> *vij[4];
                                    for(int spin=0;spin<4;spin++)
                                        vij[spin] = pvp_nc[spin][GlobalC::GridT.trace_lo[start1]/GlobalV::NPOL + iw];

                                    int iw2_lo = GlobalC::GridT.trace_lo[start2]/GlobalV::NPOL;
                                    int iw2_end = iw2_lo + atom2->nw;

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
                                    }                                    
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
                        if(GlobalV::NSPIN!=4)
                        {
                            // pvp is symmetric, only half is calculated.
                            tmp[j] = pvp[mug][nug];
                        }
                        else
                        {
                            if(i%2==0&&j%2==0)
                            {
                                //spin = 0;
                                tmp[j] = pvp_nc[0][mug0][nug0]+pvp_nc[3][mug0][nug0];
                            }	
                            else if(i%2==1&&j%2==1)
                            {
                                //spin = 3;
                                tmp[j] = pvp_nc[0][mug0][nug0]-pvp_nc[3][mug0][nug0];
                            }
                            else if(i%2==0&&j%2==1)
                            {
                                // spin = 1;
                                if(!GlobalV::DOMAG) tmp[j] = 0;
                                else tmp[j] = pvp_nc[1][mug0][nug0] - std::complex<double>(0.0,1.0) * pvp_nc[2][mug0][nug0];
                            }
                            else if(i%2==1&&j%2==0) 
                            {
                                //spin = 2;
                                if(!GlobalV::DOMAG) tmp[j] = 0;
                                else tmp[j] = pvp_nc[1][mug0][nug0] + std::complex<double>(0.0,1.0) * pvp_nc[2][mug0][nug0];
                            }
                            else
                            {
                                ModuleBase::WARNING_QUIT("Gint_k::folding_vl_k_nc","index is wrong!");
                            }                            
                        }
                    }
                    else
                    {
                        // need to get elements from the other half.
                        // I have question on this! 2011-02-22
                        if(GlobalV::NSPIN!=4)
                        {
                            tmp[j] = conj(pvp[nug][mug]);
                        }
                        else
                        {
                            if(i%2==0&&j%2==0)
                            {
                                //spin = 0;
                                tmp[j] = conj(pvp_nc[0][nug0][mug0]+pvp_nc[3][nug0][mug0]);
                            }	
                            else if(i%2==1&&j%2==1)
                            {
                                //spin = 3;
                                tmp[j] = conj(pvp_nc[0][nug0][mug0]-pvp_nc[3][nug0][mug0]);
                            }
                            else if(i%2==1&&j%2==0)
                            {
                                // spin = 1;
                                if(!GlobalV::DOMAG) tmp[j] = 0;
                                else tmp[j] = conj(pvp_nc[1][nug0][mug0] - std::complex<double>(0.0,1.0) * pvp_nc[2][nug0][mug0]);
                            }
                            else if(i%2==0&&j%2==1) 
                            {
                                //spin = 2;
                                if(!GlobalV::DOMAG) tmp[j] = 0;
                                else tmp[j] = conj(pvp_nc[1][nug0][mug0] + std::complex<double>(0.0,1.0) * pvp_nc[2][nug0][mug0]);
                            }
                            else
                            {
                                ModuleBase::WARNING_QUIT("Gint_k::folding_vl_k_nc","index is wrong!");
                            }                           
                        }
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
            if (!LM->ParaV->in_this_processor(i,j))
            {
                continue;
            }
            // set the matrix value.
            LM->set_HSk(i,j,tmp[j],'L');
        }
        delete[] tmp;
    }

    // delete the tmp matrix.
    for(int i=0; i<GlobalC::GridT.lgd; i++)
    {
        delete[] pvp[i];
    }
    delete[] pvp;
    if(GlobalV::NSPIN==4)
    {
        for(int spin =0;spin<4;spin++)
        {
            for(int i=0; i<GlobalC::GridT.lgd; i++)
            {
                delete[] pvp_nc[spin][i];
            }
            delete[] pvp_nc[spin];
        }
        delete[] pvp_nc;
    }
    ModuleBase::timer::tick("Gint_k","Distri");

    ModuleBase::timer::tick("Gint_k","folding_vl_k");
    return;
}
