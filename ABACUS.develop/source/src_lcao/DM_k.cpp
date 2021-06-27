#include "local_orbital_charge.h"
#include "../src_pw/global.h"
#include "module_base/blas_connector.h"

#include "LCAO_nnr.h"
void Local_Orbital_Charge::allocate_DM_k(void)
{
    TITLE("Local_Orbital_Charge","allocate_k");

    this->nnrg_now = LNNR.nnrg;
    //xiaohui add 'OUT_LEVEL' line, 2015-09-16
    if(OUT_LEVEL != "m") OUT(ofs_running,"nnrg_last",nnrg_last);
    if(OUT_LEVEL != "m") OUT(ofs_running,"nnrg_now",nnrg_now);

    if(this->init_DM_R)
    {
        assert(nnrg_last > 0);
        for(int is=0; is<NSPIN; is++)
        {
            delete[] DM_R[is];
        }
        delete[] DM_R;
        init_DM_R=false;
    }

    if(nnrg_now>0)
    {
        this->DM_R = new double*[NSPIN];
        for(int is=0; is<NSPIN; is++)
        {
            this->DM_R[is] = new double[nnrg_now];
            ZEROS(DM_R[is], nnrg_now);
        }
        this->nnrg_last = nnrg_now;
        this->init_DM_R = true;
        Memory::record("LocalOrbital_Charge","Density_Matrix",NSPIN*nnrg_now,"double");
    }
    else if(nnrg_now==0)
    {
        this->init_DM_R = false;
    }
    else
    {
        WARNING_QUIT("Local_Orbital_Charge::allocate_k","check init_DM_R.");
    }

	// Peize Lin test 2019-01-16 
    wfc_dm_2d.init();

    return;
}



#include "record_adj.h"
inline void cal_DM_ATOM(const Grid_Technique &gt, const complex<double> fac, Record_adj RA,
                   const int ia1, const int iw1_lo, const int nw1, const int gstart, 
                   complex<double> *WFC_PHASE, complex<double> **DM_ATOM)
{

    const char transa='N', transb='T';  
    const complex<double> alpha=1, beta=1;

    for(int ik=0; ik<kv.nks; ik++)
    {
        complex<double> **wfc = LOWF.WFC_K[ik];
        const int ispin = kv.isk[ik];
        int atom2start=0;

        for (int ia2 = 0; ia2 < RA.na_each[ia1]; ++ia2)
        {
            complex<double> *DM=&DM_ATOM[ispin][atom2start];
            const int T2 = RA.info[ia1][ia2][3];
            const int I2 = RA.info[ia1][ia2][4];
            Atom* atom2 = &ucell.atoms[T2];
            const int start2 = ucell.itiaiw2iwt(T2,I2,0);
            const int iw2_lo=gt.trace_lo[start2];
            const int nw2=atom2->nw;
            complex<double> exp_R= exp( fac * (
                        kv.kvec_d[ik].x * RA.info[ia1][ia2][0] + 
                        kv.kvec_d[ik].y * RA.info[ia1][ia2][1] + 
                        kv.kvec_d[ik].z * RA.info[ia1][ia2][2]  
                        ) );
            
            //ZEROS(WFC_PHASE, NBANDS*nw1);
            int ibStart=0;
            int nRow=0;
            for(int ib=0; ib<NBANDS; ++ib)
            {
                const double wg_local=wf.wg(ik,ib);
                if(wg_local>0)
                {
                    if(nRow==0) ibStart=ib;
                    const int iline=nRow*nw1;
                    complex<double> phase=exp_R*wg_local;
                    for(int iw1=0; iw1<nw1; ++iw1)
					{
                        WFC_PHASE[iline+iw1]=phase*conj(wfc[ib][iw1_lo+iw1]);
					}
                    ++nRow;
                }
                else
				{
                    break;
				}
            } // ib
            zgemm_(&transa, &transb, &nw2, &nw1, &nRow, &alpha,
                &wfc[ibStart][iw2_lo], &gt.lgd, 
                WFC_PHASE, &nw1,
                &beta, DM, &nw2);           
            atom2start+=nw1*nw2;
        } // ia2
    } // ik
    return;
}

//added by zhengdy-soc, for non-collinear case
inline void cal_DM_ATOM_nc(const Grid_Technique &gt, const complex<double> fac, Record_adj RA,
                   const int ia1, const int iw1_lo, const int nw1, const int gstart, 
                   complex<double> *WFC_PHASE, complex<double> **DM_ATOM)
{
    if(NSPIN !=4 ) WARNING_QUIT("Local_Orbital_Charge","NSPIN not match!");
    const char transa='N', transb='T';  
    const complex<double> alpha=1, beta=1;
    int ispin=0;

    for(int is1=0;is1<2;is1++)
    {
        for(int is2=0;is2<2;is2++)
        {
            for(int ik=0; ik<kv.nks; ik++)
            {
                complex<double> **wfc = LOWF.WFC_K[ik];
                int atom2start=0;

                for (int ia2 = 0; ia2 < RA.na_each[ia1]; ++ia2)
                {
                    complex<double> *DM=&DM_ATOM[ispin][atom2start];
                    const int T2 = RA.info[ia1][ia2][3];
                    const int I2 = RA.info[ia1][ia2][4];
                    Atom* atom2 = &ucell.atoms[T2];
                    const int start2 = ucell.itiaiw2iwt(T2,I2,0);
                    const int iw2_lo=gt.trace_lo[start2]/NPOL + gt.lgd/NPOL*is2;
                    const int nw2=atom2->nw;
                    complex<double> exp_R= exp( fac * (
                                kv.kvec_d[ik].x * RA.info[ia1][ia2][0] + 
                                kv.kvec_d[ik].y * RA.info[ia1][ia2][1] + 
                                kv.kvec_d[ik].z * RA.info[ia1][ia2][2]  
                                ) );
            
            		//ZEROS(WFC_PHASE, NBANDS*nw1);
                    int ibStart=0;
                    int nRow=0;
                    for(int ib=0; ib<NBANDS; ++ib)
                    {
                        const double w1=wf.wg(ik,ib);
                        if(w1>0)
                        {
                            if(nRow==0) 
							{
								ibStart=ib;
							}
                            const int iline=nRow*nw1;
                            complex<double> phase=exp_R*w1;
                            for(int iw1=0; iw1<nw1; ++iw1)
							{
                                WFC_PHASE[iline+iw1]=phase*conj(wfc[ib][iw1_lo+iw1 + gt.lgd/NPOL*is1]);
							}
                            ++nRow;
                        }
                        else
                            break;
                    } // ib
                    zgemm_(&transa, &transb, &nw2, &nw1, &nRow, &alpha,
                        &wfc[ibStart][iw2_lo], &gt.lgd, 
                        WFC_PHASE, &nw1,
                        &beta, DM, &nw2);           

                    atom2start+=nw1*nw2;
                } // ia2
            } // ik
            ispin++;
        }//is2
    }//is1
    return;
}


void Local_Orbital_Charge::cal_dk_k(const Grid_Technique &gt)
{
    TITLE("Local_Orbital_Charge","cal_dk_k");
    timer::tick("LCAO_Charge","cal_dk_k",'F');  
    //int nnrg = 0;
    Vector3<double> tau1, dtau;
        
    Record_adj RA;
    RA.for_grid(gt);

    int ca = 0;
    complex<double> fac = TWO_PI * IMAG_UNIT;

    complex<double> *WFC_PHASE=new complex<double>[NLOCAL*ucell.nwmax];
    
    int DM_ATOM_SIZE=1; 
    complex<double> **DM_ATOM=new complex<double> *[NSPIN];
    for(int is=0; is<NSPIN; ++is)
    {
        DM_ATOM[is]=new complex<double>[DM_ATOM_SIZE];
        ZEROS(DM_ATOM[is], DM_ATOM_SIZE);
    }
    for(int T1=0; T1<ucell.ntype; T1++)
    {
        Atom* atom1 = &ucell.atoms[T1];
        for(int I1=0; I1<atom1->na; I1++)
        {
            const int iat = ucell.itia2iat(T1,I1);
            if(gt.in_this_processor[iat])
            {
                const int start1 = ucell.itiaiw2iwt(T1,I1,0);
                const int gstart = LNNR.nlocstartg[iat];
                const int ng = LNNR.nlocdimg[iat];
                const int iw1_lo=gt.trace_lo[start1]/NPOL;
                const int nw1=atom1->nw;

                if(DM_ATOM_SIZE<ng)
                {
                    DM_ATOM_SIZE=ng;
                    for(int is=0; is<NSPIN; ++is)
					{
                        delete[] DM_ATOM[is];
					}
                    for(int is=0; is<NSPIN; ++is)
					{
                        DM_ATOM[is]=new complex<double>[DM_ATOM_SIZE];
					}
                }
                for(int is=0; is<NSPIN; ++is)
				{
                    ZEROS(DM_ATOM[is], ng);
				}
                ZEROS(WFC_PHASE, NBANDS*nw1);
                if(NSPIN!=4)
				{
					cal_DM_ATOM(gt, fac, RA, ca, iw1_lo, nw1, gstart, WFC_PHASE, DM_ATOM);
				}
                else 
				{
					cal_DM_ATOM_nc(gt, fac, RA, ca, iw1_lo, nw1, gstart, WFC_PHASE, DM_ATOM);
				}
                ++ca;

                if(NSPIN!=4)
                {
                    for(int is=0; is<NSPIN; ++is)
                    {
                        for(int iv=0; iv<ng; ++iv)
                        {
                            this->DM_R[is][gstart+iv]=DM_ATOM[is][iv].real();
                        }
                    }
                }
                else
                {//zhengdy-soc
                    for(int iv=0; iv<ng; ++iv)
                    {
						//note: storage nondiagonal term as Re[] and Im[] respectly;
						this->DM_R[0][gstart+iv]=DM_ATOM[0][iv].real() + DM_ATOM[3][iv].real();
						if(NONCOLIN)
						{//DOMAG
							this->DM_R[1][gstart+iv]=DM_ATOM[1][iv].real() + DM_ATOM[2][iv].real();
							this->DM_R[2][gstart+iv]=DM_ATOM[1][iv].imag() - DM_ATOM[2][iv].imag();
							this->DM_R[3][gstart+iv]=DM_ATOM[0][iv].real() - DM_ATOM[3][iv].real();
						}
						else if(!NONCOLIN)//DOMAG_Z
						{
							this->DM_R[1][gstart+iv]= 0.0;
							this->DM_R[1][gstart+iv]= 0.0;
							this->DM_R[3][gstart+iv]=DM_ATOM[0][iv].real() - DM_ATOM[3][iv].real();
						}
						else//soc with no mag 
						{
							this->DM_R[1][gstart+iv]= 0.0;
							this->DM_R[2][gstart+iv]= 0.0;
							this->DM_R[3][gstart+iv]= 0.0;
						}
					}
				}
            } // if gt.in_this_processor
        }// I1
    }// T1


    //------------
    // for test
    //------------
/*  cout << setprecision(3);
    for(int i=0; i<nnrg_now; i++)

    for(int ik=0; ik<kv.nkstot; ++ik)
    {
        for(int ib=0; ib<NBANDS; ++ib)
        {
            cout << " ik=" << ik << " ib=" << ib << " occ=" << wf.wg(ik,ib) << " e=" << wf.ekb[ik][ib] << endl;
        }
    }

    for(int i=0; i<10; i++)
    {
        if(DM_R[0][i]>1.0e-8)
        {
            cout << " i=" << i << " DM_R=" << DM_R[0][i] << endl;
        }
    }
*/
    for(int i=0; i<NSPIN; ++i)
	{
        delete[] DM_ATOM[i];
	}
    delete[] DM_ATOM;
    delete[] WFC_PHASE;

    RA.delete_grid();//xiaohui add 2015-02-04

    timer::tick("LCAO_Charge","cal_dk_k",'F');  
    return;
}


