#include "gint_k.h"
#include "../src_pw/global.h"
#include "LCAO_nnr.h"
#include "global_fp.h" // mohan add 2021-01-30

void Gint_k::fvl_k_RealSpace(matrix& fvl_dphi, const double *vl)
{
	TITLE("Gint_k","cal_force");
	timer::tick("Gint_k","cal_force");

	if(!this->reduced)
	{
		WARNING_QUIT("Gint_k::cal_force_k","The force with k can only with reduced H.");
	}

	int nnrg = LNNR.nnrg;
 	//xiaohui add "OUT_LEVEL", 2015-09-16
	if(OUT_LEVEL != "m") ofs_running << " LNNR.nnrg in cal_force_k = " << LNNR.nnrg << endl;
	assert(nnrg>=0);

	// just because to make thea arrys meaningful.
	if(LNNR.nnrg == 0)
	{
		nnrg = 1;
	}

	// to store < phi | vlocal | dphi>
	double* pvdpx = new double[nnrg];
	double* pvdpy = new double[nnrg];
	double* pvdpz = new double[nnrg];
	ZEROS(pvdpx, nnrg);
	ZEROS(pvdpy, nnrg);
	ZEROS(pvdpz, nnrg); 

    const double delta_r = ORB.dr_uniform;
    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const int max_size = GridT.max_atom;
    // how many meshcells in bigcell.
    const int bxyz = GridT.bxyz;	

	double*** dr;// vectors between atom and grid: [bxyz, maxsize, 3]
	double** distance; // distance between atom and grid: [bxyz, maxsize]
	double*** psir_ylm;
	bool** cal_flag;
	double* ylma;
	double*** dphi_x;
	double*** dphi_y;
	double*** dphi_z;
    if(max_size!=0)
    {
        dr = new double**[bxyz];
        distance = new double*[bxyz];
        psir_ylm = new double**[bxyz];
        cal_flag = new bool*[bxyz];
		dphi_x = new double**[bxyz];
		dphi_y = new double**[bxyz];
		dphi_z = new double**[bxyz];

        // mohan fix bug 2011-05-02
        int nn = 0;
        for(int it=0; it<ucell.ntype; it++)
        {
            nn = max(nn, (ucell.atoms[it].nwl+1)*(ucell.atoms[it].nwl+1));
        }
        ylma = new double[nn];
        ZEROS(ylma, nn);

        for(int i=0; i<bxyz; i++)
        {
            dr[i] = new double*[max_size];
            psir_ylm[i] = new double*[max_size];
            distance[i] = new double[max_size];
            cal_flag[i] = new bool[max_size];
			dphi_x[i] = new double*[max_size];
			dphi_y[i] = new double*[max_size];
			dphi_z[i] = new double*[max_size];

            ZEROS(distance[i], max_size);
            ZEROS(cal_flag[i], max_size);

            for(int j=0; j<max_size; j++)
            {
                dr[i][j] = new double[3];
                psir_ylm[i][j] = new double[ucell.nwmax];
				dphi_x[i][j] = new double[ucell.nwmax];
				dphi_y[i][j] = new double[ucell.nwmax];
				dphi_z[i][j] = new double[ucell.nwmax];
                ZEROS(dr[i][j],3);
                ZEROS(psir_ylm[i][j],ucell.nwmax);
                ZEROS(dphi_x[i][j],ucell.nwmax);
                ZEROS(dphi_y[i][j],ucell.nwmax);
                ZEROS(dphi_z[i][j],ucell.nwmax);
            }
        }
    }

    assert(this->ncxyz!=0);
    const double dv = ucell.omega/this->ncxyz;
    int vl_index=0;
    double* vldr3 = new double[bxyz];
    ZEROS(vldr3, bxyz);

	for(int i=0; i<nbx; i++)
	{
		for(int j=0; j<nby; j++)
		{
			for(int k=nbz_start; k<nbz_start+nbz; k++)
			{
				const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
				const int size = GridT.how_many_atoms[ grid_index ];
				if(size==0) continue;

				//---------------------------------
				// get the wave functions in this
				// grid.
				//---------------------------------
				this->set_ijk_atom_force(grid_index, size, 
				psir_ylm, dr, cal_flag, 
				distance, ylma, delta_r,
				dphi_x, dphi_y, dphi_z);

				int bindex = 0;
				// z is the fastest,
				for(int ii=0; ii<pw.bx; ii++)
				{
					for(int jj=0; jj<pw.by; jj++)
					{
						for(int kk=0; kk<pw.bz; kk++)
						{
							const int iii = i*pw.bx + ii;
							const int jjj = j*pw.by + jj;
							const int kkk = k*pw.bz + kk;
							vl_index = (kkk-pw.nczp_start) + jjj*pw.nczp + iii*pw.ncy*pw.nczp;
							vldr3[bindex] = vl[ vl_index ] * dv;
							//vldr3[bindex] = dv; // for overlap test
				
							++bindex;
						}
					}
				}

				this->evaluate_vl_force(grid_index, size,i,j,k,
					psir_ylm, cal_flag, vldr3, distance,
					dphi_x, dphi_y, dphi_z,
					pvdpx, pvdpy, pvdpz, GridT);

			}// int k
		}// int j
	} // int i


	//---------------------------------------
	// Folding R here
	//---------------------------------------


	//LM.DHloc_fixedR_x
	this->folding_force(fvl_dphi,pvdpx, pvdpy, pvdpz);

	delete[] pvdpx;
	delete[] pvdpy;
	delete[] pvdpz;

    delete[] vldr3;
    if(max_size!=0)
    {
        for(int i=0; i<pw.bxyz; i++)
        {
            for(int j=0; j<max_size; j++)
            {
                delete[] dr[i][j];
                delete[] psir_ylm[i][j];
				delete[] dphi_x[i][j];
				delete[] dphi_y[i][j];
				delete[] dphi_z[i][j];
            }
            delete[] dr[i];
            delete[] distance[i];
            delete[] psir_ylm[i];
            delete[] cal_flag[i];
			delete[] dphi_x[i];
			delete[] dphi_y[i];
			delete[] dphi_z[i];
        }
        delete[] dr;
        delete[] distance;
        delete[] psir_ylm;
		delete[] dphi_x;
		delete[] dphi_y;
		delete[] dphi_z;
        delete[] cal_flag;

        delete[] ylma;
    }
	timer::tick("Gint_k","cal_force");
	return;
}

void Gint_k::svl_k_RealSpace(
	matrix& fvl_dphi, 
	matrix& svl_dphi, 
	const double *vl)
{
	TITLE("Gint_k","cal_stress");
	timer::tick("Gint_k","cal_stress");

	if(!this->reduced)
	{
		WARNING_QUIT("Gint_k::cal_stress_k","The stress with k can only with reduced H.");
	}

	int nnrg = LNNR.nnrg;

	if(OUT_LEVEL != "m") ofs_running << " LNNR.nnrg in cal_force_k = " << LNNR.nnrg << endl;
	assert(nnrg>=0);

	// just because to make thea arrys meaningful.
	if(LNNR.nnrg == 0)
	{
		nnrg = 1;
	}

	// to store < phi | vlocal | dphi>
	double* pvdpx = new double[nnrg];
	double* pvdpy = new double[nnrg];
	double* pvdpz = new double[nnrg];
	double* pvdp11 = new double[nnrg];
	double* pvdp22 = new double[nnrg];
	double* pvdp33 = new double[nnrg];
	double* pvdp12 = new double[nnrg];
	double* pvdp13 = new double[nnrg];
	double* pvdp23 = new double[nnrg];
	ZEROS(pvdpx, nnrg);
	ZEROS(pvdpy, nnrg);
	ZEROS(pvdpz, nnrg);
	ZEROS(pvdp11, nnrg);
	ZEROS(pvdp22, nnrg);
	ZEROS(pvdp33, nnrg);
	ZEROS(pvdp12, nnrg);
	ZEROS(pvdp13, nnrg);
	ZEROS(pvdp23, nnrg);


	const double delta_r = ORB.dr_uniform;
	// it's a uniform grid to save orbital values, so the delta_r is a constant.
	const int max_size = GridT.max_atom;
	// how many meshcells in bigcell.
	const int bxyz = GridT.bxyz;

	double*** dr;// vectors between atom and grid: [bxyz, maxsize, 3]
	double** distance; // distance between atom and grid: [bxyz, maxsize]
	double*** psir_ylm;
	bool** cal_flag;
	double* ylma;
	double*** dphi_x;
	double*** dphi_y;
	double*** dphi_z;

	if(max_size!=0)
	{
		dr = new double**[bxyz];
		distance = new double*[bxyz];
		psir_ylm = new double**[bxyz];
		cal_flag = new bool*[bxyz];
		dphi_x = new double**[bxyz];
		dphi_y = new double**[bxyz];
		dphi_z = new double**[bxyz];

		// mohan fix bug 2011-05-02
		int nn = 0;
		for(int it=0; it<ucell.ntype; it++)
		{
			nn = max(nn, (ucell.atoms[it].nwl+1)*(ucell.atoms[it].nwl+1));
		}
		ylma = new double[nn];
		ZEROS(ylma, nn);

		for(int i=0; i<bxyz; i++)
		{
			dr[i] = new double*[max_size];
			psir_ylm[i] = new double*[max_size];
			distance[i] = new double[max_size];
			cal_flag[i] = new bool[max_size];
			dphi_x[i] = new double*[max_size];
			dphi_y[i] = new double*[max_size];
			dphi_z[i] = new double*[max_size];

			ZEROS(distance[i], max_size);
			ZEROS(cal_flag[i], max_size);

			for(int j=0; j<max_size; j++)
			{
				dr[i][j] = new double[3];
				psir_ylm[i][j] = new double[ucell.nwmax];
				dphi_x[i][j] = new double[ucell.nwmax];
				dphi_y[i][j] = new double[ucell.nwmax];
				dphi_z[i][j] = new double[ucell.nwmax];
				ZEROS(dr[i][j],3);
				ZEROS(psir_ylm[i][j],ucell.nwmax);
				ZEROS(dphi_x[i][j],ucell.nwmax);
				ZEROS(dphi_y[i][j],ucell.nwmax);
				ZEROS(dphi_z[i][j],ucell.nwmax);
			}
		}
	}

	assert(this->ncxyz!=0);
	const double dv = ucell.omega/this->ncxyz;
	int vl_index=0;
	double* vldr3 = new double[bxyz];
	ZEROS(vldr3, bxyz);

	for(int i=0; i<nbx; i++)
	{
		for(int j=0; j<nby; j++)
		{
			for(int k=nbz_start; k<nbz_start+nbz; k++)
			{
				const int grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;
				const int size = GridT.how_many_atoms[ grid_index ];
				if(size==0) continue;

				//---------------------------------
				// get the wave functions in this
				// grid.
				//---------------------------------
				this->set_ijk_atom_force(grid_index, size,
						psir_ylm, dr, cal_flag,
						distance, ylma, delta_r,
						dphi_x, dphi_y, dphi_z);

				int bindex = 0;
				// z is the fastest,
				for(int ii=0; ii<pw.bx; ii++)
				{
					for(int jj=0; jj<pw.by; jj++)
					{
						for(int kk=0; kk<pw.bz; kk++)
						{
							const int iii = i*pw.bx + ii;
							const int jjj = j*pw.by + jj;
							const int kkk = k*pw.bz + kk;
							vl_index = (kkk-pw.nczp_start) + jjj*pw.nczp + iii*pw.ncy*pw.nczp;
							vldr3[bindex] = vl[ vl_index ] * dv;
							//        vldr3[bindex] = dv; // for overlap test

							++bindex;
						}
					}
				}
				//cout<<"loop  "<<i<<" "<<j<<" "<<k<<endl;//test

				this->evaluate_vl_stress(grid_index, size,i,j,k,
						psir_ylm, cal_flag, vldr3, distance,
						dphi_x, dphi_y, dphi_z,
						pvdpx, pvdpy, pvdpz,
						pvdp11, pvdp22, pvdp33, pvdp12, pvdp13, pvdp23, dr, GridT);
			}// int k
		}// int j
	} // int i


	//---------------------------------------
	// Folding R here
	//---------------------------------------

	//LM.DHloc_fixedR_x
	this->folding_stress(fvl_dphi, svl_dphi, pvdpx, pvdpy, pvdpz,
			pvdp11, pvdp22, pvdp33, pvdp12, pvdp13, pvdp23);

	delete[] pvdpx;
	delete[] pvdpy;
	delete[] pvdpz;
	delete[] pvdp11;
	delete[] pvdp22;
	delete[] pvdp33;
	delete[] pvdp12;
	delete[] pvdp13;
	delete[] pvdp23;

	delete[] vldr3;
	if(max_size!=0)
	{
		for(int i=0; i<pw.bxyz; i++)
		{
			for(int j=0; j<max_size; j++)
			{
				delete[] dr[i][j];
				delete[] psir_ylm[i][j];
				delete[] dphi_x[i][j];
				delete[] dphi_y[i][j];
				delete[] dphi_z[i][j];
			}
			delete[] dr[i];
			delete[] distance[i];
			delete[] psir_ylm[i];
			delete[] cal_flag[i];
			delete[] dphi_x[i];
			delete[] dphi_y[i];
			delete[] dphi_z[i];
		}
		delete[] dr;
		delete[] distance;
		delete[] psir_ylm;
		delete[] dphi_x;
		delete[] dphi_y;
		delete[] dphi_z;
		delete[] cal_flag;

		delete[] ylma;
	}
	timer::tick("Gint_k","cal_stress");
	return;
}


void Gint_k::evaluate_vl_stress(
	const int &grid_index, 
	const int &size, 
	const int &i, 
	const int &j, 
	const int &k,
	double*** psir_ylm, 
	bool** cal_flag, 
	double* vldr3, 
	double** distance,
	double*** dphi_x, 
	double*** dphi_y, 
	double*** dphi_z,
	double* pvdpx, 
	double* pvdpy, 
	double* pvdpz, 
	double* pvdp11, 
	double* pvdp22, 
	double* pvdp33, 
	double* pvdp12, 
	double* pvdp13, 
	double* pvdp23, 
	double*** dr,
	const Grid_Technique &gt)
{

	double *psi1, *psi2;
	double *iw1p, *iw2p;
	double *iw1px, *iw1py, *iw1pz;//extra pointer compared to non-force grid integration.
	double *iw2px, *iw2py, *iw2pz;//extra pointer compared to non-force grid integration.
	double *end1, *end2;
	double *pvp1, *pvp2, *pvp3;//extra pointer
        double *pvp11, *pvp22, *pvp33, *pvp12, *pvp13, *pvp23;
	int iw1_lo, iw2_lo;
	int iwi, iww;
	double vpsir1, vpsir2, vpsir3;//extra pointer
        double vpsir11, vpsir22, vpsir33, vpsir12, vpsir13, vpsir23;//extra pointer
	double *psix, *psiy, *psiz;


	bool *all_out_of_range = new bool[size];
	for(int ia=0; ia<size; ++ia) //number of atoms
	{
		all_out_of_range[ia] = true;
		for(int ib=0; ib<gt.bxyz; ++ib) //number of small box in big box
		{
			if(cal_flag[ib][ia])
			{
				all_out_of_range[ia] = false;
				//break; //mohan add 2012-07-10
			}
		}
	}


	double* dmR = LOC.DM_R[CURRENT_SPIN];

	double* dmR2;
	for (int ia1=0; ia1<size; ++ia1)
	{
		if(all_out_of_range[ia1]) continue;

        const int mcell_index1 = gt.bcell_start[grid_index] + ia1;
		const int iat = gt.which_atom[mcell_index1];
        const int T1 = ucell.iat2it[iat];
        const int I1 = ucell.iat2ia[iat];
        const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
		const int iw1_start = gt.trace_lo[start1];
        Atom *atom1 = &ucell.atoms[T1];
	
        //~~~~~~~~~~~~~~~~
        // get cell R1.
        //~~~~~~~~~~~~~~~~
        const int id1 = gt.which_unitcell[mcell_index1];
        const int R1x = gt.ucell_index2x[id1];
        const int R1y = gt.ucell_index2y[id1];
        const int R1z = gt.ucell_index2z[id1];
        const int DM_start = LNNR.nlocstartg[iat];

        // get (j,beta,R2)
        for (int ia2=0; ia2<size; ++ia2)
        {
			if(all_out_of_range[ia2]) continue;

			//---------------------------------------------
			// check if we need to calculate the big cell.
			//---------------------------------------------
			bool same_flag = false;
			for(int ib=0; ib<gt.bxyz; ++ib)
			{
				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
				{
				//	cout << " same flag is = " << ib << endl;
					same_flag = true;
					break;
				}
			} 

			if(!same_flag) continue;

            const int bcell2 = gt.bcell_start[grid_index] + ia2;
            const int T2 = ucell.iat2it[ gt.which_atom[bcell2]];
			const int iat2 = gt.which_atom[bcell2];

//---------------------------------------------------------
//            if (T2 >= T1 ) //mohan fix bug 2012-07-10
//---------------------------------------------------------
            {
                Atom *atom2 = &ucell.atoms[T2];
                const int I2 = ucell.iat2ia[iat2];
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
				const int iw2_start = gt.trace_lo[start2];

	            //---------------
                // get cell R2.
	            //---------------
                const int id2 = gt.which_unitcell[bcell2];
                const int R2x = gt.ucell_index2x[id2];
                const int R2y = gt.ucell_index2y[id2];
                const int R2z = gt.ucell_index2z[id2];

				//------------------------------------------------
				// calculate the 'offset': R2 position relative
				// to R1 atom.
				//------------------------------------------------
                const int dRx = R1x - R2x;
                const int dRy = R1y - R2y;
                const int dRz = R1z - R2z;

                double rt[3];
                Mathzone::Direct_to_Cartesian(dRx,dRy,dRz,
                                              ucell.a1.x, ucell.a1.y, ucell.a1.z,
                                              ucell.a2.x, ucell.a2.y, ucell.a2.z,
                                              ucell.a3.x, ucell.a3.y, ucell.a3.z,
                                              rt[0],rt[1],rt[2]);
	
				const int index = LNNR.cal_RindexAtom(dRx, dRy, dRz, iat2);
                int offset = -1;

				int* find_start = LNNR.find_R2[iat];
				int* findend = LNNR.find_R2[iat] + LNNR.nad[iat];
				
				// the nad should be a large expense of time.
				for(int* find=find_start; find < findend; find++)
				{
					if( find[0] == index )
					{
						offset = find - find_start;
						break;
					}
				}

				if(offset == -1 )
                {
                    WARNING_QUIT("gint_k","evaluate_vl_force wrong");
                }
                assert(offset < LNNR.nad[iat]);

				//--------------------------------------------------------------- 
				// what I do above is to get 'offset' for atom pair (iat1, iat2)
				// if I want to simplify this searching for offset,
				// I should take advantage of gt.which_unitcell.
				//--------------------------------------------------------------- 

				const int iatw = DM_start + LNNR.find_R2st[iat][offset];
				
				for(int ib=0; ib<gt.bxyz; ++ib)
				{
					if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
					{
						psi1 = psir_ylm[ib][ia1];
						psi2 = psir_ylm[ib][ia2];
				
						psix = dphi_x[ib][ia1];
						psiy = dphi_y[ib][ia1];
						psiz = dphi_z[ib][ia1];
						//psix = psi1; // for overlap test
						//psiy = psi1; // for overlap test
						//psiz = psi1; // for overlap test
						
						
						end1 = psi1 + atom1->nw;
						end2 = psi2 + atom2->nw;
						iw1_lo = iw1_start;	
						//------------------------------------
						// circle for wave functions of atom 1.
						//------------------------------------
						iwi = 0;

						iw1px = psix;
						iw1py = psiy;
						iw1pz = psiz;

						for (iw1p=psi1; iw1p < end1; ++ iw1p)
						{
							//vpsir1 = iw1p[0] * vldr3[ib];
							vpsir1 = iw1px[0] * vldr3[ib];
							vpsir2 = iw1py[0] * vldr3[ib];
							vpsir3 = iw1pz[0] * vldr3[ib];
                            vpsir11 = iw1px[0] * vldr3[ib] * dr[ib][ia1][0];
                            vpsir22 = iw1py[0] * vldr3[ib] * dr[ib][ia1][1];
                            vpsir33 = iw1pz[0] * vldr3[ib] * dr[ib][ia1][2];
                            vpsir12 = iw1px[0] * vldr3[ib] * dr[ib][ia1][1];
                            vpsir13 = iw1px[0] * vldr3[ib] * dr[ib][ia1][2];
                            vpsir23 = iw1py[0] * vldr3[ib] * dr[ib][ia1][2];
							++iw1px;
							++iw1py;
							++iw1pz;

							iw2_lo = iw2_start;
							iww = iatw + iwi;// -1 because ++iww from below.

							dmR2 = &dmR[iww]; //mohan add 2012-01-05

							iwi += atom2->nw;
					
							//---------------------------------
							// only correct for one processor
							//---------------------------------
//							pvp1 = &LM.DHloc_fixedR_x[iww];
//							pvp2 = &LM.DHloc_fixedR_y[iww];
//							pvp3 = &LM.DHloc_fixedR_z[iww];

							//--------------------------------------
							// store phi_i(r) vlocal(r) * dphi_j(r)
							//--------------------------------------
							pvp1 = &pvdpx[iww]; //mohan add 2012-1-6
							pvp2 = &pvdpy[iww];
							pvp3 = &pvdpz[iww];

                            pvp11 = &pvdp11[iww]; //zhengdy add 2017/3/28
                            pvp22 = &pvdp22[iww];
                            pvp33 = &pvdp33[iww];
                            pvp12 = &pvdp12[iww]; 
                            pvp13 = &pvdp13[iww];
                            pvp23 = &pvdp23[iww];
							//------------------------------------
							// circle for wave functions of atom 2.
							//------------------------------------
							iw2px = psix;
							iw2py = psiy;
							iw2pz = psiz;

							for(iw2p=psi2; iw2p < end2; ++iw2p,
								++iw2px, ++iw2py, ++iw2pz)
							{
								// the main difference to calculate
								// the force is that the whole 
								// matrix should be calculated!
								//if( iw1_lo > iw2_lo)
								//{
								//	++iw2_lo;
								//	++pvp1;
								//	++pvp2;
								//	++pvp3;
								//	continue;
								//}
								// mohan tmp

						// bug here!!!!!!!!!!!!!!!!!!!!!!!!!!!
						// DM(R) * psi1(r) * v(r) * psi2_R(r)
								pvp1[0] += dmR2[0] * vpsir1 * iw2p[0];
								pvp2[0] += dmR2[0] * vpsir2 * iw2p[0];
								pvp3[0] += dmR2[0] * vpsir3 * iw2p[0];
                                pvp11[0] += dmR2[0] * vpsir11 * iw2p[0] ;
                                pvp22[0] += dmR2[0] * vpsir22 * iw2p[0] ;
                                pvp33[0] += dmR2[0] * vpsir33 * iw2p[0] ;
                                pvp12[0] += dmR2[0] * vpsir12 * iw2p[0] ;
                                pvp13[0] += dmR2[0] * vpsir13 * iw2p[0] ;
                                pvp23[0] += dmR2[0] * vpsir23 * iw2p[0] ;
//								pvp1[0] += vpsir1 * iw2p[0];
//								pvp2[0] += vpsir2 * iw2p[0];
//								pvp3[0] += vpsir3 * iw2p[0];
						// bug here!!!!!!!!!!!!!!!!!!!!!!!!!!!

								//pvp1[0] -= vpsir1 * iw2px[0];
								//pvp2[0] -= vpsir1 * iw2py[0];
								//pvp3[0] -= vpsir1 * iw2pz[0];

								++iw2_lo;
								++pvp1;
								++pvp2;
								++pvp3;
                                ++pvp11;
                                ++pvp22;
                                ++pvp33;
                                ++pvp12;
                                ++pvp13;
                                ++pvp23;
						//density matrix
								++dmR2;
							}
							++iw1_lo;
						}// iw
					}//end flag
				}//end ib
            }// T
        }// ia2
	}//ia1

	delete[] all_out_of_range;
	return;
}

void Gint_k::evaluate_vl_force(const int &grid_index, const int &size, const int &i, const int &j, const int &k,
        double*** psir_ylm, bool** cal_flag, double* vldr3, double** distance,
        double*** dphi_x, double*** dphi_y, double*** dphi_z,
        double* pvdpx, double* pvdpy, double* pvdpz,
        const Grid_Technique &gt)
{                       
                                
        double *psi1, *psi2;    
        double *iw1p, *iw2p;
        double *iw1px, *iw1py, *iw1pz;//extra pointer compared to non-force grid integration.
        double *iw2px, *iw2py, *iw2pz;//extra pointer compared to non-force grid integration.
        double *end1, *end2;
        double *pvp1, *pvp2, *pvp3;//extra pointer
        int iw1_lo, iw2_lo;
        int iwi, iww;
        double vpsir1, vpsir2, vpsir3;//extra pointer
        double *psix, *psiy, *psiz;


        bool *all_out_of_range = new bool[size];
        for(int ia=0; ia<size; ++ia) //number of atoms
        {
                all_out_of_range[ia] = true;
                for(int ib=0; ib<gt.bxyz; ++ib) //number of small box in big box
                {
                        if(cal_flag[ib][ia])
                        {
                                all_out_of_range[ia] = false;
                                //break; //mohan add 2012-07-10
                        }
                }
        }


        double* dmR = LOC.DM_R[CURRENT_SPIN];
        double* dmR2;
        for (int ia1=0; ia1<size; ++ia1)
        {
                if(all_out_of_range[ia1]) continue;

        const int mcell_index1 = gt.bcell_start[grid_index] + ia1;
                const int iat = gt.which_atom[mcell_index1];
        const int T1 = ucell.iat2it[iat];
        const int I1 = ucell.iat2ia[iat];
        const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const int iw1_start = gt.trace_lo[start1];
        Atom *atom1 = &ucell.atoms[T1];

        //~~~~~~~~~~~~~~~~
        // get cell R1.
        //~~~~~~~~~~~~~~~~
        const int id1 = gt.which_unitcell[mcell_index1];
        const int R1x = gt.ucell_index2x[id1];
        const int R1y = gt.ucell_index2y[id1];
        const int R1z = gt.ucell_index2z[id1];
        const int DM_start = LNNR.nlocstartg[iat];

        // get (j,beta,R2)
        for (int ia2=0; ia2<size; ++ia2)
        {
                        if(all_out_of_range[ia2]) continue;

                        //---------------------------------------------
                        // check if we need to calculate the big cell.
                        //---------------------------------------------
                        bool same_flag = false;
                        for(int ib=0; ib<gt.bxyz; ++ib)
                        {
                                if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                                {
                                //      cout << " same flag is = " << ib << endl;
                                        same_flag = true;
                                        break;
                                }
                        }

                        if(!same_flag) continue;

            const int bcell2 = gt.bcell_start[grid_index] + ia2;
            const int T2 = ucell.iat2it[ gt.which_atom[bcell2]];
                        const int iat2 = gt.which_atom[bcell2];

//---------------------------------------------------------
//            if (T2 >= T1 ) //mohan fix bug 2012-07-10
//---------------------------------------------------------
            {
                Atom *atom2 = &ucell.atoms[T2];
                const int I2 = ucell.iat2ia[iat2];
                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
                                const int iw2_start = gt.trace_lo[start2];

                    //---------------
                // get cell R2.
                    //---------------
                const int id2 = gt.which_unitcell[bcell2];
                const int R2x = gt.ucell_index2x[id2];
                const int R2y = gt.ucell_index2y[id2];
                const int R2z = gt.ucell_index2z[id2];

                                //------------------------------------------------
                                // calculate the 'offset': R2 position relative
                                // to R1 atom.
                                //------------------------------------------------
                const int dRx = R1x - R2x;
                const int dRy = R1y - R2y;
                const int dRz = R1z - R2z;

                                const int index = LNNR.cal_RindexAtom(dRx, dRy, dRz, iat2);
                int offset = -1;

                                int* find_start = LNNR.find_R2[iat];
                                int* findend = LNNR.find_R2[iat] + LNNR.nad[iat];

                                // the nad should be a large expense of time.
                                for(int* find=find_start; find < findend; find++)
                                {
                                        if( find[0] == index )
                                        {
                                                offset = find - find_start;
                                                break;
                                        }
                                }

                                if(offset == -1 )
                {
                    WARNING_QUIT("gint_k","evaluate_vl_force wrong");
                }
                assert(offset < LNNR.nad[iat]);

                                //--------------------------------------------------------------- 
                                // what I do above is to get 'offset' for atom pair (iat1, iat2)
                                // if I want to simplify this searching for offset,
                                // I should take advantage of gt.which_unitcell.
                                //--------------------------------------------------------------- 

                                const int iatw = DM_start + LNNR.find_R2st[iat][offset];

                                for(int ib=0; ib<gt.bxyz; ++ib)
                                {
                                        if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                                        {
                                                psi1 = psir_ylm[ib][ia1];
                                                psi2 = psir_ylm[ib][ia2];

                                                psix = dphi_x[ib][ia1];
                                                psiy = dphi_y[ib][ia1];
                                                psiz = dphi_z[ib][ia1];
                                                //psix = psi1; // for overlap test
                                                //psiy = psi1; // for overlap test
                                                //psiz = psi1; // for overlap test


                                                end1 = psi1 + atom1->nw;
                                                end2 = psi2 + atom2->nw;
                                                iw1_lo = iw1_start;
                                                //------------------------------------
                                                // circle for wave functions of atom 1.
                                                //------------------------------------
                                                iwi = 0;

                                                iw1px = psix;
                                                iw1py = psiy;
                                                iw1pz = psiz;
                                                for (iw1p=psi1; iw1p < end1; ++ iw1p)
                                                {
                                                        //vpsir1 = iw1p[0] * vldr3[ib];
                                                        vpsir1 = iw1px[0] * vldr3[ib];
                                                        vpsir2 = iw1py[0] * vldr3[ib];
                                                        vpsir3 = iw1pz[0] * vldr3[ib];
                                                        ++iw1px;
                                                        ++iw1py;
                                                        ++iw1pz;

                                                        iw2_lo = iw2_start;
                                                        iww = iatw + iwi;// -1 because ++iww from below.

                                                        dmR2 = &dmR[iww]; //mohan add 2012-01-05

                                                        iwi += atom2->nw;

                                                        //---------------------------------
                                                        // only correct for one processor
                                                        //---------------------------------
//                                                      pvp1 = &LM.DHloc_fixedR_x[iww];
//                                                      pvp2 = &LM.DHloc_fixedR_y[iww];
//                                                      pvp3 = &LM.DHloc_fixedR_z[iww];

                                                        //--------------------------------------
                                                        // store phi_i(r) vlocal(r) * dphi_j(r)
                                                        //--------------------------------------
                                                        pvp1 = &pvdpx[iww]; //mohan add 2012-1-6
                                                        pvp2 = &pvdpy[iww];
                                                        pvp3 = &pvdpz[iww];

                                                        //------------------------------------
                                                        // circle for wave functions of atom 2.
                                                        //------------------------------------
                                                        iw2px = psix;
                                                        iw2py = psiy;
                                                        iw2pz = psiz;

                                                        for(iw2p=psi2; iw2p < end2; ++iw2p,
                                                                ++iw2px, ++iw2py, ++iw2pz)
                                                        {
                                                                // the main difference to calculate
                                                                // the force is that the whole 
                                                                // matrix should be calculated!
                                                                //if( iw1_lo > iw2_lo)
                                                                //{
                                                                //      ++iw2_lo;
                                                                //      ++pvp1;
                                                                //      ++pvp2;
                                                                //      ++pvp3;
                                                                //      continue;
                                                                //}
                                                                // mohan tmp

                                                // bug here!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                // DM(R) * psi1(r) * v(r) * psi2_R(r)
                                                                        pvp1[0] += dmR2[0] * vpsir1 * iw2p[0];
                                                                        pvp2[0] += dmR2[0] * vpsir2 * iw2p[0];
                                                                        pvp3[0] += dmR2[0] * vpsir3 * iw2p[0];
//                                                              pvp1[0] += vpsir1 * iw2p[0];
//                                                              pvp2[0] += vpsir2 * iw2p[0];
//                                                              pvp3[0] += vpsir3 * iw2p[0];
                                                // bug here!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                                                //pvp1[0] -= vpsir1 * iw2px[0];
                                                                //pvp2[0] -= vpsir1 * iw2py[0];
                                                                //pvp3[0] -= vpsir1 * iw2pz[0];

                                                                ++iw2_lo;
                                                                ++pvp1;
                                                                ++pvp2;
                                                                ++pvp3;
                                                //density matrix
                                                                ++dmR2;
                                                        }
                                                        ++iw1_lo;
                                                }// iw
                                        }//end flag
                                }//end ib
            }// T
        }// ia2
        }//ia1

        delete[] all_out_of_range;
        return;
}


// PLEASE be aware that 'set_ijk' subroutines should be reconstructed
// since it has been used everytime grid integral is needed
// mohan add 2021-03-28
void Gint_k::set_ijk_atom_force(
	const int &grid_index, 
	const int &size,
	double*** psir_ylm, 
	double*** dr, 
	bool** cal_flag, 
	double** distance, 
	double* ylma, 
	const double &delta_r,
	double*** dphi_x, 
	double*** dphi_y, 
	double*** dphi_z)
{
	const Numerical_Orbital_Lm* pointer;
	double mt[3];
    // Peize Lin change rly, grly 2016-08-26
    static vector<double> rly;
    static vector<vector<double>> grly;
	for (int id=0; id<size; ++id)
	{
		// (2.1) get the atom type and atom index.
		const int mcell_index = GridT.bcell_start[grid_index] + id;	
		const int imcell = GridT.which_bigcell[mcell_index];
		const int iat = GridT.which_atom[mcell_index];
		const int it = ucell.iat2it[ iat ];
		const int ia = ucell.iat2ia[ iat ];
		Atom *atom = &ucell.atoms[it];

		// (2.2) get the distance between the grid and the atom.
		mt[0] = GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0];
		mt[1] = GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1];
		mt[2] = GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2];

		for(int ib=0; ib<pw.bxyz; ++ib)
		{
			// meshcell_pos: z is the fastest
			dr[ib][id][0] = GridT.meshcell_pos[ib][0] + mt[0];
			dr[ib][id][1] = GridT.meshcell_pos[ib][1] + mt[1];
			dr[ib][id][2] = GridT.meshcell_pos[ib][2] + mt[2];

			distance[ib][id] = std::sqrt(dr[ib][id][0]*dr[ib][id][0] 
			+ dr[ib][id][1]*dr[ib][id][1] + dr[ib][id][2]*dr[ib][id][2]);

			if(distance[ib][id] < ORB.Phi[it].getRcut())
			{
				cal_flag[ib][id]=true;
			}
			else
			{
				cal_flag[ib][id]=false;
				continue;
			}

			// get the 'phi' and 'dphi'.
			Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr[ib][id][0], dr[ib][id][1], dr[ib][id][2], rly, grly);

//			Ylm::sph_harm ( ucell.atoms[it].nwl,
//					dr[ib][id][0] / distance[ib][id],
//					dr[ib][id][1] / distance[ib][id],
//					dr[ib][id][2] / distance[ib][id],
//					ylma);


            //-------------------------------------------------
            // Here we can not deal with the situation on
            // r = 0, so if r = 0,  r-->1e-9
			//-------------------------------------------------

			if (distance[ib][id] < 1e-9)    // pengfei Li add 2016-3-3
			{
				distance[ib][id] = 1e-9;
			}

			// these parameters are about interpolation
			// because once we know the distance from atom to grid point,
			// we can get the parameters we need to do interpolation and
			// store them first!! these can save a lot of effort.
			const double position = distance[ib][id] / delta_r;

			const int iq = static_cast<int>(position);
			const double x0 = position - static_cast<double>(iq);
			const double x1 = 1.0 - x0;
			const double x2 = 2.0 - x0;
			const double x3 = 3.0 - x0;
			const double x12 = x1*x2 / 6.0;
			const double x03 = x0*x3 / 2.0;

			double tmp = 0.0;//mohan fix bug 2011-05-04
			double dtmp = 0.0;
			for (int iw=0; iw< atom->nw; ++iw)
			{
				if ( atom->iw2_new[iw] )
				{
					pointer = &ORB.Phi[it].PhiLN(
							atom->iw2l[iw],
							atom->iw2n[iw]);

					if(iq >= pointer->nr_uniform-4)
					{
						tmp = dtmp = 0.0;
					}
					else
					{
						// Efficient!! to get the orbital value at this point.
						tmp = x12*(pointer->psi_uniform[iq]*x3
							+pointer->psi_uniform[iq+3]*x0)
						+ x03*(pointer->psi_uniform[iq+1]*x2
						-pointer->psi_uniform[iq+2]*x1);

						dtmp = x12*(pointer->dpsi_uniform[iq]*x3
						+pointer->dpsi_uniform[iq+3]*x0)
						+ x03*(pointer->dpsi_uniform[iq+1]*x2
						-pointer->dpsi_uniform[iq+2]*x1);

						//dtmp = x12*(pointer->psi_uniform[iq]*x3
						//+pointer->psi_uniform[iq+3]*x0)
						//+ x03*(pointer->psi_uniform[iq+1]*x2
						//-pointer->psi_uniform[iq+2]*x1);
					}
				}// new l is used.

				int ll = atom->iw2l[iw];
				int idx_lm = atom->iw2_ylm[iw];
				//special case for distance[id] -> 0
				//Problems Remained
				//You have to add this two lines
				double rr = distance[ib][id];
				/*if (rr < 1e-9)
				{
					if (ll == 0)
					{
						psir_ylm[ib][id][iw] = tmp * rly[idx_lm];
						dphi_x[ib][id][iw] = dphi_y[ib][id][iw] = dphi_z[ib][id][iw] = 0.0;
					}
					else
					{
						pointer = &ORB.Phi[it].
							PhiLN(atom->iw2l[iw], atom->iw2n[iw]);

						double Zty = pointer->zty;
						psir_ylm[ib][id][iw] = Zty * rly[idx_lm];
						dphi_x[ib][id][iw] = Zty * grly[idx_lm][0];
						dphi_y[ib][id][iw] = Zty * grly[idx_lm][1];
						dphi_z[ib][id][iw] = Zty * grly[idx_lm][2];
					}
				}*/
				//else
				//{
					double rl;
					if(ll==0)
					{
						rl = 1.0;
					}
					else if(ll==1)
					{
						rl = rr;
					}
					else
					{
						rl = pow(rr, ll);
					}

					psir_ylm[ib][id][iw] = tmp * rly[idx_lm] / rl;
					double tmpdphi_rly = (dtmp  - tmp * ll / rr) / rl * rly[idx_lm] / rr;
					double tmprl = tmp/rl;
					dphi_x[ib][id][iw] = tmpdphi_rly * dr[ib][id][0]  + tmprl * grly[idx_lm][0];
					dphi_y[ib][id][iw] = tmpdphi_rly * dr[ib][id][1]  + tmprl * grly[idx_lm][1];
					dphi_z[ib][id][iw] = tmpdphi_rly * dr[ib][id][2]  + tmprl * grly[idx_lm][2];
				//}
			}
		}//end ib
	}// int id
	return;
}


