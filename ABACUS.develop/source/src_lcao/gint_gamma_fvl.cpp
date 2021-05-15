#include "gint_gamma.h"
#include "grid_technique.h"
#include "../module_ORB/ORB_read.h"
#include "../src_pw/global.h"
#include "src_global/blas_connector.h"

#include "global_fp.h" // mohan add 2021-01-30
#include "../src_global/ylm.h"

void Gint_Gamma::cal_force(const double* vlocal_in)
{
    timer::tick("Gint_Gamma","cal_force",'H');
    this->vlocal = vlocal_in;
    this->save_atoms_on_grid(GridT);
    this->gamma_force();
    timer::tick("Gint_Gamma","cal_force",'H');
}

inline void cal_psir_ylm_dphi(
	const int na_grid, 
	const int grid_index, 
	const double delta_r,
	const int*const block_index,
	const int*const block_size, 
	bool*const*const cal_flag,
	double*const*const psir_ylm, 
	double*const*const dphix, 
	double*const*const dphiy, 
	double*const*const dphiz,
	realArray& drr)
{    
    for (int id=0; id<na_grid; id++)
    {
        const int mcell_index = GridT.bcell_start[grid_index] + id;
        const int imcell = GridT.which_bigcell[mcell_index];
        int iat = GridT.which_atom[mcell_index];
        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];
        const int start = ucell.itiaiw2iwt(it, ia, 0);
        Atom *atom = &ucell.atoms[it];

		const double mt[3]={
			GridT.meshball_positions[imcell][0] - GridT.tau_in_bigcell[iat][0],
			GridT.meshball_positions[imcell][1] - GridT.tau_in_bigcell[iat][1],
			GridT.meshball_positions[imcell][2] - GridT.tau_in_bigcell[iat][2]};

        for(int ib=0; ib<pw.bxyz; ib++)
        {
            double*const p_psir_ylm=&psir_ylm[ib][block_index[id]];
            double*const p_dphix=&dphix[ib][block_index[id]];
            double*const p_dphiy=&dphiy[ib][block_index[id]];
            double*const p_dphiz=&dphiz[ib][block_index[id]];
            
			const double dr[3]={						// vectors between atom and grid
				GridT.meshcell_pos[ib][0] + mt[0],
				GridT.meshcell_pos[ib][1] + mt[1],
				GridT.meshcell_pos[ib][2] + mt[2]};

            if(STRESS)
            {
                for(int i=0;i<3;i++) 
				{
					drr(id,ib,i) = dr[i];
				}
            }

			// distance between atom and grid
            double distance = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

            if(distance > ORB.Phi[it].getRcut())
            {
                ZEROS(p_psir_ylm, block_size[id]);
                ZEROS(p_dphix, block_size[id]);
                ZEROS(p_dphiy, block_size[id]);
                ZEROS(p_dphiz, block_size[id]);
                cal_flag[ib][id]=false;
                continue;
            }

            cal_flag[ib][id]=true;
            
            //array to store spherical harmonics and its derivatives
            vector<double> rly;
            vector<vector<double>> grly;
            // >>> the old method
            // ylma[id] = new double[nnn[it]]; // liaochen found this bug 2010/03/29
            // Ylm::get_ylm_real(ucell.atoms[it].nwl+1, this->dr[id], ylma[id]);
            // <<<
            // Ylm::rlylm(ucell.atoms[it].nwl+1, dr[id].x, dr[id].y, dr[id].z, rly, grly);
            // Ylm::rlylm(ucell.atoms[it].nwl+1, dr[id].x, dr[id].y, dr[id].z, rly, grly);
            Ylm::grad_rl_sph_harm(ucell.atoms[it].nwl, dr[0], dr[1], dr[2], rly, grly);

            // 1E-7 is necessary in case of R is just on one grid
            // the following code is about interpolation,
            // here we use Polynomia interpolation,
            // iq, x0, x1, x3, x3, x12, x03 are all the variables 
            // needed in this interpolation method.
            // It's the same as the function in mathzone.cpp:
            // Mathzone::Polynomial_Interpolation
            // BUT THE IMPORTANT THING IS,
            // in order to minimiz the total operations,
            // we save some variable in advance here,
            // such as x12 and x03, which is only calculated 
            // once and saved here.
            //
            // 'distane' variable here is the distance between
            // this grid and the atom.

            //-------------------------------------------------
            // Here we can not deal with the situation on 
            // r = 0, so if r = 0,  r-->1e-9   
            //-------------------------------------------------
            if(distance < 1e-9)                            // pengfei Li 2016-3-3
            {
               distance = 1e-9;
            }

            const double position = distance / delta_r;
                       
            const double iq = static_cast<int>(position);
            const double x0 = position - iq;
            const double x1 = 1.0 - x0;
            const double x2 = 2.0 - x0;
            const double x3 = 3.0 - x0;
            const double x12 = x1*x2 / 6;
            const double x03 = x0*x3 / 2;
            
            double tmp, dtmp;

            // circle around the localized wave functions
            // for this atom, here localized wave functions
            // has the index (n,l,m), n is the multiplicity,
            // l is the angular momentum and m is the magnetic momentum,
            // for each localized wave function 'iw', we can know
            // the n according to atom->iw2n[iw] function
            // and know the l according to atom->iw2l[iw] function.
            // also the most important thing I always emphasis,
            // is to minimize the operations!!
            // so for a set of differnt m of a specific l,
            // we don't need to do interpolation for each set of (n,l,m),
            // we only need to do interpolation for each (n,l),
            // for differemtn m, only the Spherical Harmonical functions
            // are differnt, what's why we have atom->iw2_new[iw],
            // if the localized wave funcction 'iw' reprensents a new 'l',
            // we do interpolation.
            for (int iw=0; iw< atom->nw; ++iw)
            {
                // this is a new 'l', we need 1D orbital wave
                // function from interpolation method.
                if ( atom->iw2_new[iw] )
                {
                    const Numerical_Orbital_Lm &philn = ORB.Phi[it].PhiLN(
                            atom->iw2l[iw],
                            atom->iw2n[iw]);

                    //if ( iq[id] >= philn.nr_uniform-4)
                    if ( iq >= philn.nr_uniform-4)
                    {
                        tmp = dtmp = 0.0;
                    }
                    else
                    {
                        // use Polynomia Interpolation method to get the 
                        // wave functions

                        tmp = x12*(philn.psi_uniform[iq]*x3
                                +philn.psi_uniform[iq+3]*x0)
                            + x03*(philn.psi_uniform[iq+1]*x2
                                    -philn.psi_uniform[iq+2]*x1);

                        dtmp = x12*(philn.dpsi_uniform[iq]*x3
                                +philn.dpsi_uniform[iq+3]*x0)
                                + x03*(philn.dpsi_uniform[iq+1]*x2
                                    -philn.dpsi_uniform[iq+2]*x1);
                    }
                }//new l is used.
				
				// WARNING: if(!atom->iw2_new[iw]), tmp and dtmp are undefined?

                // get the 'l' of this localized wave function
                const int ll = atom->iw2l[iw];
                const int idx_lm = atom->iw2_ylm[iw];

                //special case for distance -> 0
                //Problems Remained
                //You have to add this two lines
                /*if (distance < 1e-9)
                {
                    // if l==0 for localized orbital,
                    // here is how we choose the 3D wave function psir_ylm,
                    // and the derivative at r=0 point.
                    if (ll == 0)
                    {
                        // psir_ylm is the three dimensional localized wave functions
                        // (n,l,m), which is the multiply between 1D wave function 'tmp' and
                        // spherical harmonic function rly.
                        p_psir_ylm[iw] = tmp * rly[idx_lm];
                        // the derivative of psir_ylm with respect to atom position R,
                        // it's a vector, so it has x,y,z components.
                        p_dphix[iw] = 0.0;
                        p_dphiy[iw] = 0.0;
                        p_dphiz[iw] = 0.0;
                    }
                    // if l>0, how do we choose 3D localized wave function
                    // and the derivative.
                    else
                    {
                        const Numerical_Orbital_Lm &philn = ORB.Phi[it].PhiLN(atom->iw2l[iw], atom->iw2n[iw]);

                        double Zty = philn.zty;
                        p_psir_ylm[iw] = Zty * rly[idx_lm];
                        p_dphix[iw] = Zty * grly[idx_lm][0];
                        p_dphiy[iw] = Zty * grly[idx_lm][1];
                        p_dphiz[iw] = Zty * grly[idx_lm][2];
                    } // if (ll == 0)
                }*/ 
                // if r >0, here is how we choose the 3D wave function psir_ylm,
                // and the derivative at r=0 point.
                //else
                //{
                    const double rl = pow(distance, ll);

                    // 3D wave functions
                    p_psir_ylm[iw] = tmp * rly[idx_lm] / rl;

                    // derivative of wave functions with respect to atom positions.
                    const double tmpdphi_rly = (dtmp  - tmp * ll / distance) / rl * rly[idx_lm] / distance;
                    const double tmprl = tmp/rl;

                    p_dphix[iw] = tmpdphi_rly * dr[0]  + tmprl * grly[idx_lm][0];
                    p_dphiy[iw] = tmpdphi_rly * dr[1]  + tmprl * grly[idx_lm][1];
                    p_dphiz[iw] = tmpdphi_rly * dr[2]  + tmprl * grly[idx_lm][2];
                //}// if  (distance < 1e-9)
            } // iw            
        }// ib
    }//!id //finish loop of calc pre-info for each adjacent atom

	return;
}


inline void cal_meshball_DGridV(
	const int na_grid,
	int lgd_now,
	int LD_pool, 
	const int*const block_index,
	const int*const block_iw, 
	const int*const block_size,
	const bool*const*const cal_flag,
	const double*const vldr3, 
	const double*const*const psir_ylm,
	double*const*const psir_vlbr3,
	const double*const*const dphix, const double*const*const dphiy, const double*const*const dphiz, 
	double*const*const DGridV_x,  double*const*const DGridV_y,  double*const*const DGridV_z,
	double*const*const DGridV_11, double*const*const DGridV_12, double*const*const DGridV_13,
	double*const*const DGridV_22, double*const*const DGridV_23, double*const*const DGridV_33,
	realArray& drr)
{
    const char transa='N', transb='T';
    const double alpha=-1.0, beta=1.0;
    
    const int allnw=block_index[na_grid];
    for(int i=0; i<pw.bxyz; ++i)
    {
        for(int j=0; j<allnw; ++j)
        {
            psir_vlbr3[i][j]=psir_ylm[i][j]*vldr3[i];
        }
    }

    //OUT(ofs_running,"lgd_now", lgd_now);
    //OUT(ofs_running,"LD_pool", LD_pool);
    for(int ia1=0; ia1<na_grid; ++ia1)
    {
        const int iw1_lo=block_iw[ia1];
        const int idx1=block_index[ia1];
        const int m=block_size[ia1]; 
        // OUT(ofs_running,"ia1", ia1);
        // OUT(ofs_running,"iw1_lo", iw1_lo);
        // OUT(ofs_running,"m", m);
        for(int ia2=0; ia2<na_grid; ++ia2)
        {
            const int iw2_lo=block_iw[ia2];
            const int idx2=block_index[ia2];
            const int n=block_size[ia2];
			
           //  for(int ib=0; ib<pw.bxyz; ++ib)
           //  {
           //      if(cal_flag[ib][ia1]&&cal_flag[ib][ia2]) 
           //      {
        			// double* dphi2x = &dphix[ib][idx2];
        			// double* dphi2y = &dphiy[ib][idx2];
        			// double* dphi2z = &dphiz[ib][idx2];
        			// for(int iw1=0; iw1<m; ++iw1)
        			// {
           //              double psi=psir_vlbr3[ib][idx1+iw1];
           //              double* pDGridV_x=&DGridV_x[iw1_lo+iw1][iw2_lo];
           //              double* pDGridV_y=&DGridV_y[iw1_lo+iw1][iw2_lo];
           //              double* pDGridV_z=&DGridV_z[iw1_lo+iw1][iw2_lo];
           //  			for(int iw2=0; iw2<n; ++iw2)
           //  			{
           //      			pDGridV_x[iw2]-=psi*dphi2x[iw2];
           //      			pDGridV_y[iw2]-=psi*dphi2y[iw2];
           //      			pDGridV_z[iw2]-=psi*dphi2z[iw2];
           //  			}
        			// }
           //      }
           //  }

            
            int cal_num=0;
            for(int ib=0; ib<pw.bxyz; ++ib)
            {
                if(cal_flag[ib][ia1]&&cal_flag[ib][ia2]) ++cal_num;
            }
            //++cal_flag_true;
            //OUT(ofs_running,"cal_num:", cal_num);
            if (cal_num > pw.bxyz/2)
//            if(0)
            {
                int k=pw.bxyz;
                // OUT(ofs_running,"ia2", ia2);
                // OUT(ofs_running,"iw2_lo", iw2_lo);
                // OUT(ofs_running,"n", n);
                //std::cout<<"Start calculate DGridV_x"<<endl;
                //OUT(ofs_running,"Start calculate DGridV_x");
                dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                    &dphix[0][idx2], &LD_pool, 
                    &psir_vlbr3[0][idx1], &LD_pool,  
                    &beta, &DGridV_x[iw1_lo][iw2_lo], &lgd_now);
                //std::cout<<"Start calculate DGridV_y"<<endl;
                //OUT(ofs_running,"Start calculate DGridV_y");
                dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                    &dphiy[0][idx2], &LD_pool, 
                    &psir_vlbr3[0][idx1], &LD_pool,  
                    &beta, &DGridV_y[iw1_lo][iw2_lo], &lgd_now);
                //std::cout<<"Start calculate DGridV_z"<<endl;
                //OUT(ofs_running,"Start calculate DGridV_z");
                dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                    &dphiz[0][idx2], &LD_pool, 
                    &psir_vlbr3[0][idx1], &LD_pool,  
                    &beta, &DGridV_z[iw1_lo][iw2_lo], &lgd_now);

                if(STRESS)
                {
					k=1;
					for(int ib=0; ib<pw.bxyz; ++ib)
					{
						double stress_alpha1 = alpha * drr(ia2,ib,0);
						double stress_alpha2 = alpha * drr(ia2,ib,1);
						double stress_alpha3 = alpha * drr(ia2,ib,2);
						dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha1,
							&dphix[ib][idx2], &LD_pool,
							&psir_vlbr3[ib][idx1], &LD_pool,
							&beta, &DGridV_11[iw1_lo][iw2_lo], &lgd_now);
						dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha2,
							&dphix[ib][idx2], &LD_pool,
							&psir_vlbr3[ib][idx1], &LD_pool,
							&beta, &DGridV_12[iw1_lo][iw2_lo], &lgd_now);
						dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha3,
							&dphix[ib][idx2], &LD_pool,
							&psir_vlbr3[ib][idx1], &LD_pool,
							&beta, &DGridV_13[iw1_lo][iw2_lo], &lgd_now);
						dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha2,
							&dphiy[ib][idx2], &LD_pool,
							&psir_vlbr3[ib][idx1], &LD_pool,
							&beta, &DGridV_22[iw1_lo][iw2_lo], &lgd_now);
						dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha3,
							&dphiy[ib][idx2], &LD_pool,
							&psir_vlbr3[ib][idx1], &LD_pool,
							&beta, &DGridV_23[iw1_lo][iw2_lo], &lgd_now);
						dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha3,
							&dphiz[ib][idx2], &LD_pool,
							&psir_vlbr3[ib][idx1], &LD_pool,
							&beta, &DGridV_33[iw1_lo][iw2_lo], &lgd_now);
					}
//                DGridV_11[iw1_lo*GridT.lgd + iw2_lo] += DGridV_x[iw1_lo][iw2_lo] * drr[ia2][0][0];
//                DGridV_12[iw1_lo*GridT.lgd + iw2_lo] += DGridV_x[iw1_lo][iw2_lo] * drr[ia2][0][1];
//                DGridV_13[iw1_lo*GridT.lgd + iw2_lo] += DGridV_x[iw1_lo][iw2_lo] * drr[ia2][0][2];
//                DGridV_22[iw1_lo*GridT.lgd + iw2_lo] += DGridV_y[iw1_lo][iw2_lo] * drr[ia2][0][1];
//                DGridV_23[iw1_lo*GridT.lgd + iw2_lo] += DGridV_y[iw1_lo][iw2_lo] * drr[ia2][0][2];
//                DGridV_33[iw1_lo*GridT.lgd + iw2_lo] += DGridV_z[iw1_lo][iw2_lo] * drr[ia2][0][2];
                }
            }
            else if (cal_num > 0)
            {
                int k=1;
                for(int ib=0; ib<pw.bxyz; ++ib)
                {
                    if(cal_flag[ib][ia1]&&cal_flag[ib][ia2])
                    {
                        //OUT(ofs_running,"Start calculate DGridV_x");
                        dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                            &dphix[ib][idx2], &LD_pool, 
                            &psir_vlbr3[ib][idx1], &LD_pool,  
                            &beta, &DGridV_x[iw1_lo][iw2_lo], &lgd_now);
                        //std::cout<<"Start calculate DGridV_y"<<endl;
                        //OUT(ofs_running,"Start calculate DGridV_y");
                        dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                            &dphiy[ib][idx2], &LD_pool, 
                            &psir_vlbr3[ib][idx1], &LD_pool,  
                            &beta, &DGridV_y[iw1_lo][iw2_lo], &lgd_now);
                        //std::cout<<"Start calculate DGridV_z"<<endl;
                        //OUT(ofs_running,"Start calculate DGridV_z");
                        dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                            &dphiz[ib][idx2], &LD_pool, 
                            &psir_vlbr3[ib][idx1], &LD_pool,  
                            &beta, &DGridV_z[iw1_lo][iw2_lo], &lgd_now);

                        if(STRESS)
						{
							double stress_alpha1 = alpha * drr(ia2,ib,0);
							double stress_alpha2 = alpha * drr(ia2,ib,1);
							double stress_alpha3 = alpha * drr(ia2,ib,2);
							dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha1,
									&dphix[ib][idx2], &LD_pool,
									&psir_vlbr3[ib][idx1], &LD_pool,
									&beta, &DGridV_11[iw1_lo][iw2_lo], &lgd_now);
							dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha2,
									&dphix[ib][idx2], &LD_pool,
									&psir_vlbr3[ib][idx1], &LD_pool,
									&beta, &DGridV_12[iw1_lo][iw2_lo], &lgd_now);
							dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha3,
									&dphix[ib][idx2], &LD_pool,
									&psir_vlbr3[ib][idx1], &LD_pool,
									&beta, &DGridV_13[iw1_lo][iw2_lo], &lgd_now);
							dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha2,
									&dphiy[ib][idx2], &LD_pool,
									&psir_vlbr3[ib][idx1], &LD_pool,
									&beta, &DGridV_22[iw1_lo][iw2_lo], &lgd_now);
							dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha3,
									&dphiy[ib][idx2], &LD_pool,
									&psir_vlbr3[ib][idx1], &LD_pool,
									&beta, &DGridV_23[iw1_lo][iw2_lo], &lgd_now);
							dgemm_ (&transa, &transb, &n, &m, &k, &stress_alpha3,
									&dphiz[ib][idx2], &LD_pool,
									&psir_vlbr3[ib][idx1], &LD_pool,
									&beta, &DGridV_33[iw1_lo][iw2_lo], &lgd_now);
							// DGridV_11[iw1_lo*GridT.lgd + iw2_lo] += DGridV_x[iw1_lo][iw2_lo] * drr[ia2][ib][0];
							// DGridV_12[iw1_lo*GridT.lgd + iw2_lo] += DGridV_x[iw1_lo][iw2_lo] * drr[ia2][ib][1];
							// DGridV_13[iw1_lo*GridT.lgd + iw2_lo] += DGridV_x[iw1_lo][iw2_lo] * drr[ia2][ib][2];
							// DGridV_22[iw1_lo*GridT.lgd + iw2_lo] += DGridV_y[iw1_lo][iw2_lo] * drr[ia2][ib][1];
							// DGridV_23[iw1_lo*GridT.lgd + iw2_lo] += DGridV_y[iw1_lo][iw2_lo] * drr[ia2][ib][2];
							// DGridV_33[iw1_lo*GridT.lgd + iw2_lo] += DGridV_z[iw1_lo][iw2_lo] * drr[ia2][ib][2];
							//cout<<"DGridV: "<<iw1_lo<<" "<<iw2_lo
							//<<" "<<iw1_lo*GridT.lgd + iw2_lo<<" "
							//<<DGridV_x[iw1_lo][iw2_lo]<<" "
							//<<DGridV_y[iw1_lo][iw2_lo]<<" "
							//<<DGridV_z[iw1_lo][iw2_lo]<<" "
							//<<drr[ia2][0][0]<<" "<<drr[ia2][0][1]<<" "<<drr[ia2][0][2]<<endl;
						}
					}                    
                }
            }
        }
    }
}


void Gint_Gamma::gamma_force(void)
{
    TITLE("Grid_Integral","gamma_force");
    timer::tick("Gint_Gamma","gamma_force",'I');
    // GridT.lgd: local grid dimension (sub-FFT-mesh).
    int DGridV_Size=GridT.lgd*GridT.lgd;
    //OUT(ofs_running,"Enter gamma_force, DGridV_Size", DGridV_Size);
    double *DGridV_pool=new double[3*DGridV_Size];
    ZEROS(DGridV_pool, 3*DGridV_Size);
    
    double** DGridV_x = new double*[GridT.lgd];
    double** DGridV_y = new double*[GridT.lgd];
    double** DGridV_z = new double*[GridT.lgd];
    double* DGridV_stress_pool;
    double** DGridV_11;
    double** DGridV_12;
    double** DGridV_13;
    double** DGridV_22;
    double** DGridV_23;
    double** DGridV_33;

    if(STRESS)
    {
        DGridV_stress_pool = new double[6*DGridV_Size];
        ZEROS(DGridV_stress_pool, 6*DGridV_Size);
        DGridV_11 = new double*[GridT.lgd];
        DGridV_12 = new double*[GridT.lgd];
        DGridV_13 = new double*[GridT.lgd];
        DGridV_22 = new double*[GridT.lgd];
        DGridV_23 = new double*[GridT.lgd];
        DGridV_33 = new double*[GridT.lgd];
        for (int i=0; i<GridT.lgd; ++i)
        {
            DGridV_11[i] = &DGridV_stress_pool[i*GridT.lgd];
            DGridV_12[i] = &DGridV_stress_pool[i*GridT.lgd+DGridV_Size];
            DGridV_13[i] = &DGridV_stress_pool[i*GridT.lgd+2*DGridV_Size];
            DGridV_22[i] = &DGridV_stress_pool[i*GridT.lgd+3*DGridV_Size];
            DGridV_23[i] = &DGridV_stress_pool[i*GridT.lgd+4*DGridV_Size];
            DGridV_33[i] = &DGridV_stress_pool[i*GridT.lgd+5*DGridV_Size];
        }
        Memory::record("Gint_Gamma","DGridV_stress",6*GridT.lgd*GridT.lgd,"double");
    }
    for (int i=0; i<GridT.lgd; ++i)
    {
        DGridV_x[i] = &DGridV_pool[i*GridT.lgd];
        DGridV_y[i] = &DGridV_pool[i*GridT.lgd+DGridV_Size];
        DGridV_z[i] = &DGridV_pool[i*GridT.lgd+2*DGridV_Size];
    }
    Memory::record("Gint_Gamma","DGridV",3*GridT.lgd*GridT.lgd,"double");
    //OUT(ofs_running,"DGridV was allocated");

    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const double delta_r = ORB.dr_uniform;

    int LD_pool=max_size*ucell.nwmax;
    double* psir_ylm_pool;
    double* psir_vlbr3_pool;
    double* dphi_pool;
    
    double** psir_ylm;
    double** psir_vlbr3;
    double** dphix;
    double** dphiy;
    double** dphiz;

    bool** cal_flag;

    const int ncyz=pw.ncy*pw.nczp;

/*    if(max_size<=0 || GridT.lgd <= 0) 
    {
      //OUT(ofs_running,"max_size", max_size);
      //OUT(ofs_running,"GridT.lgd", GridT.lgd);
        goto ENDandRETURN;
    }*/
    if(max_size>0 && GridT.lgd > 0)
    {
    
        psir_ylm_pool=new double [pw.bxyz*LD_pool];
        ZEROS(psir_ylm_pool, pw.bxyz*LD_pool);
        psir_vlbr3_pool=new double[pw.bxyz*LD_pool];
        ZEROS(psir_vlbr3_pool, pw.bxyz*LD_pool);
        dphi_pool=new double [3*pw.bxyz*LD_pool];
        ZEROS(dphi_pool, 3*pw.bxyz*LD_pool);
        psir_ylm = new double*[pw.bxyz];
        psir_vlbr3=new double *[pw.bxyz];
        dphix = new double*[pw.bxyz];
        dphiy = new double*[pw.bxyz];
        dphiz = new double*[pw.bxyz];    
        
        cal_flag=new bool*[pw.bxyz];
        for(int i=0; i<pw.bxyz; i++)
        {
            psir_ylm[i] = &psir_ylm_pool[i*LD_pool];
            psir_vlbr3[i]=&psir_vlbr3_pool[i*LD_pool];
            dphix[i] = &dphi_pool[i*LD_pool];
            dphiy[i] = &dphi_pool[i*LD_pool+pw.bxyz*LD_pool];
            dphiz[i] = &dphi_pool[i*LD_pool+2*pw.bxyz*LD_pool];
            cal_flag[i] = new bool[max_size];
        }

        realArray drr;//rewrite drr form by zhengdy-2019-04-02
        if(STRESS)
        {
            drr.create(max_size, pw.bxyz, 3);
            drr.zero_out();
        }    
/*        double ***drr;//store dr for stress calculate, added by zhengdy
        if(STRESS)//added by zhengdy-stress
        {
    		drr = new double**[max_size];
    		for(int id=0; id<max_size; id++)
    		{
    			drr[id] = new double*[pw.bxyz];
    			for(int ib=0; ib<pw.bxyz; ib++)
    			{
    				drr[id][ib] = new double[3];
    				ZEROS(drr[id][ib],3);
    			}
    		}
        }*/
        //OUT(ofs_running,"Data were prepared");
        //timer::tick("Gint_Gamma","prepare",'J');
        for (int i=0; i< GridT.nbx; i++)
        {
            const int ibx = i*pw.bx; 

            for (int j=0; j< GridT.nby; j++)
            {
                const int jby = j*pw.by;

                for (int k= GridT.nbzp_start; k< GridT.nbzp_start+GridT.nbzp; k++)
                {
                    const int kbz = k*pw.bz-pw.nczp_start; 
                    const int grid_index = (k-GridT.nbzp_start) + j * GridT.nbzp + i * GridT.nby * GridT.nbzp;
                    const int na_grid = GridT.how_many_atoms[ grid_index ];
                    if(na_grid==0)continue;
                    
					//------------------------------------------------------------------
					// extract the local potentials.
					//------------------------------------------------------------------
					double *vldr3 = get_vldr3(ncyz, ibx, jby, kbz);
					
					//------------------------------------------------------
					// index of wave functions for each block
					//------------------------------------------------------
					int *block_iw = get_block_iw(na_grid, grid_index, this->max_size);
					
					int* block_index = get_colidx(na_grid, grid_index);
					
					//------------------------------------------------------
					// band size: number of columns of a band
					//------------------------------------------------------------------
					int* block_size = get_bsize(na_grid, grid_index);
    
                    cal_psir_ylm_dphi(na_grid, grid_index, delta_r, 
                            block_index, block_size, cal_flag, psir_ylm, dphix, dphiy, dphiz, drr);
    
                    cal_meshball_DGridV(na_grid, GridT.lgd, LD_pool, block_index, block_iw, block_size, cal_flag, vldr3, 
                                psir_ylm, psir_vlbr3, dphix,  dphiy, dphiz, 
                                DGridV_x, DGridV_y, DGridV_z,
                                DGridV_11, DGridV_12, DGridV_13,
                                DGridV_22, DGridV_23, DGridV_33, drr);
								
					free(vldr3);		vldr3=nullptr;
					free(block_iw);		block_iw=nullptr;
					free(block_index);	block_index=nullptr;
					free(block_size);	block_size=nullptr;
                }// k
            }// j
        }// i
    
        //OUT(ofs_running,"DGridV was calculated");
        delete[] dphix;
        delete[] dphiy;
        delete[] dphiz;
        delete[] dphi_pool;
        delete[] psir_ylm;
        delete[] psir_ylm_pool;
        delete[] psir_vlbr3;
        delete[] psir_vlbr3_pool;
        for(int ib=0; ib<pw.bxyz; ++ib)
		{
            delete[] cal_flag[ib];
		}
        delete[] cal_flag;
        //OUT(ofs_running,"temp variables were deleted");

    }//end if, replace goto line
//ENDandRETURN:
    timer::tick("Gint_Gamma","gamma_force",'I');
#ifdef __MPI
    timer::tick("Gint_Gamma","gamma_force_wait",'I');
	MPI_Barrier(MPI_COMM_WORLD);
    timer::tick("Gint_Gamma","gamma_force_wait",'I');
#endif
    timer::tick("Gint_Gamma","gamma_force2",'I');


    //OUT(ofs_running,"Start reduce DGridV");

    double* tmpx = new double[NLOCAL];
    double* tmpy = new double[NLOCAL];
    double* tmpz = new double[NLOCAL];
    double* tmp11;
    double* tmp12;
    double* tmp13;
    double* tmp22;
    double* tmp23;
    double* tmp33;

	if(STRESS)
	{
		tmp11 = new double[NLOCAL];
		tmp12 = new double[NLOCAL];
		tmp13 = new double[NLOCAL];
		tmp22 = new double[NLOCAL];
		tmp23 = new double[NLOCAL];
		tmp33 = new double[NLOCAL];
	}

    for (int i=0; i<NLOCAL; i++)
    {
        ZEROS(tmpx, NLOCAL);
        ZEROS(tmpy, NLOCAL);
		ZEROS(tmpz, NLOCAL);

		if(STRESS)
		{
			ZEROS(tmp11, NLOCAL);
			ZEROS(tmp12, NLOCAL);
			ZEROS(tmp13, NLOCAL);
			ZEROS(tmp22, NLOCAL);
			ZEROS(tmp23, NLOCAL);
			ZEROS(tmp33, NLOCAL);
		}

        const int mu = GridT.trace_lo[i];
        // mohan fix bug 2010-09-05
        // lack mu>=0 and nu>=0 in previous version.
        if(mu >=0)
        {
            for (int j=0; j<NLOCAL; j++)
            {
                const int nu = GridT.trace_lo[j];
                if(nu>=0)
                {
                    tmpx[j] = DGridV_x[mu][nu];
                    tmpy[j] = DGridV_y[mu][nu];
                    tmpz[j] = DGridV_z[mu][nu];
                    if(STRESS)
                    {
                        tmp11[j] = DGridV_11[mu][nu];
                        tmp12[j] = DGridV_12[mu][nu];
                        tmp13[j] = DGridV_13[mu][nu];
                        tmp22[j] = DGridV_22[mu][nu];
                        tmp23[j] = DGridV_23[mu][nu];
                        tmp33[j] = DGridV_33[mu][nu];
                    }
                }
            }
        }

        // There may be overlaps of tmpx,y,z between different
        // processors, however, the true value is the sum of it.
        // so it would be totally correct.
        Parallel_Reduce::reduce_double_pool( tmpx, NLOCAL );
        Parallel_Reduce::reduce_double_pool( tmpy, NLOCAL );
        Parallel_Reduce::reduce_double_pool( tmpz, NLOCAL );
		if(STRESS)
		{
			Parallel_Reduce::reduce_double_pool( tmp11, NLOCAL );
			Parallel_Reduce::reduce_double_pool( tmp12, NLOCAL );
			Parallel_Reduce::reduce_double_pool( tmp13, NLOCAL );
			Parallel_Reduce::reduce_double_pool( tmp22, NLOCAL );
			Parallel_Reduce::reduce_double_pool( tmp23, NLOCAL );
			Parallel_Reduce::reduce_double_pool( tmp33, NLOCAL );
		}

        for (int j=0; j<NLOCAL; j++)
        {
            if (!ParaO.in_this_processor(i,j))
            {
                continue;
            }
            LM.set_force (i,j,tmpx[j], tmpy[j], tmpz[j],'N');
            if(STRESS)
            {
                const int irr = ParaO.trace_loc_row[ i ];
                const int icc = ParaO.trace_loc_col[ j ];
                const int index = irr * ParaO.ncol + icc;
                LM.DHloc_fixed_11[index] += tmp11[j];
                LM.DHloc_fixed_12[index] += tmp12[j];
                LM.DHloc_fixed_13[index] += tmp13[j];
                LM.DHloc_fixed_22[index] += tmp22[j];
                LM.DHloc_fixed_23[index] += tmp23[j];
                LM.DHloc_fixed_33[index] += tmp33[j];
            }
        }
    }
    delete[] tmpx;
    delete[] tmpy;
    delete[] tmpz;
    if(STRESS)
    {
        delete[] tmp11;
        delete[] tmp12;
        delete[] tmp13;
        delete[] tmp22;
        delete[] tmp23;
        delete[] tmp33;
    }
    //OUT(ofs_running,"DGridV was reduced");

   //OUT(ofs_running,"Start reduce DGridV");
   //Parallel_Reduce::reduce_double_pool(DGridV_pool, 3*DGridV_Size );
    // double* tmp = new double[3*NLOCAL*NLOCAL];
    // ZEROS(tmp, 3*NLOCAL*NLOCAL);
    // for (int i=0; i<NLOCAL; i++)
    // {
    //     ZEROS(tmp, 3*NLOCAL);
    //     double* tmpx = &tmp[i*NLOCAL];
    //     double* tmpy = &tmp[i*NLOCAL+NLOCAL*NLOCAL];
    //     double* tmpz = &tmp[i*NLOCAL+2*NLOCAL*NLOCAL];
    //     if(DGridV_Size>0)
    //     {
    //         const int mu = GridT.trace_lo[i];
    //         if(mu >=0)
    //         {
    //             for (int j=0; j<NLOCAL; j++)
    //             {
    //                 const int nu = GridT.trace_lo[j];
    //                 if(nu>=0)
    //                 {
    //                     tmpx[j] = DGridV_x[mu][nu];
    //                     tmpy[j] = DGridV_y[mu][nu];
    //                     tmpz[j] = DGridV_z[mu][nu];
    //                 }
    //             }
    //         }
    //     }
    // }
    
    // Parallel_Reduce::reduce_double_pool(tmp, 3*NLOCAL*NLOCAL);

    // for (int i=0; i<NLOCAL; i++)
    // {
    //     for(int j=0; j<NLOCAL; ++j)
    //     {
    //         if(ParaO.in_this_processor(i,j))
    //             LM.set_force (i,j,DGridV_x[i][j], DGridV_y[i][j], DGridV_z[i][j],'N');
    //     }
    // }
    // delete[] tmp;
    timer::tick("Gint_Gamma","gamma_force2",'I');   
    
    //delete DGridV_x,y,z
    delete [] DGridV_x;
    delete [] DGridV_y;
    delete [] DGridV_z;
    if(STRESS)
    {
        delete [] DGridV_11;
        delete [] DGridV_12;
        delete [] DGridV_13;
        delete [] DGridV_22;
        delete [] DGridV_23;
        delete [] DGridV_33;
        delete [] DGridV_stress_pool;
    }
    delete [] DGridV_pool;
    return;
}
