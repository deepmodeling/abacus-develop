//=========================================================
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#include "gint_gamma.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "../module_orbital/ORB_read.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

#include "global_fp.h" // mohan add 2021-01-30
#include "../module_base/ylm.h"

// calcualte the forces related to grid
void Gint_Gamma::cal_force(const double*const vlocal)
{
    ModuleBase::timer::tick("Gint_Gamma","cal_force");
    this->save_atoms_on_grid(GlobalC::GridT);
    this->gamma_force(vlocal);
    ModuleBase::timer::tick("Gint_Gamma","cal_force");
}

inline void cal_psir_ylm_dphi(
	const int na_grid,   					// how many atoms on this (i,j,k) grid
	const int grid_index, 
	const double delta_r,
	const int*const block_index,			// block_index[na_grid+1], count total number of atomis orbitals
	const int*const block_size,  			// block_size[na_grid],	number of columns of a band
	bool*const*const cal_flag,		        // cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
	double*const*const psir_ylm,	    	// psir_ylm[GlobalC::pw.bxyz][LD_pool] 
	double*const*const dphix, 
	double*const*const dphiy, 
	double*const*const dphiz,
	ModuleBase::realArray& drr)
{    
    for (int id=0; id<na_grid; id++)
    {
        const int mcell_index = GlobalC::GridT.bcell_start[grid_index] + id;
        const int imcell = GlobalC::GridT.which_bigcell[mcell_index];
        int iat = GlobalC::GridT.which_atom[mcell_index];
        const int it = GlobalC::ucell.iat2it[iat];
        const int ia = GlobalC::ucell.iat2ia[iat];
        const int start = GlobalC::ucell.itiaiw2iwt(it, ia, 0);
        Atom *atom = &GlobalC::ucell.atoms[it];

		const double mt[3]={
			GlobalC::GridT.meshball_positions[imcell][0] - GlobalC::GridT.tau_in_bigcell[iat][0],
			GlobalC::GridT.meshball_positions[imcell][1] - GlobalC::GridT.tau_in_bigcell[iat][1],
			GlobalC::GridT.meshball_positions[imcell][2] - GlobalC::GridT.tau_in_bigcell[iat][2]};

        for(int ib=0; ib<GlobalC::pw.bxyz; ib++)
        {
            double*const p_psir_ylm=&psir_ylm[ib][block_index[id]];
            double*const p_dphix=&dphix[ib][block_index[id]];
            double*const p_dphiy=&dphiy[ib][block_index[id]];
            double*const p_dphiz=&dphiz[ib][block_index[id]];
            
			const double dr[3]={						// vectors between atom and grid
				GlobalC::GridT.meshcell_pos[ib][0] + mt[0],
				GlobalC::GridT.meshcell_pos[ib][1] + mt[1],
				GlobalC::GridT.meshcell_pos[ib][2] + mt[2]};

            if(GlobalV::CAL_STRESS)
            {
                for(int i=0;i<3;i++) 
				{
					drr(id,ib,i) = dr[i];
				}
            }

			// distance between atom and grid
            double distance = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

            if(distance > GlobalC::ORB.Phi[it].getRcut())
            {
                ModuleBase::GlobalFunc::ZEROS(p_psir_ylm, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_dphix, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_dphiy, block_size[id]);
                ModuleBase::GlobalFunc::ZEROS(p_dphiz, block_size[id]);
                cal_flag[ib][id]=false;
                continue;
            }

            cal_flag[ib][id]=true;
            
            //array to store spherical harmonics and its derivatives
            std::vector<double> rly;
            std::vector<std::vector<double>> grly;
            // Ylm::get_ylm_real(GlobalC::ucell.atoms[it].nwl+1, this->dr[id], ylma[id]);
            // <<<
            // Ylm::rlylm(GlobalC::ucell.atoms[it].nwl+1, dr[id].x, dr[id].y, dr[id].z, rly, grly);
            // Ylm::rlylm(GlobalC::ucell.atoms[it].nwl+1, dr[id].x, dr[id].y, dr[id].z, rly, grly);
            ModuleBase::Ylm::grad_rl_sph_harm(GlobalC::ucell.atoms[it].nwl, dr[0], dr[1], dr[2], rly, grly);

            // 1E-7 is necessary in case of R is just on one grid
            // the following code is about interpolation,
            // here we use Polynomia interpolation,
            // iq, x0, x1, x3, x3, x12, x03 are all the variables 
            // needed in this interpolation method.
            // It's the same as the function in mathzone.cpp:
            // ModuleBase::Mathzone::Polynomial_Interpolation
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
                    const Numerical_Orbital_Lm &philn = GlobalC::ORB.Phi[it].PhiLN(
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
                        // it's a std::vector, so it has x,y,z components.
                        p_dphix[iw] = 0.0;
                        p_dphiy[iw] = 0.0;
                        p_dphiz[iw] = 0.0;
                    }
                    // if l>0, how do we choose 3D localized wave function
                    // and the derivative.
                    else
                    {
                        const Numerical_Orbital_Lm &philn = GlobalC::ORB.Phi[it].PhiLN(atom->iw2l[iw], atom->iw2n[iw]);

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
	const int na_grid,   			    	    	// how many atoms on this (i,j,k) grid
	int lgd_now,
	int LD_pool, 
	const int*const block_index,	    	    	// block_index[na_grid+1], count total number of atomis orbitals
	const int*const block_iw, 		    	    	// block_iw[na_grid],	index of wave functions for each block
	const int*const block_size, 	    	    	// block_size[na_grid],	number of columns of a band
	const bool*const*const cal_flag,		        // cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
	const double*const vldr3,               	    // vldr3[GlobalC::pw.bxyz]
	const double*const*const psir_ylm,		        // psir_ylm[GlobalC::pw.bxyz][LD_pool]
	double*const*const psir_vlbr3,              	// psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
	const double*const*const dphix, const double*const*const dphiy, const double*const*const dphiz, 
	double*const*const DGridV_x,  double*const*const DGridV_y,  double*const*const DGridV_z,
	double*const*const DGridV_11, double*const*const DGridV_12, double*const*const DGridV_13,
	double*const*const DGridV_22, double*const*const DGridV_23, double*const*const DGridV_33,
	ModuleBase::realArray& drr)
{
    const char transa='N', transb='T';
    const double alpha=-1.0, beta=1.0;
    
    const int allnw=block_index[na_grid];
    for(int i=0; i<GlobalC::pw.bxyz; ++i)
    {
        for(int j=0; j<allnw; ++j)
        {
            psir_vlbr3[i][j]=psir_ylm[i][j]*vldr3[i];
        }
    }

    //OUT(GlobalV::ofs_running,"lgd_now", lgd_now);
    //OUT(GlobalV::ofs_running,"LD_pool", LD_pool);
    for(int ia1=0; ia1<na_grid; ++ia1)
    {
        const int iw1_lo=block_iw[ia1];
        const int idx1=block_index[ia1];
        const int m=block_size[ia1]; 
        // OUT(GlobalV::ofs_running,"ia1", ia1);
        // OUT(GlobalV::ofs_running,"iw1_lo", iw1_lo);
        // OUT(GlobalV::ofs_running,"m", m);
        for(int ia2=0; ia2<na_grid; ++ia2)
        {
            const int iw2_lo=block_iw[ia2];
            const int idx2=block_index[ia2];
            const int n=block_size[ia2];
			           
            int cal_num=0;
            for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
            {
                if(cal_flag[ib][ia1]&&cal_flag[ib][ia2]) ++cal_num;
            }
            //++cal_flag_true;
            //OUT(GlobalV::ofs_running,"cal_num:", cal_num);
            if (cal_num > GlobalC::pw.bxyz/2)
//            if(0)
            {
                int k=GlobalC::pw.bxyz;
                dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                    &dphix[0][idx2], &LD_pool, 
                    &psir_vlbr3[0][idx1], &LD_pool,  
                    &beta, &DGridV_x[iw1_lo][iw2_lo], &lgd_now);
                dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                    &dphiy[0][idx2], &LD_pool, 
                    &psir_vlbr3[0][idx1], &LD_pool,  
                    &beta, &DGridV_y[iw1_lo][iw2_lo], &lgd_now);
                dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                    &dphiz[0][idx2], &LD_pool, 
                    &psir_vlbr3[0][idx1], &LD_pool,  
                    &beta, &DGridV_z[iw1_lo][iw2_lo], &lgd_now);

                if(GlobalV::CAL_STRESS)
                {
					k=1;
					for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
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
                }
            }
            else if (cal_num > 0)
            {
                int k=1;
                for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
                {
                    if(cal_flag[ib][ia1]&&cal_flag[ib][ia2])
                    {
                        //OUT(GlobalV::ofs_running,"Start calculate DGridV_x");
                        dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                            &dphix[ib][idx2], &LD_pool, 
                            &psir_vlbr3[ib][idx1], &LD_pool,  
                            &beta, &DGridV_x[iw1_lo][iw2_lo], &lgd_now);
                        //std::cout<<"Start calculate DGridV_y"<<std::endl;
                        //OUT(GlobalV::ofs_running,"Start calculate DGridV_y");
                        dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                            &dphiy[ib][idx2], &LD_pool, 
                            &psir_vlbr3[ib][idx1], &LD_pool,  
                            &beta, &DGridV_y[iw1_lo][iw2_lo], &lgd_now);
                        //std::cout<<"Start calculate DGridV_z"<<std::endl;
                        //OUT(GlobalV::ofs_running,"Start calculate DGridV_z");
                        dgemm_ (&transa, &transb, &n, &m, &k, &alpha,
                            &dphiz[ib][idx2], &LD_pool, 
                            &psir_vlbr3[ib][idx1], &LD_pool,  
                            &beta, &DGridV_z[iw1_lo][iw2_lo], &lgd_now);

                        if(GlobalV::CAL_STRESS)
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
					}                    
                }
            }
        }
    }
}


void Gint_Gamma::gamma_force(const double*const vlocal) const
{
    ModuleBase::TITLE("Grid_Integral","gamma_force");
    ModuleBase::timer::tick("Gint_Gamma","gamma_force");
    // GlobalC::GridT.lgd: local grid dimension (sub-FFT-mesh).
    int DGridV_Size=GlobalC::GridT.lgd*GlobalC::GridT.lgd;
    //OUT(GlobalV::ofs_running,"Enter gamma_force, DGridV_Size", DGridV_Size);
    double *DGridV_pool=new double[3*DGridV_Size];
    ModuleBase::GlobalFunc::ZEROS(DGridV_pool, 3*DGridV_Size);
    
    double** DGridV_x = new double*[GlobalC::GridT.lgd];
    double** DGridV_y = new double*[GlobalC::GridT.lgd];
    double** DGridV_z = new double*[GlobalC::GridT.lgd];
    double* DGridV_stress_pool;
    double** DGridV_11;
    double** DGridV_12;
    double** DGridV_13;
    double** DGridV_22;
    double** DGridV_23;
    double** DGridV_33;

    if(GlobalV::CAL_STRESS)
    {
        DGridV_stress_pool = new double[6*DGridV_Size];
        ModuleBase::GlobalFunc::ZEROS(DGridV_stress_pool, 6*DGridV_Size);
        DGridV_11 = new double*[GlobalC::GridT.lgd];
        DGridV_12 = new double*[GlobalC::GridT.lgd];
        DGridV_13 = new double*[GlobalC::GridT.lgd];
        DGridV_22 = new double*[GlobalC::GridT.lgd];
        DGridV_23 = new double*[GlobalC::GridT.lgd];
        DGridV_33 = new double*[GlobalC::GridT.lgd];
        for (int i=0; i<GlobalC::GridT.lgd; ++i)
        {
            DGridV_11[i] = &DGridV_stress_pool[i*GlobalC::GridT.lgd];
            DGridV_12[i] = &DGridV_stress_pool[i*GlobalC::GridT.lgd+DGridV_Size];
            DGridV_13[i] = &DGridV_stress_pool[i*GlobalC::GridT.lgd+2*DGridV_Size];
            DGridV_22[i] = &DGridV_stress_pool[i*GlobalC::GridT.lgd+3*DGridV_Size];
            DGridV_23[i] = &DGridV_stress_pool[i*GlobalC::GridT.lgd+4*DGridV_Size];
            DGridV_33[i] = &DGridV_stress_pool[i*GlobalC::GridT.lgd+5*DGridV_Size];
        }
        ModuleBase::Memory::record("Gint_Gamma","DGridV_stress",6*GlobalC::GridT.lgd*GlobalC::GridT.lgd,"double");
    }
    for (int i=0; i<GlobalC::GridT.lgd; ++i)
    {
        DGridV_x[i] = &DGridV_pool[i*GlobalC::GridT.lgd];
        DGridV_y[i] = &DGridV_pool[i*GlobalC::GridT.lgd+DGridV_Size];
        DGridV_z[i] = &DGridV_pool[i*GlobalC::GridT.lgd+2*DGridV_Size];
    }
    ModuleBase::Memory::record("Gint_Gamma","DGridV",3*GlobalC::GridT.lgd*GlobalC::GridT.lgd,"double");
#ifdef _OPENMP
    #pragma omp parallel
    {
        double *DGridV_pool_thread=new double[3*DGridV_Size];
        ModuleBase::GlobalFunc::ZEROS(DGridV_pool_thread, 3*DGridV_Size);

        double** DGridV_x_thread = new double*[GlobalC::GridT.lgd];
        double** DGridV_y_thread = new double*[GlobalC::GridT.lgd];
        double** DGridV_z_thread = new double*[GlobalC::GridT.lgd];
        double* DGridV_stress_pool_thread;
        double** DGridV_11_thread;
        double** DGridV_12_thread;
        double** DGridV_13_thread;
        double** DGridV_22_thread;
        double** DGridV_23_thread;
        double** DGridV_33_thread;

        if(GlobalV::CAL_STRESS)
        {
            DGridV_stress_pool_thread = new double[6*DGridV_Size];
            ModuleBase::GlobalFunc::ZEROS(DGridV_stress_pool_thread, 6*DGridV_Size);
            DGridV_11_thread = new double*[GlobalC::GridT.lgd];
            DGridV_12_thread = new double*[GlobalC::GridT.lgd];
            DGridV_13_thread = new double*[GlobalC::GridT.lgd];
            DGridV_22_thread = new double*[GlobalC::GridT.lgd];
            DGridV_23_thread = new double*[GlobalC::GridT.lgd];
            DGridV_33_thread = new double*[GlobalC::GridT.lgd];

            for(int i=0; i<GlobalC::GridT.lgd; ++i)
            {
                DGridV_11_thread[i] = &DGridV_stress_pool_thread[i*GlobalC::GridT.lgd];
                DGridV_12_thread[i] = &DGridV_stress_pool_thread[i*GlobalC::GridT.lgd+DGridV_Size];
                DGridV_13_thread[i] = &DGridV_stress_pool_thread[i*GlobalC::GridT.lgd+2*DGridV_Size];
                DGridV_22_thread[i] = &DGridV_stress_pool_thread[i*GlobalC::GridT.lgd+3*DGridV_Size];
                DGridV_23_thread[i] = &DGridV_stress_pool_thread[i*GlobalC::GridT.lgd+4*DGridV_Size];
                DGridV_33_thread[i] = &DGridV_stress_pool_thread[i*GlobalC::GridT.lgd+5*DGridV_Size];
            }
            ModuleBase::Memory::record("Gint_Gamma","DGridV_stress_thread",6*GlobalC::GridT.lgd*GlobalC::GridT.lgd,"double");
        }

        for(int i=0; i<GlobalC::GridT.lgd; ++i)
        {
            DGridV_x_thread[i] = &DGridV_pool_thread[i*GlobalC::GridT.lgd];
            DGridV_y_thread[i] = &DGridV_pool_thread[i*GlobalC::GridT.lgd+DGridV_Size];
            DGridV_z_thread[i] = &DGridV_pool_thread[i*GlobalC::GridT.lgd+2*DGridV_Size];
        }
        ModuleBase::Memory::record("Gint_Gamma","DGridV_thread",3*GlobalC::GridT.lgd*GlobalC::GridT.lgd,"double");
#endif
        // it's a uniform grid to save orbital values, so the delta_r is a constant.
        const double delta_r = GlobalC::ORB.dr_uniform;

        int LD_pool=max_size*GlobalC::ucell.nwmax;
        double* dphi_pool;

        double** dphix;
        double** dphiy;
        double** dphiz;

        bool** cal_flag;
        const int ncyz=GlobalC::pw.ncy*GlobalC::pw.nczp;

        if(max_size>0 && GlobalC::GridT.lgd > 0)
        {
            dphi_pool=new double [3*GlobalC::pw.bxyz*LD_pool];
            ModuleBase::GlobalFunc::ZEROS(dphi_pool, 3*GlobalC::pw.bxyz*LD_pool);
            dphix = new double*[GlobalC::pw.bxyz];
            dphiy = new double*[GlobalC::pw.bxyz];
            dphiz = new double*[GlobalC::pw.bxyz];

            cal_flag=new bool*[GlobalC::pw.bxyz];
            for(int i=0; i<GlobalC::pw.bxyz; i++)
            {
                dphix[i] = &dphi_pool[i*LD_pool];
                dphiy[i] = &dphi_pool[i*LD_pool+GlobalC::pw.bxyz*LD_pool];
                dphiz[i] = &dphi_pool[i*LD_pool+2*GlobalC::pw.bxyz*LD_pool];
                cal_flag[i] = new bool[max_size];
            }

            ModuleBase::realArray drr;//rewrite drr form by zhengdy-2019-04-02
            if(GlobalV::CAL_STRESS)
            {
                drr.create(max_size, GlobalC::pw.bxyz, 3);
                drr.zero_out();
            }
#ifdef _OPENMP
            //#pragma omp for schedule(dynamic)
            #pragma omp for
#endif
            for(int i=0; i< GlobalC::GridT.nbx; i++)
            {
                const int ibx = i*GlobalC::pw.bx;
                for(int j=0; j< GlobalC::GridT.nby; j++)
                {
                    const int jby = j*GlobalC::pw.by;
                    for(int k= GlobalC::GridT.nbzp_start; k< GlobalC::GridT.nbzp_start+GlobalC::GridT.nbzp; k++)
                    {
                        const int kbz = k*GlobalC::pw.bz-GlobalC::pw.nczp_start;
                        const int grid_index = (k-GlobalC::GridT.nbzp_start) + j * GlobalC::GridT.nbzp + i * GlobalC::GridT.nby * GlobalC::GridT.nbzp;
                        const int na_grid = GlobalC::GridT.how_many_atoms[ grid_index ];
                        if(na_grid==0) continue;

                        //------------------------------------------------------------------
                        // extract the local potentials.
                        //------------------------------------------------------------------
                        double *vldr3 = get_vldr3(vlocal, ncyz, ibx, jby, kbz);

                        int * block_iw, * block_index, * block_size;
                        Gint_Tools::get_block_info(na_grid, grid_index, block_iw, block_index, block_size);

                        Gint_Tools::Array_Pool<double> psir_vlbr3(GlobalC::pw.bxyz, LD_pool);
                        Gint_Tools::Array_Pool<double> psir_ylm(GlobalC::pw.bxyz, LD_pool);

                        cal_psir_ylm_dphi(na_grid, grid_index, delta_r,
                            block_index, block_size, cal_flag, psir_ylm.ptr_2D, dphix, dphiy, dphiz, drr);

#ifdef _OPENMP
                        cal_meshball_DGridV(na_grid, GlobalC::GridT.lgd, LD_pool, block_index, block_iw, block_size, cal_flag, vldr3,
                            psir_ylm.ptr_2D, psir_vlbr3.ptr_2D, dphix,  dphiy, dphiz,
                            DGridV_x_thread, DGridV_y_thread, DGridV_z_thread,
                            DGridV_11_thread, DGridV_12_thread, DGridV_13_thread,
                            DGridV_22_thread, DGridV_23_thread, DGridV_33_thread, drr);
#else
                        cal_meshball_DGridV(na_grid, GlobalC::GridT.lgd, LD_pool, block_index, block_iw, block_size, cal_flag, vldr3,
                            psir_ylm.ptr_2D, psir_vlbr3.ptr_2D, dphix,  dphiy, dphiz,
                            DGridV_x, DGridV_y, DGridV_z,
                            DGridV_11, DGridV_12, DGridV_13,
                            DGridV_22, DGridV_23, DGridV_33, drr);
#endif

                        free(vldr3);        vldr3=nullptr;
                        delete[] block_iw;
                        delete[] block_index;
                        delete[] block_size;
                    }// k
                }// j
            }// i

            delete[] dphix;
            delete[] dphiy;
            delete[] dphiz;
            delete[] dphi_pool;
            for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
            {
                delete[] cal_flag[ib];
            }
            delete[] cal_flag;
        }
#ifdef _OPENMP
        #pragma omp critical(cal_fvl)
        {
            for(int i=0; i<GlobalC::GridT.lgd; i++)
            {
                for(int j=0; j<GlobalC::GridT.lgd; j++)
                {
                    DGridV_x[i][j] += DGridV_x_thread[i][j];
                    DGridV_y[i][j] += DGridV_y_thread[i][j];
                    DGridV_z[i][j] += DGridV_z_thread[i][j];

                    if(GlobalV::CAL_STRESS)
                    {
                        DGridV_11[i][j] += DGridV_11_thread[i][j];
                        DGridV_12[i][j] += DGridV_12_thread[i][j];
                        DGridV_13[i][j] += DGridV_13_thread[i][j];
                        DGridV_22[i][j] += DGridV_22_thread[i][j];
                        DGridV_23[i][j] += DGridV_23_thread[i][j];
                        DGridV_33[i][j] += DGridV_33_thread[i][j];
                    }
                }
            }

            delete [] DGridV_x_thread;
            delete [] DGridV_y_thread;
            delete [] DGridV_z_thread;
            if(GlobalV::CAL_STRESS)
            {
                delete [] DGridV_11_thread;
                delete [] DGridV_12_thread;
                delete [] DGridV_13_thread;
                delete [] DGridV_22_thread;
                delete [] DGridV_23_thread;
                delete [] DGridV_33_thread;
                delete [] DGridV_stress_pool_thread;
            }
            delete [] DGridV_pool_thread;
        }
#endif
#ifdef _OPENMP
    }
#endif

    ModuleBase::timer::tick("Gint_Gamma","gamma_force");
#ifdef __MPI
    ModuleBase::timer::tick("Gint_Gamma","gamma_force_wait");
    MPI_Barrier(MPI_COMM_WORLD);
    ModuleBase::timer::tick("Gint_Gamma","gamma_force_wait");
#endif
    ModuleBase::timer::tick("Gint_Gamma","gamma_force2");

    //OUT(GlobalV::ofs_running,"Start reduce DGridV");

    double* tmpx = new double[GlobalV::NLOCAL];
    double* tmpy = new double[GlobalV::NLOCAL];
    double* tmpz = new double[GlobalV::NLOCAL];
    double* tmp11;
    double* tmp12;
    double* tmp13;
    double* tmp22;
    double* tmp23;
    double* tmp33;

	if(GlobalV::CAL_STRESS)
	{
		tmp11 = new double[GlobalV::NLOCAL];
		tmp12 = new double[GlobalV::NLOCAL];
		tmp13 = new double[GlobalV::NLOCAL];
		tmp22 = new double[GlobalV::NLOCAL];
		tmp23 = new double[GlobalV::NLOCAL];
		tmp33 = new double[GlobalV::NLOCAL];
	}

    for (int i=0; i<GlobalV::NLOCAL; i++)
    {
        ModuleBase::GlobalFunc::ZEROS(tmpx, GlobalV::NLOCAL);
        ModuleBase::GlobalFunc::ZEROS(tmpy, GlobalV::NLOCAL);
		ModuleBase::GlobalFunc::ZEROS(tmpz, GlobalV::NLOCAL);

		if(GlobalV::CAL_STRESS)
		{
			ModuleBase::GlobalFunc::ZEROS(tmp11, GlobalV::NLOCAL);
			ModuleBase::GlobalFunc::ZEROS(tmp12, GlobalV::NLOCAL);
			ModuleBase::GlobalFunc::ZEROS(tmp13, GlobalV::NLOCAL);
			ModuleBase::GlobalFunc::ZEROS(tmp22, GlobalV::NLOCAL);
			ModuleBase::GlobalFunc::ZEROS(tmp23, GlobalV::NLOCAL);
			ModuleBase::GlobalFunc::ZEROS(tmp33, GlobalV::NLOCAL);
		}

        const int mu = GlobalC::GridT.trace_lo[i];
        // mohan fix bug 2010-09-05
        // lack mu>=0 and nu>=0 in previous version.
        if(mu >=0)
        {
            for (int j=0; j<GlobalV::NLOCAL; j++)
            {
                const int nu = GlobalC::GridT.trace_lo[j];
                if(nu>=0)
                {
                    tmpx[j] = DGridV_x[mu][nu];
                    tmpy[j] = DGridV_y[mu][nu];
                    tmpz[j] = DGridV_z[mu][nu];
                    if(GlobalV::CAL_STRESS)
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
        Parallel_Reduce::reduce_double_pool( tmpx, GlobalV::NLOCAL );
        Parallel_Reduce::reduce_double_pool( tmpy, GlobalV::NLOCAL );
        Parallel_Reduce::reduce_double_pool( tmpz, GlobalV::NLOCAL );
		if(GlobalV::CAL_STRESS)
		{
			Parallel_Reduce::reduce_double_pool( tmp11, GlobalV::NLOCAL );
			Parallel_Reduce::reduce_double_pool( tmp12, GlobalV::NLOCAL );
			Parallel_Reduce::reduce_double_pool( tmp13, GlobalV::NLOCAL );
			Parallel_Reduce::reduce_double_pool( tmp22, GlobalV::NLOCAL );
			Parallel_Reduce::reduce_double_pool( tmp23, GlobalV::NLOCAL );
			Parallel_Reduce::reduce_double_pool( tmp33, GlobalV::NLOCAL );
        }

        const Parallel_Orbitals* pv = this->LM->ParaV;

        for (int j=0; j<GlobalV::NLOCAL; j++)
        {
            if (!pv->in_this_processor(i,j))
            {
                continue;
            }
            this->LM->set_force (i,j,tmpx[j], tmpy[j], tmpz[j],'N');
            if(GlobalV::CAL_STRESS)
            {
                const int irr = pv->trace_loc_row[ i ];
                const int icc = pv->trace_loc_col[ j ];
                const int index = irr * pv->ncol + icc;
                this->LM->DHloc_fixed_11[index] += tmp11[j];
                this->LM->DHloc_fixed_12[index] += tmp12[j];
                this->LM->DHloc_fixed_13[index] += tmp13[j];
                this->LM->DHloc_fixed_22[index] += tmp22[j];
                this->LM->DHloc_fixed_23[index] += tmp23[j];
                this->LM->DHloc_fixed_33[index] += tmp33[j];
            }
        }
    }
    delete[] tmpx;
    delete[] tmpy;
    delete[] tmpz;
    if(GlobalV::CAL_STRESS)
    {
        delete[] tmp11;
        delete[] tmp12;
        delete[] tmp13;
        delete[] tmp22;
        delete[] tmp23;
        delete[] tmp33;
    }

    //delete DGridV_x,y,z
    delete [] DGridV_x;
    delete [] DGridV_y;
    delete [] DGridV_z;
    if(GlobalV::CAL_STRESS)
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

    ModuleBase::timer::tick("Gint_Gamma","gamma_force2");   
    return;
}
