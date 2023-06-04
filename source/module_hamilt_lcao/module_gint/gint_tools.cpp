//=========================================================
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#include "gint_tools.h"

#include <cmath>

#include "module_base/timer.h"
#include "module_base/ylm.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace Gint_Tools
{
    int* get_vindex(
        const int bxyz,
        const int bx,
        const int by,
        const int bz,
        const int nplane,
        const int start_ind,
		const int ncyz)
	{
		int *vindex = new int[bxyz];
		int bindex = 0;

		for(int ii=0; ii<bx; ii++)
		{
			const int ipart = ii*ncyz;
			for(int jj=0; jj<by; jj++)
			{
				const int jpart = jj*nplane + ipart;
				for(int kk=0; kk<bz; kk++)
				{
					vindex[bindex] = start_ind + kk + jpart;
					++bindex;
				}
			}
		}
		return vindex;
	}

	// here vindex refers to local potentials
    int* get_vindex(
        const int bxyz,
        const int bx,
        const int by,
        const int bz,
        const int nplane,
        const int ncyz,
		const int ibx,
		const int jby,
		const int kbz)
	{
		int *vindex = new int[bxyz];
		int bindex=0;
		// z is the fastest,
		// ipart can be obtained by using a previously stored array
		for(int ii=0; ii<bx; ii++)
		{
			const int ipart = (ibx+ii)*ncyz;
			for(int jj=0; jj<by; jj++)
			{
				// jpart can be obtained by using a previously stored array
				const int jpart = (jby+jj)*nplane + ipart;
				for(int kk=0; kk<bz; kk++)
				{
					vindex[bindex] = kbz+kk + jpart;
					++bindex;
				}
			}
		}
		return vindex;
	}

	// extract the local potentials.
	double* get_vldr3(
        const double* const vlocal,		// vlocal[ir]
        const int bxyz,
        const int bx,
        const int by,
        const int bz,
        const int nplane,
        const int ncyz,
		const int ibx,
		const int jby,
		const int kbz,
		const double dv)
	{
		// set the index for obtaining local potentials
		int* vindex = Gint_Tools::get_vindex(bxyz, bx, by, bz, nplane, ncyz, ibx, jby, kbz);
		double *vldr3 = new double[bxyz];
		for(int ib=0; ib<bxyz; ib++)
		{
			vldr3[ib]=vlocal[vindex[ib]] * dv;
		}
		delete[] vindex;
		return vldr3;
	}

	double* get_vldr3(
        const double* const vlocal,		// vlocal[ir]
        const int bxyz,
        const int bx,
        const int by,
        const int bz,
        const int nplane,
        const int start_ind,
		const int ncyz,
		const double dv)
	{
		// set the index for obtaining local potentials
		int* vindex = Gint_Tools::get_vindex(bxyz, bx, by, bz, nplane, start_ind, ncyz);
		double *vldr3 = new double[bxyz];
		for(int ib=0; ib<bxyz; ib++)
		{
			vldr3[ib]=vlocal[vindex[ib]] * dv;
		}
		delete[] vindex;
		return vldr3;
	}

    void get_block_info(
        const Grid_Technique& gt, 
        const int bxyz,
        const int na_grid,
		const int grid_index,
		int * &block_iw,
		int * &block_index,
		int * &block_size,
		bool** &cal_flag
	)
	{
		block_iw = new int[na_grid];
		block_index = new int[na_grid+1];
		block_size = new int[na_grid];
		cal_flag = new bool* [bxyz];
		for(int ib=0; ib<bxyz; ib++)
		{
			cal_flag[ib] = new bool[na_grid];
		}

		block_index[0] = 0;
		for (int id=0; id<na_grid; id++)
		{
			const int mcell_index=gt.bcell_start[grid_index] + id;
			const int iat=gt.which_atom[mcell_index]; // index of atom
			const int it=GlobalC::ucell.iat2it[ iat ]; // index of atom type
			const int ia=GlobalC::ucell.iat2ia[ iat ]; // index of atoms within each type
			const int start=GlobalC::ucell.itiaiw2iwt(it, ia, 0); // the index of the first wave function for atom (it,ia)
			block_iw[id]=gt.trace_lo[start];
			block_index[id+1] = block_index[id]+GlobalC::ucell.atoms[it].nw;
			block_size[id]=GlobalC::ucell.atoms[it].nw;

			const int imcell=gt.which_bigcell[mcell_index];
			const double mt[3] = {
				gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
				gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
				gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};

			for(int ib=0; ib<bxyz; ib++)
			{
				// meshcell_pos: z is the fastest
				const double dr[3] = {
					gt.meshcell_pos[ib][0] + mt[0],
					gt.meshcell_pos[ib][1] + mt[1],
					gt.meshcell_pos[ib][2] + mt[2]};
				const double distance = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);	// distance between atom and grid

				if(distance > GlobalC::ORB.Phi[it].getRcut() - 1.0e-10)
					cal_flag[ib][id]=false;
				else
					cal_flag[ib][id]=true;
			}// end ib
		}
	}

    void cal_psir_ylm(
        const Grid_Technique& gt, 
        const int bxyz,
        const int na_grid, 					// number of atoms on this grid
		const int grid_index, 				// 1d index of FFT index (i,j,k)
		const double delta_r, 				// delta_r of the uniform FFT grid
		const int*const block_index,  		// block_index[na_grid+1], count total number of atomis orbitals
		const int*const block_size, 		// block_size[na_grid],	number of columns of a band
		const bool*const*const cal_flag,
		double*const*const psir_ylm) 	// cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
    {
		ModuleBase::timer::tick("Gint_Tools", "cal_psir_ylm");
        std::vector<double> ylma;
        for (int id=0; id<na_grid; id++)
		{
			// there are two parameters we want to know here:
			// in which bigcell of the meshball the atom is in?
			// what's the cartesian coordinate of the bigcell?
			const int mcell_index=gt.bcell_start[grid_index] + id;

			const int iat=gt.which_atom[mcell_index]; // index of atom
			const int it=GlobalC::ucell.iat2it[iat]; // index of atom type
			const Atom*const atom=&GlobalC::ucell.atoms[it];
			auto &OrbPhi = GlobalC::ORB.Phi[it];
			std::vector<const double*> it_psi_uniform(atom->nw);
			std::vector<const double*> it_dpsi_uniform(atom->nw);
			// preprocess index
			for (int iw=0; iw< atom->nw; ++iw)
			{
				if ( atom->iw2_new[iw] )
				{
					auto philn = &OrbPhi.PhiLN(atom->iw2l[iw], atom->iw2n[iw]);
					it_psi_uniform[iw] = &philn->psi_uniform[0];
					it_dpsi_uniform[iw] = &philn->dpsi_uniform[0];
				}
			}

			// meshball_positions should be the bigcell position in meshball
			// to the center of meshball.
			// calculated in cartesian coordinates
			// the std::vector from the grid which is now being operated to the atom position.
			// in meshball language, is the std::vector from imcell to the center cel, plus
			// tau_in_bigcell.
			const int imcell=gt.which_bigcell[mcell_index];
			const double mt[3] = {
				gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
				gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
				gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};

			// number of grids in each big cell (bxyz)
			for(int ib=0; ib<bxyz; ib++)
			{
				double *p=&psir_ylm[ib][block_index[id]];
				if(!cal_flag[ib][id])
				{
					ModuleBase::GlobalFunc::ZEROS(p, block_size[id]);
				}
				else
				{
					// meshcell_pos: z is the fastest
					const double dr[3] = {
						gt.meshcell_pos[ib][0] + mt[0],
						gt.meshcell_pos[ib][1] + mt[1],
						gt.meshcell_pos[ib][2] + mt[2]};
					double distance = std::sqrt( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );	// distance between atom and grid
					//if(distance[id] > gt.orbital_rmax) continue;
					if (distance < 1.0E-9) distance += 1.0E-9;

					//------------------------------------------------------
					// spherical harmonic functions Ylm
					//------------------------------------------------------
					//	Ylm::get_ylm_real(this->nnn[it], this->dr[id], ylma);
					ModuleBase::Ylm::sph_harm ( GlobalC::ucell.atoms[it].nwl,
							dr[0] / distance,
							dr[1] / distance,
							dr[2] / distance,
							ylma);
					// these parameters are related to interpolation
					// because once the distance from atom to grid point is known,
					// we can obtain the parameters for interpolation and
					// store them first! these operations can save lots of efforts.
					const double position = distance / delta_r;
					const int ip = static_cast<int>(position);
					const double dx = position - ip;
					const double dx2 = dx * dx;
					const double dx3 = dx2 * dx;

					const double c3 = 3.0*dx2-2.0*dx3;
					const double c1 = 1.0-c3;
					const double c2 = (dx-2.0*dx2+dx3)*delta_r;
					const double c4 = (dx3-dx2)*delta_r;

					double phi=0;
					for (int iw=0; iw< atom->nw; ++iw)
					{
						if ( atom->iw2_new[iw] )
						{
							auto psi_uniform = it_psi_uniform[iw];
							auto dpsi_uniform = it_dpsi_uniform[iw];
							phi = c1*psi_uniform[ip] + c2*dpsi_uniform[ip]			 // radial wave functions
								+ c3*psi_uniform[ip+1] + c4*dpsi_uniform[ip+1];
						}
						p[iw]=phi * ylma[atom->iw2_ylm[iw]];
					} // end iw
				}// end distance<=(GlobalC::ORB.Phi[it].getRcut()-1.0e-15)
			}// end ib
		}// end id
		ModuleBase::timer::tick("Gint_Tools", "cal_psir_ylm");
		return;
	}

    void cal_dpsir_ylm(
        const Grid_Technique& gt, 
        const int bxyz,
        const int na_grid, 					// number of atoms on this grid
		const int grid_index, 				// 1d index of FFT index (i,j,k)
		const double delta_r, 				// delta_r of the uniform FFT grid
		const int*const block_index,  		// block_index[na_grid+1], count total number of atomis orbitals
		const int*const block_size, 		// block_size[na_grid],	number of columns of a band
		const bool*const*const cal_flag,    // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		double*const*const psir_ylm,
		double*const*const dpsir_ylm_x,
		double*const*const dpsir_ylm_y,
		double*const*const dpsir_ylm_z)
	{
		ModuleBase::timer::tick("Gint_Tools", "cal_dpsir_ylm");
		for (int id=0; id<na_grid; id++)
		{
			const int mcell_index = gt.bcell_start[grid_index] + id;
			const int imcell = gt.which_bigcell[mcell_index];
			int iat = gt.which_atom[mcell_index];
			const int it = GlobalC::ucell.iat2it[iat];
			const int ia = GlobalC::ucell.iat2ia[iat];
			Atom *atom = &GlobalC::ucell.atoms[it];

			const double mt[3]={
				gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
				gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
				gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};

			for(int ib=0; ib<bxyz; ib++)
			{
				double*const p_psi=&psir_ylm[ib][block_index[id]];
				double*const p_dpsi_x=&dpsir_ylm_x[ib][block_index[id]];
				double*const p_dpsi_y=&dpsir_ylm_y[ib][block_index[id]];
				double*const p_dpsi_z=&dpsir_ylm_z[ib][block_index[id]];
				if(!cal_flag[ib][id])
				{
					ModuleBase::GlobalFunc::ZEROS(p_psi, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_x, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_y, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_z, block_size[id]);
				}
				else
				{
					const double dr[3]={						// vectors between atom and grid
						gt.meshcell_pos[ib][0] + mt[0],
						gt.meshcell_pos[ib][1] + mt[1],
						gt.meshcell_pos[ib][2] + mt[2]};
					double distance = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

					//array to store spherical harmonics and its derivatives
					std::vector<double> rly;
					std::vector<std::vector<double>> grly;
					ModuleBase::Ylm::grad_rl_sph_harm(GlobalC::ucell.atoms[it].nwl, dr[0], dr[1], dr[2], rly, grly);
					if(distance < 1e-9)  distance = 1e-9;

					const double position = distance / delta_r;

					const double iq = static_cast<int>(position);
					const double x0 = position - iq;
					const double x1 = 1.0 - x0;
					const double x2 = 2.0 - x0;
					const double x3 = 3.0 - x0;
					const double x12 = x1*x2 / 6;
					const double x03 = x0*x3 / 2;

					double tmp, dtmp;

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

						// get the 'l' of this localized wave function
						const int ll = atom->iw2l[iw];
						const int idx_lm = atom->iw2_ylm[iw];

						const double rl = pow(distance, ll);

						// 3D wave functions
						p_psi[iw] = tmp * rly[idx_lm] / rl;

						// derivative of wave functions with respect to atom positions.
						const double tmpdphi_rly = (dtmp  - tmp * ll / distance) / rl * rly[idx_lm] / distance;
						const double tmprl = tmp/rl;

						p_dpsi_x[iw] = tmpdphi_rly * dr[0]  + tmprl * grly[idx_lm][0];
						p_dpsi_y[iw] = tmpdphi_rly * dr[1]  + tmprl * grly[idx_lm][1];
						p_dpsi_z[iw] = tmpdphi_rly * dr[2]  + tmprl * grly[idx_lm][2];
					}//iw
				}//else
			}
		}
		ModuleBase::timer::tick("Gint_Tools", "cal_dpsir_ylm");
		return;
	}

    void cal_ddpsir_ylm(
        const Grid_Technique& gt, 
        const int bxyz,
        const int na_grid, 					// number of atoms on this grid
		const int grid_index, 				// 1d index of FFT index (i,j,k)
		const double delta_r, 				// delta_r of the uniform FFT grid
		const int*const block_index,  		// block_index[na_grid+1], count total number of atomis orbitals
		const int*const block_size, 		// block_size[na_grid],	number of columns of a band
		const bool*const*const cal_flag,    // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		double*const*const ddpsir_ylm_xx,
		double*const*const ddpsir_ylm_xy,
		double*const*const ddpsir_ylm_xz,
		double*const*const ddpsir_ylm_yy,
		double*const*const ddpsir_ylm_yz,
		double*const*const ddpsir_ylm_zz)
	{
		ModuleBase::timer::tick("Gint_Tools", "cal_ddpsir_ylm");
		for (int id=0; id<na_grid; id++)
		{
			const int mcell_index = gt.bcell_start[grid_index] + id;
			const int imcell = gt.which_bigcell[mcell_index];
			int iat = gt.which_atom[mcell_index];
			const int it = GlobalC::ucell.iat2it[iat];
			const int ia = GlobalC::ucell.iat2ia[iat];
			Atom *atom = &GlobalC::ucell.atoms[it];

			const double mt[3]={
				gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
				gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
				gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};

			for(int ib=0; ib<bxyz; ib++)
			{
				double*const p_ddpsi_xx=&ddpsir_ylm_xx[ib][block_index[id]];
				double*const p_ddpsi_xy=&ddpsir_ylm_xy[ib][block_index[id]];
				double*const p_ddpsi_xz=&ddpsir_ylm_xz[ib][block_index[id]];
				double*const p_ddpsi_yy=&ddpsir_ylm_yy[ib][block_index[id]];
				double*const p_ddpsi_yz=&ddpsir_ylm_yz[ib][block_index[id]];
				double*const p_ddpsi_zz=&ddpsir_ylm_zz[ib][block_index[id]];
				if(!cal_flag[ib][id])
				{
					ModuleBase::GlobalFunc::ZEROS(p_ddpsi_xx, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_ddpsi_xy, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_ddpsi_xz, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_ddpsi_yy, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_ddpsi_yz, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_ddpsi_zz, block_size[id]);
				}
				else
				{
					const double dr[3]={						// vectors between atom and grid
						gt.meshcell_pos[ib][0] + mt[0],
						gt.meshcell_pos[ib][1] + mt[1],
						gt.meshcell_pos[ib][2] + mt[2]};
					double distance = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);

					// for some unknown reason, the finite difference between dpsi and ddpsi
					// using analytical expression is always wrong; as a result,
					// I switch to explicit finite difference method for evaluating
					// the second derivatives of the orbitals
					if(/*distance < 1e-9*/true)
					{
						double*** dpsi = new double** [atom->nw];
						for(int i=0;i<atom->nw;i++)
						{
							dpsi[i] = new double*[6];
							for(int j=0;j<6;j++)
							{
								dpsi[i][j] = new double[3];
								ModuleBase::GlobalFunc::ZEROS(dpsi[i][j], 3);
							}
						}

						double* dr1 = new double[3];

						double** displ = new double* [6];
						for(int i=0;i<6;i++)
						{
							displ[i] = new double[3];
							ModuleBase::GlobalFunc::ZEROS(displ[i], 3);
						}
						displ[0][0] = 0.0001; // in x direction
						displ[1][0] = -0.0001;
						displ[2][1] = 0.0001; // in y direction
						displ[3][1] = -0.0001;
						displ[4][2] = 0.0001; // in z direction
						displ[5][2] = -0.0001;

						for(int i=0;i<6;i++)
						{
							dr1[0] = dr[0] + displ[i][0];
							dr1[1] = dr[1] + displ[i][1];
							dr1[2] = dr[2] + displ[i][2];

							//array to store spherical harmonics and its derivatives
							std::vector<double> rly;
							std::vector<std::vector<double>> grly;
							ModuleBase::Ylm::grad_rl_sph_harm(GlobalC::ucell.atoms[it].nwl, dr1[0], dr1[1], dr1[2], rly, grly);

							double distance1 = std::sqrt(dr1[0]*dr1[0] + dr1[1]*dr1[1] + dr1[2]*dr1[2]);
							if(distance1 < 1e-9)  distance1 = 1e-9;

							const double position = distance1 / delta_r;

							const double iq = static_cast<int>(position);
							const double x0 = position - iq;
							const double x1 = 1.0 - x0;
							const double x2 = 2.0 - x0;
							const double x3 = 3.0 - x0;
							const double x12 = x1*x2 / 6;
							const double x03 = x0*x3 / 2;

							double tmp, dtmp;

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

								// get the 'l' of this localized wave function
								const int ll = atom->iw2l[iw];
								const int idx_lm = atom->iw2_ylm[iw];

								const double rl = pow(distance1, ll);

								// derivative of wave functions with respect to atom positions.
								const double tmpdphi_rly = (dtmp  - tmp * ll / distance1) / rl * rly[idx_lm] / distance1;
								const double tmprl = tmp/rl;

								dpsi[iw][i][0] = tmpdphi_rly * dr1[0]  + tmprl * grly[idx_lm][0];
								dpsi[iw][i][1] = tmpdphi_rly * dr1[1]  + tmprl * grly[idx_lm][1];
								dpsi[iw][i][2] = tmpdphi_rly * dr1[2]  + tmprl * grly[idx_lm][2];
							} // end iw
						}//end i = 0-6

						for(int iw=0;iw<atom->nw;iw++)
						{
							p_ddpsi_xx[iw] = (dpsi[iw][0][0] - dpsi[iw][1][0]) / 0.0002;
							p_ddpsi_xy[iw] = ((dpsi[iw][2][0] - dpsi[iw][3][0])+(dpsi[iw][0][1] - dpsi[iw][1][1])) / 0.0004;
							p_ddpsi_xz[iw] = ((dpsi[iw][4][0] - dpsi[iw][5][0])+(dpsi[iw][0][2] - dpsi[iw][1][2])) / 0.0004;
							p_ddpsi_yy[iw] = (dpsi[iw][2][1] - dpsi[iw][3][1]) / 0.0002;
							p_ddpsi_yz[iw] = ((dpsi[iw][4][1] - dpsi[iw][5][1])+(dpsi[iw][2][2] - dpsi[iw][3][2])) / 0.0004;
							p_ddpsi_zz[iw] = (dpsi[iw][4][2] - dpsi[iw][5][2]) / 0.0002;
						}

						for(int i=0;i<atom->nw;i++)
						{
							for(int j=0;j<6;j++)
							{
								delete[] dpsi[i][j];
							}
							delete dpsi[i];
						}
						delete[] dpsi;

						delete[] dr1;
						for(int i=0;i<6;i++)
						{
							delete[] displ[i];
						}
						delete[] displ;
					}
					else
					// the analytical method for evaluating 2nd derivatives
					// it is not used currently
					{
						//array to store spherical harmonics and its derivatives
						std::vector<double> rly;
						std::vector<std::vector<double>> grly;
						std::vector<std::vector<double>> hrly;
						ModuleBase::Ylm::grad_rl_sph_harm(GlobalC::ucell.atoms[it].nwl, dr[0], dr[1], dr[2], rly, grly);
						ModuleBase::Ylm::hes_rl_sph_harm(GlobalC::ucell.atoms[it].nwl, dr[0], dr[1], dr[2], hrly);
						const double position = distance / delta_r;

						const double iq = static_cast<int>(position);
						const double x0 = position - iq;
						const double x1 = 1.0 - x0;
						const double x2 = 2.0 - x0;
						const double x3 = 3.0 - x0;
						const double x12 = x1*x2 / 6;
						const double x03 = x0*x3 / 2;

						double tmp, dtmp, ddtmp;

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
									tmp = dtmp = ddtmp = 0.0;
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

								    ddtmp = x12*(philn.ddpsi_uniform[iq]*x3
											+philn.ddpsi_uniform[iq+3]*x0)
											+ x03*(philn.ddpsi_uniform[iq+1]*x2
												-philn.ddpsi_uniform[iq+2]*x1);
								}
							}//new l is used.

							// get the 'l' of this localized wave function
							const int ll = atom->iw2l[iw];
							const int idx_lm = atom->iw2_ylm[iw];

							const double rl = pow(distance, ll);
							const double r_lp2 = pow(distance, ll+2);

							// d/dr (R_l / r^l)
							const double tmpdphi = (dtmp - tmp * ll / distance) / rl;
							const double term1 = ddtmp / r_lp2;
							const double term2 = (2 * ll + 1) * dtmp / r_lp2 / distance;
							const double term3 = ll * (ll + 2) * tmp / r_lp2 / distance / distance;
							const double term4 = tmpdphi / distance;
							const double term5 = term1 - term2 + term3;

							// hessian of (R_l / r^l)
							const double term_xx = term4 + dr[0]*dr[0]*term5;
							const double term_xy = dr[0]*dr[1]*term5;
							const double term_xz = dr[0]*dr[2]*term5;
							const double term_yy = term4 + dr[1]*dr[1]*term5;
							const double term_yz = dr[1]*dr[2]*term5;
							const double term_zz = term4 + dr[2]*dr[2]*term5;

							// d/dr (R_l / r^l) * alpha / r
							const double term_1x = dr[0] * term4;
							const double term_1y = dr[1] * term4;
							const double term_1z = dr[2] * term4;

							p_ddpsi_xx[iw] = term_xx * rly[idx_lm] + 2.0*term_1x*grly[idx_lm][0] + tmp/rl*hrly[idx_lm][0];
							p_ddpsi_xy[iw] = term_xy * rly[idx_lm] + term_1x*grly[idx_lm][1] + term_1y*grly[idx_lm][0] + tmp/rl*hrly[idx_lm][1];
							p_ddpsi_xz[iw] = term_xz * rly[idx_lm] + term_1x*grly[idx_lm][2] + term_1z*grly[idx_lm][0] + tmp/rl*hrly[idx_lm][2];
							p_ddpsi_yy[iw] = term_yy * rly[idx_lm] + 2.0*term_1y*grly[idx_lm][1] + tmp/rl*hrly[idx_lm][3];
							p_ddpsi_yz[iw] = term_yz * rly[idx_lm] + term_1y*grly[idx_lm][2] + term_1z*grly[idx_lm][1] + tmp/rl*hrly[idx_lm][4];
							p_ddpsi_zz[iw] = term_zz * rly[idx_lm] + 2.0*term_1z*grly[idx_lm][2] + tmp/rl*hrly[idx_lm][5];

						}//iw
					} // end if
				}//else
			}//end ib
		}//end id(atom)
		ModuleBase::timer::tick("Gint_Tools", "cal_ddpsir_ylm");
		return;
	}

    void cal_dpsirr_ylm(
        const Grid_Technique& gt,
        const int bxyz,
        const int na_grid, 					// number of atoms on this grid
		const int grid_index, 				// 1d index of FFT index (i,j,k)
		const int*const block_index,  		// block_index[na_grid+1], count total number of atomis orbitals
		const int*const block_size, 		// block_size[na_grid],	number of columns of a band
		const bool*const*const cal_flag,    // cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		double*const*const dpsir_ylm_x,
		double*const*const dpsir_ylm_y,
		double*const*const dpsir_ylm_z,
		double*const*const dpsir_ylm_xx,
		double*const*const dpsir_ylm_xy,
		double*const*const dpsir_ylm_xz,
		double*const*const dpsir_ylm_yy,
		double*const*const dpsir_ylm_yz,
		double*const*const dpsir_ylm_zz)
	{
		ModuleBase::timer::tick("Gint_Tools", "cal_dpsirr_ylm");
		for (int id=0; id<na_grid; id++)
		{
			const int mcell_index = gt.bcell_start[grid_index] + id;
			const int imcell = gt.which_bigcell[mcell_index];
			int iat = gt.which_atom[mcell_index];
			const int it = GlobalC::ucell.iat2it[iat];
			Atom *atom = &GlobalC::ucell.atoms[it];

			const double mt[3]={
				gt.meshball_positions[imcell][0] - gt.tau_in_bigcell[iat][0],
				gt.meshball_positions[imcell][1] - gt.tau_in_bigcell[iat][1],
				gt.meshball_positions[imcell][2] - gt.tau_in_bigcell[iat][2]};

			for(int ib=0; ib<bxyz; ib++)
			{
				double*const p_dpsi_x=&dpsir_ylm_x[ib][block_index[id]];
				double*const p_dpsi_y=&dpsir_ylm_y[ib][block_index[id]];
				double*const p_dpsi_z=&dpsir_ylm_z[ib][block_index[id]];
				double*const p_dpsi_xx=&dpsir_ylm_xx[ib][block_index[id]];
				double*const p_dpsi_xy=&dpsir_ylm_xy[ib][block_index[id]];
				double*const p_dpsi_xz=&dpsir_ylm_xz[ib][block_index[id]];
				double*const p_dpsi_yy=&dpsir_ylm_yy[ib][block_index[id]];
				double*const p_dpsi_yz=&dpsir_ylm_yz[ib][block_index[id]];
				double*const p_dpsi_zz=&dpsir_ylm_zz[ib][block_index[id]];
				if(!cal_flag[ib][id])
				{
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_xx, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_xy, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_xz, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_yy, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_yz, block_size[id]);
					ModuleBase::GlobalFunc::ZEROS(p_dpsi_zz, block_size[id]);
				}
				else
				{
					const double dr[3]={						// vectors between atom and grid
						gt.meshcell_pos[ib][0] + mt[0],
						gt.meshcell_pos[ib][1] + mt[1],
						gt.meshcell_pos[ib][2] + mt[2]};

					for (int iw=0; iw< atom->nw; ++iw)
					{

						p_dpsi_xx[iw] = p_dpsi_x[iw]*dr[0];
						p_dpsi_xy[iw] = p_dpsi_x[iw]*dr[1];
						p_dpsi_xz[iw] = p_dpsi_x[iw]*dr[2];
						p_dpsi_yy[iw] = p_dpsi_y[iw]*dr[1];
						p_dpsi_yz[iw] = p_dpsi_y[iw]*dr[2];
						p_dpsi_zz[iw] = p_dpsi_z[iw]*dr[2];

					}//iw
				}//else
			}
		}
		ModuleBase::timer::tick("Gint_Tools", "cal_dpsirr_ylm");
		return;
	}

	// atomic basis sets
	// psir_vlbr3[bxyz][LD_pool]
    Gint_Tools::Array_Pool<double> get_psir_vlbr3(
        const int bxyz,
        const int na_grid,  					    // how many atoms on this (i,j,k) grid
		const int LD_pool,
		const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
		const bool*const*const cal_flag,	    	// cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const vldr3,			    	// vldr3[bxyz]
		const double*const*const psir_ylm)		    // psir_ylm[bxyz][LD_pool]
	{
		Gint_Tools::Array_Pool<double> psir_vlbr3(bxyz, LD_pool);
		for(int ib=0; ib<bxyz; ++ib)
		{
			for(int ia=0; ia<na_grid; ++ia)
			{
				if(cal_flag[ib][ia])
				{
					for(int i=block_index[ia]; i<block_index[ia+1]; ++i)
					{
						psir_vlbr3.ptr_2D[ib][i]=psir_ylm[ib][i]*vldr3[ib];
					}
				}
				else
				{
					for(int i=block_index[ia]; i<block_index[ia+1]; ++i)
					{
						psir_vlbr3.ptr_2D[ib][i]=0;
					}
				}

			}
		}
		return psir_vlbr3;
	}

    void mult_psi_DM(
        const Grid_Technique& gt, 
        const int bxyz,
        const int na_grid,  					    // how many atoms on this (i,j,k) grid
		const int LD_pool,
		const int*const block_iw,				    // block_iw[na_grid],	index of wave functions for each block
		const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
		const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
		const bool*const*const cal_flag,	    	// cal_flag[bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const*const psi,	    // psir_vlbr3[bxyz][LD_pool]
		double ** psi_DM,
		const double*const*const DM,
		const int job) // 1: density, 2: force
	{
		constexpr char side='L', uplo='U';
		constexpr char transa='N', transb='N';
		constexpr double alpha_symm=1, beta=1;
		constexpr int inc=1;
		double alpha_gemm;

		switch(job)
		{
			case 1:
				alpha_gemm=2.0;
				break;
			case 2:
				alpha_gemm=1.0;
				break;
			default:
				ModuleBase::WARNING_QUIT("psir_dm","job can only be 1 or 2");
		}

		for (int ia1=0; ia1<na_grid; ia1++)
		{
			const int iw1_lo=block_iw[ia1];
			if(job==1)//density
			{
            	//ia1==ia2, diagonal part
				// find the first ib and last ib for non-zeros cal_flag
				int first_ib=0, last_ib=0;
				for(int ib=0; ib<bxyz; ++ib)
				{
					if(cal_flag[ib][ia1])
					{
						first_ib=ib;
						break;
					}
				}
				for(int ib=bxyz-1; ib>=0; --ib)
				{
					if(cal_flag[ib][ia1])
					{
						last_ib=ib+1;
						break;
					}
				}
				const int ib_length=last_ib-first_ib;
				if(ib_length<=0) continue;

				int cal_num=0;
				for(int ib=first_ib; ib<last_ib; ++ib)
				{
					cal_num += cal_flag[ib][ia1];
				}
				// if enough cal_flag is nonzero
				if(cal_num>ib_length/4)
				{
					dsymm_(&side, &uplo, &block_size[ia1], &ib_length,
						&alpha_symm, &DM[iw1_lo][iw1_lo], &gt.lgd,
						&psi[first_ib][block_index[ia1]], &LD_pool,
						&beta, &psi_DM[first_ib][block_index[ia1]], &LD_pool);
				}
				else
				{
					// int k=1;
					for(int ib=first_ib; ib<last_ib; ++ib)
					{
						if(cal_flag[ib][ia1])
						{
							dsymv_(&uplo, &block_size[ia1],
								&alpha_symm, &DM[iw1_lo][iw1_lo], &gt.lgd,
								&psi[ib][block_index[ia1]], &inc,
								&beta, &psi_DM[ib][block_index[ia1]], &inc);
						}
					}
				}
			}

			int start;
			switch(job)
			{
				case 1:
					start=ia1+1;
					break;
				case 2:
					start=0;
					break;
				default:
					ModuleBase::WARNING_QUIT("psi_dm","job can only be 1 or 2");
			}

			for (int ia2=start; ia2<na_grid; ia2++)
			{
				int first_ib=0, last_ib=0;
				for(int ib=0; ib<bxyz; ++ib)
				{
					if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
					{
						first_ib=ib;
						break;
					}
				}
				for(int ib=bxyz-1; ib>=0; --ib)
				{
					if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
					{
						last_ib=ib+1;
						break;
					}
				}
				const int ib_length=last_ib-first_ib;
				if(ib_length<=0) continue;

				int cal_pair_num=0;
				for(int ib=first_ib; ib<last_ib; ++ib)
				{
					cal_pair_num += cal_flag[ib][ia1] && cal_flag[ib][ia2];
				}
				const int iw2_lo=block_iw[ia2];
				if(cal_pair_num>ib_length/4)
				{
                    dgemm_(&transa, &transb, &block_size[ia2], &ib_length, &block_size[ia1],
                        &alpha_gemm, &DM[iw1_lo][iw2_lo], &gt.lgd,
                        &psi[first_ib][block_index[ia1]], &LD_pool,
                        &beta, &psi_DM[first_ib][block_index[ia2]], &LD_pool);
				}
                else
                {
                    for(int ib=first_ib; ib<last_ib; ++ib)
                    {
                        if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                        {
                            dgemv_(&transa, &block_size[ia2], &block_size[ia1],
                                &alpha_gemm, &DM[iw1_lo][iw2_lo], &gt.lgd,
                                &psi[ib][block_index[ia1]], &inc,
                                &beta, &psi_DM[ib][block_index[ia2]], &inc);
                        }
                    }
                }
			}// ia2
		} // ia1
	}

//calculating (psi_DMR)_mu = sum_nu DMR_mu,nu psi_nu
//note : there is a difference between rho and force
//in calculating rho, due to symmetry, the summation over mu,nu
//can be done as sum_mu,mu + 2 sum_mu<nu, saving some time
//but for force, we cannot exchange the index
    void mult_psi_DMR(
        const Grid_Technique& gt, 
        const int bxyz,
        const int& grid_index,
		const int &na_grid,
		const int*const block_index,
		const int*const block_size,
		bool** cal_flag,
		double** psi,
		double ** psi_DMR,
		double* DMR,
		const int job)
	{
		double *psi2, *psi2_dmr;
		int iwi, iww;
		const int LD_pool = gt.max_atom*GlobalC::ucell.nwmax;

		bool *all_out_of_range = new bool[na_grid];
		for(int ia=0; ia<na_grid; ++ia) //number of atoms
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

		//parameters for lapack subroutiens
		const char trans='N';
		const double alpha=1.0, beta=1.0;
		const int inc=1;
		double alpha1;
		switch(job)
		{
			case 1:
				alpha1=2.0;
				break;
			case 2:
				alpha1=1.0;
				break;
			default:
				ModuleBase::WARNING_QUIT("psir_dmr","job can only be 1 or 2");
		}

		for (int ia1=0; ia1<na_grid; ia1++)
		{
			if(all_out_of_range[ia1]) continue;

			const int mcell_index1 = gt.bcell_start[grid_index] + ia1;
			const int iat = gt.which_atom[mcell_index1];
			const int T1 = GlobalC::ucell.iat2it[iat];
			const int I1 = GlobalC::ucell.iat2ia[iat];
			Atom *atom1 = &GlobalC::ucell.atoms[T1];

			//~~~~~~~~~~~~~~~~
			// get cell R1.
			//~~~~~~~~~~~~~~~~
			const int id1 = gt.which_unitcell[mcell_index1];
			const int R1x = gt.ucell_index2x[id1];
			const int R1y = gt.ucell_index2y[id1];
			const int R1z = gt.ucell_index2z[id1];

			if(job==1) //density
			{
				const int idx1=block_index[ia1];
				int* find_start = gt.find_R2[iat];
				int* find_end = gt.find_R2[iat] + gt.nad[iat];
				//ia2==ia1
				int cal_num=0;
				for(int ib=0; ib<bxyz; ++ib)
				{
					if(cal_flag[ib][ia1])
					{
						++cal_num;
					}
				}

				int offset;
				if(cal_num>0)
				{
					//find offset
					const int index = gt.cal_RindexAtom(0, 0, 0, iat);
					offset = -1;
					for(int* find=find_start; find < find_end; find++)
					{
						//--------------------------------------------------------------
						// start positions of adjacent atom of 'iat'
						//--------------------------------------------------------------
						if( find[0] == index )
						{
							offset = find - find_start; // start positions of adjacent atom of 'iat'
							break;
						}
					}

					assert(offset!=-1);
					assert(offset < gt.nad[iat]);
				}

				if(cal_num>bxyz/4)
				{
					const int DM_start = gt.nlocstartg[iat]+ gt.find_R2st[iat][offset];
					dgemm_(&trans, &trans, &block_size[ia1], &bxyz, &block_size[ia1], &alpha,
						&DMR[DM_start], &block_size[ia1],
						&psi[0][idx1], &LD_pool,
						&beta, &psi_DMR[0][idx1], &LD_pool);
				}
				else if(cal_num>0)
				{
					const int DM_start = gt.nlocstartg[iat]+ gt.find_R2st[iat][offset];
					for(int ib=0; ib<bxyz; ++ib					)
					{
						if(cal_flag[ib][ia1])
						{
							dgemv_(&trans, &block_size[ia1], &block_size[ia1], &alpha,
									&DMR[DM_start], &block_size[ia1],
									&psi[ib][idx1], &inc,
									&beta, &psi_DMR[ib][idx1], &inc);
						}
					}
				}
			}

			// get (j,beta,R2)
			int start;
			switch(job)
			{
				case 1:
					start=ia1+1;
					break;
				case 2:
					start=0;
					break;
				default:
					ModuleBase::WARNING_QUIT("psi_dmr","job can only be 1 or 2");
			}

			for (int ia2=start; ia2<na_grid; ia2++)
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
						same_flag = true;
						break;
					}
				}

				if(!same_flag) continue;

				const int bcell2 = gt.bcell_start[grid_index] + ia2;
				const int T2 = GlobalC::ucell.iat2it[ gt.which_atom[bcell2]];
				const int iat2 = gt.which_atom[bcell2];

				Atom *atom2 = &GlobalC::ucell.atoms[T2];

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

				const int index = gt.cal_RindexAtom(dRx, dRy, dRz, iat2);
				int offset = -1;

				int* find_start = gt.find_R2[iat];
				int* findend = gt.find_R2[iat] + gt.nad[iat];

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
					ModuleBase::WARNING_QUIT("gint_k","mult_psi_DMR wrong");
				}
				assert(offset < gt.nad[iat]);

				//---------------------------------------------------------------
				// what I do above is to get 'offset' for atom std::pair (iat1, iat2)
				// if I want to simplify this searching for offset,
				// I should take advantage of gt.which_unitcell.
				//---------------------------------------------------------------

				int cal_num=0;
   				for(int ib=0; ib<bxyz; ++ib)
    			{
        			if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
        			    ++cal_num;
    			}

				if(cal_num>bxyz/4)
				{
					const int idx1=block_index[ia1];
			        const int idx2=block_index[ia2];
    				const int DM_start = gt.nlocstartg[iat]+ gt.find_R2st[iat][offset];
    				dgemm_(&trans, &trans, &block_size[ia2], &bxyz, &block_size[ia1], &alpha1,
    					&DMR[DM_start], &block_size[ia2],
    					&psi[0][idx1], &LD_pool,
    					&beta, &psi_DMR[0][idx2], &LD_pool);
				}
				else if(cal_num>0)
				{
					const int idx1=block_index[ia1];
					const int idx2=block_index[ia2];
    				const int DM_start = gt.nlocstartg[iat]+ gt.find_R2st[iat][offset];

    				for(int ib=0; ib<bxyz; ++ib)
    				{
        				if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
        				{
            				dgemv_(&trans, &block_size[ia2], &block_size[ia1], &alpha1,
            					&DMR[DM_start], &block_size[ia2],
            					&psi[ib][idx1], &inc,
            					&beta, &psi_DMR[ib][idx2], &inc);
        				}
    				}
				} // cal_num
			}// ia2
		}//ia1

		delete[] all_out_of_range;

	}
}