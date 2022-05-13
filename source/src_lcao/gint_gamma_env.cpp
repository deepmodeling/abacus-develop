#include "gint_gamma.h"
#include "grid_technique.h"
#include "../module_orbital/ORB_read.h"
#include "../src_pw/global.h"
#include "../module_base/ylm.h"
#include "../module_base/timer.h"

void Gint_Gamma::cal_env(const double* wfc, double* rho)
{
    ModuleBase::TITLE("Gint_Gamma","cal_env");
    ModuleBase::timer::tick("Gint_Gamma","cal_env");

	this->max_size = GlobalC::GridT.max_atom;
    this->gamma_envelope(wfc, rho);

    ModuleBase::timer::tick("Gint_Gamma","cal_env");
    return;
}


void Gint_Gamma::gamma_envelope(const double* wfc, double* rho)
{
    ModuleBase::TITLE("Grid_Integral","gamma_charge");

    // it's a uniform grid to save orbital values, so the delta_r is a constant.
    const double delta_r = GlobalC::ORB.dr_uniform;
	const Numerical_Orbital_Lm* pointer;

	double*** dr; // vectors between atom and grid: [bxyz, maxsize, 3]
	double** distance; // distance between atom and grid: [bxyz, maxsize]
	double*** psir_ylm;	
	bool** cal_flag;
	if(max_size!=0) 
	{
		dr = new double**[GlobalC::pw.bxyz];
		distance = new double*[GlobalC::pw.bxyz];
		psir_ylm = new double**[GlobalC::pw.bxyz];
		cal_flag = new bool*[GlobalC::pw.bxyz];

		for(int i=0; i<GlobalC::pw.bxyz; i++)
		{
			dr[i] = new double*[max_size];
			distance[i] = new double[max_size];
			psir_ylm[i] = new double*[max_size];
			cal_flag[i] = new bool[max_size];

			ModuleBase::GlobalFunc::ZEROS(distance[i], max_size);
			ModuleBase::GlobalFunc::ZEROS(cal_flag[i], max_size);

			for(int j=0; j<max_size; j++) 
			{
				dr[i][j] = new double[3];
				psir_ylm[i][j] = new double[GlobalC::ucell.nwmax];
				ModuleBase::GlobalFunc::ZEROS(dr[i][j],3);
				ModuleBase::GlobalFunc::ZEROS(psir_ylm[i][j],GlobalC::ucell.nwmax);
			}
		}
	}

	double mt[3]={0,0,0};
	double *vldr3 = new double[GlobalC::pw.bxyz];
	double v1 = 0.0;
	int* vindex=new int[GlobalC::pw.bxyz];
	ModuleBase::GlobalFunc::ZEROS(vldr3, GlobalC::pw.bxyz);
	ModuleBase::GlobalFunc::ZEROS(vindex, GlobalC::pw.bxyz);
	double phi=0.0;

	const int nbx = GlobalC::GridT.nbx;
	const int nby = GlobalC::GridT.nby;
	const int nbz_start = GlobalC::GridT.nbzp_start;
	const int nbz = GlobalC::GridT.nbzp;

    for (int i=0; i<nbx; i++)
    {
        for (int j=0; j<nby; j++)
        {
            for (int k=nbz_start; k<nbz_start+nbz; k++) // FFT grid
            {
                this->grid_index = (k-nbz_start) + j * nbz + i * nby * nbz;

                // get the value: how many atoms has orbital value on this grid.
                const int size = GlobalC::GridT.how_many_atoms[ this->grid_index ];
				if(size==0) continue;

				// (1) initialized the phi * Ylm.
                for (int id=0; id<size; id++)
                {
                    // there are two parameters we want to know here:
                    // in which bigcell of the meshball the atom in?
                    // what's the cartesian coordinate of the bigcell?
                    const int mcell_index = GlobalC::GridT.bcell_start[grid_index] + id;
                    const int imcell = GlobalC::GridT.which_bigcell[mcell_index];

                    int iat = GlobalC::GridT.which_atom[mcell_index];

                    const int it = GlobalC::ucell.iat2it[ iat ];
                    const int ia = GlobalC::ucell.iat2ia[ iat ];

                    // meshball_positions should be the bigcell position in meshball
                    // to the center of meshball.
                    // calculated in cartesian coordinates
                    // the std::vector from the grid which is now being operated to the atom position.
                    // in meshball language, is the std::vector from imcell to the center cel, plus
                    // tau_in_bigcell.
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
						+ dr[ib][id][1]*dr[ib][id][1] 
						+ dr[ib][id][2]*dr[ib][id][2]);

						if(distance[ib][id] <= GlobalC::ORB.Phi[it].getRcut()) 
						{
							cal_flag[ib][id]=true;
						}
						else 
						{
							cal_flag[ib][id]=false;
							continue;
						}
						
						std::vector<double> ylma;
						//if(distance[id] > GlobalC::GridT.orbital_rmax) continue;
						if (distance[ib][id] < 1.0E-9) distance[ib][id] += 1.0E-9;
						
						ModuleBase::Ylm::sph_harm (	GlobalC::ucell.atoms[it].nwl,
								dr[ib][id][0] / distance[ib][id],
								dr[ib][id][1] / distance[ib][id],
								dr[ib][id][2] / distance[ib][id],
								ylma);
						// these parameters are about interpolation
						// because once we know the distance from atom to grid point,
						// we can get the parameters we need to do interpolation and
						// store them first!! these can save a lot of effort.
						const double position = distance[ib][id] / delta_r;

						int ip;
						double dx, dx2, dx3;
						double c1, c2, c3, c4;

						ip = static_cast<int>(position);
						dx = position - ip;
						dx2 = dx * dx;
						dx3 = dx2 * dx;

						c3 = 3.0*dx2-2.0*dx3;
						c1 = 1.0-c3;
						c2 = (dx-2.0*dx2+dx3)*delta_r;
						c4 = (dx3-dx2)*delta_r;

						//		int ip = this->iq[id];
						//		double A = ip+1.0-position/delta_r;
						//		double B = 1.0-A;
						//		double coef1 = (A*A*A-A)/6.0*delta_r*delta_r;
						//		double coef2 = (B*B*B-B)/6.0*delta_r*delta_r;

						Atom* atom1 = &GlobalC::ucell.atoms[it];
						for (int iw=0; iw< atom1->nw; iw++)
						{
							if ( atom1->iw2_new[iw] )
							{
								pointer = &GlobalC::ORB.Phi[it].PhiLN(
										atom1->iw2l[iw],
										atom1->iw2n[iw]);
								phi = c1*pointer->psi_uniform[ip]+c2*pointer->dpsi_uniform[ip]
									+ c3*pointer->psi_uniform[ip+1] + c4*pointer->dpsi_uniform[ip+1];
							}
							psir_ylm[ib][id][iw] = phi * ylma[atom1->iw2_ylm[iw]];
							//psir_ylm[ib][id][iw] = 1;//for test
						}
					}// end ib
                }// end id


				int bindex = 0;
				// z is the fastest, 
				for(int ii=0; ii<GlobalC::pw.bx; ii++)
				{
					const int ipart = (i*GlobalC::pw.bx + ii) * GlobalC::pw.ncy*GlobalC::pw.nczp;
					for(int jj=0; jj<GlobalC::pw.by; jj++)
					{
						const int jpart = (j*GlobalC::pw.by + jj) * GlobalC::pw.nczp;
						for(int kk=0; kk<GlobalC::pw.bz; kk++)
						{
							vindex[bindex] = (k*GlobalC::pw.bz + kk-GlobalC::pw.nczp_start) + jpart + ipart; 
						//	assert(vindex[bindex] < GlobalC::pw.nrxx);
							++bindex;
						}
					}
				}

                for (int ia1=0; ia1<size; ia1++)
                {
                    const int mcell_index1 = GlobalC::GridT.bcell_start[grid_index] + ia1;
					const int iat = GlobalC::GridT.which_atom[mcell_index1];
					const int T1 = GlobalC::ucell.iat2it[iat];
					Atom *atom1 = &GlobalC::ucell.atoms[T1];
					const int I1 = GlobalC::ucell.iat2ia[iat];
					// get the start index of local orbitals.
					const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
					for (int ib=0; ib<GlobalC::pw.bxyz; ib++)
					{
						if(cal_flag[ib][ia1])
						{
							int iw1_lo = GlobalC::GridT.trace_lo[start1];
							double* psi1 = psir_ylm[ib][ia1];
							double tmp = 0.0;
							for (int iw=0; iw< atom1->nw; ++iw, ++iw1_lo)
							{
								tmp += psi1[iw] * wfc[iw1_lo];
							}//iw
							rho[ vindex[ib] ] += tmp;
						}// cal_flag
					}//ib
				}// ia1

            }// k
        }// j
    }// i

	delete[] vindex;
	delete[] vldr3;
	
	if(max_size!=0) 
	{
		for(int i=0; i<GlobalC::pw.bxyz; i++)
		{
			for(int j=0; j<max_size; j++) 
			{
				delete[] dr[i][j];
				delete[] psir_ylm[i][j];
			}
			delete[] dr[i];
			delete[] distance[i];
			delete[] psir_ylm[i];
			delete[] cal_flag[i];
		}
		delete[] dr;
		delete[] distance;
		delete[] psir_ylm;
		delete[] cal_flag;
	}

    return;
}

