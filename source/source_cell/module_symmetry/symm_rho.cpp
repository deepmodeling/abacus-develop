#include "symmetry.h"
using namespace ModuleSymmetry;

#include "source_base/libm/libm.h"
#include "module_parameter/parameter.h"

void Symmetry::rho_symmetry( double *rho,
                             const int &nr1, const int &nr2, const int &nr3)
{
    ModuleBase::timer::tick("Symmetry","rho_symmetry");

	// allocate flag for each FFT grid.
    bool* symflag = new bool[nr1 * nr2 * nr3];
    for (int i=0; i<nr1*nr2*nr3; i++)
    {
        symflag[i] = false;
    }

    assert(nrotk >0 );
    assert(nrotk <=48 );
    int *ri = new int[nrotk];
    int *rj = new int[nrotk];
    int *rk = new int[nrotk];

    int ci = 0;
    for (int i = 0; i< nr1; ++i)
    {
        for (int j = 0; j< nr2; ++j)
        {
            for (int k = 0; k< nr3; ++k)
            {
                if (!symflag[i * nr2 * nr3 + j * nr3 + k])
                {
                    double sum = 0;

                    for (int isym = 0; isym < nrotk; ++isym)
                    {
                        this->rotate(gmatrix[isym], gtrans[isym], i, j, k, nr1, nr2, nr3, ri[isym], rj[isym], rk[isym]);
                        const int index = ri[isym] * nr2 * nr3 + rj[isym] * nr3 + rk[isym];
                        sum += rho[ index ];
                    }
                    sum /= nrotk;

                    for (int isym = 0; isym < nrotk; ++isym)
                    {
                        const int index = ri[isym] * nr2 * nr3 + rj[isym] * nr3 + rk[isym];
                        rho[index] = sum;
                        symflag[index] = true;
                    }
                }
            }
        }
    }

    delete[] symflag;
    delete[] ri;
    delete[] rj;
    delete[] rk;
    ModuleBase::timer::tick("Symmetry","rho_symmetry");
}

void Symmetry::rhog_symmetry(std::complex<double> *rhogtot, 
    int* ixyz2ipw, const int &nx, const int &ny, const int &nz, 
    const int &fftnx, const int &fftny, const int &fftnz)
{
	ModuleBase::timer::tick("Symmetry","rhog_symmetry");
	// ----------------------------------------------------------------------
	// the current way is to cluster the FFT grid points into groups in advance.
	// and use OpenMP to realize parallel calculation, one thread works in one group.
	// ----------------------------------------------------------------------

    const int nxyz = fftnx*fftny*fftnz;
    assert(nxyz>0);

	// allocate flag for each FFT grid.
    // which group the grid belongs to
	int* symflag = new int[nxyz];

    // which rotration operation the grid corresponds to
	int(*isymflag)[48] = new int[nxyz][48];

    // group information
	int(*table_xyz)[48] = new int[nxyz][48];

    // how many symmetry operations have been covered
	int* count_xyz = new int[nxyz];

	for (int i = 0; i < nxyz; i++)
	{
		symflag[i] = -1;
	}
	int group_index = 0;

	assert(nrotk >0 );
	assert(nrotk <=48 );

	//map the gmatrix to inv
    std::vector<int>invmap(this->nrotk, -1);
    this->gmatrix_invmap(kgmatrix, nrotk, invmap.data());

	// ------------------------------------------------------------------------
	// This code defines a lambda function called "rotate_recip" that takes 
	// a 3x3 matrix and a 3D vector as input. It performs a rotation operation 
	// on the vector using the matrix and returns the rotated vector. 
	// Specifically, it calculates the new coordinates of the vector after 
	// the rotation and applies periodic boundary conditions to ensure that 
	// the coordinates are within the FFT-grid dimensions. 
	// The rotated vector is returned by modifying the input vector.
	// ------------------------------------------------------------------------
    //rotate function (different from real space, without scaling gmatrix)
    auto rotate_recip = [&] (ModuleBase::Matrix3& g, ModuleBase::Vector3<int>& g0, int& ii, int& jj, int& kk) 
    {
        ii = int(g.e11 * g0.x + g.e21 * g0.y + g.e31 * g0.z) ;
        if (ii < 0)
        {
            ii += 10 * nx;
        }
        ii = ii%nx;
        jj = int(g.e12 * g0.x + g.e22 * g0.y + g.e32 * g0.z) ;
        if (jj < 0)
        {
            jj += 10 * ny;
        }
        jj = jj%ny;
        kk = int(g.e13 * g0.x + g.e23 * g0.y + g.e33 * g0.z);
        if (kk < 0)
        {
            kk += 10 * nz;
        }
        kk = kk%nz;
        return;
    };

	// ------------------------------------------------------------------------
    // Trying to group fft grids first.
    // It iterates over each FFT-grid point and checks if it is within the 
    // PW-sphere. If it is, put all the FFT-grid points connected by the 
    // rotation operation into one group( the index is stored in int(*table_xyz)).
    // The code marks the point as processed to avoid redundant calculations
    // by using int* symflag.
	// ------------------------------------------------------------------------

    ModuleBase::timer::tick("Symmetry","group_fft_grids");
    for (int i = 0; i< fftnx; ++i)
    {
        //tmp variable
        ModuleBase::Vector3<int> tmp_gdirect0(0, 0, 0);
        tmp_gdirect0.x=(i>int(nx/2)+1)?(i-nx):i;
        for (int j = 0; j< fftny; ++j)
        {
            tmp_gdirect0.y=(j>int(ny/2)+1)?(j-ny):j;
            for (int k = 0; k< fftnz; ++k)
            {
                int ixyz0=(i*fftny+j)*fftnz+k;
                if (symflag[ixyz0] == -1)
                {
                    int ipw0=ixyz2ipw[ixyz0];
                    //if a fft-grid is not in pw-sphere, just do not consider it.
                    if (ipw0 == -1) {
                        continue;
                    }
                    tmp_gdirect0.z=(k>int(nz/2)+1)?(k-nz):k;
                    int rot_count=0;
                    for (int isym = 0; isym < nrotk; ++isym)
                    {
                        if (invmap[isym] < 0 || invmap[isym] > nrotk) { continue; }
                        //tmp variables  
                        int ii, jj, kk=0;
                        rotate_recip(kgmatrix[invmap[isym]], tmp_gdirect0, ii, jj, kk);
                        if(ii>=fftnx || jj>=fftny || kk>= fftnz)
                        {
                            if(!PARAM.globalv.gamma_only_pw)
                            {
                                std::cout << " ROTATE OUT OF FFT-GRID IN RHOG_SYMMETRY !" << std::endl;
		                        ModuleBase::QUIT();
                            }
                            // for gamma_only_pw, just do not consider this rotation.
                            continue;
                        }
                        int ixyz=(ii*fftny+jj)*fftnz+kk;
                        //fft-grid index to (ip, ig)
                        int ipw=ixyz2ipw[ixyz];
                        if(ipw==-1) //not in pw-sphere
                        {
                            continue;   //else, just skip it
                        }
                        symflag[ixyz] = group_index;
                        isymflag[group_index][rot_count] = invmap[isym];
                        table_xyz[group_index][rot_count] = ixyz;
                        ++rot_count;
                        assert(rot_count <= nrotk);
                        count_xyz[group_index] = rot_count;
                    }
                group_index++;
                }
            }
        }
    }
    ModuleBase::timer::tick("Symmetry","group_fft_grids");

	// -------------------------------------------------------------------
	//  This code performs symmetry operations on the reciprocal space 
	//	charge density using FFT-grids. It iterates over each FFT-grid 
	//	point in a particular group, applies a phase factor and sums the 
	//	charge density over the symmetry operations, and then divides by 
	//	the number of symmetry operations. Finally, it updates the charge
	//	density for each FFT-grid point using the calculated sum.
	// -------------------------------------------------------------------

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int g_index = 0; g_index < group_index; g_index++)
	{
		// record the index and gphase but not the final gdirect for each symm-opt
		int *ipw_record = new int[nrotk];
		int *ixyz_record = new int[nrotk];
		std::complex<double>* gphase_record = new std::complex<double> [nrotk];
		std::complex<double> sum(0, 0);
		int rot_count=0;

		for (int c_index = 0; c_index < count_xyz[g_index]; ++c_index)
		{
			int ixyz0 = table_xyz[g_index][c_index];
			int ipw0 = ixyz2ipw[ixyz0];

			if (symflag[ixyz0] == g_index)
			{
				// note : do not use PBC after rotation. 
				// we need a real gdirect to get the correspoding rhogtot.
				int k = ixyz0%fftnz;
				int j = ((ixyz0-k)/fftnz)%fftny;
				int i = ((ixyz0-k)/fftnz-j)/fftny;

				//fft-grid index to gdirect
				ModuleBase::Vector3<double> tmp_gdirect_double(0.0, 0.0, 0.0);
				tmp_gdirect_double.x=static_cast<double>((i>int(nx/2)+1)?(i-nx):i);
				tmp_gdirect_double.y=static_cast<double>((j>int(ny/2)+1)?(j-ny):j);
				tmp_gdirect_double.z=static_cast<double>((k>int(nz/2)+1)?(k-nz):k);

				//calculate phase factor
				tmp_gdirect_double = tmp_gdirect_double * ModuleBase::TWO_PI;

				double cos_arg = 0.0, sin_arg = 0.0;
				double arg_gtrans = tmp_gdirect_double * gtrans[isymflag[g_index][c_index]];

				std::complex<double> phase_gtrans (ModuleBase::libm::cos(arg_gtrans), 
						ModuleBase::libm::sin(arg_gtrans));

				// for each pricell in supercell:
				for (int ipt = 0;ipt < ((ModuleSymmetry::Symmetry::pricell_loop) ? this->ncell : 1);++ipt)
				{
					double arg = tmp_gdirect_double * ptrans[ipt];
					double tmp_cos = 0.0, tmp_sin = 0.0;
					ModuleBase::libm::sincos(arg, &tmp_sin, &tmp_cos);
					cos_arg += tmp_cos;
					sin_arg += tmp_sin;
				}

				// add nothing to sum, so don't consider this isym into rot_count
				cos_arg/=static_cast<double>(ncell);
				sin_arg/=static_cast<double>(ncell);

				//deal with double-zero
				if (equal(cos_arg, 0.0) && equal(sin_arg, 0.0)) 
				{
					continue;
				}

				std::complex<double> gphase(cos_arg, sin_arg);
				gphase = phase_gtrans * gphase;

				//deal with small difference from 1
				if (equal(gphase.real(), 1.0) && equal(gphase.imag(), 0)) 
				{
					gphase = std::complex<double>(1.0, 0.0);
				}

				gphase_record[rot_count]=gphase;
				sum += rhogtot[ipw0]*gphase;
				//record
				ipw_record[rot_count]=ipw0;
				ixyz_record[rot_count]=ixyz0;
				++rot_count;
				//assert(rot_count<=nrotk);
			}//end if section
		}//end c_index loop
		sum /= rot_count;
		for (int isym = 0; isym < rot_count; ++isym)
		{
			rhogtot[ipw_record[isym]] = sum/gphase_record[isym];
		}
		
		//Clean the records variables for each fft grid point
		delete[] ipw_record;
		delete[] ixyz_record;
		delete[] gphase_record;
	}//end g_index loop

	delete[] symflag;
	delete[] isymflag;
	delete[] table_xyz;
	delete[] count_xyz;
	ModuleBase::timer::tick("Symmetry","rhog_symmetry");
}
