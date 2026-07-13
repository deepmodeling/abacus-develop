#include "symmetry.h"
using namespace ModuleSymmetry;

#include "source_base/libm/libm.h"

namespace
{
	// ------------------------------------------------------------------------
	// Rotating reciprocal-space FFT-grid vector (with PBC)
	// The rotated vector is returned via ii, jj, kk.
	// ------------------------------------------------------------------------
	//rotate function (different from real space, without scaling gmatrix)
	static inline void rotate_recip(ModuleBase::Matrix3& g, ModuleBase::Vector3<int>& g0, int& ii, int& jj, int& kk,
		const int& nx, const int& ny, const int& nz)
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
	}

	// ------------------------------------------------------------------------
	// Trying to group fft grids first.
	// It iterates over each FFT-grid point and checks if it is within the
	// PW-sphere. If it is, put all the FFT-grid points connected by the
	// rotation operation into one group( the index is stored in int(*table_xyz)).
	// The code marks the point as processed to avoid redundant calculations
	// by using int* symflag.
	// This grouping is purely spatial (depends only on kgmatrix/invmap and the
	// FFT-grid geometry), so it is shared between rhog_symmetry and rhog_symmetry_soc;
	// the differing spin/phase accumulation happens after this call.
	// ------------------------------------------------------------------------
	static void group_fft_grids(const int& nrotk, ModuleBase::Matrix3* kgmatrix, const std::vector<int>& invmap,
		int* ixyz2ipw, const int& nx, const int& ny, const int& nz,
		const int& fftnx, const int& fftny, const int& fftnz,
		int* symflag, int (*isymflag)[48], int (*table_xyz)[48], int* count_xyz, int& group_index)
	{
	    ModuleBase::timer::start("Symmetry","group_fft_grids");
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
	                        rotate_recip(kgmatrix[invmap[isym]], tmp_gdirect0, ii, jj, kk, nx, ny, nz);
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
	    ModuleBase::timer::end("Symmetry","group_fft_grids");
	}
} // namespace

void Symmetry::rho_symmetry( double *rho,
                             const int &nr1, const int &nr2, const int &nr3)
{
    ModuleBase::timer::start("Symmetry","rho_symmetry");

    assert(nr1>0);
    assert(nr2>0);
    assert(nr3>0);

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
    ModuleBase::timer::end("Symmetry","rho_symmetry");
}

void Symmetry::rhog_symmetry(std::complex<double> *rhogtot, 
    int* ixyz2ipw, const int &nx, const int &ny, const int &nz, 
    const int &fftnx, const int &fftny, const int &fftnz,
    const bool gamma_only_pw)
{
    ModuleBase::timer::start("Symmetry","rhog_symmetry");
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
    // Group the FFT grids connected by symmetry (spatial grouping only).
    // gamma_only_pw handling is preserved inside the shared group_fft_grids helper
    // (it reads PARAM.globalv.gamma_only_pw).
    // ------------------------------------------------------------------------
    group_fft_grids(nrotk, kgmatrix, invmap, ixyz2ipw, nx, ny, nz, fftnx, fftny, fftnz,
        symflag, isymflag, table_xyz, count_xyz, group_index);

    // -------------------------------------------------------------------
    //  This code performs symmetry operations on the reciprocal space 
    //    charge density using FFT-grids. It iterates over each FFT-grid 
    //    point in a particular group, applies a phase factor and sums the 
    //    charge density over the symmetry operations, and then divides by 
    //    the number of symmetry operations. Finally, it updates the charge
    //    density for each FFT-grid point using the calculated sum.
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
        if (rot_count!=0) sum/= rot_count;
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
    ModuleBase::timer::end("Symmetry","rhog_symmetry");
}

void Symmetry::rhog_symmetry_soc(std::complex<double>* rhogtot_x, std::complex<double>* rhogtot_y,
    std::complex<double>* rhogtot_z, const ModuleBase::Matrix3* wspin,
    int* ixyz2ipw, const int &nx, const int &ny, const int &nz,
    const int & fftnx, const int &fftny, const int &fftnz)
{
	ModuleBase::timer::start("Symmetry","rhog_symmetry_soc");
	// The grouping of FFT grid points into symmetry-connected orbits is purely spatial and
	// therefore identical to rhog_symmetry. Only the accumulation/write-back is changed:
	// the three spin components are mixed by W(g) (rotated to the orbit-representative frame
	// with W(g)^T on the way in, and back with W(g) on the way out), exactly as the scalar
	// version uses the phase factor gphase.

    const int nxyz = fftnx*fftny*fftnz;
    assert(nxyz>0);

	int* symflag = new int[nxyz];
	int(*isymflag)[48] = new int[nxyz][48];
	int(*table_xyz)[48] = new int[nxyz][48];
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

    // Group the FFT grids connected by symmetry (spatial grouping only, shared
    // with rhog_symmetry); see the shared helpers rotate_recip/group_fft_grids
    // defined above.
    group_fft_grids(nrotk, kgmatrix, invmap, ixyz2ipw, nx, ny, nz, fftnx, fftny, fftnz,
        symflag, isymflag, table_xyz, count_xyz, group_index);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int g_index = 0; g_index < group_index; g_index++)
	{
		int *ipw_record = new int[nrotk];
		int *ixyz_record = new int[nrotk];
		int *sym_record = new int[nrotk];
		std::complex<double>* gphase_record = new std::complex<double> [nrotk];
		// orbit-representative-frame spin vector accumulated over the symmetry operations
		std::complex<double> sum_x(0, 0), sum_y(0, 0), sum_z(0, 0);
		int rot_count=0;

		for (int c_index = 0; c_index < count_xyz[g_index]; ++c_index)
		{
			int ixyz0 = table_xyz[g_index][c_index];
			int ipw0 = ixyz2ipw[ixyz0];

			if (symflag[ixyz0] == g_index)
			{
				int k = ixyz0%fftnz;
				int j = ((ixyz0-k)/fftnz)%fftny;
				int i = ((ixyz0-k)/fftnz-j)/fftny;

				ModuleBase::Vector3<double> tmp_gdirect_double(0.0, 0.0, 0.0);
				tmp_gdirect_double.x=static_cast<double>((i>int(nx/2)+1)?(i-nx):i);
				tmp_gdirect_double.y=static_cast<double>((j>int(ny/2)+1)?(j-ny):j);
				tmp_gdirect_double.z=static_cast<double>((k>int(nz/2)+1)?(k-nz):k);

				tmp_gdirect_double = tmp_gdirect_double * ModuleBase::TWO_PI;

				double cos_arg = 0.0, sin_arg = 0.0;
				double arg_gtrans = tmp_gdirect_double * gtrans[isymflag[g_index][c_index]];

				std::complex<double> phase_gtrans (ModuleBase::libm::cos(arg_gtrans),
						ModuleBase::libm::sin(arg_gtrans));

				for (int ipt = 0;ipt < ((ModuleSymmetry::Symmetry::pricell_loop) ? this->ncell : 1);++ipt)
				{
					double arg = tmp_gdirect_double * ptrans[ipt];
					double tmp_cos = 0.0, tmp_sin = 0.0;
					ModuleBase::libm::sincos(arg, &tmp_sin, &tmp_cos);
					cos_arg += tmp_cos;
					sin_arg += tmp_sin;
				}

				cos_arg/=static_cast<double>(ncell);
				sin_arg/=static_cast<double>(ncell);

				if (equal(cos_arg, 0.0) && equal(sin_arg, 0.0))
				{
					continue;
				}

				std::complex<double> gphase(cos_arg, sin_arg);
				gphase = phase_gtrans * gphase;

				if (equal(gphase.real(), 1.0) && equal(gphase.imag(), 0))
				{
					gphase = std::complex<double>(1.0, 0.0);
				}

				// pull this orbit member back to the representative frame: multiply by gphase
				// (removes the translation/phase, as in the scalar version) then by W(g)^T.
				const int isym = isymflag[g_index][c_index];
				const ModuleBase::Matrix3& W = wspin[isym];
				const std::complex<double> vx = rhogtot_x[ipw0] * gphase;
				const std::complex<double> vy = rhogtot_y[ipw0] * gphase;
				const std::complex<double> vz = rhogtot_z[ipw0] * gphase;
				sum_x += W.e11 * vx + W.e21 * vy + W.e31 * vz; // (W^T v)_x
				sum_y += W.e12 * vx + W.e22 * vy + W.e32 * vz; // (W^T v)_y
				sum_z += W.e13 * vx + W.e23 * vy + W.e33 * vz; // (W^T v)_z

				gphase_record[rot_count]=gphase;
				ipw_record[rot_count]=ipw0;
				ixyz_record[rot_count]=ixyz0;
				sym_record[rot_count]=isym;
				++rot_count;
			}//end if section
		}//end c_index loop
		if (rot_count!=0)
		{
			sum_x/= rot_count;
			sum_y/= rot_count;
			sum_z/= rot_count;
		}
		for (int ir = 0; ir < rot_count; ++ir)
		{
			// push the representative-frame value back out to this member: W(g) * S / gphase.
			const ModuleBase::Matrix3& W = wspin[sym_record[ir]];
			const std::complex<double> inv_gphase = 1.0 / gphase_record[ir];
			rhogtot_x[ipw_record[ir]] = (W.e11 * sum_x + W.e12 * sum_y + W.e13 * sum_z) * inv_gphase;
			rhogtot_y[ipw_record[ir]] = (W.e21 * sum_x + W.e22 * sum_y + W.e23 * sum_z) * inv_gphase;
			rhogtot_z[ipw_record[ir]] = (W.e31 * sum_x + W.e32 * sum_y + W.e33 * sum_z) * inv_gphase;
		}

		delete[] ipw_record;
		delete[] ixyz_record;
		delete[] sym_record;
		delete[] gphase_record;
	}//end g_index loop

	delete[] symflag;
	delete[] isymflag;
	delete[] table_xyz;
	delete[] count_xyz;
	ModuleBase::timer::end("Symmetry","rhog_symmetry_soc");
}
