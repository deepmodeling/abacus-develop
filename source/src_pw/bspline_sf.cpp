#include "pw_basis.h"
#include "global.h"
#include "unistd.h"
#include "../module_base/math_bspline.h"

//
//DESCRIPTION:
//    Calculate structure factor with Cardinal B-spline interpolation
//    Ref: J. Chem. Phys. 103, 8577 (1995)
//    qianrui create 2021-9-17
//INPUT LIST:
//    norder: the order of Cardinal B-spline base functions
//FURTHER OPTIMIZATION:
//    1. Use "r2c" fft
//    2. Add parallel algorithm for fftw or na loop
//
void PW_Basis::bspline_sf(const int norder)
{
    double *r = new double [ncxyz]; 
    double *tmpr = new double[nrxx];
    double *zpiece = new double[nx*ny];
    complex<double> *b1 = new complex<double> [nx];
    complex<double> *b2 = new complex<double> [ny];
    complex<double> *b3 = new complex<double> [nz];
    GlobalC::UFFT.allocate();
    

    for (int it=0; it<Ucell->ntype; it++)
    {
		const int na = Ucell->atoms[it].na;
		const ModuleBase::Vector3<double> * const taud = Ucell->atoms[it].taud;
        ModuleBase::GlobalFunc::ZEROS(r,ncxyz);

        //A parallel algorithm can be added in the future.
        for(int ia = 0 ; ia < na ; ++ia)
        {
            double gridx = taud[ia].x * nx;
            double gridy = taud[ia].y * ny;
            double gridz = taud[ia].z * nz;
            double dx = gridx - floor(gridx);
            double dy = gridy - floor(gridy);
            double dz = gridz - floor(gridz);
            //I'm not sure if there is a mod function for double data

            ModuleBase::Bspline bsx, bsy, bsz;
            bsx.init(norder, 1, 0);
            bsy.init(norder, 1, 0);
            bsz.init(norder, 1, 0);
            bsx.getbspline(dx);
            bsy.getbspline(dy);
            bsz.getbspline(dz);

            for(int iz = 0 ; iz <= norder ; ++iz)
            {
                int icz = int(nz*10-iz+floor(gridz))%nz;
                for(int iy = 0 ; iy <= norder ; ++iy)
                {
                    int icy = int(ny*10-iy+floor(gridy))%ny;
                    for(int ix = 0 ; ix <= norder ; ++ix )
                    {
                        int icx = int(nx*10-ix+floor(gridx))%nx;
                        r[icz*ny*nx + icx*ny + icy] += bsz.bezier_ele(iz) 
                                                 * bsy.bezier_ele(iy) 
                                                 * bsx.bezier_ele(ix); 
                    }
                }
            }
        }
        
        //distribute data to different processors for UFFT
        //---------------------------------------------------
        for(int iz = 0; iz < nz; iz++)
	    {
	    	ModuleBase::GlobalFunc::ZEROS(zpiece, nx*ny);
	    	if(GlobalV::MY_RANK==0)
	    	{
	    		for(int ir = 0; ir < ny * nx; ir++)
	    		{
	    			zpiece[ir] = r[iz*ny*nx + ir];
	    		}
	    	}
        
        #ifdef __MPI
	    	GlobalC::Pgrid.zpiece_to_all(zpiece, iz, tmpr);
        #endif
        
	    }
        //---------------------------------------------------

        //It should be optimized with r2c
        GlobalC::UFFT.ToReciSpace(tmpr, &strucFac(it,0));
        this->bsplinecoef(b1,b2,b3,norder);
        for(int ig = 0 ; ig < ngmc ; ++ig)
        {
           int idx = int(gdirect[ig].x+0.1+this->nx)%this->nx;
           int idy = int(gdirect[ig].y+0.1+this->ny)%this->ny;
           int idz = int(gdirect[ig].z+0.1+this->nz)%this->nz;
           strucFac(it,ig) *= ( b1[idx] * b2[idy] * b3[idz] * double(this->ncxyz) );
        }
    }   
    delete[] r;
    delete[] tmpr;
    delete[] zpiece; 
    delete[] b1;
    delete[] b2;
    delete[] b3;

    return;
}

void PW_Basis:: bsplinecoef(complex<double> *b1, complex<double> *b2, complex<double> *b3, const int norder)
{
    const std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
    ModuleBase::Bspline bsp;
    bsp.init(norder, 1, 0);
    bsp.getbspline(1.0);
    for(int ix = 0 ; ix < this->nx ; ++ix)
    {
        complex<double> fracx=0;
        for(int io = 0 ; io < norder - 1 ; ++io)
        {
            fracx += bsp.bezier_ele(io)*exp(ci_tpi*double(ix)/double(nx)*double(io));
        }
        b1[ix] = exp(ci_tpi*double(norder*ix)/double(nx))/fracx;
    }
    for(int iy = 0 ; iy < this->ny ; ++iy)
    {
        complex<double> fracy=0;
        for(int io = 0 ; io < norder - 1 ; ++io)
        {
            fracy += bsp.bezier_ele(io)*exp(ci_tpi*double(iy)/double(ny)*double(io));
        }
        b2[iy] = exp(ci_tpi*double(norder*iy)/double(ny))/fracy;
    }
    for(int iz = 0 ; iz < this->nz ; ++iz)
    {
        complex<double> fracz=0;
        for(int io = 0 ; io < norder - 1 ; ++io)
        {
            fracz += bsp.bezier_ele(io)*exp(ci_tpi*double(iz)/double(nz)*double(io));
        }
        b3[iz] = exp(ci_tpi*double(norder*iz)/double(nz))/fracz;
    }
}
