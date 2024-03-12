#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "structure_factor.h"
#include "module_base/constants.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/math_bspline.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/libm/libm.h"


#ifdef _OPENMP
#include <omp.h>
#endif

Structure_Factor::Structure_Factor()
{

}

Structure_Factor::~Structure_Factor()
{
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            delmem_cd_op()(gpu_ctx, this->c_eigts1);
            delmem_cd_op()(gpu_ctx, this->c_eigts2);
            delmem_cd_op()(gpu_ctx, this->c_eigts3);
        }
        delmem_zd_op()(gpu_ctx, this->z_eigts1);
        delmem_zd_op()(gpu_ctx, this->z_eigts2);
        delmem_zd_op()(gpu_ctx, this->z_eigts3);
    }
    else {
        if (GlobalV::precision_flag == "single") {
            delmem_ch_op()(cpu_ctx, this->c_eigts1);
            delmem_ch_op()(cpu_ctx, this->c_eigts2);
            delmem_ch_op()(cpu_ctx, this->c_eigts3);
        }
        // There's no need to delete double precision pointers while in a CPU environment.
    }
}

// called in input.cpp
void Structure_Factor::set(const ModulePW::PW_Basis* rho_basis_in, const int& nbspline_in)
{
    ModuleBase::TITLE("PW_Basis","set");
    this->rho_basis = rho_basis_in;
    this->nbspline = nbspline_in;
    return;
}

// Peize Lin optimize and add OpenMP 2021.04.01
//  Calculate structure factor
void Structure_Factor::setup_structure_factor(UnitCell* Ucell, const ModulePW::PW_Basis* rho_basis)
{
    ModuleBase::TITLE("PW_Basis","setup_structure_factor");
    ModuleBase::timer::tick("PW_Basis","setup_struc_factor");
    const std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;

    this->strucFac.create(Ucell->ntype, rho_basis->npw);
    ModuleBase::Memory::record("SF::strucFac", sizeof(std::complex<double>) * Ucell->ntype*rho_basis->npw);

//	std::string outstr;
//	outstr = GlobalV::global_out_dir + "strucFac.dat"; 
//	std::ofstream ofs( outstr.c_str() ) ;
    bool usebspline;
    if(nbspline > 0)   usebspline = true;
    else    usebspline = false;
    
    if(usebspline)
    {
        nbspline = int((nbspline+1)/2)*2; // nbspline must be a positive even number.
        this->bspline_sf(nbspline,Ucell, rho_basis);
    }
    else
    {
        for (int it=0; it<Ucell->ntype; it++)
        {
	    	const int na = Ucell->atoms[it].na;
	    	const ModuleBase::Vector3<double> * const tau = Ucell->atoms[it].tau;
#ifdef _OPENMP
		    #pragma omp parallel for
#endif
            for (int ig=0; ig<rho_basis->npw; ig++)
            {
		    	const ModuleBase::Vector3<double> gcar_ig = rho_basis->gcar[ig];
                std::complex<double> sum_phase = ModuleBase::ZERO;
                for (int ia=0; ia<na; ia++)
                {
                    // e^{-i G*tau}
                    sum_phase += ModuleBase::libm::exp( ci_tpi * (gcar_ig * tau[ia]) );
                }
                this->strucFac(it,ig) = sum_phase;
            }
        }
    }

//	ofs.close();

    int i,j; //ng;
    this->eigts1.create(Ucell->nat, 2*rho_basis->nx + 1);
    this->eigts2.create(Ucell->nat, 2*rho_basis->ny + 1);
    this->eigts3.create(Ucell->nat, 2*rho_basis->nz + 1);

    ModuleBase::Memory::record("SF::eigts123",sizeof(std::complex<double>) * (Ucell->nat*2 * (rho_basis->nx + rho_basis->ny + rho_basis->nz) + 3));

    ModuleBase::Vector3<double> gtau;
    int inat = 0;
    for (i = 0; i < Ucell->ntype; i++)
    {
        if (GlobalV::test_pw > 1)
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"eigts",i);
        }
        for (j = 0; j < Ucell->atoms[i].na;j++)
        {
            gtau = Ucell->G * Ucell->atoms[i].tau[j];  //HLX: fixed on 10/13/2006
#ifdef _OPENMP
#pragma omp parallel
{
		    #pragma omp for schedule(static, 16)
#endif
            for (int n1 = -rho_basis->nx; n1 <= rho_basis->nx;n1++)
            {
                double arg = n1 * gtau.x;
                this->eigts1(inat, n1 + rho_basis->nx) = ModuleBase::libm::exp( ci_tpi*arg  );
            }
#ifdef _OPENMP
		    #pragma omp for schedule(static, 16)
#endif
            for (int n2 = -rho_basis->ny; n2 <= rho_basis->ny;n2++)
            {
                double arg = n2 * gtau.y;
                this->eigts2(inat, n2 + rho_basis->ny) = ModuleBase::libm::exp( ci_tpi*arg );
            }
#ifdef _OPENMP
		    #pragma omp for schedule(static, 16)
#endif
            for (int n3 = -rho_basis->nz; n3 <= rho_basis->nz;n3++)
            {
                double arg = n3 * gtau.z;
                this->eigts3(inat, n3 + rho_basis->nz) = ModuleBase::libm::exp( ci_tpi*arg );
            }
#ifdef _OPENMP
}
#endif
            inat++;
        }
    }
    if (GlobalV::device_flag == "gpu") {
        if (GlobalV::precision_flag == "single") {
            resmem_cd_op()(gpu_ctx, this->c_eigts1, Ucell->nat * (2 * rho_basis->nx + 1));
            resmem_cd_op()(gpu_ctx, this->c_eigts2, Ucell->nat * (2 * rho_basis->ny + 1));
            resmem_cd_op()(gpu_ctx, this->c_eigts3, Ucell->nat * (2 * rho_basis->nz + 1));
            castmem_z2c_h2d_op()(gpu_ctx, cpu_ctx, this->c_eigts1, this->eigts1.c, Ucell->nat * (2 * rho_basis->nx + 1));
            castmem_z2c_h2d_op()(gpu_ctx, cpu_ctx, this->c_eigts2, this->eigts2.c, Ucell->nat * (2 * rho_basis->ny + 1));
            castmem_z2c_h2d_op()(gpu_ctx, cpu_ctx, this->c_eigts3, this->eigts3.c, Ucell->nat * (2 * rho_basis->nz + 1));
        }
        resmem_zd_op()(gpu_ctx, this->z_eigts1, Ucell->nat * (2 * rho_basis->nx + 1));
        resmem_zd_op()(gpu_ctx, this->z_eigts2, Ucell->nat * (2 * rho_basis->ny + 1));
        resmem_zd_op()(gpu_ctx, this->z_eigts3, Ucell->nat * (2 * rho_basis->nz + 1));
        syncmem_z2z_h2d_op()(gpu_ctx, cpu_ctx, this->z_eigts1, this->eigts1.c, Ucell->nat * (2 * rho_basis->nx + 1));
        syncmem_z2z_h2d_op()(gpu_ctx, cpu_ctx, this->z_eigts2, this->eigts2.c, Ucell->nat * (2 * rho_basis->ny + 1));
        syncmem_z2z_h2d_op()(gpu_ctx, cpu_ctx, this->z_eigts3, this->eigts3.c, Ucell->nat * (2 * rho_basis->nz + 1));
    }
    else {
        if (GlobalV::precision_flag == "single") {
            resmem_ch_op()(cpu_ctx, this->c_eigts1, Ucell->nat * (2 * rho_basis->nx + 1));
            resmem_ch_op()(cpu_ctx, this->c_eigts2, Ucell->nat * (2 * rho_basis->ny + 1));
            resmem_ch_op()(cpu_ctx, this->c_eigts3, Ucell->nat * (2 * rho_basis->nz + 1));
            castmem_z2c_h2h_op()(cpu_ctx, cpu_ctx, this->c_eigts1, this->eigts1.c, Ucell->nat * (2 * rho_basis->nx + 1));
            castmem_z2c_h2h_op()(cpu_ctx, cpu_ctx, this->c_eigts2, this->eigts2.c, Ucell->nat * (2 * rho_basis->ny + 1));
            castmem_z2c_h2h_op()(cpu_ctx, cpu_ctx, this->c_eigts3, this->eigts3.c, Ucell->nat * (2 * rho_basis->nz + 1));
        }
        this->z_eigts1 = this->eigts1.c;
        this->z_eigts2 = this->eigts2.c;
        this->z_eigts3 = this->eigts3.c;
        // There's no need to delete double precision pointers while in a CPU environment.
    }
    ModuleBase::timer::tick("PW_Basis","setup_struc_factor"); 
    return;
}

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
void Structure_Factor::bspline_sf(const int norder, UnitCell* Ucell, const ModulePW::PW_Basis* rho_basis)
{
    double *r = new double [rho_basis->nxyz]; 
    double *tmpr = new double[rho_basis->nrxx];
    double *zpiece = new double[rho_basis->nxy];
    std::complex<double> *b1 = new std::complex<double> [rho_basis->nx];
    std::complex<double> *b2 = new std::complex<double> [rho_basis->ny];
    std::complex<double> *b3 = new std::complex<double> [rho_basis->nz];

    for (int it=0; it<Ucell->ntype; it++)
    {
		const int na = Ucell->atoms[it].na;
		const ModuleBase::Vector3<double> * const taud = Ucell->atoms[it].taud;
        ModuleBase::GlobalFunc::ZEROS(r,rho_basis->nxyz);

        //A parallel algorithm can be added in the future.
#ifdef _OPENMP
		#pragma omp parallel for
#endif
        for(int ia = 0 ; ia < na ; ++ia)
        {
            double gridx = taud[ia].x * rho_basis->nx;
            double gridy = taud[ia].y * rho_basis->ny;
            double gridz = taud[ia].z * rho_basis->nz;
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
                int icz = int(rho_basis->nz*10-iz+floor(gridz))%rho_basis->nz;
                for(int iy = 0 ; iy <= norder ; ++iy)
                {
                    int icy = int(rho_basis->ny*10-iy+floor(gridy))%rho_basis->ny;
                    for(int ix = 0 ; ix <= norder ; ++ix )
                    {
                        int icx = int(rho_basis->nx*10-ix+floor(gridx))%rho_basis->nx;
#ifdef _OPENMP
		                #pragma omp atomic
#endif
                        r[icz*rho_basis->ny*rho_basis->nx + icx*rho_basis->ny + icy] += bsz.bezier_ele(iz) 
                                                 * bsy.bezier_ele(iy) 
                                                 * bsx.bezier_ele(ix); 
                    }
                }
            }
        }
        
        //distribute data to different processors for UFFT
        //---------------------------------------------------
        for(int iz = 0; iz < rho_basis->nz; iz++)
	    {
	    	if(GlobalV::MY_RANK==0)
	    	{
#ifdef _OPENMP
		    #pragma omp parallel for schedule(static, 512)
#endif
	    		for(int ir = 0; ir < rho_basis->nxy; ir++)
	    		{
	    			zpiece[ir] = r[iz*rho_basis->nxy + ir];
	    		}
	    	}
        
        #ifdef __MPI
	    	GlobalC::Pgrid.zpiece_to_all(zpiece, iz, tmpr);
        #endif
        
	    }
        //---------------------------------------------------

        //It should be optimized with r2c
        rho_basis->real2recip(tmpr, &strucFac(it,0));
        this->bsplinecoef(b1,b2,b3,rho_basis->nx, rho_basis->ny, rho_basis->nz, norder);
#ifdef _OPENMP
		#pragma omp parallel for schedule(static, 128)
#endif
        for(int ig = 0 ; ig < rho_basis->npw ; ++ig)
        {
           int idx = int(rho_basis->gdirect[ig].x+0.1+rho_basis->nx)%rho_basis->nx;
           int idy = int(rho_basis->gdirect[ig].y+0.1+rho_basis->ny)%rho_basis->ny;
           int idz = int(rho_basis->gdirect[ig].z+0.1+rho_basis->nz)%rho_basis->nz;
           strucFac(it,ig) *= ( b1[idx] * b2[idy] * b3[idz] * double(rho_basis->nxyz) );
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

void Structure_Factor:: bsplinecoef(std::complex<double> *b1, std::complex<double> *b2, std::complex<double> *b3, 
                        const int nx, const int ny, const int nz, const int norder)
{
    const std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
    ModuleBase::Bspline bsp;
    bsp.init(norder, 1, 0);
    bsp.getbspline(1.0);
#ifdef _OPENMP
#pragma omp parallel
{
	#pragma omp for schedule(static, 16)
#endif
    for(int ix = 0 ; ix < nx ; ++ix)
    {
        std::complex<double> fracx=0;
        for(int io = 0 ; io < norder - 1 ; ++io)
        {
            fracx += bsp.bezier_ele(io)*ModuleBase::libm::exp(ci_tpi*double(ix)/double(nx)*double(io));
        }
        b1[ix] = ModuleBase::libm::exp(ci_tpi*double(norder*ix)/double(nx))/fracx;
    }
#ifdef _OPENMP
	#pragma omp for schedule(static, 16)
#endif
    for(int iy = 0 ; iy < ny ; ++iy)
    {
        std::complex<double> fracy=0;
        for(int io = 0 ; io < norder - 1 ; ++io)
        {
            fracy += bsp.bezier_ele(io)*ModuleBase::libm::exp(ci_tpi*double(iy)/double(ny)*double(io));
        }
        b2[iy] = ModuleBase::libm::exp(ci_tpi*double(norder*iy)/double(ny))/fracy;
    }
#ifdef _OPENMP
	#pragma omp for schedule(static, 16)
#endif
    for(int iz = 0 ; iz < nz ; ++iz)
    {
        std::complex<double> fracz=0;
        for(int io = 0 ; io < norder - 1 ; ++io)
        {
            fracz += bsp.bezier_ele(io)*ModuleBase::libm::exp(ci_tpi*double(iz)/double(nz)*double(io));
        }
        b3[iz] = ModuleBase::libm::exp(ci_tpi*double(norder*iz)/double(nz))/fracz;
    }
#ifdef _OPENMP
}
#endif
}

template <>
std::complex<float> * Structure_Factor::get_eigts1_data() const
{
    return this->c_eigts1;
}
template <>
std::complex<double> * Structure_Factor::get_eigts1_data() const
{
    return this->z_eigts1;
}

template <>
std::complex<float> * Structure_Factor::get_eigts2_data() const
{
    return this->c_eigts2;
}
template <>
std::complex<double> * Structure_Factor::get_eigts2_data() const
{
    return this->z_eigts2;
}

template <>
std::complex<float> * Structure_Factor::get_eigts3_data() const
{
    return this->c_eigts3;
}
template <>
std::complex<double> * Structure_Factor::get_eigts3_data() const
{
    return this->z_eigts3;
}