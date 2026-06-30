#include "source_base/global_function.h"
#include "source_io/module_parameter/parameter.h"
#include "structure_factor.h"
#include "source_base/constants.h"
#include "source_base/math_bspline.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_base/libm/libm.h"

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

Structure_Factor::Structure_Factor()
{
    // LCAO basis doesn't support GPU acceleration on this function currently.
    if(PARAM.inp.basis_type == "pw")
    {
        this->device = PARAM.inp.device;
    }
}

Structure_Factor::~Structure_Factor()
{
    if (device == "gpu")
    {
        delmem_cd_op()(this->c_eigts1);
        delmem_cd_op()(this->c_eigts2);
        delmem_cd_op()(this->c_eigts3);
        delmem_zd_op()(this->z_eigts1);
        delmem_zd_op()(this->z_eigts2);
        delmem_zd_op()(this->z_eigts3);
    }
    else
    {
        delmem_ch_op()(this->c_eigts1);
        delmem_ch_op()(this->c_eigts2);
        delmem_ch_op()(this->c_eigts3);
        // There's no need to delete double precision pointers while in a CPU environment.
    }
}

// called in input.cpp
void Structure_Factor::set(const ModulePW::PW_Basis* rho_basis_in, const int& nbspline_in)
{
    ModuleBase::TITLE("Structure_Factor","set");
    this->rho_basis = rho_basis_in;
    this->nbspline = nbspline_in;
    return;
}

// Peize Lin optimize and add OpenMP 2021.04.01
//  Calculate structure factor
void Structure_Factor::setup(const UnitCell* Ucell, const Parallel_Grid& pgrid, const ModulePW::PW_Basis* rho_basis)
{
    ModuleBase::TITLE("Structure_Factor","setup");
    ModuleBase::timer::start("Structure_Factor","setup");

    const std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
    this->ucell = Ucell;
    this->strucFac.create(Ucell->ntype, rho_basis->npw);
    ModuleBase::Memory::record("SF::strucFac", sizeof(std::complex<double>) * Ucell->ntype*rho_basis->npw);

//	std::string outstr;
//	outstr = PARAM.globalv.global_out_dir + "strucFac.dat"; 
//	std::ofstream ofs( outstr.c_str() ) ;
	bool usebspline;
	if(nbspline > 0) 
	{   
		usebspline = true;
	} 
	else 
	{    
		usebspline = false;
	}
    
    if(usebspline)
    {
        nbspline = int((nbspline+1)/2)*2; // nbspline must be a positive even number.
        this->bspline_sf(nbspline, Ucell, pgrid, rho_basis);
    }
    else
    {
        for (int it=0; it<Ucell->ntype; it++)
        {
            const int na = Ucell->atoms[it].na;
            const ModuleBase::Vector3<double> * const tau = Ucell->atoms[it].tau.data();
            // Data race fix: cache shared data to local const variables before OpenMP parallel region
            // TSan detected race condition when accessing rho_basis->npw and rho_basis->gcar directly
            // in parallel loop, even though they are logically read-only
            const int npw = rho_basis->npw;
            const ModuleBase::Vector3<double> * const gcar = rho_basis->gcar;
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for (int ig=0; ig<npw; ig++)
            {
                const ModuleBase::Vector3<double> gcar_ig = gcar[ig];
                std::complex<double> sum_phase = ModuleBase::ZERO;
                for (int ia=0; ia<na; ia++)
                {
                    sum_phase += ModuleBase::libm::exp( ci_tpi * (gcar_ig * tau[ia]) );
                }
                this->strucFac(it,ig) = sum_phase;
            }
        }
    }

//	ofs.close();

    int i=0;
    int j=0;
 
    this->eigts1.create(Ucell->nat, 2*rho_basis->nx + 1);
    this->eigts2.create(Ucell->nat, 2*rho_basis->ny + 1);
    this->eigts3.create(Ucell->nat, 2*rho_basis->nz + 1);

    ModuleBase::Memory::record("SF::eigts123",sizeof(std::complex<double>) 
    * (Ucell->nat*2 * (rho_basis->nx + rho_basis->ny + rho_basis->nz) + 3));

    ModuleBase::Vector3<double> gtau;
    int inat = 0;
    for (i = 0; i < Ucell->ntype; i++)
    {
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
    
    if (device == "gpu") {
        if (PARAM.globalv.has_float_data) {
            resmem_cd_op()(this->c_eigts1, Ucell->nat * (2 * rho_basis->nx + 1));
            resmem_cd_op()(this->c_eigts2, Ucell->nat * (2 * rho_basis->ny + 1));
            resmem_cd_op()(this->c_eigts3, Ucell->nat * (2 * rho_basis->nz + 1));
            castmem_z2c_h2d_op()(this->c_eigts1, this->eigts1.c, Ucell->nat * (2 * rho_basis->nx + 1));
            castmem_z2c_h2d_op()(this->c_eigts2, this->eigts2.c, Ucell->nat * (2 * rho_basis->ny + 1));
            castmem_z2c_h2d_op()(this->c_eigts3, this->eigts3.c, Ucell->nat * (2 * rho_basis->nz + 1));
        }
        resmem_zd_op()(this->z_eigts1, Ucell->nat * (2 * rho_basis->nx + 1));
        resmem_zd_op()(this->z_eigts2, Ucell->nat * (2 * rho_basis->ny + 1));
        resmem_zd_op()(this->z_eigts3, Ucell->nat * (2 * rho_basis->nz + 1));
        syncmem_z2z_h2d_op()(this->z_eigts1, this->eigts1.c, Ucell->nat * (2 * rho_basis->nx + 1));
        syncmem_z2z_h2d_op()(this->z_eigts2, this->eigts2.c, Ucell->nat * (2 * rho_basis->ny + 1));
        syncmem_z2z_h2d_op()(this->z_eigts3, this->eigts3.c, Ucell->nat * (2 * rho_basis->nz + 1));
    }
    else {
        if (PARAM.globalv.has_float_data) {
            resmem_ch_op()(this->c_eigts1, Ucell->nat * (2 * rho_basis->nx + 1));
            resmem_ch_op()(this->c_eigts2, Ucell->nat * (2 * rho_basis->ny + 1));
            resmem_ch_op()(this->c_eigts3, Ucell->nat * (2 * rho_basis->nz + 1));
            castmem_z2c_h2h_op()(this->c_eigts1, this->eigts1.c, Ucell->nat * (2 * rho_basis->nx + 1));
            castmem_z2c_h2h_op()(this->c_eigts2, this->eigts2.c, Ucell->nat * (2 * rho_basis->ny + 1));
            castmem_z2c_h2h_op()(this->c_eigts3, this->eigts3.c, Ucell->nat * (2 * rho_basis->nz + 1));
        }
        this->z_eigts1 = this->eigts1.c;
        this->z_eigts2 = this->eigts2.c;
        this->z_eigts3 = this->eigts3.c;
        // There's no need to delete double precision pointers while in a CPU environment.
    }
    ModuleBase::timer::end("Structure_Factor","setup");
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
//
void Structure_Factor::bspline_sf(const int norder,
                                  const UnitCell* Ucell,
                                  const Parallel_Grid& pgrid,
                                  const ModulePW::PW_Basis* rho_basis)
{
    (void)pgrid;
    std::vector<double> tmpr(rho_basis->nrxx);
    std::vector<std::complex<double>> b1(rho_basis->nx);
    std::vector<std::complex<double>> b2(rho_basis->ny);
    std::vector<std::complex<double>> b3(rho_basis->nz);
    const int nplane = rho_basis->nplane;
    const int startz = rho_basis->startz_current;

    // Each rank owns the same atoms; populate only its local FFT z slab.
    for (int it = 0; it < Ucell->ntype; it++)
    {
        const int na = Ucell->atoms[it].na;
        const ModuleBase::Vector3<double>* const taud = Ucell->atoms[it].taud.data();
        ModuleBase::GlobalFunc::ZEROS(tmpr.data(), rho_basis->nrxx);

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ia = 0; ia < na; ++ia)
        {
            const double gridx = taud[ia].x * rho_basis->nx;
            const double gridy = taud[ia].y * rho_basis->ny;
            const double gridz = taud[ia].z * rho_basis->nz;
            const double dx = gridx - floor(gridx);
            const double dy = gridy - floor(gridy);
            const double dz = gridz - floor(gridz);
            //I'm not sure if there is a mod function for double data

            ModuleBase::Bspline bsx, bsy, bsz;
            bsx.init(norder, 1, 0);
            bsy.init(norder, 1, 0);
            bsz.init(norder, 1, 0);
            bsx.getbspline(dx);
            bsy.getbspline(dy);
            bsz.getbspline(dz);

            for (int iz = 0; iz <= norder; ++iz)
            {
                const int icz = int(rho_basis->nz * 10 - iz + floor(gridz)) % rho_basis->nz;
                if (icz < startz || icz >= startz + nplane)
                {
                    continue;
                }
                const int local_z = icz - startz;
                for (int iy = 0; iy <= norder; ++iy)
                {
                    const int icy = int(rho_basis->ny * 10 - iy + floor(gridy)) % rho_basis->ny;
                    for (int ix = 0; ix <= norder; ++ix)
                    {
                        const int icx = int(rho_basis->nx * 10 - ix + floor(gridx)) % rho_basis->nx;
#ifdef _OPENMP
#pragma omp atomic
#endif
                        tmpr[(icx * rho_basis->ny + icy) * nplane + local_z]
                            += bsz.bezier_ele(iz) * bsy.bezier_ele(iy) * bsx.bezier_ele(ix);
                    }
                }
            }
        }

        //It should be optimized with r2c
        rho_basis->real2recip(tmpr.data(), &strucFac(it, 0));
        this->bsplinecoef(b1.data(),
                          b2.data(),
                          b3.data(),
                          rho_basis->nx,
                          rho_basis->ny,
                          rho_basis->nz,
                          norder);
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 128)
#endif
        for (int ig = 0; ig < rho_basis->npw; ++ig)
        {
            const int idx = int(rho_basis->gdirect[ig].x + 0.1 + rho_basis->nx) % rho_basis->nx;
            const int idy = int(rho_basis->gdirect[ig].y + 0.1 + rho_basis->ny) % rho_basis->ny;
            const int idz = int(rho_basis->gdirect[ig].z + 0.1 + rho_basis->nz) % rho_basis->nz;
            strucFac(it, ig) *= (b1[idx] * b2[idy] * b3[idz] * double(rho_basis->nxyz));
        }
    }

    return;
}

void Structure_Factor::bsplinecoef(std::complex<double> *b1, std::complex<double> *b2, std::complex<double> *b3, 
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
