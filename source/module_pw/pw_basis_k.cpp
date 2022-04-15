#include "pw_basis_k.h"
#include "../module_base/constants.h"
namespace ModulePW
{

PW_Basis_K::PW_Basis_K()
{
    nks = 1;
    kvec_d = NULL;
    kvec_c = NULL;
    npwk = NULL;
    igl2isz_k = NULL;
}
PW_Basis_K::~PW_Basis_K()
{
    if(kvec_d != NULL) delete[] kvec_d;
    if(kvec_c != NULL) delete[] kvec_c;
    if(npwk != NULL) delete[] npwk;
    if(igl2isz_k != NULL) delete[] igl2isz_k;
}

void PW_Basis_K:: initparameters(
    bool gamma_only_in,
    double gk_ecut_in,
    int nks_in, //number of k points in this pool
    ModuleBase::Vector3<double> *kvec_d_in, // Direct coordinates of k points
    int poolnproc_in, // Number of processors in this pool
    int poolrank_in, // Rank in this pool
    int distribution_type_in
)
{
    this->nks = nks_in;
    this->kvec_d = new ModuleBase::Vector3<double> [nks];
    this->kvec_c = new ModuleBase::Vector3<double> [nks];

    double kmaxmod = 0;
    for(int ik = 0 ; ik < this->nks ; ++ik)
    {
        this->kvec_d[ik] = kvec_d_in[ik];
        this->kvec_c[ik] = this->kvec_d[ik] * this->G;
        double kmod = sqrt(this->kvec_c[ik] * this->kvec_c[ik]);
        if(kmod > kmaxmod)  kmaxmod = kmod;
    }
    // MPI_Allreduce(MPI_IN_PLACE, &kmaxmod, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / this->lat0 / this->lat0;
    this->gk_ecut = gk_ecut_in/tpiba2;
    this->ggecut = pow(sqrt(this->gk_ecut) + kmaxmod, 2);

    this->gamma_only = gamma_only_in;
    if(kmaxmod > 0)     this->gamma_only = false; //if it is not the gamma point, we do not use gamma_only
    if (this->gamma_only)   this->ny = int(this->bigny / 2) + 1;
    else                    this->ny = bigny;
    this->nxy = this->nx * this->ny;
    this->nxyz = this->nxy * this->nz;

    this->poolnproc = poolnproc_in;
    this->poolrank = poolrank_in;
    this->distribution_type = distribution_type_in;
    return;
}

void PW_Basis_K::setupIndGk()
{
    //count npwk
    this->npwk_max = 0;
    this->npwk = new int [this->nks];
    for (int ik = 0; ik < this->nks; ik++)
    {
        int ng = 0;
        for (int ig = 0; ig < npw ; ig++)
        {
            const double gk2 = this->get_GPlusK_cartesian(ik, ig).norm2();       
            if (gk2 <= this->gk_ecut)
            {
                ++ng;
            }
        }
        this->npwk[ik] = ng;
        if ( this->npwk_max < ng)
        {
            this->npwk_max = ng;
        }
    }

    //get igl2isz_k
    this->igl2isz_k = new int [this->nks * this->npwk_max];
    for (int ik = 0; ik < this->nks; ik++)
    {
        int igl = 0;
        for (int ig = 0; ig < npw ; ig++)
        {
            const double gk2 = this->get_GPlusK_cartesian(ik, ig).norm2();       
            if (gk2 <= this->gk_ecut)
            {
                this->igl2isz_k[ik*npwk_max + igl] = this->ig2isz[ig];
                ++igl;
            }
        }
    }

    delete[] this->ig2isz;
    this->ig2isz = NULL;

    return;
}
void PW_Basis_K::setuptransform()
{
    this->distribute_r();
    this->distribute_g();
    this->getstartgr();
    this->setupIndGk();
    this->ft.initfft(this->nx,this->bigny,this->nz,this->liy,this->riy,this->nst,this->nplane,this->poolnproc,this->gamma_only);
    this->ft.setupFFT();
}

void PW_Basis_K::collect_local_pw()
{
    if(gg != NULL) delete[] gg;
    if(gdirect != NULL) delete[] gdirect;
    if(gcar != NULL) delete[] gcar;
    this->gg = new double[this->npwk_max * this->nks];
    this->gdirect = new ModuleBase::Vector3<double>[this->npwk_max * this->nks];
    this->gcar = new ModuleBase::Vector3<double>[this->npwk_max * this->nks];

    ModuleBase::Vector3<double> f;
    for(int ik = 0 ; ik < this->nks ; ++ik)
    {
        for(int igl = 0 ; igl < this-> npwk[ik] ; ++igl)
        {
            int isz = this->igl2isz_k[ik * npwk_max + igl];
            int iz = isz % this->nz;
            int is = isz / this->nz;
            int ixy = this->is2ixy[is];
            int ix = ixy / this->ny;
            int iy = ixy % this->ny;
            if (ix >= int(this->nx/2) + 1) ix -= this->nx;
            if (iy >= int(this->bigny/2) + 1) iy -= this->bigny;
            if (iz >= int(this->nz/2) + 1) iz -= this->nz;
            f.x = ix;
            f.y = iy;
            f.z = iz;
            this->gg[ik * npwk_max + igl] = f * (this->GGT * f);
            this->gdirect[ik * npwk_max + igl] = f;
            this->gcar[ik * npwk_max + igl] = f * this->G;
        }
    }
}





}