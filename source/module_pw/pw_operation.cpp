#include "pw_basis_k.h"
#include <assert.h>
namespace ModulePW
{

ModuleBase::Vector3<double> PW_Basis_K:: get_GPlusK_cartesian(const int ik, const int ig) const {
    assert(ig>=0 && ig<this->npw && ik>=0 && ik<this->nks);
    int isz = this->ig2isz[ig];
    int iz = isz % this->nz;
    int is = isz / this->nz;
    int ix = this->is2ixy[is] / this->ny;
    int iy = this->is2ixy[is] % this->ny;
    if (ix >= int(this->nx/2) + 1) ix -= this->nx;
    if (iy >= int(this->bigny/2) + 1) iy -= this->bigny;
    if (iz >= int(this->nz/2) + 1) iz -= this->nz;
    ModuleBase::Vector3<double> f;
    f.x = ix;
    f.y = iy;
    f.z = iz;
    f = f * this->G;
    ModuleBase::Vector3<double> g_temp_ = this->kvec_c[ik] + f;
    return g_temp_;
}


double PW_Basis_K::get_GPlusK_cartesian_projection(const int ik, const int ig, const int axis) const
{
    assert(ig>=0 && ig<this->npw && ik>=0 && ik<this->nks);
    int isz = this->ig2isz[ig];
    int iz = isz % this->nz;
    int is = isz / this->nz;
    int ix = this->is2ixy[is] / this->ny;
    int iy = this->is2ixy[is] % this->ny;
    if (ix >= int(this->nx/2) + 1) ix -= this->nx;
    if (iy >= int(this->bigny/2) + 1) iy -= this->bigny;
    if (iz >= int(this->nz/2) + 1) iz -= this->nz;
    ModuleBase::Vector3<double> f;
    f.x = ix;
    f.y = iy;
    f.z = iz;
    f = f * this->G;
    ModuleBase::Vector3<double> g_temp_ = this->kvec_c[ik] + f;
    if (axis == 0)
    {
        return g_temp_.x;
    }
    else if (axis == 1)
    {
        return g_temp_.y;
    }
    else if (axis == 2)
    {
        return g_temp_.z;
    }
    return 0.0;
}

double PW_Basis_K::get_SquareGPlusK_cartesian(const int ik, const int ig) const 
{
    assert(ig>=0 && ig<this->npw && ik>=0 && ik<this->nks);
    int isz = this->ig2isz[ig];
    int iz = isz % this->nz;
    int is = isz / this->nz;
    int ix = this->is2ixy[is] / this->ny;
    int iy = this->is2ixy[is] % this->ny;
    if (ix >= int(this->nx/2) + 1) ix -= this->nx;
    if (iy >= int(this->bigny/2) + 1) iy -= this->bigny;
    if (iz >= int(this->nz/2) + 1) iz -= this->nz;
    ModuleBase::Vector3<double> f;
    f.x = ix;
    f.y = iy;
    f.z = iz;
    f = f * this->G;
    ModuleBase::Vector3<double> g_temp_ = this->kvec_c[ik] + f;
    return (g_temp_ * g_temp_);
}

}