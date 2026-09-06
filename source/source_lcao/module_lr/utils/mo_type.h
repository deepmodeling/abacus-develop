#pragma once
#include <stdexcept>
namespace LR_Util
{
#ifndef MO_TYPE_H
#define MO_TYPE_H
    enum MO_TYPE { OO, VO, OV, VV, ALL };
#endif
/* the first index is contiguous in memory
MO_TYPE: OO   VO    OV    VV    ALL
nmo1     nocc nocc  nvirt nvirt nocc+nvirt
nmo2     nocc nvirt nocc  nvirt nocc+nvirt
imo1     0    0     nocc  nocc  0
imo2     0    nocc  0     nocc  0
*/
inline void set_dim(const MO_TYPE type, const int& nocc, const int& nvirt,
    int& nmo1, int& nmo2, int& imo1, int& imo2)
{
    switch(type)
    {
    case MO_TYPE::OO:
        nmo1 = nocc; nmo2 = nocc; imo1 = 0; imo2 = 0;
        break;
    case MO_TYPE::VO:
        nmo1 = nocc; nmo2 = nvirt; imo1 = 0; imo2 = nocc;
        break;
    case MO_TYPE::OV:
        nmo1 = nvirt; nmo2 = nocc; imo1 = nocc; imo2 = 0;
        break;
    case MO_TYPE::VV:
        nmo1 = nvirt; nmo2 = nvirt; imo1 = nocc; imo2 = nocc;
        break;
    case MO_TYPE::ALL:
        nmo1 = nocc + nvirt;
        nmo2 = nocc + nvirt;
        imo1 = 0;
        imo2 = 0;
        break;
    default:
        throw std::runtime_error("Error in LR::set_dim: unknown MO_TYPE");
    }
}
}