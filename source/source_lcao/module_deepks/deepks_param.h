#ifndef LCAO_DEEPKS_PARAM
#define LCAO_DEEPKS_PARAM

#include <vector>
namespace ModuleBase
{
struct IntArray;
}

struct DeePKS_Param
{
    int lmaxd = 0;
    int nmaxd = 0;
    int inlmax = 0;
    int n_descriptor = 0;
    int des_per_atom = 0;
    std::vector<int> inl2l;
    ModuleBase::IntArray* inl_index = nullptr;
};

#endif