#include "grid_meshk.h"

#include "module_hamilt_pw/hamilt_pwdft/global.h"

Grid_MeshK::Grid_MeshK()
{
    ucell_index2x = nullptr;
    ucell_index2y = nullptr;
    ucell_index2z = nullptr;
}

Grid_MeshK::~Grid_MeshK()
{
    delete[] ucell_index2x;
    delete[] ucell_index2y;
    delete[] ucell_index2z;
}

int Grid_MeshK::cal_Rindex(const int& u1, const int& u2, const int& u3) const
{
    const int x1 = u1 - this->minu1;
    const int x2 = u2 - this->minu2;
    const int x3 = u3 - this->minu3;

    if (x1 < 0 || x2 < 0 || x3 < 0)
    {
        std::cout << " u1=" << u1 << " minu1=" << minu1 << std::endl;
        std::cout << " u2=" << u2 << " minu2=" << minu2 << std::endl;
        std::cout << " u3=" << u3 << " minu3=" << minu3 << std::endl;
        ModuleBase::WARNING_QUIT("Grid_MeshK::cal_Rindex",
                                 "x1<0 || x2<0 || x3<0 !");
    }

    assert(x1 >= 0);
    assert(x2 >= 0);
    assert(x3 >= 0);
    return (x3 + x2 * this->nu3 + x1 * this->nu2 * this->nu3);
}
