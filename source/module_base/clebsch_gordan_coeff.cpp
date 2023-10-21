#include "clebsch_gordan_coeff.h"

#include "module_base/constants.h"
#include "module_base/inverse_matrix.h"
#include "module_base/math_ylmreal.h"

namespace ModuleBase
{

Clebsch_Gordan::Clebsch_Gordan()
{
}
Clebsch_Gordan::~Clebsch_Gordan()
{
}

void Clebsch_Gordan::clebsch_gordan(const int& lli,
                                    ModuleBase::realArray& ap,
                                    ModuleBase::IntArray& lpx,
                                    ModuleBase::IntArray& lpl)
{
    if (lli < 0)
    {
        std::cout << "Clebsch_Gordan: lmaxkb + 1 < 0" << std::endl;
        exit(1);
    }

    const int llx = (2 * lli - 1) * (2 * lli - 1);
    ModuleBase::Vector3<double>* r = new ModuleBase::Vector3<double>[llx];
    ModuleBase::matrix ylm(llx, llx);
    ModuleBase::matrix mly(llx, llx);

    // generate an array of random vectors (uniform deviate on unitary sphere)
    gen_rndm_r(llx, r);

    // generate the real spherical harmonics for the array: ylm(ir,lm)
    ModuleBase::YlmReal::Ylm_Real(llx, llx, r, ylm);

    // store the inverse of ylm(ir,lm) in mly(lm,ir)
    ModuleBase::Inverse_Matrix_Real(llx, ylm.c, mly.c);

    // for each li,lj compute ap(l,li,lj) and the indices, lpx and lpl
    lpx.create(lli * lli, lli * lli);
    lpl.create(lli * lli, lli * lli, llx);
    ap.create(llx, lli * lli, lli * lli);
    for (int li = 0; li < lli * lli; li++)
    {
        for (int lj = 0; lj < lli * lli; lj++)
        {
            lpx(li, lj) = 0;
            for (int L = 0; L < llx; L++)
            {
                ap(L, li, lj) = compute_ap(L, li, lj, llx, ylm, mly);
                if (std::abs(ap(L, li, lj)) > 1.0e-3)
                {
                    lpl(li, lj, lpx(li, lj)) = L;
                    lpx(li, lj)++;
                }
            }
        }
    }

    delete[] r;
}

void Clebsch_Gordan::gen_rndm_r(const int& llx, ModuleBase::Vector3<double>* r)
{
    for (int i = 0; i < llx; i++)
    {
        double costheta = 2.0 * static_cast<double>(std::rand()) / RAND_MAX - 1.0;
        double sintheta = std::sqrt(1.0 - costheta * costheta);
        double phi = ModuleBase::TWO_PI * static_cast<double>(std::rand()) / RAND_MAX;
        r[i].x = sintheta * std::cos(phi);
        r[i].y = sintheta * std::sin(phi);
        r[i].z = costheta;
    }
}

double Clebsch_Gordan::compute_ap(const int& L,
                                  const int& li,
                                  const int& lj,
                                  const int& llx,
                                  const ModuleBase::matrix& ylm,
                                  const ModuleBase::matrix& mly)
{
    double compute_ap = 0.0;
    for (int ir = 0; ir < llx; ir++)
    {
        compute_ap += mly(ir, L) * ylm(li, ir) * ylm(lj, ir);
    }
    return compute_ap;
}

} // namespace ModuleBase