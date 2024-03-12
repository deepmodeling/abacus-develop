//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
// UPDATE : 2019-04-26
//==========================================================

#ifndef VDWD2_H
#define VDWD2_H

#include "vdw.h"

namespace vdw
{

class Vdwd2 : public Vdw
{

  public:
    Vdwd2(const UnitCell &unit_in) : Vdw(unit_in) {}

    ~Vdwd2() = default;

    Vdwd2Parameters &parameter() { return para_; }
    const Vdwd2Parameters &parameter() const { return para_; }

  private:
    Vdwd2Parameters para_;

    void cal_energy() override;
    void cal_force() override;
    void cal_stress() override;

    template <typename F> void index_loops(F &&f)
    {
        int xidx = para_.period().x / 2;
        int yidx = para_.period().y / 2;
        int zidx = para_.period().z / 2;

        for (int it1 = 0; it1 != ucell_.ntype; ++it1)
        {
            for (int it2 = 0; it2 != ucell_.ntype; ++it2)
            {
                const double C6_product
                    = sqrt(para_.C6().at(ucell_.atoms[it1].ncpp.psd) * para_.C6().at(ucell_.atoms[it2].ncpp.psd))
                      / pow(ucell_.lat0, 6);
                const double R0_sum
                    = (para_.R0().at(ucell_.atoms[it1].ncpp.psd) + para_.R0().at(ucell_.atoms[it2].ncpp.psd)) / ucell_.lat0;
                if (!R0_sum)
                {
                    ModuleBase::WARNING_QUIT("Input", "R0_sum can not be 0");
                }
                for (int ia1 = 0; ia1 != ucell_.atoms[it1].na; ++ia1)
                {
                    const ModuleBase::Vector3<double> tau1 = ucell_.atoms[it1].tau[ia1];
                    for (int ia2 = 0; ia2 != ucell_.atoms[it2].na; ++ia2)
                    {
                        ModuleBase::Vector3<int> ilat_loop;
                        for (ilat_loop.x = -xidx; ilat_loop.x <= xidx; ++ilat_loop.x)
                        {
                            for (ilat_loop.y = -yidx; ilat_loop.y <= yidx; ++ilat_loop.y)
                            {
                                for (ilat_loop.z = -zidx; ilat_loop.z <= zidx; ++ilat_loop.z)
                                {
                                    if ((!(ilat_loop.x || ilat_loop.y || ilat_loop.z)) && (it1 == it2) && (ia1 == ia2))
                                        continue;
                                    const ModuleBase::Vector3<double> tau2
                                        = ucell_.atoms[it2].tau[ia2] + ilat_loop * ucell_.latvec;
                                    const double r_sqr = (tau1 - tau2).norm2();
                                    const double r = sqrt(r_sqr);
                                    // calculations happen in f
                                    f(r, R0_sum, C6_product, r_sqr, it1, ia1, tau1, tau2);
                                }
                            }
                        } // end for ilat_loop
                    } // end for ia2
                } // end for ia1
            } // end for it2
        } // end for it1
    }
};

} // namespace vdw

#endif // VDWD2_H
