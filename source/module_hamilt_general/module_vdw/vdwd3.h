//==========================================================
// AUTHOR : Yuyang Ji
// DATE : 2019-04-22
// UPDATE : 2021-4-19
//==========================================================

#ifndef VDWD3_H
#define VDWD3_H

#include "vdw.h"

namespace vdw
{

class Vdwd3 : public Vdw
{

  public:
    Vdwd3(const UnitCell &unit_in) : Vdw(unit_in) { }

    ~Vdwd3() = default;

    Vdwd3Parameters &parameter() { return para_; }
    const Vdwd3Parameters &parameter() const { return para_; }

  private:
    Vdwd3Parameters para_;

    std::vector<ModuleBase::Vector3<double>> lat_;
    std::vector<int> iz_;
    std::vector<ModuleBase::Vector3<double>> xyz_;
    std::vector<int> rep_vdw_;
    std::vector<int> rep_cn_;

    void cal_energy() override;
    void cal_force() override;
    void cal_stress() override;

    void init();

    void set_criteria(double rthr, const std::vector<ModuleBase::Vector3<double>> &lat, std::vector<double> &tau_max);

    std::vector<double> atom_kind();

    void get_c6(int iat, int jat, double nci, double ncj, double &c6);

    void pbc_ncoord(std::vector<double> &cn);

    void pbc_three_body(const std::vector<int> &iz,
                        const std::vector<ModuleBase::Vector3<double>> &lat,
                        const std::vector<ModuleBase::Vector3<double>> &xyz,
                        const std::vector<int> &rep_cn,
                        const std::vector<double> &cc6ab,
                        double &eabc);

    void pbc_gdisp(std::vector<ModuleBase::Vector3<double>> &g, ModuleBase::matrix &smearing_sigma);

    void get_dc6_dcnij(int mxci, int mxcj, double cni, double cnj, int izi, int izj, int iat, int jat,
                       double &c6check, double &dc6i, double &dc6j);

    int lin(int i1, int i2)
    {
        int idum1 = std::max(i1 + 1, i2 + 1);
        int idum2 = std::min(i1 + 1, i2 + 1);
        int res = idum2 + idum1 * (idum1 - 1) / 2 - 1;
        return res;
    }
};

} // namespace vdw

#endif // VDWD3_H
