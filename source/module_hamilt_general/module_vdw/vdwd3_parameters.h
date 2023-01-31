//==========================================================
// AUTHOR : Yuyang Ji
// DATE : 2019-04-22
// UPDATE : 2021-4-19
//==========================================================

#ifndef VDWD3_PARAMETERS_H
#define VDWD3_PARAMETERS_H

#include "module_io/input.h"
#include "vdw_parameters.h"

namespace vdw
{

class Vdwd3Parameters : public VdwParameters
{

  public:
    Vdwd3Parameters() : VdwParameters() {};

    ~Vdwd3Parameters() = default;

    void initial_parameters(const Input &input);

    inline const std::string &version() const { return version_; }

    inline bool abc() const { return abc_; }
    inline double rthr2() const { return rthr2_; }
    inline double cn_thr2() const { return cn_thr2_; }
    inline double s6() const { return s6_; }
    inline double rs6() const { return rs6_; }
    inline double s18() const { return s18_; }
    inline double rs18() const { return rs18_; }

    inline const std::vector<int> &mxc() const { return mxc_; }
    inline const std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> &c6ab() const { return c6ab_; }
    inline const std::vector<double> &r2r4() const { return r2r4_; }
    inline const std::vector<double> &rcov() { return rcov_; }
    inline const std::vector<std::vector<double>> &r0ab() { return r0ab_; }

    inline double k1() const { return k1_; }
    inline double k2() const { return k2_; }
    inline double k3() const { return k3_; }
    inline double alp6() const { return alp6_; }
    inline double alp8() const { return alp8_; }
    inline double alp10() const { return alp10_; }

  private:
    std::string version_;

    bool abc_; // third-order term?
    double rthr2_; // R^2 distance neglect threshold (important for speed in case of large systems) (a.u.)
    double cn_thr2_; // R^2 distance to cutoff for CN_calculation (a.u.)
    double s6_;
    double rs6_;
    double s18_;
    double rs18_;

    static constexpr size_t max_elem_ = 94;
    static constexpr double k1_ = 16.0, k2_ = 4.0 / 3.0, k3_ = -4.0;
    static constexpr double alp6_ = 14.0, alp8_ = alp6_ + 2, alp10_ = alp8_ + 2;

    std::vector<int> mxc_;
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> c6ab_;
    std::vector<double> r2r4_;
    std::vector<double> rcov_;
    std::vector<std::vector<double>> r0ab_;

    void init_C6();
    void init_r2r4();
    void init_rcov();
    void init_r0ab();

    int limit(int &i);
};

} // namespace vdw

#endif // VDWD3_PARAMETERS_H
