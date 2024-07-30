#include <sstream>
#include <cassert>
#include <cmath>
#include <numeric>
#include "module_hamilt_pw/hamilt_pwdft/radial_projection.h"
#include "module_base/constants.h"

void RadialProjection::RadialProjector::sbtfft(const int nr, 
                                               const double* r, 
                                               const double* in, 
                                               const int l, 
                                               std::complex<double>* out) const
{
    /* user should take care of the memory allocation by own, the size required
    would be (2*l+1)*npw */
    const int npw = qnorm_.size();
    std::vector<double> jlq(npw);
    sbt_->direct(l, nr, r, in, npw, qnorm_.data(), jlq.data());
    for(int m = -l; m <= l; m++)
    {
        int lm = l*l + m;
        for(int iq = 0; iq < npw; iq++)
        {
            out[m+l+lm*npw] = ModuleBase::FOUR_PI/std::sqrt(omega_) * std::pow(ModuleBase::IMAG_UNIT, l) * jlq[iq] * ylm_(lm, iq);
        }
    }
}

void RadialProjection::RadialProjector::sbtfft(const std::vector<double>& r,
                                               const std::vector<double>& in,
                                               const int l,
                                               std::vector<std::complex<double>>& out) const
{
    const int nr = r.size();
    const int npw = qnorm_.size();
    out.resize((2*l+1)*npw);
    sbtfft(nr, r.data(), in.data(), l, out.data()); 
}

void RadialProjection::_radial_indexing(const int ntype,
                                        const std::vector<int>& lmax,
                                        const std::vector<std::vector<int>>& nzeta,
                                        std::map<std::tuple<int, int, int>, int>& map,
                                        std::vector<std::tuple<int, int, int>>& rmap)
{
    int iproj = 0;
    for(int it = 0; it < ntype; it++)
    {
        for(int l = 0; l <= lmax[it]; l++)
        {
            for(int zeta = 0; zeta < nzeta[it][l]; zeta++)
            {
                map[std::make_tuple(it, l, zeta)] = iproj;
                rmap.push_back(std::make_tuple(it, l, zeta));
                iproj++;
            }
        }
    }
}

void RadialProjection::_mask_func(std::vector<double>& mask)
{
    /* mask function is hard coded here, eta = 15 */
    mask.resize(201);
    std::string src;
    src += " 0.10000000E+01 0.99948662E+00 0.99863154E+00 0.99743557E+00";
    src += " 0.99589985E+00 0.99402586E+00 0.99181538E+00 0.98927052E+00";
    src += " 0.98639370E+00 0.98318766E+00 0.97965544E+00 0.97580040E+00";
    src += " 0.97162618E+00 0.96713671E+00 0.96233623E+00 0.95722924E+00";
    src += " 0.95182053E+00 0.94611516E+00 0.94011842E+00 0.93383589E+00";
    src += " 0.92727338E+00 0.92043693E+00 0.91333282E+00 0.90596753E+00";
    src += " 0.89834777E+00 0.89048044E+00 0.88237263E+00 0.87403161E+00";
    src += " 0.86546483E+00 0.85667987E+00 0.84768450E+00 0.83848659E+00";
    src += " 0.82909416E+00 0.81951535E+00 0.80975838E+00 0.79983160E+00";
    src += " 0.78974340E+00 0.77950227E+00 0.76911677E+00 0.75859548E+00";
    src += " 0.74794703E+00 0.73718009E+00 0.72630334E+00 0.71532544E+00";
    src += " 0.70425508E+00 0.69310092E+00 0.68187158E+00 0.67057566E+00";
    src += " 0.65922170E+00 0.64781819E+00 0.63637355E+00 0.62489612E+00";
    src += " 0.61339415E+00 0.60187581E+00 0.59034914E+00 0.57882208E+00";
    src += " 0.56730245E+00 0.55579794E+00 0.54431609E+00 0.53286431E+00";
    src += " 0.52144984E+00 0.51007978E+00 0.49876105E+00 0.48750040E+00";
    src += " 0.47630440E+00 0.46517945E+00 0.45413176E+00 0.44316732E+00";
    src += " 0.43229196E+00 0.42151128E+00 0.41083069E+00 0.40025539E+00";
    src += " 0.38979038E+00 0.37944042E+00 0.36921008E+00 0.35910371E+00";
    src += " 0.34912542E+00 0.33927912E+00 0.32956851E+00 0.31999705E+00";
    src += " 0.31056799E+00 0.30128436E+00 0.29214897E+00 0.28316441E+00";
    src += " 0.27433307E+00 0.26565709E+00 0.25713844E+00 0.24877886E+00";
    src += " 0.24057988E+00 0.23254283E+00 0.22466884E+00 0.21695884E+00";
    src += " 0.20941357E+00 0.20203357E+00 0.19481920E+00 0.18777065E+00";
    src += " 0.18088790E+00 0.17417080E+00 0.16761900E+00 0.16123200E+00";
    src += " 0.15500913E+00 0.14894959E+00 0.14305240E+00 0.13731647E+00";
    src += " 0.13174055E+00 0.12632327E+00 0.12106315E+00 0.11595855E+00";
    src += " 0.11100775E+00 0.10620891E+00 0.10156010E+00 0.97059268E-01";
    src += " 0.92704295E-01 0.88492966E-01 0.84422989E-01 0.80492001E-01";
    src += " 0.76697569E-01 0.73037197E-01 0.69508335E-01 0.66108380E-01";
    src += " 0.62834685E-01 0.59684561E-01 0.56655284E-01 0.53744102E-01";
    src += " 0.50948236E-01 0.48264886E-01 0.45691239E-01 0.43224469E-01";
    src += " 0.40861744E-01 0.38600231E-01 0.36437098E-01 0.34369520E-01";
    src += " 0.32394681E-01 0.30509780E-01 0.28712032E-01 0.26998673E-01";
    src += " 0.25366964E-01 0.23814193E-01 0.22337676E-01 0.20934765E-01";
    src += " 0.19602844E-01 0.18339338E-01 0.17141711E-01 0.16007467E-01";
    src += " 0.14934157E-01 0.13919377E-01 0.12960772E-01 0.12056034E-01";
    src += " 0.11202905E-01 0.10399183E-01 0.96427132E-02 0.89313983E-02";
    src += " 0.82631938E-02 0.76361106E-02 0.70482151E-02 0.64976294E-02";
    src += " 0.59825322E-02 0.55011581E-02 0.50517982E-02 0.46327998E-02";
    src += " 0.42425662E-02 0.38795566E-02 0.35422853E-02 0.32293218E-02";
    src += " 0.29392897E-02 0.26708663E-02 0.24227820E-02 0.21938194E-02";
    src += " 0.19828122E-02 0.17886449E-02 0.16102512E-02 0.14466132E-02";
    src += " 0.12967606E-02 0.11597692E-02 0.10347601E-02 0.92089812E-03";
    src += " 0.81739110E-03 0.72348823E-03 0.63847906E-03 0.56169212E-03";
    src += " 0.49249371E-03 0.43028657E-03 0.37450862E-03 0.32463165E-03";
    src += " 0.28016004E-03 0.24062948E-03 0.20560566E-03 0.17468305E-03";
    src += " 0.14748362E-03 0.12365560E-03 0.10287226E-03 0.84830727E-04";
    src += " 0.69250769E-04 0.55873673E-04 0.44461100E-04 0.34793983E-04";
    src += " 0.26671449E-04 0.19909778E-04 0.14341381E-04 0.98138215E-05";
    std::stringstream ss(src);
    for(int i = 0; i < mask.size(); i++)
    {
        ss >> mask[i];
    }
}

void RadialProjection::_do_mask_on_radial(const int nr1,
                                          const double* r,
                                          const double* in,
                                          const int nr2,
                                          const double* mask,
                                          double* out)
{
    /* the key here is to avoid any float-point overflow */
}