#include "NCLibxc.h"
#include "interface_to_libxc.h"
#include <vector>
#include <xc.h>

namespace NCXC{

///////////////////////////////////////////////////////////////////////////////////lda

// post-processing of libxc. get partial derivatives from libxc and integrate them to get the 0th,1st,2nd derivatives of the functional
//Note that e is the exchange and correlation energy per electron per volume. You need to multiply by \rho before the integration.
/*
v1=0, v2=1, f1=0,0, f2=0,1, f3=1,1
*/

void NCLibxc::postlibxc_lda(int xc_id, const std::vector<double>& rho_up, const std::vector<double>& rho_down, 
                   std::vector<double>& e, std::vector<double>& v1, std::vector<double>& v2, 
                   std::vector<double>& f1, std::vector<double>& f2, std::vector<double>& f3)
{
    LibxcInterface libxc(xc_id, true); // xc_id now passed from the caller

    std::vector<double> exc = libxc.lda_exc(rho_up, rho_down);
    std::vector<double> vrho_1(rho_up.size()), vrho_2(rho_down.size());
    libxc.lda_vxc(rho_up, rho_down, vrho_1, vrho_2);

    std::vector<double> v2rho2_1(rho_up.size()), v2rho2_2(rho_down.size()), v2rho2_3(rho_up.size());
    libxc.lda_fxc(rho_up, rho_down, v2rho2_1, v2rho2_2, v2rho2_3);

    
    e.resize(rho_up.size());
    v1.resize(rho_up.size());
    v2.resize(rho_down.size());
    f1.resize(rho_up.size());
    f2.resize(rho_down.size());
    f3.resize(rho_up.size());

    for (size_t i = 0; i < rho_up.size(); ++i) {
        e[i] = exc[i] ; 
        v1[i] = vrho_1[i] ; 
        v2[i] = vrho_2[i] ; 
        f1[i] = v2rho2_1[i] ; 
        f2[i] = v2rho2_2[i] ; 
        f3[i] = v2rho2_3[i] ; 
    }
}


}// end of namespace NCXC