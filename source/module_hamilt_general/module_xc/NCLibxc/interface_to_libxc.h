// Programmed by Xiaoyu Zhang at Peking University, Beijing, China 2024/08/26

#ifndef INTERFACE_TO_LIBXC_H
#define INTERFACE_TO_LIBXC_H

#include <xc.h>
#include <vector>
#include <stdexcept>

class LibxcInterface
{
private:
    xc_func_type func;

public:
    // Constructor with spin polarization
    LibxcInterface(int xc_id, bool spin_polarized = true);

    // Destructor
    ~LibxcInterface();

    std::vector<double> lda_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down);
    void lda_vxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &vrho_1, std::vector<double> &vrho_2);
    void lda_fxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v2rho2_1, std::vector<double> &v2rho2_2, std::vector<double> &v2rho2_3);
    void lda_kxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v3rho3_1, std::vector<double> &v3rho3_2, std::vector<double> &v3rho3_3, std::vector<double> &v3rho3_4);
    void lda_lxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v4rho4_1, std::vector<double> &v4rho4_2, std::vector<double> &v4rho4_3, std::vector<double> &v4rho4_4, std::vector<double> &v4rho4_5);

    std::vector<double> gga_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down, const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3);
    void gga_vxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3, std::vector<double> &vrho_1, std::vector<double> &vrho_2,
                 std::vector<double> &vsigma_1, std::vector<double> &vsigma_2, std::vector<double> &vsigma_3);
    void gga_fxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
                 std::vector<double> &v2rho2_1, std::vector<double> &v2rho2_2, std::vector<double> &v2rho2_3,
                 std::vector<double> &v2rhosigma_1, std::vector<double> &v2rhosigma_2, std::vector<double> &v2rhosigma_3, std::vector<double> &v2rhosigma_4, std::vector<double> &v2rhosigma_5, std::vector<double> &v2rhosigma_6, std::vector<double> &v2sigma2_1, std::vector<double> &v2sigma2_2, std::vector<double> &v2sigma2_3, std::vector<double> &v2sigma2_4, std::vector<double> &v2sigma2_5, std::vector<double> &v2sigma2_6);
    void gga_kxc(
        const std::vector<double> &rho_up, const std::vector<double> &rho_down,
        const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
        std::vector<double> &v3rho3_1, std::vector<double> &v3rho3_2, std::vector<double> &v3rho3_3,std::vector<double> &v3rho3_4,
        std::vector<double> &v3rho2sigma_1, std::vector<double> &v3rho2sigma_2, std::vector<double> &v3rho2sigma_3,std::vector<double> &v3rho2sigma_4, std::vector<double> &v3rho2sigma_5, std::vector<double> &v3rho2sigma_6, std::vector<double> &v3rho2sigma_7, std::vector<double> &v3rho2sigma_8, std::vector<double> &v3rho2sigma_9,
        std::vector<double> &v3rhosigma2_1, std::vector<double> &v3rhosigma2_2, std::vector<double> &v3rhosigma2_3,
        std::vector<double> &v3rhosigma2_4, std::vector<double> &v3rhosigma2_5, std::vector<double> &v3rhosigma2_6,
        std::vector<double> &v3rhosigma2_7, std::vector<double> &v3rhosigma2_8, std::vector<double> &v3rhosigma2_9, std::vector<double> &v3rhosigma2_10, std::vector<double> &v3rhosigma2_11, std::vector<double> &v3rhosigma2_12,
        std::vector<double> &v3sigma3_1, std::vector<double> &v3sigma3_2, std::vector<double> &v3sigma3_3,
        std::vector<double> &v3sigma3_4, std::vector<double> &v3sigma3_5, std::vector<double> &v3sigma3_6,
        std::vector<double> &v3sigma3_7, std::vector<double> &v3sigma3_8, std::vector<double> &v3sigma3_9, std::vector<double> &v3sigma3_10
    );
    void gga_lxc(
        const std::vector<double> &rho_up, const std::vector<double> &rho_down,
        const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
        std::vector<double> &v4rho4_1, std::vector<double> &v4rho4_2, std::vector<double> &v4rho4_3,
        std::vector<double> &v4rho4_4, std::vector<double> &v4rho4_5, std::vector<double> &v4rho3sigma_1,
        std::vector<double> &v4rho3sigma_2, std::vector<double> &v4rho3sigma_3, std::vector<double> &v4rho3sigma_4, std::vector<double> &v4rho3sigma_5,std::vector<double> &v4rho3sigma_6,
        std::vector<double> &v4rho3sigma_7, std::vector<double> &v4rho3sigma_8, std::vector<double> &v4rho3sigma_9, std::vector<double> &v4rho3sigma_10, std::vector<double> &v4rho3sigma_11,
        std::vector<double> &v4rho3sigma_12,
        std::vector<double> &v4rho2sigma2_1, std::vector<double> &v4rho2sigma2_2,
        std::vector<double> &v4rho2sigma2_3, std::vector<double> &v4rho2sigma2_4, std::vector<double> &v4rho2sigma2_5,std::vector<double> &v4rho2sigma2_6,
        std::vector<double> &v4rho2sigma2_7, std::vector<double> &v4rho2sigma2_8, std::vector<double> &v4rho2sigma2_9, std::vector<double> &v4rho2sigma2_10,
        std::vector<double> &v4rho2sigma2_11, std::vector<double> &v4rho2sigma2_12,
        std::vector<double> &v4rho2sigma2_13, std::vector<double> &v4rho2sigma2_14, std::vector<double> &v4rho2sigma2_15,
        std::vector<double> &v4rho2sigma2_16, std::vector<double> &v4rho2sigma2_17, std::vector<double> &v4rho2sigma2_18,  
        std::vector<double> &v4rhosigma3_1, std::vector<double> &v4rhosigma3_2, std::vector<double> &v4rhosigma3_3,
        std::vector<double> &v4rhosigma3_4, std::vector<double> &v4rhosigma3_5, std::vector<double> &v4rhosigma3_6,
        std::vector<double> &v4rhosigma3_7, std::vector<double> &v4rhosigma3_8, std::vector<double> &v4rhosigma3_9,
        std::vector<double> &v4rhosigma3_10, std::vector<double> &v4rhosigma3_11, std::vector<double> &v4rhosigma3_12,
        std::vector<double> &v4rhosigma3_13, std::vector<double> &v4rhosigma3_14, std::vector<double> &v4rhosigma3_15,
        std::vector<double> &v4rhosigma3_16, std::vector<double> &v4rhosigma3_17, std::vector<double> &v4rhosigma3_18,
        std::vector<double> &v4rhosigma3_19, std::vector<double> &v4rhosigma3_20,
        std::vector<double> &v4sigma4_1, std::vector<double> &v4sigma4_2, std::vector<double> &v4sigma4_3,
        std::vector<double> &v4sigma4_4, std::vector<double> &v4sigma4_5, std::vector<double> &v4sigma4_6,
        std::vector<double> &v4sigma4_7, std::vector<double> &v4sigma4_8, std::vector<double> &v4sigma4_9,
        std::vector<double> &v4sigma4_10, std::vector<double> &v4sigma4_11, std::vector<double> &v4sigma4_12,
        std::vector<double> &v4sigma4_13, std::vector<double> &v4sigma4_14, std::vector<double> &v4sigma4_15
    );
    std::vector<double> mgga_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down, 
                                             const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, 
                                             const std::vector<double> &sigma_3, const std::vector<double> &lapl_up, 
                                             const std::vector<double> &lapl_down, const std::vector<double> &tau_up, 
                                             const std::vector<double> &tau_down);
    void mgga_vxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                              const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
                              const std::vector<double> &lapl_up, const std::vector<double> &lapl_down,
                              const std::vector<double> &tau_up, const std::vector<double> &tau_down,
                              std::vector<double> &vrho_1, std::vector<double> &vrho_2,
                              std::vector<double> &vsigma_1, std::vector<double> &vsigma_2, std::vector<double> &vsigma_3,
                              std::vector<double> &vlapl_1, std::vector<double> &vlapl_2,
                              std::vector<double> &vtau_1, std::vector<double> &vtau_2);
    void mgga_fxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                              const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
                              const std::vector<double> &lapl_up, const std::vector<double> &lapl_down,
                              const std::vector<double> &tau_up, const std::vector<double> &tau_down,
                              std::vector<double> &v2rho2_1, std::vector<double> &v2rho2_2, std::vector<double> &v2rho2_3,
                              std::vector<double> &v2rhosigma_1, std::vector<double> &v2rhosigma_2, std::vector<double> &v2rhosigma_3,
                              std::vector<double> &v2rhosigma_4, std::vector<double> &v2rhosigma_5, std::vector<double> &v2rhosigma_6,
                              std::vector<double> &v2rholapl_1, std::vector<double> &v2rholapl_2, std::vector<double> &v2rholapl_3, std::vector<double> &v2rholapl_4,
                              std::vector<double> &v2rhotau_1, std::vector<double> &v2rhotau_2, std::vector<double> &v2rhotau_3,std::vector<double> &v2rhotau_4,
                              std::vector<double> &v2sigma2_1, std::vector<double> &v2sigma2_2, std::vector<double> &v2sigma2_3,
                              std::vector<double> &v2sigma2_4, std::vector<double> &v2sigma2_5, std::vector<double> &v2sigma2_6,
                              std::vector<double> &v2sigmalapl_1, std::vector<double> &v2sigmalapl_2, std::vector<double> &v2sigmalapl_3,
                              std::vector<double> &v2sigmalapl_4, std::vector<double> &v2sigmalapl_5, std::vector<double> &v2sigmalapl_6,
                              std::vector<double> &v2sigmatau_1, std::vector<double> &v2sigmatau_2, std::vector<double> &v2sigmatau_3,
                              std::vector<double> &v2sigmatau_4, std::vector<double> &v2sigmatau_5, std::vector<double> &v2sigmatau_6,
                              std::vector<double> &v2lapl2_1, std::vector<double> &v2lapl2_2, std::vector<double> &v2lapl2_3,
                              std::vector<double> &v2lapltau_1, std::vector<double> &v2lapltau_2,  std::vector<double> &v2lapltau_3,  std::vector<double> &v2lapltau_4,
                              std::vector<double> &v2tau2_1, std::vector<double> &v2tau2_2, std::vector<double> &v2tau2_3);
    

    void example_lda_spin();
    void example_gga_spin();
    void example_mgga_spin();

    //void test_ggas_with_lxc();

};

#endif // INTERFACE_TO_LIBXC_H