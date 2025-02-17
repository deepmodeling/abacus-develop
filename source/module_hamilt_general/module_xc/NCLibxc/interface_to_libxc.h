// Programmed by Xiaoyu Zhang at Peking University, Beijing, China 2024/08/26

#ifndef INTERFACE_TO_LIBXC_H
#define INTERFACE_TO_LIBXC_H

#include <xc.h>
#include <vector>
#include <stdexcept>

namespace NCXC{

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

    void mgga_kxc(
    const std::vector<double> &rho_up, const std::vector<double> &rho_down,
    const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
    const std::vector<double> &lapl_up, const std::vector<double> &lapl_down,
    const std::vector<double> &tau_up, const std::vector<double> &tau_down,
    std::vector<double> &v3rho3_1, std::vector<double> &v3rho3_2, std::vector<double> &v3rho3_3, std::vector<double> &v3rho3_4,
    std::vector<double> &v3rho2sigma_1, std::vector<double> &v3rho2sigma_2, std::vector<double> &v3rho2sigma_3,
    std::vector<double> &v3rho2sigma_4, std::vector<double> &v3rho2sigma_5, std::vector<double> &v3rho2sigma_6,
    std::vector<double> &v3rho2sigma_7, std::vector<double> &v3rho2sigma_8, std::vector<double> &v3rho2sigma_9,
    std::vector<double> &v3rho2lapl_1, std::vector<double> &v3rho2lapl_2, std::vector<double> &v3rho2lapl_3,
    std::vector<double> &v3rho2lapl_4, std::vector<double> &v3rho2lapl_5, std::vector<double> &v3rho2lapl_6,
    std::vector<double> &v3rho2tau_1, std::vector<double> &v3rho2tau_2, std::vector<double> &v3rho2tau_3,
    std::vector<double> &v3rho2tau_4, std::vector<double> &v3rho2tau_5, std::vector<double> &v3rho2tau_6,
    std::vector<double> &v3rhosigma2_1, std::vector<double> &v3rhosigma2_2, std::vector<double> &v3rhosigma2_3,
    std::vector<double> &v3rhosigma2_4, std::vector<double> &v3rhosigma2_5, std::vector<double> &v3rhosigma2_6,
    std::vector<double> &v3rhosigma2_7, std::vector<double> &v3rhosigma2_8, std::vector<double> &v3rhosigma2_9,
    std::vector<double> &v3rhosigma2_10, std::vector<double> &v3rhosigma2_11, std::vector<double> &v3rhosigma2_12,
    std::vector<double> &v3rhosigmalapl_1, std::vector<double> &v3rhosigmalapl_2, std::vector<double> &v3rhosigmalapl_3,
    std::vector<double> &v3rhosigmalapl_4, std::vector<double> &v3rhosigmalapl_5, std::vector<double> &v3rhosigmalapl_6,
    std::vector<double> &v3rhosigmalapl_7, std::vector<double> &v3rhosigmalapl_8, std::vector<double> &v3rhosigmalapl_9,
    std::vector<double> &v3rhosigmalapl_10, std::vector<double> &v3rhosigmalapl_11, std::vector<double> &v3rhosigmalapl_12,
    std::vector<double> &v3rhosigmatau_1, std::vector<double> &v3rhosigmatau_2, std::vector<double> &v3rhosigmatau_3,
    std::vector<double> &v3rhosigmatau_4, std::vector<double> &v3rhosigmatau_5, std::vector<double> &v3rhosigmatau_6,
    std::vector<double> &v3rhosigmatau_7, std::vector<double> &v3rhosigmatau_8, std::vector<double> &v3rhosigmatau_9,
    std::vector<double> &v3rhosigmatau_10, std::vector<double> &v3rhosigmatau_11, std::vector<double> &v3rhosigmatau_12,
    std::vector<double> &v3rholapl2_1, std::vector<double> &v3rholapl2_2, std::vector<double> &v3rholapl2_3,
    std::vector<double> &v3rholapl2_4, std::vector<double> &v3rholapl2_5, std::vector<double> &v3rholapl2_6,
    std::vector<double> &v3rholapltau_1, std::vector<double> &v3rholapltau_2, std::vector<double> &v3rholapltau_3,
    std::vector<double> &v3rholapltau_4, std::vector<double> &v3rholapltau_5, std::vector<double> &v3rholapltau_6,
    std::vector<double> &v3rholapltau_7, std::vector<double> &v3rholapltau_8,
    std::vector<double> &v3rhotau2_1, std::vector<double> &v3rhotau2_2, std::vector<double> &v3rhotau2_3,
    std::vector<double> &v3rhotau2_4, std::vector<double> &v3rhotau2_5, std::vector<double> &v3rhotau2_6,
    std::vector<double> &v3sigma3_1, std::vector<double> &v3sigma3_2, std::vector<double> &v3sigma3_3,
    std::vector<double> &v3sigma3_4, std::vector<double> &v3sigma3_5, std::vector<double> &v3sigma3_6,
    std::vector<double> &v3sigma3_7, std::vector<double> &v3sigma3_8, std::vector<double> &v3sigma3_9,
    std::vector<double> &v3sigma3_10,
    std::vector<double> &v3sigma2lapl_1, std::vector<double> &v3sigma2lapl_2, std::vector<double> &v3sigma2lapl_3,
    std::vector<double> &v3sigma2lapl_4, std::vector<double> &v3sigma2lapl_5, std::vector<double> &v3sigma2lapl_6,
    std::vector<double> &v3sigma2lapl_7, std::vector<double> &v3sigma2lapl_8, std::vector<double> &v3sigma2lapl_9,
    std::vector<double> &v3sigma2lapl_10, std::vector<double> &v3sigma2lapl_11, std::vector<double> &v3sigma2lapl_12,
    std::vector<double> &v3sigma2tau_1, std::vector<double> &v3sigma2tau_2, std::vector<double> &v3sigma2tau_3,
    std::vector<double> &v3sigma2tau_4, std::vector<double> &v3sigma2tau_5, std::vector<double> &v3sigma2tau_6,
    std::vector<double> &v3sigma2tau_7, std::vector<double> &v3sigma2tau_8, std::vector<double> &v3sigma2tau_9,
    std::vector<double> &v3sigma2tau_10, std::vector<double> &v3sigma2tau_11, std::vector<double> &v3sigma2tau_12,
    std::vector<double> &v3sigmalapl2_1, std::vector<double> &v3sigmalapl2_2, std::vector<double> &v3sigmalapl2_3,
    std::vector<double> &v3sigmalapl2_4, std::vector<double> &v3sigmalapl2_5, std::vector<double> &v3sigmalapl2_6,
    std::vector<double> &v3sigmalapl2_7, std::vector<double> &v3sigmalapl2_8, std::vector<double> &v3sigmalapl2_9,
    std::vector<double> &v3sigmalapltau_1, std::vector<double> &v3sigmalapltau_2, std::vector<double> &v3sigmalapltau_3,
    std::vector<double> &v3sigmalapltau_4, std::vector<double> &v3sigmalapltau_5, std::vector<double> &v3sigmalapltau_6,
    std::vector<double> &v3sigmalapltau_7, std::vector<double> &v3sigmalapltau_8, std::vector<double> &v3sigmalapltau_9,
    std::vector<double> &v3sigmalapltau_10, std::vector<double> &v3sigmalapltau_11, std::vector<double> &v3sigmalapltau_12,
    std::vector<double> &v3sigmatau2_1, std::vector<double> &v3sigmatau2_2, std::vector<double> &v3sigmatau2_3,
    std::vector<double> &v3sigmatau2_4, std::vector<double> &v3sigmatau2_5, std::vector<double> &v3sigmatau2_6,
    std::vector<double> &v3sigmatau2_7, std::vector<double> &v3sigmatau2_8, std::vector<double> &v3sigmatau2_9,
    std::vector<double> &v3lapl3_1, std::vector<double> &v3lapl3_2, std::vector<double> &v3lapl3_3, std::vector<double> &v3lapl3_4,
    std::vector<double> &v3lapl2tau_1, std::vector<double> &v3lapl2tau_2, std::vector<double> &v3lapl2tau_3,
    std::vector<double> &v3lapl2tau_4, std::vector<double> &v3lapl2tau_5, std::vector<double> &v3lapl2tau_6,
    std::vector<double> &v3lapltau2_1, std::vector<double> &v3lapltau2_2, std::vector<double> &v3lapltau2_3,
    std::vector<double> &v3lapltau2_4, std::vector<double> &v3lapltau2_5, std::vector<double> &v3lapltau2_6,
    std::vector<double> &v3tau3_1, std::vector<double> &v3tau3_2, std::vector<double> &v3tau3_3, std::vector<double> &v3tau3_4);
    
    void mgga_lxc(
    const std::vector<double> &rho_up, const std::vector<double> &rho_down,
    const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
    const std::vector<double> &lapl_up, const std::vector<double> &lapl_down,
    const std::vector<double> &tau_up, const std::vector<double> &tau_down,
    // Output arrays:
    // v4rho4 => 5*np
    std::vector<double> &v4rho4_1, std::vector<double> &v4rho4_2, std::vector<double> &v4rho4_3,
    std::vector<double> &v4rho4_4, std::vector<double> &v4rho4_5,
    // v4rho3sigma => 12*np
    std::vector<double> &v4rho3sigma_1, std::vector<double> &v4rho3sigma_2, std::vector<double> &v4rho3sigma_3,
    std::vector<double> &v4rho3sigma_4, std::vector<double> &v4rho3sigma_5, std::vector<double> &v4rho3sigma_6,
    std::vector<double> &v4rho3sigma_7, std::vector<double> &v4rho3sigma_8, std::vector<double> &v4rho3sigma_9,
    std::vector<double> &v4rho3sigma_10, std::vector<double> &v4rho3sigma_11, std::vector<double> &v4rho3sigma_12,
    // v4rho3lapl => 8*np
    std::vector<double> &v4rho3lapl_1, std::vector<double> &v4rho3lapl_2, std::vector<double> &v4rho3lapl_3,
    std::vector<double> &v4rho3lapl_4, std::vector<double> &v4rho3lapl_5, std::vector<double> &v4rho3lapl_6,
    std::vector<double> &v4rho3lapl_7, std::vector<double> &v4rho3lapl_8,
    // v4rho3tau => 8*np
    std::vector<double> &v4rho3tau_1, std::vector<double> &v4rho3tau_2, std::vector<double> &v4rho3tau_3,
    std::vector<double> &v4rho3tau_4, std::vector<double> &v4rho3tau_5, std::vector<double> &v4rho3tau_6,
    std::vector<double> &v4rho3tau_7, std::vector<double> &v4rho3tau_8,
    // v4rho2sigma2 => 18*np
    std::vector<double> &v4rho2sigma2_1, std::vector<double> &v4rho2sigma2_2, std::vector<double> &v4rho2sigma2_3,
    std::vector<double> &v4rho2sigma2_4, std::vector<double> &v4rho2sigma2_5, std::vector<double> &v4rho2sigma2_6,
    std::vector<double> &v4rho2sigma2_7, std::vector<double> &v4rho2sigma2_8, std::vector<double> &v4rho2sigma2_9,
    std::vector<double> &v4rho2sigma2_10, std::vector<double> &v4rho2sigma2_11, std::vector<double> &v4rho2sigma2_12,
    std::vector<double> &v4rho2sigma2_13, std::vector<double> &v4rho2sigma2_14, std::vector<double> &v4rho2sigma2_15,
    std::vector<double> &v4rho2sigma2_16, std::vector<double> &v4rho2sigma2_17, std::vector<double> &v4rho2sigma2_18,
    // v4rho2sigmalapl => 18*np
    std::vector<double> &v4rho2sigmalapl_1, std::vector<double> &v4rho2sigmalapl_2, std::vector<double> &v4rho2sigmalapl_3,
    std::vector<double> &v4rho2sigmalapl_4, std::vector<double> &v4rho2sigmalapl_5, std::vector<double> &v4rho2sigmalapl_6,
    std::vector<double> &v4rho2sigmalapl_7, std::vector<double> &v4rho2sigmalapl_8, std::vector<double> &v4rho2sigmalapl_9,
    std::vector<double> &v4rho2sigmalapl_10, std::vector<double> &v4rho2sigmalapl_11, std::vector<double> &v4rho2sigmalapl_12,
    std::vector<double> &v4rho2sigmalapl_13, std::vector<double> &v4rho2sigmalapl_14, std::vector<double> &v4rho2sigmalapl_15,
    std::vector<double> &v4rho2sigmalapl_16, std::vector<double> &v4rho2sigmalapl_17, std::vector<double> &v4rho2sigmalapl_18,
    // v4rho2sigmatau => 18*np
    std::vector<double> &v4rho2sigmatau_1, std::vector<double> &v4rho2sigmatau_2, std::vector<double> &v4rho2sigmatau_3,
    std::vector<double> &v4rho2sigmatau_4, std::vector<double> &v4rho2sigmatau_5, std::vector<double> &v4rho2sigmatau_6,
    std::vector<double> &v4rho2sigmatau_7, std::vector<double> &v4rho2sigmatau_8, std::vector<double> &v4rho2sigmatau_9,
    std::vector<double> &v4rho2sigmatau_10, std::vector<double> &v4rho2sigmatau_11, std::vector<double> &v4rho2sigmatau_12,
    std::vector<double> &v4rho2sigmatau_13, std::vector<double> &v4rho2sigmatau_14, std::vector<double> &v4rho2sigmatau_15,
    std::vector<double> &v4rho2sigmatau_16, std::vector<double> &v4rho2sigmatau_17, std::vector<double> &v4rho2sigmatau_18,
    // v4rho2lapl2 => 9*np
    std::vector<double> &v4rho2lapl2_1, std::vector<double> &v4rho2lapl2_2, std::vector<double> &v4rho2lapl2_3,
    std::vector<double> &v4rho2lapl2_4, std::vector<double> &v4rho2lapl2_5, std::vector<double> &v4rho2lapl2_6,
    std::vector<double> &v4rho2lapl2_7, std::vector<double> &v4rho2lapl2_8, std::vector<double> &v4rho2lapl2_9,
    // v4rho2lapltau => 12*np
    std::vector<double> &v4rho2lapltau_1, std::vector<double> &v4rho2lapltau_2, std::vector<double> &v4rho2lapltau_3,
    std::vector<double> &v4rho2lapltau_4, std::vector<double> &v4rho2lapltau_5, std::vector<double> &v4rho2lapltau_6,
    std::vector<double> &v4rho2lapltau_7, std::vector<double> &v4rho2lapltau_8, std::vector<double> &v4rho2lapltau_9,
    std::vector<double> &v4rho2lapltau_10, std::vector<double> &v4rho2lapltau_11, std::vector<double> &v4rho2lapltau_12,
    // v4rho2tau2 => 9*np
    std::vector<double> &v4rho2tau2_1, std::vector<double> &v4rho2tau2_2, std::vector<double> &v4rho2tau2_3,
    std::vector<double> &v4rho2tau2_4, std::vector<double> &v4rho2tau2_5, std::vector<double> &v4rho2tau2_6,
    std::vector<double> &v4rho2tau2_7, std::vector<double> &v4rho2tau2_8, std::vector<double> &v4rho2tau2_9,
    // v4rhosigma3 => 20*np
    std::vector<double> &v4rhosigma3_1, std::vector<double> &v4rhosigma3_2, std::vector<double> &v4rhosigma3_3,
    std::vector<double> &v4rhosigma3_4, std::vector<double> &v4rhosigma3_5, std::vector<double> &v4rhosigma3_6,
    std::vector<double> &v4rhosigma3_7, std::vector<double> &v4rhosigma3_8, std::vector<double> &v4rhosigma3_9,
    std::vector<double> &v4rhosigma3_10, std::vector<double> &v4rhosigma3_11, std::vector<double> &v4rhosigma3_12,
    std::vector<double> &v4rhosigma3_13, std::vector<double> &v4rhosigma3_14, std::vector<double> &v4rhosigma3_15,
    std::vector<double> &v4rhosigma3_16, std::vector<double> &v4rhosigma3_17, std::vector<double> &v4rhosigma3_18,
    std::vector<double> &v4rhosigma3_19, std::vector<double> &v4rhosigma3_20,
    // v4rhosigma2lapl => 24*np
    std::vector<double> &v4rhosigma2lapl_1, std::vector<double> &v4rhosigma2lapl_2, std::vector<double> &v4rhosigma2lapl_3,
    std::vector<double> &v4rhosigma2lapl_4, std::vector<double> &v4rhosigma2lapl_5, std::vector<double> &v4rhosigma2lapl_6,
    std::vector<double> &v4rhosigma2lapl_7, std::vector<double> &v4rhosigma2lapl_8, std::vector<double> &v4rhosigma2lapl_9,
    std::vector<double> &v4rhosigma2lapl_10, std::vector<double> &v4rhosigma2lapl_11, std::vector<double> &v4rhosigma2lapl_12,
    std::vector<double> &v4rhosigma2lapl_13, std::vector<double> &v4rhosigma2lapl_14, std::vector<double> &v4rhosigma2lapl_15,
    std::vector<double> &v4rhosigma2lapl_16, std::vector<double> &v4rhosigma2lapl_17, std::vector<double> &v4rhosigma2lapl_18,
    std::vector<double> &v4rhosigma2lapl_19, std::vector<double> &v4rhosigma2lapl_20, std::vector<double> &v4rhosigma2lapl_21,
    std::vector<double> &v4rhosigma2lapl_22, std::vector<double> &v4rhosigma2lapl_23, std::vector<double> &v4rhosigma2lapl_24,
    // v4rhosigma2tau => 24*np
    std::vector<double> &v4rhosigma2tau_1, std::vector<double> &v4rhosigma2tau_2, std::vector<double> &v4rhosigma2tau_3,
    std::vector<double> &v4rhosigma2tau_4, std::vector<double> &v4rhosigma2tau_5, std::vector<double> &v4rhosigma2tau_6,
    std::vector<double> &v4rhosigma2tau_7, std::vector<double> &v4rhosigma2tau_8, std::vector<double> &v4rhosigma2tau_9,
    std::vector<double> &v4rhosigma2tau_10, std::vector<double> &v4rhosigma2tau_11, std::vector<double> &v4rhosigma2tau_12,
    std::vector<double> &v4rhosigma2tau_13, std::vector<double> &v4rhosigma2tau_14, std::vector<double> &v4rhosigma2tau_15,
    std::vector<double> &v4rhosigma2tau_16, std::vector<double> &v4rhosigma2tau_17, std::vector<double> &v4rhosigma2tau_18,
    std::vector<double> &v4rhosigma2tau_19, std::vector<double> &v4rhosigma2tau_20, std::vector<double> &v4rhosigma2tau_21,
    std::vector<double> &v4rhosigma2tau_22, std::vector<double> &v4rhosigma2tau_23, std::vector<double> &v4rhosigma2tau_24,
    // v4rhosigmalapl2 => 18*np
    std::vector<double> &v4rhosigmalapl2_1, std::vector<double> &v4rhosigmalapl2_2, std::vector<double> &v4rhosigmalapl2_3,
    std::vector<double> &v4rhosigmalapl2_4, std::vector<double> &v4rhosigmalapl2_5, std::vector<double> &v4rhosigmalapl2_6,
    std::vector<double> &v4rhosigmalapl2_7, std::vector<double> &v4rhosigmalapl2_8, std::vector<double> &v4rhosigmalapl2_9,
    std::vector<double> &v4rhosigmalapl2_10, std::vector<double> &v4rhosigmalapl2_11, std::vector<double> &v4rhosigmalapl2_12,
    std::vector<double> &v4rhosigmalapl2_13, std::vector<double> &v4rhosigmalapl2_14, std::vector<double> &v4rhosigmalapl2_15,
    std::vector<double> &v4rhosigmalapl2_16, std::vector<double> &v4rhosigmalapl2_17, std::vector<double> &v4rhosigmalapl2_18,
    // v4rhosigmalapltau => 24*np
    std::vector<double> &v4rhosigmalapltau_1, std::vector<double> &v4rhosigmalapltau_2, std::vector<double> &v4rhosigmalapltau_3,
    std::vector<double> &v4rhosigmalapltau_4, std::vector<double> &v4rhosigmalapltau_5, std::vector<double> &v4rhosigmalapltau_6,
    std::vector<double> &v4rhosigmalapltau_7, std::vector<double> &v4rhosigmalapltau_8, std::vector<double> &v4rhosigmalapltau_9,
    std::vector<double> &v4rhosigmalapltau_10, std::vector<double> &v4rhosigmalapltau_11, std::vector<double> &v4rhosigmalapltau_12,
    std::vector<double> &v4rhosigmalapltau_13, std::vector<double> &v4rhosigmalapltau_14, std::vector<double> &v4rhosigmalapltau_15,
    std::vector<double> &v4rhosigmalapltau_16, std::vector<double> &v4rhosigmalapltau_17, std::vector<double> &v4rhosigmalapltau_18,
    std::vector<double> &v4rhosigmalapltau_19, std::vector<double> &v4rhosigmalapltau_20, std::vector<double> &v4rhosigmalapltau_21,
    std::vector<double> &v4rhosigmalapltau_22, std::vector<double> &v4rhosigmalapltau_23, std::vector<double> &v4rhosigmalapltau_24,
    // v4rhosigmatau2 => 18*np
    std::vector<double> &v4rhosigmatau2_1, std::vector<double> &v4rhosigmatau2_2, std::vector<double> &v4rhosigmatau2_3,
    std::vector<double> &v4rhosigmatau2_4, std::vector<double> &v4rhosigmatau2_5, std::vector<double> &v4rhosigmatau2_6,
    std::vector<double> &v4rhosigmatau2_7, std::vector<double> &v4rhosigmatau2_8, std::vector<double> &v4rhosigmatau2_9,
    std::vector<double> &v4rhosigmatau2_10, std::vector<double> &v4rhosigmatau2_11, std::vector<double> &v4rhosigmatau2_12,
    std::vector<double> &v4rhosigmatau2_13, std::vector<double> &v4rhosigmatau2_14, std::vector<double> &v4rhosigmatau2_15,
    std::vector<double> &v4rhosigmatau2_16, std::vector<double> &v4rhosigmatau2_17, std::vector<double> &v4rhosigmatau2_18,
    // v4rholapl3 => 8*np
std::vector<double> &v4rholapl3_1, std::vector<double> &v4rholapl3_2, std::vector<double> &v4rholapl3_3,
std::vector<double> &v4rholapl3_4, std::vector<double> &v4rholapl3_5, std::vector<double> &v4rholapl3_6,
std::vector<double> &v4rholapl3_7, std::vector<double> &v4rholapl3_8,

// v4rholapl2tau => 12*np
std::vector<double> &v4rholapl2tau_1, std::vector<double> &v4rholapl2tau_2, std::vector<double> &v4rholapl2tau_3,
std::vector<double> &v4rholapl2tau_4, std::vector<double> &v4rholapl2tau_5, std::vector<double> &v4rholapl2tau_6,
std::vector<double> &v4rholapl2tau_7, std::vector<double> &v4rholapl2tau_8, std::vector<double> &v4rholapl2tau_9,
std::vector<double> &v4rholapl2tau_10, std::vector<double> &v4rholapl2tau_11, std::vector<double> &v4rholapl2tau_12,

// v4rholapltau2 => 12*np
std::vector<double> &v4rholapltau2_1, std::vector<double> &v4rholapltau2_2, std::vector<double> &v4rholapltau2_3,
std::vector<double> &v4rholapltau2_4, std::vector<double> &v4rholapltau2_5, std::vector<double> &v4rholapltau2_6,
std::vector<double> &v4rholapltau2_7, std::vector<double> &v4rholapltau2_8, std::vector<double> &v4rholapltau2_9,
std::vector<double> &v4rholapltau2_10, std::vector<double> &v4rholapltau2_11, std::vector<double> &v4rholapltau2_12,

// v4rhotau3 => 8*np
std::vector<double> &v4rhotau3_1, std::vector<double> &v4rhotau3_2, std::vector<double> &v4rhotau3_3,
std::vector<double> &v4rhotau3_4, std::vector<double> &v4rhotau3_5, std::vector<double> &v4rhotau3_6,
std::vector<double> &v4rhotau3_7, std::vector<double> &v4rhotau3_8,

// v4sigma4 => 15*np
std::vector<double> &v4sigma4_1, std::vector<double> &v4sigma4_2, std::vector<double> &v4sigma4_3,
std::vector<double> &v4sigma4_4, std::vector<double> &v4sigma4_5, std::vector<double> &v4sigma4_6,
std::vector<double> &v4sigma4_7, std::vector<double> &v4sigma4_8, std::vector<double> &v4sigma4_9,
std::vector<double> &v4sigma4_10, std::vector<double> &v4sigma4_11, std::vector<double> &v4sigma4_12,
std::vector<double> &v4sigma4_13, std::vector<double> &v4sigma4_14, std::vector<double> &v4sigma4_15,

// v4sigma3lapl => 20*np
std::vector<double> &v4sigma3lapl_1, std::vector<double> &v4sigma3lapl_2, std::vector<double> &v4sigma3lapl_3,
std::vector<double> &v4sigma3lapl_4, std::vector<double> &v4sigma3lapl_5, std::vector<double> &v4sigma3lapl_6,
std::vector<double> &v4sigma3lapl_7, std::vector<double> &v4sigma3lapl_8, std::vector<double> &v4sigma3lapl_9,
std::vector<double> &v4sigma3lapl_10, std::vector<double> &v4sigma3lapl_11, std::vector<double> &v4sigma3lapl_12,
std::vector<double> &v4sigma3lapl_13, std::vector<double> &v4sigma3lapl_14, std::vector<double> &v4sigma3lapl_15,
std::vector<double> &v4sigma3lapl_16, std::vector<double> &v4sigma3lapl_17, std::vector<double> &v4sigma3lapl_18,
std::vector<double> &v4sigma3lapl_19, std::vector<double> &v4sigma3lapl_20,

// v4sigma3tau => 20*np
std::vector<double> &v4sigma3tau_1, std::vector<double> &v4sigma3tau_2, std::vector<double> &v4sigma3tau_3,
std::vector<double> &v4sigma3tau_4, std::vector<double> &v4sigma3tau_5, std::vector<double> &v4sigma3tau_6,
std::vector<double> &v4sigma3tau_7, std::vector<double> &v4sigma3tau_8, std::vector<double> &v4sigma3tau_9,
std::vector<double> &v4sigma3tau_10, std::vector<double> &v4sigma3tau_11, std::vector<double> &v4sigma3tau_12,
std::vector<double> &v4sigma3tau_13, std::vector<double> &v4sigma3tau_14, std::vector<double> &v4sigma3tau_15,
std::vector<double> &v4sigma3tau_16, std::vector<double> &v4sigma3tau_17, std::vector<double> &v4sigma3tau_18,
std::vector<double> &v4sigma3tau_19, std::vector<double> &v4sigma3tau_20,

// v4sigma2lapl2 => 18*np
std::vector<double> &v4sigma2lapl2_1, std::vector<double> &v4sigma2lapl2_2, std::vector<double> &v4sigma2lapl2_3,
std::vector<double> &v4sigma2lapl2_4, std::vector<double> &v4sigma2lapl2_5, std::vector<double> &v4sigma2lapl2_6,
std::vector<double> &v4sigma2lapl2_7, std::vector<double> &v4sigma2lapl2_8, std::vector<double> &v4sigma2lapl2_9,
std::vector<double> &v4sigma2lapl2_10, std::vector<double> &v4sigma2lapl2_11, std::vector<double> &v4sigma2lapl2_12,
std::vector<double> &v4sigma2lapl2_13, std::vector<double> &v4sigma2lapl2_14, std::vector<double> &v4sigma2lapl2_15,
std::vector<double> &v4sigma2lapl2_16, std::vector<double> &v4sigma2lapl2_17, std::vector<double> &v4sigma2lapl2_18,

// v4sigma2lapltau => 24*np
std::vector<double> &v4sigma2lapltau_1, std::vector<double> &v4sigma2lapltau_2, std::vector<double> &v4sigma2lapltau_3,
std::vector<double> &v4sigma2lapltau_4, std::vector<double> &v4sigma2lapltau_5, std::vector<double> &v4sigma2lapltau_6,
std::vector<double> &v4sigma2lapltau_7, std::vector<double> &v4sigma2lapltau_8, std::vector<double> &v4sigma2lapltau_9,
std::vector<double> &v4sigma2lapltau_10, std::vector<double> &v4sigma2lapltau_11, std::vector<double> &v4sigma2lapltau_12,
std::vector<double> &v4sigma2lapltau_13, std::vector<double> &v4sigma2lapltau_14, std::vector<double> &v4sigma2lapltau_15,
std::vector<double> &v4sigma2lapltau_16, std::vector<double> &v4sigma2lapltau_17, std::vector<double> &v4sigma2lapltau_18,
std::vector<double> &v4sigma2lapltau_19, std::vector<double> &v4sigma2lapltau_20, std::vector<double> &v4sigma2lapltau_21,
std::vector<double> &v4sigma2lapltau_22, std::vector<double> &v4sigma2lapltau_23, std::vector<double> &v4sigma2lapltau_24,

// v4sigma2tau2 => 18*np
std::vector<double> &v4sigma2tau2_1, std::vector<double> &v4sigma2tau2_2, std::vector<double> &v4sigma2tau2_3,
std::vector<double> &v4sigma2tau2_4, std::vector<double> &v4sigma2tau2_5, std::vector<double> &v4sigma2tau2_6,
std::vector<double> &v4sigma2tau2_7, std::vector<double> &v4sigma2tau2_8, std::vector<double> &v4sigma2tau2_9,
std::vector<double> &v4sigma2tau2_10, std::vector<double> &v4sigma2tau2_11, std::vector<double> &v4sigma2tau2_12,
std::vector<double> &v4sigma2tau2_13, std::vector<double> &v4sigma2tau2_14, std::vector<double> &v4sigma2tau2_15,
std::vector<double> &v4sigma2tau2_16, std::vector<double> &v4sigma2tau2_17, std::vector<double> &v4sigma2tau2_18,

// v4sigmalapl3 => 12*np
std::vector<double> &v4sigmalapl3_1, std::vector<double> &v4sigmalapl3_2, std::vector<double> &v4sigmalapl3_3,
std::vector<double> &v4sigmalapl3_4, std::vector<double> &v4sigmalapl3_5, std::vector<double> &v4sigmalapl3_6,
std::vector<double> &v4sigmalapl3_7, std::vector<double> &v4sigmalapl3_8, std::vector<double> &v4sigmalapl3_9,
std::vector<double> &v4sigmalapl3_10, std::vector<double> &v4sigmalapl3_11, std::vector<double> &v4sigmalapl3_12,

// v4sigmalapl2tau => 18*np
std::vector<double> &v4sigmalapl2tau_1, std::vector<double> &v4sigmalapl2tau_2, std::vector<double> &v4sigmalapl2tau_3,
std::vector<double> &v4sigmalapl2tau_4, std::vector<double> &v4sigmalapl2tau_5, std::vector<double> &v4sigmalapl2tau_6,
std::vector<double> &v4sigmalapl2tau_7, std::vector<double> &v4sigmalapl2tau_8, std::vector<double> &v4sigmalapl2tau_9,
std::vector<double> &v4sigmalapl2tau_10, std::vector<double> &v4sigmalapl2tau_11, std::vector<double> &v4sigmalapl2tau_12,
std::vector<double> &v4sigmalapl2tau_13, std::vector<double> &v4sigmalapl2tau_14, std::vector<double> &v4sigmalapl2tau_15,
std::vector<double> &v4sigmalapl2tau_16, std::vector<double> &v4sigmalapl2tau_17, std::vector<double> &v4sigmalapl2tau_18,

// v4sigmalapltau2 => 18*np
std::vector<double> &v4sigmalapltau2_1, std::vector<double> &v4sigmalapltau2_2, std::vector<double> &v4sigmalapltau2_3,
std::vector<double> &v4sigmalapltau2_4, std::vector<double> &v4sigmalapltau2_5, std::vector<double> &v4sigmalapltau2_6,
std::vector<double> &v4sigmalapltau2_7, std::vector<double> &v4sigmalapltau2_8, std::vector<double> &v4sigmalapltau2_9,
std::vector<double> &v4sigmalapltau2_10, std::vector<double> &v4sigmalapltau2_11, std::vector<double> &v4sigmalapltau2_12,
std::vector<double> &v4sigmalapltau2_13, std::vector<double> &v4sigmalapltau2_14, std::vector<double> &v4sigmalapltau2_15,
std::vector<double> &v4sigmalapltau2_16, std::vector<double> &v4sigmalapltau2_17, std::vector<double> &v4sigmalapltau2_18,

// v4sigmatau3 => 12*np
std::vector<double> &v4sigmatau3_1, std::vector<double> &v4sigmatau3_2, std::vector<double> &v4sigmatau3_3,
std::vector<double> &v4sigmatau3_4, std::vector<double> &v4sigmatau3_5, std::vector<double> &v4sigmatau3_6,
std::vector<double> &v4sigmatau3_7, std::vector<double> &v4sigmatau3_8, std::vector<double> &v4sigmatau3_9,
std::vector<double> &v4sigmatau3_10, std::vector<double> &v4sigmatau3_11, std::vector<double> &v4sigmatau3_12,

// v4lapl4 => 5*np
std::vector<double> &v4lapl4_1, std::vector<double> &v4lapl4_2, std::vector<double> &v4lapl4_3,
std::vector<double> &v4lapl4_4, std::vector<double> &v4lapl4_5,

// v4lapl3tau => 8*np
std::vector<double> &v4lapl3tau_1, std::vector<double> &v4lapl3tau_2, std::vector<double> &v4lapl3tau_3,
std::vector<double> &v4lapl3tau_4, std::vector<double> &v4lapl3tau_5, std::vector<double> &v4lapl3tau_6,
std::vector<double> &v4lapl3tau_7, std::vector<double> &v4lapl3tau_8,

// v4lapl2tau2 => 9*np
std::vector<double> &v4lapl2tau2_1, std::vector<double> &v4lapl2tau2_2, std::vector<double> &v4lapl2tau2_3,
std::vector<double> &v4lapl2tau2_4, std::vector<double> &v4lapl2tau2_5, std::vector<double> &v4lapl2tau2_6,
std::vector<double> &v4lapl2tau2_7, std::vector<double> &v4lapl2tau2_8, std::vector<double> &v4lapl2tau2_9,

// v4lapltau3 => 8*np
std::vector<double> &v4lapltau3_1, std::vector<double> &v4lapltau3_2, std::vector<double> &v4lapltau3_3,
std::vector<double> &v4lapltau3_4, std::vector<double> &v4lapltau3_5, std::vector<double> &v4lapltau3_6,
std::vector<double> &v4lapltau3_7, std::vector<double> &v4lapltau3_8,
    std::vector<double> &v4tau4_1, std::vector<double> &v4tau4_2, std::vector<double> &v4tau4_3,
    std::vector<double> &v4tau4_4, std::vector<double> &v4tau4_5
);


    void example_lda_spin();
    void example_gga_spin();
    void example_mgga_spin();

    //void test_ggas_with_lxc();

};
}
#endif // INTERFACE_TO_LIBXC_H