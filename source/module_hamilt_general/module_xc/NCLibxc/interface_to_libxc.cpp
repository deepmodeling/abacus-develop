
// This program is used to interface with libxc which is a library for exchange-correlation functionals in density functional theory.
/*
How to compile as an independent program:
source /public1/home/scg0213/software-scg0213/libxc-7.0.0/install/libxc.sh
g++ -std=c++11 -o my_program interface_to_libxc.cpp -I/public1/home/scg0213/software-scg0213/libxc-7.0.0/install/include -L/public1/home/scg0213/software-scg0213/libxc-7.0.0/install/lib -lxc
./my_program
*/
// func_id can be found in the libxc website https://libxc.gitlab.io/functionals/
#include "interface_to_libxc.h"
#include <iostream>
#include <cstdio>

LibxcInterface::LibxcInterface(int xc_id, bool spin_polarized)
{
    int spin_option = spin_polarized ? XC_POLARIZED : XC_UNPOLARIZED;
    if (xc_func_init(&func, xc_id, spin_option) != 0)
    {
        std::cerr << "Failed to initialize functional." << std::endl;
        throw std::runtime_error("Libxc initialization error");
    }
}

LibxcInterface::~LibxcInterface()
{
    xc_func_end(&func);
}

    /////////////////////////////////////////////////////////////
    // LDA START, up to the fourth derivative
    // LDA Energy Density for spin-polarized systems
    std::vector<double> LibxcInterface::lda_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
        }
        std::vector<double> exc(np);
        xc_lda_exc(&func, np, rho.data(), exc.data());
        return exc;
    }

    // LDA Potential for spin-polarized systems
    void LibxcInterface::lda_vxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &vrho_1, std::vector<double> &vrho_2)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> vrho(2 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
        }
        xc_lda_vxc(&func, np, rho.data(), vrho.data());
        for (int i = 0; i < np; ++i)
        {
            vrho_1[i] = vrho[2 * i];
            vrho_2[i] = vrho[2 * i + 1];
        }
    }

    // LDA the second derivative for spin-polarized systems
    void LibxcInterface::lda_fxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v2rho2_1, std::vector<double> &v2rho2_2, std::vector<double> &v2rho2_3)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> v2rho2(3 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
        }
        xc_lda_fxc(&func, np, rho.data(), v2rho2.data());
        for (int i = 0; i < np; ++i)
        {
            v2rho2_1[i] = v2rho2[3 * i];
            v2rho2_2[i] = v2rho2[3 * i + 1];
            v2rho2_3[i] = v2rho2[3 * i + 2];
        }
    }

    // LDA the third derivative for spin-polarized systems
    void LibxcInterface::lda_kxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v3rho3_1, std::vector<double> &v3rho3_2, std::vector<double> &v3rho3_3, std::vector<double> &v3rho3_4)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> v3rho3(4 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
        }
        xc_lda_kxc(&func, np, rho.data(), v3rho3.data());
        for (int i = 0; i < np; ++i)
        {
            v3rho3_1[i] = v3rho3[4 * i];
            v3rho3_2[i] = v3rho3[4 * i + 1];
            v3rho3_3[i] = v3rho3[4 * i + 2];
            v3rho3_4[i] = v3rho3[4 * i + 3];
        }
    }

    // LDA the fourth derivative for spin-polarized systems
    void LibxcInterface::lda_lxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 std::vector<double> &v4rho4_1, std::vector<double> &v4rho4_2, std::vector<double> &v4rho4_3, std::vector<double> &v4rho4_4, std::vector<double> &v4rho4_5)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> v4rho4(5 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
        }
        xc_lda_lxc(&func, np, rho.data(), v4rho4.data());
        for (int i = 0; i < np; ++i)
        {
            v4rho4_1[i] = v4rho4[5 * i];
            v4rho4_2[i] = v4rho4[5 * i + 1];
            v4rho4_3[i] = v4rho4[5 * i + 2];
            v4rho4_4[i] = v4rho4[5 * i + 3];
            v4rho4_5[i] = v4rho4[5 * i + 4];
        }
    }
    //##################  LDA END ###############################################


    /////////////////////////////////////////////////////////////
    //GGA START, up to the fourth derivative
    // GGA Energy Density for spin-polarized systems
    std::vector<double> LibxcInterface::gga_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down, const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
            sigma[3 * i] = sigma_1[i];
            sigma[3 * i + 1] = sigma_2[i];
            sigma[3 * i + 2] = sigma_3[i];
        }
        std::vector<double> exc(np);
        xc_gga_exc(&func, np, rho.data(), sigma.data(), exc.data());
        return exc;
    }

    // GGA Potential for spin-polarized systems
    void LibxcInterface::gga_vxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3, std::vector<double> &vrho_1, std::vector<double> &vrho_2,
                 std::vector<double> &vsigma_1, std::vector<double> &vsigma_2, std::vector<double> &vsigma_3)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np);
        std::vector<double> vrho(2 * np);
        std::vector<double> vsigma(3 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
            sigma[3 * i] = sigma_1[i];
            sigma[3 * i + 1] = sigma_2[i];
            sigma[3 * i + 2] = sigma_3[i];
        }
        xc_gga_vxc(&func, np, rho.data(), sigma.data(), vrho.data(), vsigma.data());
        for (int i = 0; i < np; ++i)
        {
            vrho_1[i] = vrho[2 * i];
            vrho_2[i] = vrho[2 * i + 1];
            vsigma_1[i] = vsigma[3 * i];
            vsigma_2[i] = vsigma[3 * i + 1];
            vsigma_3[i] = vsigma[3 * i + 2];
        }
    }

    // GGA the second derivative for spin-polarized systems
    void LibxcInterface::gga_fxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                 const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
                 std::vector<double> &v2rho2_1, std::vector<double> &v2rho2_2, std::vector<double> &v2rho2_3,
                 std::vector<double> &v2rhosigma_1, std::vector<double> &v2rhosigma_2, std::vector<double> &v2rhosigma_3, std::vector<double> &v2rhosigma_4, std::vector<double> &v2rhosigma_5, std::vector<double> &v2rhosigma_6, std::vector<double> &v2sigma2_1, std::vector<double> &v2sigma2_2, std::vector<double> &v2sigma2_3, std::vector<double> &v2sigma2_4, std::vector<double> &v2sigma2_5, std::vector<double> &v2sigma2_6)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np);
        std::vector<double> v2rho2(3 * np);
        std::vector<double> v2rhosigma(6 * np);
        std::vector<double> v2sigma2(6 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
            sigma[3 * i] = sigma_1[i];
            sigma[3 * i + 1] = sigma_2[i];
            sigma[3 * i + 2] = sigma_3[i];
        }
        xc_gga_fxc(&func, np, rho.data(), sigma.data(), v2rho2.data(), v2rhosigma.data(), v2sigma2.data());
        for (int i = 0; i < np; ++i)
        {
            v2rho2_1[i] = v2rho2[3 * i];
            v2rho2_2[i] = v2rho2[3 * i + 1];
            v2rho2_3[i] = v2rho2[3 * i + 2];
            v2rhosigma_1[i] = v2rhosigma[6 * i];
            v2rhosigma_2[i] = v2rhosigma[6 * i + 1];
            v2rhosigma_3[i] = v2rhosigma[6 * i + 2];
            v2rhosigma_4[i] = v2rhosigma[6 * i + 3];
            v2rhosigma_5[i] = v2rhosigma[6 * i + 4];
            v2rhosigma_6[i] = v2rhosigma[6 * i + 5];
            v2sigma2_1[i] = v2sigma2[6 * i];
            v2sigma2_2[i] = v2sigma2[6 * i + 1];
            v2sigma2_3[i] = v2sigma2[6 * i + 2];
            v2sigma2_4[i] = v2sigma2[6 * i + 3];
            v2sigma2_5[i] = v2sigma2[6 * i + 4];
            v2sigma2_6[i] = v2sigma2[6 * i + 5];
        }
    }

    // GGA the third derivative for spin-polarized systems
    void LibxcInterface::gga_kxc(
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
    )
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np), v3rho3(4 * np), v3rho2sigma(9 * np), v3rhosigma2(12 * np), v3sigma3(10 * np);
    
        for(int i = 0; i < np; ++i) {
            rho[2*i]     = rho_up[i];
            rho[2*i + 1] = rho_down[i];
            sigma[3*i]   = sigma_1[i];
            sigma[3*i+1] = sigma_2[i];
            sigma[3*i+2] = sigma_3[i];
        }
    
        xc_gga_kxc(&func, np, rho.data(), sigma.data(),
                   v3rho3.data(), v3rho2sigma.data(), v3rhosigma2.data(), v3sigma3.data());

    
        // Assign results back
        for(int i = 0; i < np; ++i) {
            v3rho3_1[i]         = v3rho3[4 * i + 0];
            v3rho3_2[i]         = v3rho3[4 * i + 1];
            v3rho3_3[i]         = v3rho3[4 * i + 2];
            v3rho3_4[i]         = v3rho3[4 * i + 3];
            v3rho2sigma_1[i]    = v3rho2sigma[9 * i + 0];
            v3rho2sigma_2[i]    = v3rho2sigma[9 * i + 1];
            v3rho2sigma_3[i]    = v3rho2sigma[9 * i + 2];
            v3rho2sigma_4[i]    = v3rho2sigma[9 * i + 3];
            v3rho2sigma_5[i]    = v3rho2sigma[9 * i + 4];
            v3rho2sigma_6[i]    = v3rho2sigma[9 * i + 5];
            v3rho2sigma_7[i]    = v3rho2sigma[9 * i + 6];
            v3rho2sigma_8[i]    = v3rho2sigma[9 * i + 7];
            v3rho2sigma_9[i]    = v3rho2sigma[9 * i + 8];
            v3rhosigma2_1[i]    = v3rhosigma2[12 * i + 0];
            v3rhosigma2_2[i]    = v3rhosigma2[12 * i + 1];
            v3rhosigma2_3[i]    = v3rhosigma2[12 * i + 2];
            v3rhosigma2_4[i]    = v3rhosigma2[12 * i + 3];
            v3rhosigma2_5[i]    = v3rhosigma2[12 * i + 4];
            v3rhosigma2_6[i]    = v3rhosigma2[12 * i + 5];
            v3rhosigma2_7[i]    = v3rhosigma2[12 * i + 6];
            v3rhosigma2_8[i]    = v3rhosigma2[12 * i + 7];
            v3rhosigma2_9[i]    = v3rhosigma2[12 * i + 8];
            v3rhosigma2_10[i]   = v3rhosigma2[12 * i + 9];
            v3rhosigma2_11[i]   = v3rhosigma2[12 * i + 10];
            v3rhosigma2_12[i]   = v3rhosigma2[12 * i + 11];
            v3sigma3_1[i]       = v3sigma3[10 * i + 0];
            v3sigma3_2[i]       = v3sigma3[10 * i + 1];
            v3sigma3_3[i]       = v3sigma3[10 * i + 2];
            v3sigma3_4[i]       = v3sigma3[10 * i + 3];
            v3sigma3_5[i]       = v3sigma3[10 * i + 4];
            v3sigma3_6[i]       = v3sigma3[10 * i + 5];
            v3sigma3_7[i]       = v3sigma3[10 * i + 6];
            v3sigma3_8[i]       = v3sigma3[10 * i + 7];
            v3sigma3_9[i]       = v3sigma3[10 * i + 8];
            v3sigma3_10[i]      = v3sigma3[10 * i + 9];
        }
    }
    
    // GGA the fourth derivative for spin-polarized systems
    void LibxcInterface::gga_lxc(
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
    )
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np), v4rho4(5 * np), v4rho3sigma(12 * np), v4rho2sigma2(18 * np),v4rhosigma3(20 * np), v4sigma4(15 * np);
    
        for(int i = 0; i < np; ++i) {
            rho[2*i]     = rho_up[i];
            rho[2*i + 1] = rho_down[i];
            sigma[3*i]   = sigma_1[i];
            sigma[3*i+1] = sigma_2[i];
            sigma[3*i+2] = sigma_3[i];
        }
    
        xc_gga_lxc(&func, np, rho.data(), sigma.data(),
                   v4rho4.data(), v4rho3sigma.data(), v4rho2sigma2.data(), v4rhosigma3.data(), v4sigma4.data());

    
        // Assign results back
        for(int i = 0; i < np; ++i) {
            v4rho4_1[i]         = v4rho4[5*i + 0];
            v4rho4_2[i]         = v4rho4[5*i + 1];
            v4rho4_3[i]         = v4rho4[5*i + 2];
            v4rho4_4[i]         = v4rho4[5*i + 3];
            v4rho4_5[i]         = v4rho4[5*i + 4];
            v4rho3sigma_1[i]    = v4rho3sigma[12*i + 0];
            v4rho3sigma_2[i]    = v4rho3sigma[12*i + 1];
            v4rho3sigma_3[i]    = v4rho3sigma[12*i + 2];
            v4rho3sigma_4[i]    = v4rho3sigma[12*i + 3];
            v4rho3sigma_5[i]    = v4rho3sigma[12*i + 4];
            v4rho3sigma_6[i]    = v4rho3sigma[12*i + 5];
            v4rho3sigma_7[i]    = v4rho3sigma[12*i + 6];
            v4rho3sigma_8[i]    = v4rho3sigma[12*i + 7];
            v4rho3sigma_9[i]    = v4rho3sigma[12*i + 8];
            v4rho3sigma_10[i]   = v4rho3sigma[12*i + 9];
            v4rho3sigma_11[i]   = v4rho3sigma[12*i + 10];
            v4rho3sigma_12[i]   = v4rho3sigma[12*i + 11];
            v4rho2sigma2_1[i]   = v4rho2sigma2[18*i + 0];
            v4rho2sigma2_2[i]   = v4rho2sigma2[18*i + 1];
            v4rho2sigma2_3[i]   = v4rho2sigma2[18*i + 2];
            v4rho2sigma2_4[i]   = v4rho2sigma2[18*i + 3];
            v4rho2sigma2_5[i]   = v4rho2sigma2[18*i + 4];
            v4rho2sigma2_6[i]   = v4rho2sigma2[18*i + 5];
            v4rho2sigma2_7[i]   = v4rho2sigma2[18*i + 6];
            v4rho2sigma2_8[i]   = v4rho2sigma2[18*i + 7];
            v4rho2sigma2_9[i]   = v4rho2sigma2[18*i + 8];
            v4rho2sigma2_10[i]  = v4rho2sigma2[18*i + 9];
            v4rho2sigma2_11[i]  = v4rho2sigma2[18*i + 10];
            v4rho2sigma2_12[i]  = v4rho2sigma2[18*i + 11];
            v4rho2sigma2_13[i]  = v4rho2sigma2[18*i + 12];
            v4rho2sigma2_14[i]  = v4rho2sigma2[18*i + 13];
            v4rho2sigma2_15[i]  = v4rho2sigma2[18*i + 14];
            v4rho2sigma2_16[i]  = v4rho2sigma2[18*i + 15];
            v4rho2sigma2_17[i]  = v4rho2sigma2[18*i + 16];
            v4rho2sigma2_18[i]  = v4rho2sigma2[18*i + 17];
            v4rhosigma3_1[i]    = v4rhosigma3[20*i + 0];
            v4rhosigma3_2[i]    = v4rhosigma3[20*i + 1];
            v4rhosigma3_3[i]    = v4rhosigma3[20*i + 2];
            v4rhosigma3_4[i]    = v4rhosigma3[20*i + 3];
            v4rhosigma3_5[i]    = v4rhosigma3[20*i + 4];
            v4rhosigma3_6[i]    = v4rhosigma3[20*i + 5];
            v4rhosigma3_7[i]    = v4rhosigma3[20*i + 6];
            v4rhosigma3_8[i]    = v4rhosigma3[20*i + 7];
            v4rhosigma3_9[i]    = v4rhosigma3[20*i + 8];
            v4rhosigma3_10[i]   = v4rhosigma3[20*i + 9];
            v4rhosigma3_11[i]   = v4rhosigma3[20*i + 10];
            v4rhosigma3_12[i]   = v4rhosigma3[20*i + 11];
            v4rhosigma3_13[i]   = v4rhosigma3[20*i + 12];
            v4rhosigma3_14[i]   = v4rhosigma3[20*i + 13];
            v4rhosigma3_15[i]   = v4rhosigma3[20*i + 14];
            v4rhosigma3_16[i]   = v4rhosigma3[20*i + 15];
            v4rhosigma3_17[i]   = v4rhosigma3[20*i + 16];
            v4rhosigma3_18[i]   = v4rhosigma3[20*i + 17];
            v4rhosigma3_19[i]   = v4rhosigma3[20*i + 18];
            v4rhosigma3_20[i]   = v4rhosigma3[20*i + 19];      
            v4sigma4_1[i]       = v4sigma4[15*i + 0];
            v4sigma4_2[i]       = v4sigma4[15*i + 1];
            v4sigma4_3[i]       = v4sigma4[15*i + 2];
            v4sigma4_4[i]       = v4sigma4[15*i + 3];
            v4sigma4_5[i]       = v4sigma4[15*i + 4];
            v4sigma4_6[i]       = v4sigma4[15*i + 5];
            v4sigma4_7[i]       = v4sigma4[15*i + 6];
            v4sigma4_8[i]       = v4sigma4[15*i + 7];
            v4sigma4_9[i]       = v4sigma4[15*i + 8];
            v4sigma4_10[i]      = v4sigma4[15*i + 9];
            v4sigma4_11[i]      = v4sigma4[15*i + 10];
            v4sigma4_12[i]      = v4sigma4[15*i + 11];
            v4sigma4_13[i]      = v4sigma4[15*i + 12];
            v4sigma4_14[i]      = v4sigma4[15*i + 13];
            v4sigma4_15[i]      = v4sigma4[15*i + 14];
        }
    }
    

    // GGA END
    /////////////////////////////////////////////////////////////
    
    // meta-GGA START, up to the second derivative
    // meta-GGA Energy Density for spin-polarized systems
    std::vector<double> LibxcInterface::mgga_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down, 
                                             const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, 
                                             const std::vector<double> &sigma_3, const std::vector<double> &lapl_up, 
                                             const std::vector<double> &lapl_down, const std::vector<double> &tau_up, 
                                             const std::vector<double> &tau_down)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np);
        std::vector<double> lapl(2 * np);
        std::vector<double> tau(2 * np);

        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
            sigma[3 * i] = sigma_1[i];
            sigma[3 * i + 1] = sigma_2[i];
            sigma[3 * i + 2] = sigma_3[i];
            lapl[2 * i] = lapl_up[i];
            lapl[2 * i + 1] = lapl_down[i];
            tau[2 * i] = tau_up[i];
            tau[2 * i + 1] = tau_down[i];
        }

    std::vector<double> exc(np);
    xc_mgga_exc(&func, np, rho.data(), sigma.data(), lapl.data(), tau.data(), exc.data());

    return exc;
    }
    // meta-GGA potential for spin-polarized systems
    // MGGA Potential for spin-polarized systems
    void LibxcInterface::mgga_vxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
                              const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3,
                              const std::vector<double> &lapl_up, const std::vector<double> &lapl_down,
                              const std::vector<double> &tau_up, const std::vector<double> &tau_down,
                              std::vector<double> &vrho_1, std::vector<double> &vrho_2,
                              std::vector<double> &vsigma_1, std::vector<double> &vsigma_2, std::vector<double> &vsigma_3,
                              std::vector<double> &vlapl_1, std::vector<double> &vlapl_2,
                              std::vector<double> &vtau_1, std::vector<double> &vtau_2)
{
    int np = rho_up.size();
    std::vector<double> rho(2 * np);
    std::vector<double> sigma(3 * np);
    std::vector<double> lapl(2 * np);
    std::vector<double> tau(2 * np);
    std::vector<double> vrho(2 * np);
    std::vector<double> vsigma(3 * np);
    std::vector<double> vlapl(2 * np);
    std::vector<double> vtau(2 * np);

    for (int i = 0; i < np; ++i)
    {
        rho[2 * i] = rho_up[i];
        rho[2 * i + 1] = rho_down[i];
        sigma[3 * i] = sigma_1[i];
        sigma[3 * i + 1] = sigma_2[i];
        sigma[3 * i + 2] = sigma_3[i];
        lapl[2 * i] = lapl_up[i];
        lapl[2 * i + 1] = lapl_down[i];
        tau[2 * i] = tau_up[i];
        tau[2 * i + 1] = tau_down[i];
    }

    xc_mgga_vxc(&func, np, rho.data(), sigma.data(), lapl.data(), tau.data(), vrho.data(), vsigma.data(), vlapl.data(), vtau.data());

    for (int i = 0; i < np; ++i)
    {
        vrho_1[i] = vrho[2 * i];
        vrho_2[i] = vrho[2 * i + 1];
        vsigma_1[i] = vsigma[3 * i];
        vsigma_2[i] = vsigma[3 * i + 1];
        vsigma_3[i] = vsigma[3 * i + 2];
        vlapl_1[i] = vlapl[2 * i];
        vlapl_2[i] = vlapl[2 * i + 1];
        vtau_1[i] = vtau[2 * i];
        vtau_2[i] = vtau[2 * i + 1];
    }
}
    // MGGA the second derivative for spin-polarized systems
    void LibxcInterface::mgga_fxc(const std::vector<double> &rho_up, const std::vector<double> &rho_down,
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
                              std::vector<double> &v2tau2_1, std::vector<double> &v2tau2_2, std::vector<double> &v2tau2_3)
{
    int np = rho_up.size();
    std::vector<double> rho(2 * np);
    std::vector<double> sigma(3 * np);
    std::vector<double> lapl(2 * np);
    std::vector<double> tau(2 * np);
    std::vector<double> v2rho2(3 * np);
    std::vector<double> v2rhosigma(6 * np);
    std::vector<double> v2rholapl(4 * np);
    std::vector<double> v2rhotau(4 * np);
    std::vector<double> v2sigma2(6 * np);
    std::vector<double> v2sigmalapl(6 * np);
    std::vector<double> v2sigmatau(6 * np);
    std::vector<double> v2lapl2(3 * np);
    std::vector<double> v2lapltau(4 * np);
    std::vector<double> v2tau2(3 * np);

    for (int i = 0; i < np; ++i)
    {
        rho[2 * i] = rho_up[i];
        rho[2 * i + 1] = rho_down[i];
        sigma[3 * i] = sigma_1[i];
        sigma[3 * i + 1] = sigma_2[i];
        sigma[3 * i + 2] = sigma_3[i];
        lapl[2 * i] = lapl_up[i];
        lapl[2 * i + 1] = lapl_down[i];
        tau[2 * i] = tau_up[i];
        tau[2 * i + 1] = tau_down[i];
    }

    xc_mgga_fxc(&func, np, rho.data(), sigma.data(), lapl.data(), tau.data(), v2rho2.data(), v2rhosigma.data(), v2rholapl.data(), v2rhotau.data(), v2sigma2.data(), v2sigmalapl.data(), v2sigmatau.data(), v2lapl2.data(), v2lapltau.data(), v2tau2.data());

    for (int i = 0; i < np; ++i)
    {
        v2rho2_1[i] = v2rho2[3 * i];
        v2rho2_2[i] = v2rho2[3 * i + 1];
        v2rho2_3[i] = v2rho2[3 * i + 2];
        v2rhosigma_1[i] = v2rhosigma[6 * i];
        v2rhosigma_2[i] = v2rhosigma[6 * i + 1];
        v2rhosigma_3[i] = v2rhosigma[6 * i + 2];
        v2rhosigma_4[i] = v2rhosigma[6 * i + 3];
        v2rhosigma_5[i] = v2rhosigma[6 * i + 4];
        v2rhosigma_6[i] = v2rhosigma[6 * i + 5];
        v2rholapl_1[i] = v2rholapl[4 * i];
        v2rholapl_2[i] = v2rholapl[4 * i + 1];
        v2rholapl_3[i] = v2rholapl[4 * i + 2];
        v2rholapl_4[i] = v2rholapl[4 * i + 3];
        v2rhotau_1[i] = v2rhotau[4 * i];
        v2rhotau_2[i] = v2rhotau[4 * i + 1];
        v2rhotau_3[i] = v2rhotau[4 * i + 2];
        v2rhotau_4[i] = v2rhotau[4 * i + 3];
        v2sigma2_1[i] = v2sigma2[6 * i];
        v2sigma2_2[i] = v2sigma2[6 * i + 1];
        v2sigma2_3[i] = v2sigma2[6 * i + 2];
        v2sigma2_4[i] = v2sigma2[6 * i + 3];
        v2sigma2_5[i] = v2sigma2[6 * i + 4];
        v2sigma2_6[i] = v2sigma2[6 * i + 5];
        v2sigmalapl_1[i] = v2sigmalapl[6 * i];
        v2sigmalapl_2[i] = v2sigmalapl[6 * i + 1];
        v2sigmalapl_3[i] = v2sigmalapl[6 * i + 2];
        v2sigmalapl_4[i] = v2sigmalapl[6 * i + 3];
        v2sigmalapl_5[i] = v2sigmalapl[6 * i + 4];
        v2sigmalapl_6[i] = v2sigmalapl[6 * i + 5];
        v2sigmatau_1[i] = v2sigmatau[6 * i];
        v2sigmatau_2[i] = v2sigmatau[6 * i + 1];
        v2sigmatau_3[i] = v2sigmatau[6 * i + 2];
        v2sigmatau_4[i] = v2sigmatau[6 * i + 3];
        v2sigmatau_5[i] = v2sigmatau[6 * i + 4];
        v2sigmatau_6[i] = v2sigmatau[6 * i + 5];
        v2lapl2_1[i] = v2lapl2[3 * i];
        v2lapl2_2[i] = v2lapl2[3 * i + 1];
        v2lapl2_3[i] = v2lapl2[3 * i + 2];
        v2lapltau_1[i] = v2lapltau[4 * i];
        v2lapltau_2[i] = v2lapltau[4 * i + 1];
        v2lapltau_3[i] = v2lapltau[4 * i + 2];
        v2lapltau_4[i] = v2lapltau[4 * i + 3];
        v2tau2_1[i] = v2tau2[3 * i];
        v2tau2_2[i] = v2tau2[3 * i + 1];
        v2tau2_3[i] = v2tau2[3 * i + 2];
    }
}
// meta-GGA END
//#########################################################

    // Example usage for spin-polarized LDA
    // The input rho_up and rho_down are the densities for spin-up and spin-down electrons, respectively.
    // rhoup = (n+s)/2, rhodown = (n-s)/2
    //output:
    //exc[]: the energy density per unit particle
    //vrho[]: first derivative of the energy per unit volume
    // v2rho2[]: second derivative of the energy per unit volume
    //v3rho3[]: third derivative of the energy per unit volume
    // v4rho4[]: fourth derivative of the energy per unit volume
    // detailed expression can be found in https://libxc.gitlab.io/manual/libxc-5.1.x/
    //vrho [(0),(1)]
    //v2rho2 [(0,0),(0,1),(1,1)]
    //v3rho3 [(0,0,0),(0,0,1),(0,1,1),(1,1,1)]
    //v4rho4 [(0,0,0,0),(0,0,0,1),(0,0,1,1),(0,1,1,1),(1,1,1,1)]
    // 0:up 1:down (spin)

    void LibxcInterface::example_lda_spin()
    {
        std::vector<double> rho_up = {0.1, 0.2, 0.3};
        std::vector<double> rho_down = {0.2, 0.1, 0.3};

        auto exc = lda_exc(rho_up, rho_down);
        std::vector<double> vrho_1(rho_up.size()), vrho_2(rho_down.size());
        lda_vxc(rho_up, rho_down, vrho_1, vrho_2);

        std::vector<double> v2rho2_1(rho_up.size()), v2rho2_2(rho_down.size()), v2rho2_3(rho_up.size());
        lda_fxc(rho_up, rho_down, v2rho2_1, v2rho2_2, v2rho2_3);

        std::vector<double> v3rho3_1(rho_up.size()), v3rho3_2(rho_down.size()), v3rho3_3(rho_up.size()), v3rho3_4(rho_down.size());
        lda_kxc(rho_up, rho_down, v3rho3_1, v3rho3_2, v3rho3_3, v3rho3_4);

        //        std::vector<double> v4rho4_1(rho_up.size()), v4rho4_2(rho_down.size()), v4rho4_3(rho_up.size()), v4rho4_4(rho_down.size()), v4rho4_5(rho_up.size());
        //        lda_lxc(rho_up, rho_down, v4rho4_1, v4rho4_2, v4rho4_3, v4rho4_4, v4rho4_5);

        std::cout << "Exc: ";
        for (auto e : exc)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "VRho_1: ";
        for (auto v : vrho_1)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "VRho_2: ";
        for (auto v : vrho_2)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "V2Rho2_1: ";
        for (auto f : v2rho2_1)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "V2Rho2_2: ";
        for (auto f : v2rho2_2)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "V2Rho2_3: ";
        for (auto f : v2rho2_3)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "V3Rho3_1: ";
        for (auto k : v3rho3_1) std::cout << k << " ";
        std::cout << std::endl;

        //        std::cout << "V3Rho3_2: ";
        //        for (auto k : v3rho3_2) std::cout << k << " ";
        //        std::cout << std::endl;

        //        std::cout << "V3Rho3_3: ";
        //        for (auto k : v3rho3_3) std::cout << k << " ";
        //        std::cout << std::endl;

        //        std::cout << "V3Rho3_4: ";
        //        for (auto k : v3rho3_4) std::cout << k << " ";
        //        std::cout << std::endl;

        //        std::cout << "V4Rho4_1: ";
        //        for (auto l : v4rho4_1) std::cout << l << " ";
        //        std::cout << std::endl;

        //        std::cout << "V4Rho4_2: ";
        //        for (auto l : v4rho4_2) std::cout << l << " ";
        //        std::cout << std::endl;

        //        std::cout << "V4Rho4_3: ";
        //        for (auto l : v4rho4_3) std::cout << l << " ";
        //        std::cout << std::endl;

        //        std::cout << "V4Rho4_4: ";
        //        for (auto l : v4rho4_4) std::cout << l << " ";
        //        std::cout << std::endl;

        //        std::cout << "V4Rho4_5: ";
        //        for (auto l : v4rho4_5) std::cout << l << " ";
        //        std::cout << std::endl;
    }

    // // Example usage for spin-polarized GGA
    //    input:
    //      np: number of points
    //       rho[]: the density
    //      sigma[]: contracted gradients of the density
    //     output:
    //      exc[]: the energy per unit particle
    //      vrho[]: first partial derivative of the energy per unit volume in terms of the density
    //      vsigma[]: first partial derivative of the energy per unit volume in terms of sigma
    //      v2rho2[]: second partial derivative of the energy per unit volume in terms of the density
    //      v2rhosigma[]: second partial derivative of the energy per unit volume in terms of the density and sigma
    //      v2sigma2[]: second partial derivative of the energy per unit volume in terms of sigma
    //    // detailed expression can be found in https://libxc.gitlab.io/manual/libxc-5.1.x/

    void LibxcInterface::example_gga_spin()
    {
        // 定义自旋极化系统的示例密度和梯度数据
        std::vector<double> rho_up = {0.1, 0.1, 0.3};
        std::vector<double> rho_down = {0.1, 0.2, 0.3};
        std::vector<double> sigma_1 = {0.01, 0.02, 0.03};
        std::vector<double> sigma_2 = {0.02, 0.01, 0.03};
        std::vector<double> sigma_3 = {0.01, 0.02, 0.03};

        // 计算交换-相关能量密度
        auto exc = gga_exc(rho_up, rho_down, sigma_1, sigma_2, sigma_3);

        std::cout << "GGA Exc: ";
        for (auto e : exc)
            std::cout << e << " ";
        std::cout << std::endl;

        // 计算交换-相关势
        std::vector<double> vrho_1(rho_up.size()), vrho_2(rho_down.size());
        std::vector<double> vsigma_1(sigma_1.size()), vsigma_2(sigma_2.size()), vsigma_3(sigma_3.size());
        gga_vxc(rho_up, rho_down, sigma_1, sigma_2, sigma_3, vrho_1, vrho_2, vsigma_1, vsigma_2, vsigma_3);

        std::cout << "GGA VRho_1: ";
        for (auto v : vrho_1)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "GGA VRho_2: ";
        for (auto v : vrho_2)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "GGA VSigma_1: ";
        for (auto v : vsigma_1)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "GGA VSigma_2: ";
        for (auto v : vsigma_2)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "GGA VSigma_3: ";
        for (auto v : vsigma_3)
            std::cout << v << " ";
        std::cout << std::endl;

        // 计算交换-相关势的二阶导数
        std::vector<double> v2rho2_1(rho_up.size()), v2rho2_2(rho_down.size()), v2rho2_3(rho_up.size());
        std::vector<double> v2rhosigma_1(rho_up.size()), v2rhosigma_2(rho_up.size()), v2rhosigma_3(rho_up.size()), v2rhosigma_4(rho_up.size()), v2rhosigma_5(rho_up.size()), v2rhosigma_6(rho_up.size());
        std::vector<double> v2sigma2_1(sigma_1.size()), v2sigma2_2(sigma_1.size()), v2sigma2_3(sigma_1.size()), v2sigma2_4(sigma_1.size()), v2sigma2_5(sigma_1.size()), v2sigma2_6(sigma_1.size());
        gga_fxc(rho_up, rho_down, sigma_1, sigma_2, sigma_3, v2rho2_1, v2rho2_2, v2rho2_3, v2rhosigma_1, v2rhosigma_2, v2rhosigma_3, v2rhosigma_4, v2rhosigma_5, v2rhosigma_6, v2sigma2_1, v2sigma2_2, v2sigma2_3, v2sigma2_4, v2sigma2_5, v2sigma2_6);

        std::cout << "GGA V2Rho2_1: ";
        for (auto f : v2rho2_1)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Rho2_2: ";
        for (auto f : v2rho2_2)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Rho2_3: ";
        for (auto f : v2rho2_3)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_1: ";
        for (auto f : v2rhosigma_1)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_2: ";
        for (auto f : v2rhosigma_2)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_3: ";
        for (auto f : v2rhosigma_3)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_4: ";
        for (auto f : v2rhosigma_4)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_5: ";
        for (auto f : v2rhosigma_5)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2RhoSigma_6: ";
        for (auto f : v2rhosigma_6)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_1: ";
        for (auto f : v2sigma2_1)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_2: ";
        for (auto f : v2sigma2_2)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_3: ";
        for (auto f : v2sigma2_3)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_4: ";
        for (auto f : v2sigma2_4)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_5: ";
        for (auto f : v2sigma2_5)
            std::cout << f << " ";
        std::cout << std::endl;

        std::cout << "GGA V2Sigma2_6: ";
        for (auto f : v2sigma2_6)
            std::cout << f << " ";
        std::cout << std::endl;

        // calculate the third derivative 
        std::vector<double> v3rho3_1(rho_up.size()), v3rho3_2(rho_up.size()), v3rho3_3(rho_up.size()), v3rho3_4(rho_up.size());
        std::vector<double> v3rho2sigma_1(rho_up.size()), v3rho2sigma_2(rho_up.size()), v3rho2sigma_3(rho_up.size()), v3rho2sigma_4(rho_up.size());
        std::vector<double> v3rho2sigma_5(rho_up.size()), v3rho2sigma_6(rho_up.size()), v3rho2sigma_7(rho_up.size()), v3rho2sigma_8(rho_up.size()), v3rho2sigma_9(rho_up.size());
        std::vector<double> v3rhosigma2_1(rho_up.size()), v3rhosigma2_2(rho_up.size()), v3rhosigma2_3(rho_up.size()), v3rhosigma2_4(rho_up.size()), v3rhosigma2_5(rho_up.size()), v3rhosigma2_6(rho_up.size()), v3rhosigma2_7(rho_up.size()), v3rhosigma2_8(rho_up.size()), v3rhosigma2_9(rho_up.size()), v3rhosigma2_10(rho_up.size()), v3rhosigma2_11(rho_up.size()), v3rhosigma2_12(rho_up.size());
        std::vector<double> v3sigma3_1(rho_up.size()), v3sigma3_2(rho_up.size()), v3sigma3_3(rho_up.size()), v3sigma3_4(rho_up.size()), v3sigma3_5(rho_up.size()), v3sigma3_6(rho_up.size()), v3sigma3_7(rho_up.size()), v3sigma3_8(rho_up.size()), v3sigma3_9(rho_up.size()), v3sigma3_10(rho_up.size());

        gga_kxc(rho_up, rho_down, sigma_1, sigma_2, sigma_3,
            v3rho3_1, v3rho3_2, v3rho3_3, v3rho3_4,
            v3rho2sigma_1, v3rho2sigma_2, v3rho2sigma_3, v3rho2sigma_4, v3rho2sigma_5, v3rho2sigma_6, v3rho2sigma_7, v3rho2sigma_8, v3rho2sigma_9,
            v3rhosigma2_1, v3rhosigma2_2, v3rhosigma2_3, v3rhosigma2_4, v3rhosigma2_5, v3rhosigma2_6, v3rhosigma2_7, v3rhosigma2_8, v3rhosigma2_9, v3rhosigma2_10, v3rhosigma2_11, v3rhosigma2_12,
            v3sigma3_1, v3sigma3_2, v3sigma3_3, v3sigma3_4, v3sigma3_5, v3sigma3_6, v3sigma3_7, v3sigma3_8, v3sigma3_9, v3sigma3_10);

        std::cout << "GGA 3rd derivative (sample) v3rho3_1: ";
        for (auto k : v3rho3_1) std::cout << k << " ";
        std::cout << std::endl;

        std::cout << "GGA 3rd derivative (sample) v3rho3_4: ";
        for (auto k : v3rho3_4) std::cout << k << " ";
        std::cout << std::endl;


        // calculate the fourth derivative
        std::vector<double> v4rho4_1(rho_up.size()), v4rho4_2(rho_up.size()), v4rho4_3(rho_up.size()), v4rho4_4(rho_up.size()), v4rho4_5(rho_up.size());
        std::vector<double> v4rho3sigma_1(rho_up.size()), v4rho3sigma_2(rho_up.size()), v4rho3sigma_3(rho_up.size()), v4rho3sigma_4(rho_up.size()), v4rho3sigma_5(rho_up.size()), v4rho3sigma_6(rho_up.size()), v4rho3sigma_7(rho_up.size()), v4rho3sigma_8(rho_up.size()), v4rho3sigma_9(rho_up.size()), v4rho3sigma_10(rho_up.size()), v4rho3sigma_11(rho_up.size()), v4rho3sigma_12(rho_up.size());

        std::vector<double> v4rho2sigma2_1(rho_up.size()), v4rho2sigma2_2(rho_up.size()), v4rho2sigma2_3(rho_up.size()), v4rho2sigma2_4(rho_up.size()), v4rho2sigma2_5(rho_up.size()), v4rho2sigma2_6(rho_up.size()), v4rho2sigma2_7(rho_up.size()), v4rho2sigma2_8(rho_up.size()), v4rho2sigma2_9(rho_up.size()), v4rho2sigma2_10(rho_up.size()), v4rho2sigma2_11(rho_up.size()), v4rho2sigma2_12(rho_up.size()), v4rho2sigma2_13(rho_up.size()), v4rho2sigma2_14(rho_up.size()), v4rho2sigma2_15(rho_up.size()),v4rho2sigma2_16(rho_up.size()), v4rho2sigma2_17(rho_up.size()), v4rho2sigma2_18(rho_up.size());

        std::vector<double> v4rhosigma3_1(rho_up.size()), v4rhosigma3_2(rho_up.size()), v4rhosigma3_3(rho_up.size()), v4rhosigma3_4(rho_up.size()), v4rhosigma3_5(rho_up.size()), v4rhosigma3_6(rho_up.size()), v4rhosigma3_7(rho_up.size()), v4rhosigma3_8(rho_up.size()), v4rhosigma3_9(rho_up.size()), v4rhosigma3_10(rho_up.size()), v4rhosigma3_11(rho_up.size()), v4rhosigma3_12(rho_up.size()), v4rhosigma3_13(rho_up.size()), v4rhosigma3_14(rho_up.size()), v4rhosigma3_15(rho_up.size()), v4rhosigma3_16(rho_up.size()), v4rhosigma3_17(rho_up.size()), v4rhosigma3_18(rho_up.size()), v4rhosigma3_19(rho_up.size()), v4rhosigma3_20(rho_up.size());

        std::vector<double> v4sigma4_1(rho_up.size()), v4sigma4_2(rho_up.size()), v4sigma4_3(rho_up.size()), v4sigma4_4(rho_up.size()), v4sigma4_5(rho_up.size()), v4sigma4_6(rho_up.size()), v4sigma4_7(rho_up.size()), v4sigma4_8(rho_up.size()), v4sigma4_9(rho_up.size()), v4sigma4_10(rho_up.size()), v4sigma4_11(rho_up.size()), v4sigma4_12(rho_up.size()), v4sigma4_13(rho_up.size()), v4sigma4_14(rho_up.size()), v4sigma4_15(rho_up.size());

        gga_lxc(rho_up, rho_down, sigma_1, sigma_2, sigma_3,
            v4rho4_1, v4rho4_2, v4rho4_3, v4rho4_4, v4rho4_5,
            v4rho3sigma_1, v4rho3sigma_2, v4rho3sigma_3, v4rho3sigma_4, v4rho3sigma_5, v4rho3sigma_6, v4rho3sigma_7, v4rho3sigma_8, v4rho3sigma_9, v4rho3sigma_10, v4rho3sigma_11, v4rho3sigma_12,
            v4rho2sigma2_1, v4rho2sigma2_2, v4rho2sigma2_3, v4rho2sigma2_4, v4rho2sigma2_5, v4rho2sigma2_6, v4rho2sigma2_7, v4rho2sigma2_8, v4rho2sigma2_9, v4rho2sigma2_10, v4rho2sigma2_11, v4rho2sigma2_12, v4rho2sigma2_13, v4rho2sigma2_14, v4rho2sigma2_15,v4rho2sigma2_16, v4rho2sigma2_17, v4rho2sigma2_18,
            v4rhosigma3_1, v4rhosigma3_2, v4rhosigma3_3, v4rhosigma3_4, v4rhosigma3_5, v4rhosigma3_6, v4rhosigma3_7, v4rhosigma3_8, v4rhosigma3_9, v4rhosigma3_10, v4rhosigma3_11, v4rhosigma3_12, v4rhosigma3_13, v4rhosigma3_14, v4rhosigma3_15, v4rhosigma3_16, v4rhosigma3_17, v4rhosigma3_18, v4rhosigma3_19, v4rhosigma3_20,
            v4sigma4_1, v4sigma4_2, v4sigma4_3, v4sigma4_4, v4sigma4_5, v4sigma4_6, v4sigma4_7, v4sigma4_8, v4sigma4_9, v4sigma4_10, v4sigma4_11, v4sigma4_12, v4sigma4_13, v4sigma4_14, v4sigma4_15);

        std::cout << "GGA 4th derivative (sample) v4rhosigma3_1: ";
        for (auto l : v4rhosigma3_1) std::cout << l << " ";
        std::cout << std::endl;

        std::cout << "GGA 4th derivative (sample) v4rhosigma3_20: ";
        for (auto l : v4rhosigma3_20) std::cout << l << " ";
        std::cout << std::endl;
    }

    void LibxcInterface::example_mgga_spin()
{
    // Define example density and gradient data for a spin-polarized system
    std::vector<double> rho_up = {0.1, 0.2, 0.3};
    std::vector<double> rho_down = {0.1, 0.2, 0.3};
    std::vector<double> sigma_1 = {0.01, 0.02, 0.03};
    std::vector<double> sigma_2 = {0.02, 0.01, 0.03};
    std::vector<double> sigma_3 = {0.01, 0.02, 0.03};
    std::vector<double> lapl_up = {0.001, 0.002, 0.003};
    std::vector<double> lapl_down = {0.001, 0.002, 0.003};
    std::vector<double> tau_up = {0.0001, 0.0002, 0.0003};
    std::vector<double> tau_down = {0.0001, 0.0002, 0.0003};

    // Compute exchange-correlation energy density
    auto exc = mgga_exc(rho_up, rho_down, sigma_1, sigma_2, sigma_3, lapl_up, lapl_down, tau_up, tau_down);

    // Compute exchange-correlation potential
    std::vector<double> vrho_1(rho_up.size()), vrho_2(rho_down.size());
    std::vector<double> vsigma_1(sigma_1.size()), vsigma_2(sigma_2.size()), vsigma_3(sigma_3.size());
    std::vector<double> vlapl_1(lapl_up.size()), vlapl_2(lapl_down.size());
    std::vector<double> vtau_1(tau_up.size()), vtau_2(tau_down.size());
    mgga_vxc(rho_up, rho_down, sigma_1, sigma_2, sigma_3, lapl_up, lapl_down, tau_up, tau_down, vrho_1, vrho_2, vsigma_1, vsigma_2, vsigma_3, vlapl_1, vlapl_2, vtau_1, vtau_2);

    // Compute second derivatives of the exchange-correlation potential
    std::vector<double> v2rho2_1(rho_up.size()), v2rho2_2(rho_down.size()), v2rho2_3(rho_up.size());
    std::vector<double> v2rhosigma_1(rho_up.size()), v2rhosigma_2(rho_up.size()), v2rhosigma_3(rho_up.size()), v2rhosigma_4(rho_up.size()), v2rhosigma_5(rho_up.size()), v2rhosigma_6(rho_up.size());
    std::vector<double> v2rholapl_1(rho_up.size()), v2rholapl_2(rho_up.size()), v2rholapl_3(rho_up.size()), v2rholapl_4(rho_up.size());
    std::vector<double> v2rhotau_1(rho_up.size()), v2rhotau_2(rho_up.size()), v2rhotau_3(rho_up.size()), v2rhotau_4(rho_up.size());
    std::vector<double> v2sigma2_1(sigma_1.size()), v2sigma2_2(sigma_1.size()), v2sigma2_3(sigma_1.size()), v2sigma2_4(sigma_1.size()), v2sigma2_5(sigma_1.size()), v2sigma2_6(sigma_1.size());
    std::vector<double> v2sigmalapl_1(sigma_1.size()), v2sigmalapl_2(sigma_1.size()), v2sigmalapl_3(sigma_1.size()), v2sigmalapl_4(sigma_1.size()), v2sigmalapl_5(sigma_1.size()), v2sigmalapl_6(sigma_1.size());
    std::vector<double> v2sigmatau_1(sigma_1.size()), v2sigmatau_2(sigma_1.size()), v2sigmatau_3(sigma_1.size()), v2sigmatau_4(sigma_1.size()), v2sigmatau_5(sigma_1.size()), v2sigmatau_6(sigma_1.size());
    std::vector<double> v2lapl2_1(lapl_up.size()), v2lapl2_2(lapl_up.size()), v2lapl2_3(lapl_up.size());
    std::vector<double> v2lapltau_1(lapl_up.size()), v2lapltau_2(lapl_up.size()), v2lapltau_3(lapl_up.size()), v2lapltau_4(lapl_up.size());
    std::vector<double> v2tau2_1(tau_up.size()), v2tau2_2(tau_up.size()), v2tau2_3(tau_up.size());
    mgga_fxc(rho_up, rho_down, sigma_1, sigma_2, sigma_3, lapl_up, lapl_down, tau_up, tau_down, v2rho2_1, v2rho2_2, v2rho2_3, v2rhosigma_1, v2rhosigma_2, v2rhosigma_3, v2rhosigma_4, v2rhosigma_5, v2rhosigma_6, v2rholapl_1, v2rholapl_2, v2rholapl_3, v2rholapl_4, v2rhotau_1, v2rhotau_2, v2rhotau_3, v2rhotau_4, v2sigma2_1, v2sigma2_2, v2sigma2_3, v2sigma2_4, v2sigma2_5, v2sigma2_6, v2sigmalapl_1, v2sigmalapl_2, v2sigmalapl_3, v2sigmalapl_4, v2sigmalapl_5, v2sigmalapl_6, v2sigmatau_1, v2sigmatau_2, v2sigmatau_3, v2sigmatau_4, v2sigmatau_5, v2sigmatau_6, v2lapl2_1, v2lapl2_2, v2lapl2_3, v2lapltau_1, v2lapltau_2, v2lapltau_3, v2lapltau_4, v2tau2_1, v2tau2_2, v2tau2_3);

    // Output results
    std::cout << "MGGA Exc: ";
    for (auto e : exc)
        std::cout << e << " ";
    std::cout << std::endl;

    std::cout << "MGGA VRho_1: ";
    for (auto v : vrho_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA VRho_2: ";
    for (auto v : vrho_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA VSigma_1: ";
    for (auto v : vsigma_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA VSigma_2: ";
    for (auto v : vsigma_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA VSigma_3: ";
    for (auto v : vsigma_3)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA VLapl_1: ";
    for (auto v : vlapl_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA VLapl_2: ";
    for (auto v : vlapl_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA VTau_1: ";
    for (auto v : vtau_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA VTau_2: ";
    for (auto v : vtau_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Rho2_1: ";
    for (auto v : v2rho2_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Rho2_2: ";
    for (auto v : v2rho2_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Rho2_3: ";
    for (auto v : v2rho2_3)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoSigma_1: ";
    for (auto v : v2rhosigma_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoSigma_2: ";
    for (auto v : v2rhosigma_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoSigma_3: ";
    for (auto v : v2rhosigma_3)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoSigma_4: ";
    for (auto v : v2rhosigma_4)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoSigma_5: ";
    for (auto v : v2rhosigma_5)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoSigma_6: ";
    for (auto v : v2rhosigma_6)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoLapl_1: ";
    for (auto v : v2rholapl_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoLapl_2: ";
    for (auto v : v2rholapl_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoLapl_3: ";
    for (auto v : v2rholapl_3)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoLapl_4: ";
    for (auto v : v2rholapl_4)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoTau_1: ";
    for (auto v : v2rhotau_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoTau_2: ";
    for (auto v : v2rhotau_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoTau_3: ";
    for (auto v : v2rhotau_3)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2RhoTau_4: ";
    for (auto v : v2rhotau_4)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Sigma2_1: ";
    for (auto v : v2sigma2_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Sigma2_2: ";
    for (auto v : v2sigma2_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Sigma2_3: ";
    for (auto v : v2sigma2_3)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Sigma2_4: ";
    for (auto v : v2sigma2_4)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Sigma2_5: ";
    for (auto v : v2sigma2_5)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Sigma2_6: ";
    for (auto v : v2sigma2_6)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaLapl_1: ";
    for (auto v : v2sigmalapl_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaLapl_2: ";
    for (auto v : v2sigmalapl_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaLapl_3: ";
    for (auto v : v2sigmalapl_3)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaLapl_4: ";
    for (auto v : v2sigmalapl_4)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaLapl_5: ";
    for (auto v : v2sigmalapl_5)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaLapl_6: ";
    for (auto v : v2sigmalapl_6)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaTau_1: ";
    for (auto v : v2sigmatau_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaTau_2: ";
    for (auto v : v2sigmatau_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaTau_3: ";
    for (auto v : v2sigmatau_3)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaTau_4: ";
    for (auto v : v2sigmatau_4)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaTau_5: ";
    for (auto v : v2sigmatau_5)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2SigmaTau_6: ";
    for (auto v : v2sigmatau_6)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Lapl2_1: ";
    for (auto v : v2lapl2_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Lapl2_2: ";
    for (auto v : v2lapl2_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Lapl2_3: ";
    for (auto v : v2lapl2_3)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2LaplTau_1: ";
    for (auto v : v2lapltau_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2LaplTau_2: ";
    for (auto v : v2lapltau_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2LaplTau_3: ";
    for (auto v : v2lapltau_3)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2LaplTau_4: ";
    for (auto v : v2lapltau_4)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Tau2_1: ";
    for (auto v : v2tau2_1)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Tau2_2: ";
    for (auto v : v2tau2_2)
        std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "MGGA V2Tau2_3: ";
    for (auto v : v2tau2_3)
        std::cout << v << " ";
}

  


    /////////////////////////////////////////////////////////////
    //meta-GGA START, up to the second derivative
    // meta-GGA Energy Density for spin-polarized systems
    /*
    std::vector<double> gga_exc(const std::vector<double> &rho_up, const std::vector<double> &rho_down, const std::vector<double> &sigma_1, const std::vector<double> &sigma_2, const std::vector<double> &sigma_3)
    {
        int np = rho_up.size();
        std::vector<double> rho(2 * np);
        std::vector<double> sigma(3 * np);
        for (int i = 0; i < np; ++i)
        {
            rho[2 * i] = rho_up[i];
            rho[2 * i + 1] = rho_down[i];
            sigma[3 * i] = sigma_1[i];
            sigma[3 * i + 1] = sigma_2[i];
            sigma[3 * i + 2] = sigma_3[i];
        }
        std::vector<double> exc(np);
        xc_gga_exc(&func, np, rho.data(), sigma.data(), exc.data());
        return exc;
    }
*/
/*

void LibxcInterface::test_ggas_with_lxc()
{
    std::cout << "Testing GGA functionals that support up to 4th derivatives (lxc)...\n";

    // 遍历所有可能的 xc_id
    for(int id = 1; id < 800; ++id)
    {
        xc_func_type localFunc;
        if(xc_func_init(&localFunc, id, XC_POLARIZED) == 0)
        {
            const xc_func_info_type* info = xc_func_get_info(&localFunc);
            if(info && info->family == XC_FAMILY_GGA)
            {
                // 检查是否支持 lxc，假设 n_deriv >= 4 表示支持
                if(info->n_deriv >= 4)
                {
                    // 进一步验证具体的高阶导数函数是否存在
                    bool has_kxc = (localFunc.kxc != NULL);
                    bool has_lxc = (localFunc.lxc != NULL);

                    if(has_kxc && has_lxc)
                    {
                        std::cout << "Functional ID: " << id
                                  << " (" << info->name << ")"
                                  << " supports up to " << info->n_deriv
                                  << " derivatives (has kxc and lxc)." << std::endl;
                    }
                    else
                    {
                        std::cout << "Functional ID: " << id
                                  << " (" << info->name << ")"
                                  << " claims n_deriv >=4 but lacks ";
                        if(!has_kxc) std::cout << "kxc ";
                        if(!has_lxc) std::cout << "lxc ";
                        std::cout << std::endl;
                    }
                }
            }
            xc_func_end(&localFunc);
        }
    }

    std::cout << "Finished testing GGA functionals.\n";
}

*/





//////////////////////////////////////////////////////////////////////////////////////////////////

// MAIN FOR INDENPENDENT TEST

/*

int main() 
{
    int vmajor, vminor, vmicro;
    //Get the libxc version 
    xc_version(&vmajor, &vminor, &vmicro);
    printf("Libxc version: %d.%d.%d\n", vmajor, vminor, vmicro);

    // //LDA  test 
   try 
   {
       LibxcInterface libxc(1, true); // Example with spin-polarized LDA (xc_id = 1)
       libxc.example_lda_spin();
   } 
   catch (const std::exception& ex) 
   {
       std::cerr << "Error: " << ex.what() << std::endl;
   }
    std::cout << "####################################################\n" << std::endl; 

    // //GGA  test
    try
    {
        LibxcInterface libxc(101, true); // Example with spin-polarized GGA (xc_id = 101)
        libxc.example_gga_spin();
        //libxc.test_ggas_with_lxc();
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
    std::cout << "####################################################\n" << std::endl;

 

    // //meta-GGA test
    try 
   {
        LibxcInterface libxc(202, true); // Example with spin-polarized meta-GGA (xc_id = 202)
        libxc.example_mgga_spin();
    } 
    catch (const std::exception& ex) 
    {
        std::cerr << "Error: " << ex.what() << std::endl;
    }

    return 0;
}

*/

//////////////////////////////////////////////////////////////////////////////////////////////////

