
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

namespace NCXC{

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
    
    // meta-GGA START, up to the fourth derivative
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

void LibxcInterface::mgga_kxc(
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
    std::vector<double> &v3tau3_1, std::vector<double> &v3tau3_2, std::vector<double> &v3tau3_3, std::vector<double> &v3tau3_4)
{
    int np = rho_up.size();
    std::vector<double> rho(2 * np), sigma(3 * np), lapl(2 * np), tau(2 * np);
    std::vector<double> v3rho3(4 * np), v3rho2sigma(9 * np), v3rho2lapl(6 * np), v3rho2tau(6 * np);
    std::vector<double> v3rhosigma2(12 * np), v3rhosigmalapl(12 * np), v3rhosigmatau(12 * np);
    std::vector<double> v3rholapl2(6 * np), v3rholapltau(8 * np), v3rhotau2(6 * np);
    std::vector<double> v3sigma3(10 * np), v3sigma2lapl(12 * np), v3sigma2tau(12 * np);
    std::vector<double> v3sigmalapl2(9 * np), v3sigmalapltau(12 * np), v3sigmatau2(9 * np);
    std::vector<double> v3lapl3(4 * np), v3lapl2tau(6 * np), v3lapltau2(6 * np), v3tau3(4 * np);

    for (int i = 0; i < np; ++i)
    {
        rho[2*i] = rho_up[i];
        rho[2*i+1] = rho_down[i];
        sigma[3*i] = sigma_1[i];
        sigma[3*i+1] = sigma_2[i];
        sigma[3*i+2] = sigma_3[i];
        lapl[2*i] = lapl_up[i];
        lapl[2*i+1] = lapl_down[i];
        tau[2*i] = tau_up[i];
        tau[2*i+1] = tau_down[i];
    }

    xc_mgga_kxc(&func, np, rho.data(), sigma.data(), lapl.data(), tau.data(),
                v3rho3.data(), v3rho2sigma.data(), v3rho2lapl.data(), v3rho2tau.data(),
                v3rhosigma2.data(), v3rhosigmalapl.data(), v3rhosigmatau.data(),
                v3rholapl2.data(), v3rholapltau.data(), v3rhotau2.data(),
                v3sigma3.data(), v3sigma2lapl.data(), v3sigma2tau.data(),
                v3sigmalapl2.data(), v3sigmalapltau.data(), v3sigmatau2.data(),
                v3lapl3.data(), v3lapl2tau.data(), v3lapltau2.data(), v3tau3.data());

    for (int i = 0; i < np; ++i)
{
    // 1) v3rho3 (4)
    v3rho3_1[i] = v3rho3[4*i + 0];
    v3rho3_2[i] = v3rho3[4*i + 1];
    v3rho3_3[i] = v3rho3[4*i + 2];
    v3rho3_4[i] = v3rho3[4*i + 3];

    // 2) v3rho2sigma (9)
    v3rho2sigma_1[i] = v3rho2sigma[9*i + 0];
    v3rho2sigma_2[i] = v3rho2sigma[9*i + 1];
    v3rho2sigma_3[i] = v3rho2sigma[9*i + 2];
    v3rho2sigma_4[i] = v3rho2sigma[9*i + 3];
    v3rho2sigma_5[i] = v3rho2sigma[9*i + 4];
    v3rho2sigma_6[i] = v3rho2sigma[9*i + 5];
    v3rho2sigma_7[i] = v3rho2sigma[9*i + 6];
    v3rho2sigma_8[i] = v3rho2sigma[9*i + 7];
    v3rho2sigma_9[i] = v3rho2sigma[9*i + 8];

    // 3) v3rho2lapl (6)
    v3rho2lapl_1[i] = v3rho2lapl[6*i + 0];
    v3rho2lapl_2[i] = v3rho2lapl[6*i + 1];
    v3rho2lapl_3[i] = v3rho2lapl[6*i + 2];
    v3rho2lapl_4[i] = v3rho2lapl[6*i + 3];
    v3rho2lapl_5[i] = v3rho2lapl[6*i + 4];
    v3rho2lapl_6[i] = v3rho2lapl[6*i + 5];

    // 4) v3rho2tau (6)
    v3rho2tau_1[i] = v3rho2tau[6*i + 0];
    v3rho2tau_2[i] = v3rho2tau[6*i + 1];
    v3rho2tau_3[i] = v3rho2tau[6*i + 2];
    v3rho2tau_4[i] = v3rho2tau[6*i + 3];
    v3rho2tau_5[i] = v3rho2tau[6*i + 4];
    v3rho2tau_6[i] = v3rho2tau[6*i + 5];

    // 5) v3rhosigma2 (12)
    v3rhosigma2_1[i]  = v3rhosigma2[12*i + 0];
    v3rhosigma2_2[i]  = v3rhosigma2[12*i + 1];
    v3rhosigma2_3[i]  = v3rhosigma2[12*i + 2];
    v3rhosigma2_4[i]  = v3rhosigma2[12*i + 3];
    v3rhosigma2_5[i]  = v3rhosigma2[12*i + 4];
    v3rhosigma2_6[i]  = v3rhosigma2[12*i + 5];
    v3rhosigma2_7[i]  = v3rhosigma2[12*i + 6];
    v3rhosigma2_8[i]  = v3rhosigma2[12*i + 7];
    v3rhosigma2_9[i]  = v3rhosigma2[12*i + 8];
    v3rhosigma2_10[i] = v3rhosigma2[12*i + 9];
    v3rhosigma2_11[i] = v3rhosigma2[12*i + 10];
    v3rhosigma2_12[i] = v3rhosigma2[12*i + 11];

    // 6) v3rhosigmalapl (12)
    v3rhosigmalapl_1[i]  = v3rhosigmalapl[12*i + 0];
    v3rhosigmalapl_2[i]  = v3rhosigmalapl[12*i + 1];
    v3rhosigmalapl_3[i]  = v3rhosigmalapl[12*i + 2];
    v3rhosigmalapl_4[i]  = v3rhosigmalapl[12*i + 3];
    v3rhosigmalapl_5[i]  = v3rhosigmalapl[12*i + 4];
    v3rhosigmalapl_6[i]  = v3rhosigmalapl[12*i + 5];
    v3rhosigmalapl_7[i]  = v3rhosigmalapl[12*i + 6];
    v3rhosigmalapl_8[i]  = v3rhosigmalapl[12*i + 7];
    v3rhosigmalapl_9[i]  = v3rhosigmalapl[12*i + 8];
    v3rhosigmalapl_10[i] = v3rhosigmalapl[12*i + 9];
    v3rhosigmalapl_11[i] = v3rhosigmalapl[12*i + 10];
    v3rhosigmalapl_12[i] = v3rhosigmalapl[12*i + 11];

    // 7) v3rhosigmatau (12)
    v3rhosigmatau_1[i]  = v3rhosigmatau[12*i + 0];
    v3rhosigmatau_2[i]  = v3rhosigmatau[12*i + 1];
    v3rhosigmatau_3[i]  = v3rhosigmatau[12*i + 2];
    v3rhosigmatau_4[i]  = v3rhosigmatau[12*i + 3];
    v3rhosigmatau_5[i]  = v3rhosigmatau[12*i + 4];
    v3rhosigmatau_6[i]  = v3rhosigmatau[12*i + 5];
    v3rhosigmatau_7[i]  = v3rhosigmatau[12*i + 6];
    v3rhosigmatau_8[i]  = v3rhosigmatau[12*i + 7];
    v3rhosigmatau_9[i]  = v3rhosigmatau[12*i + 8];
    v3rhosigmatau_10[i] = v3rhosigmatau[12*i + 9];
    v3rhosigmatau_11[i] = v3rhosigmatau[12*i + 10];
    v3rhosigmatau_12[i] = v3rhosigmatau[12*i + 11];

    // 8) v3rholapl2 (6)
    v3rholapl2_1[i] = v3rholapl2[6*i + 0];
    v3rholapl2_2[i] = v3rholapl2[6*i + 1];
    v3rholapl2_3[i] = v3rholapl2[6*i + 2];
    v3rholapl2_4[i] = v3rholapl2[6*i + 3];
    v3rholapl2_5[i] = v3rholapl2[6*i + 4];
    v3rholapl2_6[i] = v3rholapl2[6*i + 5];

    // 9) v3rholapltau (8)
    v3rholapltau_1[i] = v3rholapltau[8*i + 0];
    v3rholapltau_2[i] = v3rholapltau[8*i + 1];
    v3rholapltau_3[i] = v3rholapltau[8*i + 2];
    v3rholapltau_4[i] = v3rholapltau[8*i + 3];
    v3rholapltau_5[i] = v3rholapltau[8*i + 4];
    v3rholapltau_6[i] = v3rholapltau[8*i + 5];
    v3rholapltau_7[i] = v3rholapltau[8*i + 6];
    v3rholapltau_8[i] = v3rholapltau[8*i + 7];

    // 10) v3rhotau2 (6)
    v3rhotau2_1[i] = v3rhotau2[6*i + 0];
    v3rhotau2_2[i] = v3rhotau2[6*i + 1];
    v3rhotau2_3[i] = v3rhotau2[6*i + 2];
    v3rhotau2_4[i] = v3rhotau2[6*i + 3];
    v3rhotau2_5[i] = v3rhotau2[6*i + 4];
    v3rhotau2_6[i] = v3rhotau2[6*i + 5];

    // 11) v3sigma3 (10)
    v3sigma3_1[i] = v3sigma3[10*i + 0];
    v3sigma3_2[i] = v3sigma3[10*i + 1];
    v3sigma3_3[i] = v3sigma3[10*i + 2];
    v3sigma3_4[i] = v3sigma3[10*i + 3];
    v3sigma3_5[i] = v3sigma3[10*i + 4];
    v3sigma3_6[i] = v3sigma3[10*i + 5];
    v3sigma3_7[i] = v3sigma3[10*i + 6];
    v3sigma3_8[i] = v3sigma3[10*i + 7];
    v3sigma3_9[i] = v3sigma3[10*i + 8];
    v3sigma3_10[i] = v3sigma3[10*i + 9];

    // 12) v3sigma2lapl (12)
    v3sigma2lapl_1[i]  = v3sigma2lapl[12*i + 0];
    v3sigma2lapl_2[i]  = v3sigma2lapl[12*i + 1];
    v3sigma2lapl_3[i]  = v3sigma2lapl[12*i + 2];
    v3sigma2lapl_4[i]  = v3sigma2lapl[12*i + 3];
    v3sigma2lapl_5[i]  = v3sigma2lapl[12*i + 4];
    v3sigma2lapl_6[i]  = v3sigma2lapl[12*i + 5];
    v3sigma2lapl_7[i]  = v3sigma2lapl[12*i + 6];
    v3sigma2lapl_8[i]  = v3sigma2lapl[12*i + 7];
    v3sigma2lapl_9[i]  = v3sigma2lapl[12*i + 8];
    v3sigma2lapl_10[i] = v3sigma2lapl[12*i + 9];
    v3sigma2lapl_11[i] = v3sigma2lapl[12*i + 10];
    v3sigma2lapl_12[i] = v3sigma2lapl[12*i + 11];

    // 13) v3sigma2tau (12)
    v3sigma2tau_1[i]  = v3sigma2tau[12*i + 0];
    v3sigma2tau_2[i]  = v3sigma2tau[12*i + 1];
    v3sigma2tau_3[i]  = v3sigma2tau[12*i + 2];
    v3sigma2tau_4[i]  = v3sigma2tau[12*i + 3];
    v3sigma2tau_5[i]  = v3sigma2tau[12*i + 4];
    v3sigma2tau_6[i]  = v3sigma2tau[12*i + 5];
    v3sigma2tau_7[i]  = v3sigma2tau[12*i + 6];
    v3sigma2tau_8[i]  = v3sigma2tau[12*i + 7];
    v3sigma2tau_9[i]  = v3sigma2tau[12*i + 8];
    v3sigma2tau_10[i] = v3sigma2tau[12*i + 9];
    v3sigma2tau_11[i] = v3sigma2tau[12*i + 10];
    v3sigma2tau_12[i] = v3sigma2tau[12*i + 11];

    // 14) v3sigmalapl2 (9)
    v3sigmalapl2_1[i] = v3sigmalapl2[9*i + 0];
    v3sigmalapl2_2[i] = v3sigmalapl2[9*i + 1];
    v3sigmalapl2_3[i] = v3sigmalapl2[9*i + 2];
    v3sigmalapl2_4[i] = v3sigmalapl2[9*i + 3];
    v3sigmalapl2_5[i] = v3sigmalapl2[9*i + 4];
    v3sigmalapl2_6[i] = v3sigmalapl2[9*i + 5];
    v3sigmalapl2_7[i] = v3sigmalapl2[9*i + 6];
    v3sigmalapl2_8[i] = v3sigmalapl2[9*i + 7];
    v3sigmalapl2_9[i] = v3sigmalapl2[9*i + 8];

    // 15) v3sigmalapltau (12)
    v3sigmalapltau_1[i]  = v3sigmalapltau[12*i + 0];
    v3sigmalapltau_2[i]  = v3sigmalapltau[12*i + 1];
    v3sigmalapltau_3[i]  = v3sigmalapltau[12*i + 2];
    v3sigmalapltau_4[i]  = v3sigmalapltau[12*i + 3];
    v3sigmalapltau_5[i]  = v3sigmalapltau[12*i + 4];
    v3sigmalapltau_6[i]  = v3sigmalapltau[12*i + 5];
    v3sigmalapltau_7[i]  = v3sigmalapltau[12*i + 6];
    v3sigmalapltau_8[i]  = v3sigmalapltau[12*i + 7];
    v3sigmalapltau_9[i]  = v3sigmalapltau[12*i + 8];
    v3sigmalapltau_10[i] = v3sigmalapltau[12*i + 9];
    v3sigmalapltau_11[i] = v3sigmalapltau[12*i + 10];
    v3sigmalapltau_12[i] = v3sigmalapltau[12*i + 11];

    // 16) v3sigmatau2 (9)
    v3sigmatau2_1[i] = v3sigmatau2[9*i + 0];
    v3sigmatau2_2[i] = v3sigmatau2[9*i + 1];
    v3sigmatau2_3[i] = v3sigmatau2[9*i + 2];
    v3sigmatau2_4[i] = v3sigmatau2[9*i + 3];
    v3sigmatau2_5[i] = v3sigmatau2[9*i + 4];
    v3sigmatau2_6[i] = v3sigmatau2[9*i + 5];
    v3sigmatau2_7[i] = v3sigmatau2[9*i + 6];
    v3sigmatau2_8[i] = v3sigmatau2[9*i + 7];
    v3sigmatau2_9[i] = v3sigmatau2[9*i + 8];

    // 17) v3lapl3 (4)
    v3lapl3_1[i] = v3lapl3[4*i + 0];
    v3lapl3_2[i] = v3lapl3[4*i + 1];
    v3lapl3_3[i] = v3lapl3[4*i + 2];
    v3lapl3_4[i] = v3lapl3[4*i + 3];

    // 18) v3lapl2tau (6)
    v3lapl2tau_1[i] = v3lapl2tau[6*i + 0];
    v3lapl2tau_2[i] = v3lapl2tau[6*i + 1];
    v3lapl2tau_3[i] = v3lapl2tau[6*i + 2];
    v3lapl2tau_4[i] = v3lapl2tau[6*i + 3];
    v3lapl2tau_5[i] = v3lapl2tau[6*i + 4];
    v3lapl2tau_6[i] = v3lapl2tau[6*i + 5];

    // 19) v3lapltau2 (6)
    v3lapltau2_1[i] = v3lapltau2[6*i + 0];
    v3lapltau2_2[i] = v3lapltau2[6*i + 1];
    v3lapltau2_3[i] = v3lapltau2[6*i + 2];
    v3lapltau2_4[i] = v3lapltau2[6*i + 3];
    v3lapltau2_5[i] = v3lapltau2[6*i + 4];
    v3lapltau2_6[i] = v3lapltau2[6*i + 5];

    // 20) v3tau3 (4)
    v3tau3_1[i] = v3tau3[4*i + 0];
    v3tau3_2[i] = v3tau3[4*i + 1];
    v3tau3_3[i] = v3tau3[4*i + 2];
    v3tau3_4[i] = v3tau3[4*i + 3];
}
}



void LibxcInterface::mgga_lxc(
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
)
{
    int np = static_cast<int>(rho_up.size());
    // 准备输入数组
    std::vector<double> rho(2 * np), sigma(3 * np), lapl(2 * np), tau(2 * np);
    for(int i = 0; i < np; ++i)
    {
        rho[2*i]      = rho_up[i];
        rho[2*i + 1]  = rho_down[i];
        sigma[3*i]    = sigma_1[i];
        sigma[3*i+1]  = sigma_2[i];
        sigma[3*i+2]  = sigma_3[i];
        lapl[2*i]     = lapl_up[i];
        lapl[2*i + 1] = lapl_down[i];
        tau[2*i]      = tau_up[i];
        tau[2*i + 1]  = tau_down[i];
    }

    std::vector<double> v4rho4(5 * np), v4rho3sigma(12 * np),
        v4rho3lapl(8 * np), v4rho3tau(8 * np), v4rho2sigma2(18 * np),
        v4rho2sigmalapl(18 * np), v4rho2sigmatau(18 * np), v4rho2lapl2(9 * np),
        v4rho2lapltau(12 * np), v4rho2tau2(9 * np), v4rhosigma3(20 * np),
        v4rhosigma2lapl(24 * np), v4rhosigma2tau(24 * np),
        v4rhosigmalapl2(18 * np), v4rhosigmalapltau(24 * np),
        v4rhosigmatau2(18 * np), v4rholapl3(8 * np), v4rholapl2tau(12 * np),
        v4rholapltau2(12 * np), v4rhotau3(8 * np), v4sigma4(15 * np),
        v4sigma3lapl(20 * np), v4sigma3tau(20 * np),
        v4sigma2lapl2(18 * np), v4sigma2lapltau(24 * np),
        v4sigma2tau2(18 * np), v4sigmalapl3(12 * np),
        v4sigmalapl2tau(18 * np), v4sigmalapltau2(18 * np),
        v4sigmatau3(12 * np), v4lapl4(5 * np), v4lapl3tau(8 * np),
        v4lapl2tau2(9 * np), v4lapltau3(8 * np), v4tau4(5 * np);

    xc_mgga_lxc(&func, np, rho.data(), sigma.data(), lapl.data(), tau.data(),
                v4rho4.data(), v4rho3sigma.data(), v4rho3lapl.data(), v4rho3tau.data(),
                v4rho2sigma2.data(), v4rho2sigmalapl.data(), v4rho2sigmatau.data(),
                v4rho2lapl2.data(), v4rho2lapltau.data(), v4rho2tau2.data(),
                v4rhosigma3.data(), v4rhosigma2lapl.data(), v4rhosigma2tau.data(),
                v4rhosigmalapl2.data(), v4rhosigmalapltau.data(), v4rhosigmatau2.data(),
                v4rholapl3.data(), v4rholapl2tau.data(), v4rholapltau2.data(),
                v4rhotau3.data(), v4sigma4.data(), v4sigma3lapl.data(), v4sigma3tau.data(),
                v4sigma2lapl2.data(), v4sigma2lapltau.data(), v4sigma2tau2.data(),
                v4sigmalapl3.data(), v4sigmalapl2tau.data(), v4sigmalapltau2.data(),
                v4sigmatau3.data(), v4lapl4.data(), v4lapl3tau.data(), v4lapl2tau2.data(),
                v4lapltau3.data(), v4tau4.data());

    for(int i = 0; i < np; ++i)
    {
        // v4rho4: 5*np
        v4rho4_1[i] = v4rho4[5 * i + 0];
        v4rho4_2[i] = v4rho4[5 * i + 1];
        v4rho4_3[i] = v4rho4[5 * i + 2];
        v4rho4_4[i] = v4rho4[5 * i + 3];
        v4rho4_5[i] = v4rho4[5 * i + 4];

        // v4rho3sigma: 12*np
        v4rho3sigma_1[i]  = v4rho3sigma[12 * i + 0];
        v4rho3sigma_2[i]  = v4rho3sigma[12 * i + 1];
        v4rho3sigma_3[i]  = v4rho3sigma[12 * i + 2];
        v4rho3sigma_4[i]  = v4rho3sigma[12 * i + 3];
        v4rho3sigma_5[i]  = v4rho3sigma[12 * i + 4];
        v4rho3sigma_6[i]  = v4rho3sigma[12 * i + 5];
        v4rho3sigma_7[i]  = v4rho3sigma[12 * i + 6];
        v4rho3sigma_8[i]  = v4rho3sigma[12 * i + 7];
        v4rho3sigma_9[i]  = v4rho3sigma[12 * i + 8];
        v4rho3sigma_10[i] = v4rho3sigma[12 * i + 9];
        v4rho3sigma_11[i] = v4rho3sigma[12 * i + 10];
        v4rho3sigma_12[i] = v4rho3sigma[12 * i + 11];

        // v4rho3lapl: 8*np
        v4rho3lapl_1[i] = v4rho3lapl[8 * i + 0];
        v4rho3lapl_2[i] = v4rho3lapl[8 * i + 1];
        v4rho3lapl_3[i] = v4rho3lapl[8 * i + 2];
        v4rho3lapl_4[i] = v4rho3lapl[8 * i + 3];
        v4rho3lapl_5[i] = v4rho3lapl[8 * i + 4];
        v4rho3lapl_6[i] = v4rho3lapl[8 * i + 5];
        v4rho3lapl_7[i] = v4rho3lapl[8 * i + 6];
        v4rho3lapl_8[i] = v4rho3lapl[8 * i + 7];

        // v4rho3tau: 8*np
        v4rho3tau_1[i] = v4rho3tau[8 * i + 0];
        v4rho3tau_2[i] = v4rho3tau[8 * i + 1];
        v4rho3tau_3[i] = v4rho3tau[8 * i + 2];
        v4rho3tau_4[i] = v4rho3tau[8 * i + 3];
        v4rho3tau_5[i] = v4rho3tau[8 * i + 4];
        v4rho3tau_6[i] = v4rho3tau[8 * i + 5];
        v4rho3tau_7[i] = v4rho3tau[8 * i + 6];
        v4rho3tau_8[i] = v4rho3tau[8 * i + 7];

        // v4rho2sigma2: 18*np
        v4rho2sigma2_1[i]  = v4rho2sigma2[18 * i + 0];
        v4rho2sigma2_2[i]  = v4rho2sigma2[18 * i + 1];
        v4rho2sigma2_3[i]  = v4rho2sigma2[18 * i + 2];
        v4rho2sigma2_4[i]  = v4rho2sigma2[18 * i + 3];
        v4rho2sigma2_5[i]  = v4rho2sigma2[18 * i + 4];
        v4rho2sigma2_6[i]  = v4rho2sigma2[18 * i + 5];
        v4rho2sigma2_7[i]  = v4rho2sigma2[18 * i + 6];
        v4rho2sigma2_8[i]  = v4rho2sigma2[18 * i + 7];
        v4rho2sigma2_9[i]  = v4rho2sigma2[18 * i + 8];
        v4rho2sigma2_10[i] = v4rho2sigma2[18 * i + 9];
        v4rho2sigma2_11[i] = v4rho2sigma2[18 * i + 10];
        v4rho2sigma2_12[i] = v4rho2sigma2[18 * i + 11];
        v4rho2sigma2_13[i] = v4rho2sigma2[18 * i + 12];
        v4rho2sigma2_14[i] = v4rho2sigma2[18 * i + 13];
        v4rho2sigma2_15[i] = v4rho2sigma2[18 * i + 14];
        v4rho2sigma2_16[i] = v4rho2sigma2[18 * i + 15];
        v4rho2sigma2_17[i] = v4rho2sigma2[18 * i + 16];
        v4rho2sigma2_18[i] = v4rho2sigma2[18 * i + 17];

        // v4rho2sigmalapl: 18*np
        v4rho2sigmalapl_1[i]  = v4rho2sigmalapl[18 * i + 0];
        v4rho2sigmalapl_2[i]  = v4rho2sigmalapl[18 * i + 1];
        v4rho2sigmalapl_3[i]  = v4rho2sigmalapl[18 * i + 2];
        v4rho2sigmalapl_4[i]  = v4rho2sigmalapl[18 * i + 3];
        v4rho2sigmalapl_5[i]  = v4rho2sigmalapl[18 * i + 4];
        v4rho2sigmalapl_6[i]  = v4rho2sigmalapl[18 * i + 5];
        v4rho2sigmalapl_7[i]  = v4rho2sigmalapl[18 * i + 6];
        v4rho2sigmalapl_8[i]  = v4rho2sigmalapl[18 * i + 7];
        v4rho2sigmalapl_9[i]  = v4rho2sigmalapl[18 * i + 8];
        v4rho2sigmalapl_10[i] = v4rho2sigmalapl[18 * i + 9];
        v4rho2sigmalapl_11[i] = v4rho2sigmalapl[18 * i + 10];
        v4rho2sigmalapl_12[i] = v4rho2sigmalapl[18 * i + 11];
        v4rho2sigmalapl_13[i] = v4rho2sigmalapl[18 * i + 12];
        v4rho2sigmalapl_14[i] = v4rho2sigmalapl[18 * i + 13];
        v4rho2sigmalapl_15[i] = v4rho2sigmalapl[18 * i + 14];
        v4rho2sigmalapl_16[i] = v4rho2sigmalapl[18 * i + 15];
        v4rho2sigmalapl_17[i] = v4rho2sigmalapl[18 * i + 16];
        v4rho2sigmalapl_18[i] = v4rho2sigmalapl[18 * i + 17];

        // v4rho2sigmatau: 18*np
        v4rho2sigmatau_1[i]  = v4rho2sigmatau[18 * i + 0];
        v4rho2sigmatau_2[i]  = v4rho2sigmatau[18 * i + 1];
        v4rho2sigmatau_3[i]  = v4rho2sigmatau[18 * i + 2];
        v4rho2sigmatau_4[i]  = v4rho2sigmatau[18 * i + 3];
        v4rho2sigmatau_5[i]  = v4rho2sigmatau[18 * i + 4];
        v4rho2sigmatau_6[i]  = v4rho2sigmatau[18 * i + 5];
        v4rho2sigmatau_7[i]  = v4rho2sigmatau[18 * i + 6];
        v4rho2sigmatau_8[i]  = v4rho2sigmatau[18 * i + 7];
        v4rho2sigmatau_9[i]  = v4rho2sigmatau[18 * i + 8];
        v4rho2sigmatau_10[i] = v4rho2sigmatau[18 * i + 9];
        v4rho2sigmatau_11[i] = v4rho2sigmatau[18 * i + 10];
        v4rho2sigmatau_12[i] = v4rho2sigmatau[18 * i + 11];
        v4rho2sigmatau_13[i] = v4rho2sigmatau[18 * i + 12];
        v4rho2sigmatau_14[i] = v4rho2sigmatau[18 * i + 13];
        v4rho2sigmatau_15[i] = v4rho2sigmatau[18 * i + 14];
        v4rho2sigmatau_16[i] = v4rho2sigmatau[18 * i + 15];
        v4rho2sigmatau_17[i] = v4rho2sigmatau[18 * i + 16];
        v4rho2sigmatau_18[i] = v4rho2sigmatau[18 * i + 17];

        // v4rho2lapl2: 9*np
        v4rho2lapl2_1[i] = v4rho2lapl2[9 * i + 0];
        v4rho2lapl2_2[i] = v4rho2lapl2[9 * i + 1];
        v4rho2lapl2_3[i] = v4rho2lapl2[9 * i + 2];
        v4rho2lapl2_4[i] = v4rho2lapl2[9 * i + 3];
        v4rho2lapl2_5[i] = v4rho2lapl2[9 * i + 4];
        v4rho2lapl2_6[i] = v4rho2lapl2[9 * i + 5];
        v4rho2lapl2_7[i] = v4rho2lapl2[9 * i + 6];
        v4rho2lapl2_8[i] = v4rho2lapl2[9 * i + 7];
        v4rho2lapl2_9[i] = v4rho2lapl2[9 * i + 8];

        // v4rho2lapltau: 12*np
        v4rho2lapltau_1[i]  = v4rho2lapltau[12 * i + 0];
        v4rho2lapltau_2[i]  = v4rho2lapltau[12 * i + 1];
        v4rho2lapltau_3[i]  = v4rho2lapltau[12 * i + 2];
        v4rho2lapltau_4[i]  = v4rho2lapltau[12 * i + 3];
        v4rho2lapltau_5[i]  = v4rho2lapltau[12 * i + 4];
        v4rho2lapltau_6[i]  = v4rho2lapltau[12 * i + 5];
        v4rho2lapltau_7[i]  = v4rho2lapltau[12 * i + 6];
        v4rho2lapltau_8[i]  = v4rho2lapltau[12 * i + 7];
        v4rho2lapltau_9[i]  = v4rho2lapltau[12 * i + 8];
        v4rho2lapltau_10[i] = v4rho2lapltau[12 * i + 9];
        v4rho2lapltau_11[i] = v4rho2lapltau[12 * i + 10];
        v4rho2lapltau_12[i] = v4rho2lapltau[12 * i + 11];

        // v4rho2tau2: 9*np
        v4rho2tau2_1[i] = v4rho2tau2[9 * i + 0];
        v4rho2tau2_2[i] = v4rho2tau2[9 * i + 1];
        v4rho2tau2_3[i] = v4rho2tau2[9 * i + 2];
        v4rho2tau2_4[i] = v4rho2tau2[9 * i + 3];
        v4rho2tau2_5[i] = v4rho2tau2[9 * i + 4];
        v4rho2tau2_6[i] = v4rho2tau2[9 * i + 5];
        v4rho2tau2_7[i] = v4rho2tau2[9 * i + 6];
        v4rho2tau2_8[i] = v4rho2tau2[9 * i + 7];
        v4rho2tau2_9[i] = v4rho2tau2[9 * i + 8];

        // v4rhosigma3: 20*np
        v4rhosigma3_1[i]  = v4rhosigma3[20 * i + 0];
        v4rhosigma3_2[i]  = v4rhosigma3[20 * i + 1];
        v4rhosigma3_3[i]  = v4rhosigma3[20 * i + 2];
        v4rhosigma3_4[i]  = v4rhosigma3[20 * i + 3];
        v4rhosigma3_5[i]  = v4rhosigma3[20 * i + 4];
        v4rhosigma3_6[i]  = v4rhosigma3[20 * i + 5];
        v4rhosigma3_7[i]  = v4rhosigma3[20 * i + 6];
        v4rhosigma3_8[i]  = v4rhosigma3[20 * i + 7];
        v4rhosigma3_9[i]  = v4rhosigma3[20 * i + 8];
        v4rhosigma3_10[i] = v4rhosigma3[20 * i + 9];
        v4rhosigma3_11[i] = v4rhosigma3[20 * i + 10];
        v4rhosigma3_12[i] = v4rhosigma3[20 * i + 11];
        v4rhosigma3_13[i] = v4rhosigma3[20 * i + 12];
        v4rhosigma3_14[i] = v4rhosigma3[20 * i + 13];
        v4rhosigma3_15[i] = v4rhosigma3[20 * i + 14];
        v4rhosigma3_16[i] = v4rhosigma3[20 * i + 15];
        v4rhosigma3_17[i] = v4rhosigma3[20 * i + 16];
        v4rhosigma3_18[i] = v4rhosigma3[20 * i + 17];
        v4rhosigma3_19[i] = v4rhosigma3[20 * i + 18];
        v4rhosigma3_20[i] = v4rhosigma3[20 * i + 19];

        // v4rhosigma2lapl: 24*np
        v4rhosigma2lapl_1[i]  = v4rhosigma2lapl[24 * i + 0];
        v4rhosigma2lapl_2[i]  = v4rhosigma2lapl[24 * i + 1];
        v4rhosigma2lapl_3[i]  = v4rhosigma2lapl[24 * i + 2];
        v4rhosigma2lapl_4[i]  = v4rhosigma2lapl[24 * i + 3];
        v4rhosigma2lapl_5[i]  = v4rhosigma2lapl[24 * i + 4];
        v4rhosigma2lapl_6[i]  = v4rhosigma2lapl[24 * i + 5];
        v4rhosigma2lapl_7[i]  = v4rhosigma2lapl[24 * i + 6];
        v4rhosigma2lapl_8[i]  = v4rhosigma2lapl[24 * i + 7];
        v4rhosigma2lapl_9[i]  = v4rhosigma2lapl[24 * i + 8];
        v4rhosigma2lapl_10[i] = v4rhosigma2lapl[24 * i + 9];
        v4rhosigma2lapl_11[i] = v4rhosigma2lapl[24 * i + 10];
        v4rhosigma2lapl_12[i] = v4rhosigma2lapl[24 * i + 11];
        v4rhosigma2lapl_13[i] = v4rhosigma2lapl[24 * i + 12];
        v4rhosigma2lapl_14[i] = v4rhosigma2lapl[24 * i + 13];
        v4rhosigma2lapl_15[i] = v4rhosigma2lapl[24 * i + 14];
        v4rhosigma2lapl_16[i] = v4rhosigma2lapl[24 * i + 15];
        v4rhosigma2lapl_17[i] = v4rhosigma2lapl[24 * i + 16];
        v4rhosigma2lapl_18[i] = v4rhosigma2lapl[24 * i + 17];
        v4rhosigma2lapl_19[i] = v4rhosigma2lapl[24 * i + 18];
        v4rhosigma2lapl_20[i] = v4rhosigma2lapl[24 * i + 19];
        v4rhosigma2lapl_21[i] = v4rhosigma2lapl[24 * i + 20];
        v4rhosigma2lapl_22[i] = v4rhosigma2lapl[24 * i + 21];
        v4rhosigma2lapl_23[i] = v4rhosigma2lapl[24 * i + 22];
        v4rhosigma2lapl_24[i] = v4rhosigma2lapl[24 * i + 23];

        // v4rhosigma2tau: 24*np
        v4rhosigma2tau_1[i]  = v4rhosigma2tau[24 * i + 0];
        v4rhosigma2tau_2[i]  = v4rhosigma2tau[24 * i + 1];
        v4rhosigma2tau_3[i]  = v4rhosigma2tau[24 * i + 2];
        v4rhosigma2tau_4[i]  = v4rhosigma2tau[24 * i + 3];
        v4rhosigma2tau_5[i]  = v4rhosigma2tau[24 * i + 4];
        v4rhosigma2tau_6[i]  = v4rhosigma2tau[24 * i + 5];
        v4rhosigma2tau_7[i]  = v4rhosigma2tau[24 * i + 6];
        v4rhosigma2tau_8[i]  = v4rhosigma2tau[24 * i + 7];
        v4rhosigma2tau_9[i]  = v4rhosigma2tau[24 * i + 8];
    v4rhosigma2tau_10[i] = v4rhosigma2tau[24 * i + 9];
    v4rhosigma2tau_11[i] = v4rhosigma2tau[24 * i + 10];
    v4rhosigma2tau_12[i] = v4rhosigma2tau[24 * i + 11];
    v4rhosigma2tau_13[i] = v4rhosigma2tau[24 * i + 12];
    v4rhosigma2tau_14[i] = v4rhosigma2tau[24 * i + 13];
    v4rhosigma2tau_15[i] = v4rhosigma2tau[24 * i + 14];
    v4rhosigma2tau_16[i] = v4rhosigma2tau[24 * i + 15];
    v4rhosigma2tau_17[i] = v4rhosigma2tau[24 * i + 16];
    v4rhosigma2tau_18[i] = v4rhosigma2tau[24 * i + 17];
    v4rhosigma2tau_19[i] = v4rhosigma2tau[24 * i + 18];
    v4rhosigma2tau_20[i] = v4rhosigma2tau[24 * i + 19];
    v4rhosigma2tau_21[i] = v4rhosigma2tau[24 * i + 20];
    v4rhosigma2tau_22[i] = v4rhosigma2tau[24 * i + 21];
    v4rhosigma2tau_23[i] = v4rhosigma2tau[24 * i + 22];
    v4rhosigma2tau_24[i] = v4rhosigma2tau[24 * i + 23];

    // --- v4rhosigmalapl2: 18 * np ---
    v4rhosigmalapl2_1[i]  = v4rhosigmalapl2[18 * i + 0];
    v4rhosigmalapl2_2[i]  = v4rhosigmalapl2[18 * i + 1];
    v4rhosigmalapl2_3[i]  = v4rhosigmalapl2[18 * i + 2];
    v4rhosigmalapl2_4[i]  = v4rhosigmalapl2[18 * i + 3];
    v4rhosigmalapl2_5[i]  = v4rhosigmalapl2[18 * i + 4];
    v4rhosigmalapl2_6[i]  = v4rhosigmalapl2[18 * i + 5];
    v4rhosigmalapl2_7[i]  = v4rhosigmalapl2[18 * i + 6];
    v4rhosigmalapl2_8[i]  = v4rhosigmalapl2[18 * i + 7];
    v4rhosigmalapl2_9[i]  = v4rhosigmalapl2[18 * i + 8];
    v4rhosigmalapl2_10[i] = v4rhosigmalapl2[18 * i + 9];
    v4rhosigmalapl2_11[i] = v4rhosigmalapl2[18 * i + 10];
    v4rhosigmalapl2_12[i] = v4rhosigmalapl2[18 * i + 11];
    v4rhosigmalapl2_13[i] = v4rhosigmalapl2[18 * i + 12];
    v4rhosigmalapl2_14[i] = v4rhosigmalapl2[18 * i + 13];
    v4rhosigmalapl2_15[i] = v4rhosigmalapl2[18 * i + 14];
    v4rhosigmalapl2_16[i] = v4rhosigmalapl2[18 * i + 15];
    v4rhosigmalapl2_17[i] = v4rhosigmalapl2[18 * i + 16];
    v4rhosigmalapl2_18[i] = v4rhosigmalapl2[18 * i + 17];

    // --- v4rhosigmalapltau: 24*np ---
    v4rhosigmalapltau_1[i]  = v4rhosigmalapltau[24 * i + 0];
    v4rhosigmalapltau_2[i]  = v4rhosigmalapltau[24 * i + 1];
    v4rhosigmalapltau_3[i]  = v4rhosigmalapltau[24 * i + 2];
    v4rhosigmalapltau_4[i]  = v4rhosigmalapltau[24 * i + 3];
    v4rhosigmalapltau_5[i]  = v4rhosigmalapltau[24 * i + 4];
    v4rhosigmalapltau_6[i]  = v4rhosigmalapltau[24 * i + 5];
    v4rhosigmalapltau_7[i]  = v4rhosigmalapltau[24 * i + 6];
    v4rhosigmalapltau_8[i]  = v4rhosigmalapltau[24 * i + 7];
    v4rhosigmalapltau_9[i]  = v4rhosigmalapltau[24 * i + 8];
    v4rhosigmalapltau_10[i] = v4rhosigmalapltau[24 * i + 9];
    v4rhosigmalapltau_11[i] = v4rhosigmalapltau[24 * i + 10];
    v4rhosigmalapltau_12[i] = v4rhosigmalapltau[24 * i + 11];
    v4rhosigmalapltau_13[i] = v4rhosigmalapltau[24 * i + 12];
    v4rhosigmalapltau_14[i] = v4rhosigmalapltau[24 * i + 13];
    v4rhosigmalapltau_15[i] = v4rhosigmalapltau[24 * i + 14];
    v4rhosigmalapltau_16[i] = v4rhosigmalapltau[24 * i + 15];
    v4rhosigmalapltau_17[i] = v4rhosigmalapltau[24 * i + 16];
    v4rhosigmalapltau_18[i] = v4rhosigmalapltau[24 * i + 17];
    v4rhosigmalapltau_19[i] = v4rhosigmalapltau[24 * i + 18];
    v4rhosigmalapltau_20[i] = v4rhosigmalapltau[24 * i + 19];
    v4rhosigmalapltau_21[i] = v4rhosigmalapltau[24 * i + 20];
    v4rhosigmalapltau_22[i] = v4rhosigmalapltau[24 * i + 21];
    v4rhosigmalapltau_23[i] = v4rhosigmalapltau[24 * i + 22];
    v4rhosigmalapltau_24[i] = v4rhosigmalapltau[24 * i + 23];

    // --- v4rhosigmatau2: 18*np ---
    v4rhosigmatau2_1[i]  = v4rhosigmatau2[18 * i + 0];
    v4rhosigmatau2_2[i]  = v4rhosigmatau2[18 * i + 1];
    v4rhosigmatau2_3[i]  = v4rhosigmatau2[18 * i + 2];
    v4rhosigmatau2_4[i]  = v4rhosigmatau2[18 * i + 3];
    v4rhosigmatau2_5[i]  = v4rhosigmatau2[18 * i + 4];
    v4rhosigmatau2_6[i]  = v4rhosigmatau2[18 * i + 5];
    v4rhosigmatau2_7[i]  = v4rhosigmatau2[18 * i + 6];
    v4rhosigmatau2_8[i]  = v4rhosigmatau2[18 * i + 7];
    v4rhosigmatau2_9[i]  = v4rhosigmatau2[18 * i + 8];
    v4rhosigmatau2_10[i] = v4rhosigmatau2[18 * i + 9];
    v4rhosigmatau2_11[i] = v4rhosigmatau2[18 * i + 10];
    v4rhosigmatau2_12[i] = v4rhosigmatau2[18 * i + 11];
    v4rhosigmatau2_13[i] = v4rhosigmatau2[18 * i + 12];
    v4rhosigmatau2_14[i] = v4rhosigmatau2[18 * i + 13];
    v4rhosigmatau2_15[i] = v4rhosigmatau2[18 * i + 14];
    v4rhosigmatau2_16[i] = v4rhosigmatau2[18 * i + 15];
    v4rhosigmatau2_17[i] = v4rhosigmatau2[18 * i + 16];
    v4rhosigmatau2_18[i] = v4rhosigmatau2[18 * i + 17];

    // --- v4rholapl3: 8*np ---
    v4rholapl3_1[i] = v4rholapl3[8 * i + 0];
    v4rholapl3_2[i] = v4rholapl3[8 * i + 1];
    v4rholapl3_3[i] = v4rholapl3[8 * i + 2];
    v4rholapl3_4[i] = v4rholapl3[8 * i + 3];
    v4rholapl3_5[i] = v4rholapl3[8 * i + 4];
    v4rholapl3_6[i] = v4rholapl3[8 * i + 5];
    v4rholapl3_7[i] = v4rholapl3[8 * i + 6];
    v4rholapl3_8[i] = v4rholapl3[8 * i + 7];

    // --- v4rholapl2tau: 12*np ---
    v4rholapl2tau_1[i]  = v4rholapl2tau[12 * i + 0];
    v4rholapl2tau_2[i]  = v4rholapl2tau[12 * i + 1];
    v4rholapl2tau_3[i]  = v4rholapl2tau[12 * i + 2];
    v4rholapl2tau_4[i]  = v4rholapl2tau[12 * i + 3];
    v4rholapl2tau_5[i]  = v4rholapl2tau[12 * i + 4];
    v4rholapl2tau_6[i]  = v4rholapl2tau[12 * i + 5];
    v4rholapl2tau_7[i]  = v4rholapl2tau[12 * i + 6];
    v4rholapl2tau_8[i]  = v4rholapl2tau[12 * i + 7];
    v4rholapl2tau_9[i]  = v4rholapl2tau[12 * i + 8];
    v4rholapl2tau_10[i] = v4rholapl2tau[12 * i + 9];
    v4rholapl2tau_11[i] = v4rholapl2tau[12 * i + 10];
    v4rholapl2tau_12[i] = v4rholapl2tau[12 * i + 11];

    // --- v4rholapltau2: 12*np ---
    v4rholapltau2_1[i]  = v4rholapltau2[12 * i + 0];
    v4rholapltau2_2[i]  = v4rholapltau2[12 * i + 1];
    v4rholapltau2_3[i]  = v4rholapltau2[12 * i + 2];
    v4rholapltau2_4[i]  = v4rholapltau2[12 * i + 3];
    v4rholapltau2_5[i]  = v4rholapltau2[12 * i + 4];
    v4rholapltau2_6[i]  = v4rholapltau2[12 * i + 5];
    v4rholapltau2_7[i]  = v4rholapltau2[12 * i + 6];
    v4rholapltau2_8[i]  = v4rholapltau2[12 * i + 7];
    v4rholapltau2_9[i]  = v4rholapltau2[12 * i + 8];
    v4rholapltau2_10[i] = v4rholapltau2[12 * i + 9];
    v4rholapltau2_11[i] = v4rholapltau2[12 * i + 10];
    v4rholapltau2_12[i] = v4rholapltau2[12 * i + 11];

    // --- v4rhotau3: 8*np ---
    v4rhotau3_1[i] = v4rhotau3[8 * i + 0];
    v4rhotau3_2[i] = v4rhotau3[8 * i + 1];
    v4rhotau3_3[i] = v4rhotau3[8 * i + 2];
    v4rhotau3_4[i] = v4rhotau3[8 * i + 3];
    v4rhotau3_5[i] = v4rhotau3[8 * i + 4];
    v4rhotau3_6[i] = v4rhotau3[8 * i + 5];
    v4rhotau3_7[i] = v4rhotau3[8 * i + 6];
    v4rhotau3_8[i] = v4rhotau3[8 * i + 7];

    // --- v4sigma4: 15*np ---
    v4sigma4_1[i]  = v4sigma4[15 * i + 0];
    v4sigma4_2[i]  = v4sigma4[15 * i + 1];
    v4sigma4_3[i]  = v4sigma4[15 * i + 2];
    v4sigma4_4[i]  = v4sigma4[15 * i + 3];
    v4sigma4_5[i]  = v4sigma4[15 * i + 4];
    v4sigma4_6[i]  = v4sigma4[15 * i + 5];
    v4sigma4_7[i]  = v4sigma4[15 * i + 6];
    v4sigma4_8[i]  = v4sigma4[15 * i + 7];
    v4sigma4_9[i]  = v4sigma4[15 * i + 8];
    v4sigma4_10[i] = v4sigma4[15 * i + 9];
    v4sigma4_11[i] = v4sigma4[15 * i + 10];
    v4sigma4_12[i] = v4sigma4[15 * i + 11];
    v4sigma4_13[i] = v4sigma4[15 * i + 12];
    v4sigma4_14[i] = v4sigma4[15 * i + 13];
    v4sigma4_15[i] = v4sigma4[15 * i + 14];
       
       // --- v4sigma3lapl: 20*np ---
    v4sigma3lapl_1[i]  = v4sigma3lapl[20 * i + 0];
    v4sigma3lapl_2[i]  = v4sigma3lapl[20 * i + 1];
    v4sigma3lapl_3[i]  = v4sigma3lapl[20 * i + 2];
    v4sigma3lapl_4[i]  = v4sigma3lapl[20 * i + 3];
    v4sigma3lapl_5[i]  = v4sigma3lapl[20 * i + 4];
    v4sigma3lapl_6[i]  = v4sigma3lapl[20 * i + 5];
    v4sigma3lapl_7[i]  = v4sigma3lapl[20 * i + 6];
    v4sigma3lapl_8[i]  = v4sigma3lapl[20 * i + 7];
    v4sigma3lapl_9[i]  = v4sigma3lapl[20 * i + 8];
    v4sigma3lapl_10[i]  = v4sigma3lapl[20 * i + 9];
    v4sigma3lapl_11[i]  = v4sigma3lapl[20 * i + 10];
    v4sigma3lapl_12[i]  = v4sigma3lapl[20 * i + 11];
    v4sigma3lapl_13[i]  = v4sigma3lapl[20 * i + 12];
    v4sigma3lapl_14[i]  = v4sigma3lapl[20 * i + 13];
    v4sigma3lapl_15[i]  = v4sigma3lapl[20 * i + 14];
    v4sigma3lapl_16[i]  = v4sigma3lapl[20 * i + 15];
    v4sigma3lapl_17[i]  = v4sigma3lapl[20 * i + 16];
    v4sigma3lapl_18[i]  = v4sigma3lapl[20 * i + 17];
    v4sigma3lapl_19[i]  = v4sigma3lapl[20 * i + 18];
    v4sigma3lapl_20[i]  = v4sigma3lapl[20 * i + 19];

// --- v4sigma3tau: 20*np ---
    v4sigma3tau_1[i]  = v4sigma3tau[20 * i + 0];
    v4sigma3tau_2[i]  = v4sigma3tau[20 * i + 1];
    v4sigma3tau_3[i]  = v4sigma3tau[20 * i + 2];
    v4sigma3tau_4[i]  = v4sigma3tau[20 * i + 3];
    v4sigma3tau_5[i]  = v4sigma3tau[20 * i + 4];
    v4sigma3tau_6[i]  = v4sigma3tau[20 * i + 5];
    v4sigma3tau_7[i]  = v4sigma3tau[20 * i + 6];
    v4sigma3tau_8[i]  = v4sigma3tau[20 * i + 7];
    v4sigma3tau_9[i]  = v4sigma3tau[20 * i + 8];
    v4sigma3tau_10[i]  = v4sigma3tau[20 * i + 9];
    v4sigma3tau_11[i]  = v4sigma3tau[20 * i + 10];
    v4sigma3tau_12[i]  = v4sigma3tau[20 * i + 11];
    v4sigma3tau_13[i]  = v4sigma3tau[20 * i + 12];
    v4sigma3tau_14[i]  = v4sigma3tau[20 * i + 13];
    v4sigma3tau_15[i]  = v4sigma3tau[20 * i + 14];
    v4sigma3tau_16[i]  = v4sigma3tau[20 * i + 15];
    v4sigma3tau_17[i]  = v4sigma3tau[20 * i + 16];
    v4sigma3tau_18[i]  = v4sigma3tau[20 * i + 17];
    v4sigma3tau_19[i]  = v4sigma3tau[20 * i + 18];
    v4sigma3tau_20[i]  = v4sigma3tau[20 * i + 19];

// --- v4sigma2lapl2: 18*np ---
    v4sigma2lapl2_1[i]  = v4sigma2lapl2[18 * i + 0];
    v4sigma2lapl2_2[i]  = v4sigma2lapl2[18 * i + 1];
    v4sigma2lapl2_3[i]  = v4sigma2lapl2[18 * i + 2];
    v4sigma2lapl2_4[i]  = v4sigma2lapl2[18 * i + 3];
    v4sigma2lapl2_5[i]  = v4sigma2lapl2[18 * i + 4];
    v4sigma2lapl2_6[i]  = v4sigma2lapl2[18 * i + 5];
    v4sigma2lapl2_7[i]  = v4sigma2lapl2[18 * i + 6];
    v4sigma2lapl2_8[i]  = v4sigma2lapl2[18 * i + 7];
    v4sigma2lapl2_9[i]  = v4sigma2lapl2[18 * i + 8];
    v4sigma2lapl2_10[i]  = v4sigma2lapl2[18 * i + 9];
    v4sigma2lapl2_11[i]  = v4sigma2lapl2[18 * i + 10];
    v4sigma2lapl2_12[i]  = v4sigma2lapl2[18 * i + 11];
    v4sigma2lapl2_13[i]  = v4sigma2lapl2[18 * i + 12];
    v4sigma2lapl2_14[i]  = v4sigma2lapl2[18 * i + 13];
    v4sigma2lapl2_15[i]  = v4sigma2lapl2[18 * i + 14];
    v4sigma2lapl2_16[i]  = v4sigma2lapl2[18 * i + 15];
    v4sigma2lapl2_17[i]  = v4sigma2lapl2[18 * i + 16];
    v4sigma2lapl2_18[i]  = v4sigma2lapl2[18 * i + 17];

// --- v4sigma2lapltau: 24*np ---
    v4sigma2lapltau_1[i]  = v4sigma2lapltau[24 * i + 0];
    v4sigma2lapltau_2[i]  = v4sigma2lapltau[24 * i + 1];
    v4sigma2lapltau_3[i]  = v4sigma2lapltau[24 * i + 2];
    v4sigma2lapltau_4[i]  = v4sigma2lapltau[24 * i + 3];
    v4sigma2lapltau_5[i]  = v4sigma2lapltau[24 * i + 4];
    v4sigma2lapltau_6[i]  = v4sigma2lapltau[24 * i + 5];
    v4sigma2lapltau_7[i]  = v4sigma2lapltau[24 * i + 6];
    v4sigma2lapltau_8[i]  = v4sigma2lapltau[24 * i + 7];
    v4sigma2lapltau_9[i]  = v4sigma2lapltau[24 * i + 8];
    v4sigma2lapltau_10[i]  = v4sigma2lapltau[24 * i + 9];
    v4sigma2lapltau_11[i]  = v4sigma2lapltau[24 * i + 10];
    v4sigma2lapltau_12[i]  = v4sigma2lapltau[24 * i + 11];
    v4sigma2lapltau_13[i]  = v4sigma2lapltau[24 * i + 12];
    v4sigma2lapltau_14[i]  = v4sigma2lapltau[24 * i + 13];
    v4sigma2lapltau_15[i]  = v4sigma2lapltau[24 * i + 14];
    v4sigma2lapltau_16[i]  = v4sigma2lapltau[24 * i + 15];
    v4sigma2lapltau_17[i]  = v4sigma2lapltau[24 * i + 16];
    v4sigma2lapltau_18[i]  = v4sigma2lapltau[24 * i + 17];
    v4sigma2lapltau_19[i]  = v4sigma2lapltau[24 * i + 18];
    v4sigma2lapltau_20[i]  = v4sigma2lapltau[24 * i + 19];
    v4sigma2lapltau_21[i]  = v4sigma2lapltau[24 * i + 20];
    v4sigma2lapltau_22[i]  = v4sigma2lapltau[24 * i + 21];
    v4sigma2lapltau_23[i]  = v4sigma2lapltau[24 * i + 22];
    v4sigma2lapltau_24[i]  = v4sigma2lapltau[24 * i + 23];

// --- v4sigma2tau2: 18*np ---
    v4sigma2tau2_1[i]  = v4sigma2tau2[18 * i + 0];
    v4sigma2tau2_2[i]  = v4sigma2tau2[18 * i + 1];
    v4sigma2tau2_3[i]  = v4sigma2tau2[18 * i + 2];
    v4sigma2tau2_4[i]  = v4sigma2tau2[18 * i + 3];
    v4sigma2tau2_5[i]  = v4sigma2tau2[18 * i + 4];
    v4sigma2tau2_6[i]  = v4sigma2tau2[18 * i + 5];
    v4sigma2tau2_7[i]  = v4sigma2tau2[18 * i + 6];
    v4sigma2tau2_8[i]  = v4sigma2tau2[18 * i + 7];
    v4sigma2tau2_9[i]  = v4sigma2tau2[18 * i + 8];
    v4sigma2tau2_10[i]  = v4sigma2tau2[18 * i + 9];
    v4sigma2tau2_11[i]  = v4sigma2tau2[18 * i + 10];
    v4sigma2tau2_12[i]  = v4sigma2tau2[18 * i + 11];
    v4sigma2tau2_13[i]  = v4sigma2tau2[18 * i + 12];
    v4sigma2tau2_14[i]  = v4sigma2tau2[18 * i + 13];
    v4sigma2tau2_15[i]  = v4sigma2tau2[18 * i + 14];
    v4sigma2tau2_16[i]  = v4sigma2tau2[18 * i + 15];
    v4sigma2tau2_17[i]  = v4sigma2tau2[18 * i + 16];
    v4sigma2tau2_18[i]  = v4sigma2tau2[18 * i + 17];

// --- v4sigmalapl3: 12*np ---
    v4sigmalapl3_1[i]  = v4sigmalapl3[12 * i + 0];
    v4sigmalapl3_2[i]  = v4sigmalapl3[12 * i + 1];
    v4sigmalapl3_3[i]  = v4sigmalapl3[12 * i + 2];
    v4sigmalapl3_4[i]  = v4sigmalapl3[12 * i + 3];
    v4sigmalapl3_5[i]  = v4sigmalapl3[12 * i + 4];
    v4sigmalapl3_6[i]  = v4sigmalapl3[12 * i + 5];
    v4sigmalapl3_7[i]  = v4sigmalapl3[12 * i + 6];
    v4sigmalapl3_8[i]  = v4sigmalapl3[12 * i + 7];
    v4sigmalapl3_9[i]  = v4sigmalapl3[12 * i + 8];
    v4sigmalapl3_10[i]  = v4sigmalapl3[12 * i + 9];
    v4sigmalapl3_11[i]  = v4sigmalapl3[12 * i + 10];
    v4sigmalapl3_12[i]  = v4sigmalapl3[12 * i + 11];

// --- v4sigmalapl2tau: 18*np ---
    v4sigmalapl2tau_1[i]  = v4sigmalapl2tau[18 * i + 0];
    v4sigmalapl2tau_2[i]  = v4sigmalapl2tau[18 * i + 1];
    v4sigmalapl2tau_3[i]  = v4sigmalapl2tau[18 * i + 2];
    v4sigmalapl2tau_4[i]  = v4sigmalapl2tau[18 * i + 3];
    v4sigmalapl2tau_5[i]  = v4sigmalapl2tau[18 * i + 4];
    v4sigmalapl2tau_6[i]  = v4sigmalapl2tau[18 * i + 5];
    v4sigmalapl2tau_7[i]  = v4sigmalapl2tau[18 * i + 6];
    v4sigmalapl2tau_8[i]  = v4sigmalapl2tau[18 * i + 7];
    v4sigmalapl2tau_9[i]  = v4sigmalapl2tau[18 * i + 8];
    v4sigmalapl2tau_10[i]  = v4sigmalapl2tau[18 * i + 9];
    v4sigmalapl2tau_11[i]  = v4sigmalapl2tau[18 * i + 10];
    v4sigmalapl2tau_12[i]  = v4sigmalapl2tau[18 * i + 11];
    v4sigmalapl2tau_13[i]  = v4sigmalapl2tau[18 * i + 12];
    v4sigmalapl2tau_14[i]  = v4sigmalapl2tau[18 * i + 13];
    v4sigmalapl2tau_15[i]  = v4sigmalapl2tau[18 * i + 14];
    v4sigmalapl2tau_16[i]  = v4sigmalapl2tau[18 * i + 15];
    v4sigmalapl2tau_17[i]  = v4sigmalapl2tau[18 * i + 16];
    v4sigmalapl2tau_18[i]  = v4sigmalapl2tau[18 * i + 17];

// --- v4sigmalapltau2: 18*np ---
    v4sigmalapltau2_1[i]  = v4sigmalapltau2[18 * i + 0];
    v4sigmalapltau2_2[i]  = v4sigmalapltau2[18 * i + 1];
    v4sigmalapltau2_3[i]  = v4sigmalapltau2[18 * i + 2];
    v4sigmalapltau2_4[i]  = v4sigmalapltau2[18 * i + 3];
    v4sigmalapltau2_5[i]  = v4sigmalapltau2[18 * i + 4];
    v4sigmalapltau2_6[i]  = v4sigmalapltau2[18 * i + 5];
    v4sigmalapltau2_7[i]  = v4sigmalapltau2[18 * i + 6];
    v4sigmalapltau2_8[i]  = v4sigmalapltau2[18 * i + 7];
    v4sigmalapltau2_9[i]  = v4sigmalapltau2[18 * i + 8];
    v4sigmalapltau2_10[i]  = v4sigmalapltau2[18 * i + 9];
    v4sigmalapltau2_11[i]  = v4sigmalapltau2[18 * i + 10];
    v4sigmalapltau2_12[i]  = v4sigmalapltau2[18 * i + 11];
    v4sigmalapltau2_13[i]  = v4sigmalapltau2[18 * i + 12];
    v4sigmalapltau2_14[i]  = v4sigmalapltau2[18 * i + 13];
    v4sigmalapltau2_15[i]  = v4sigmalapltau2[18 * i + 14];
    v4sigmalapltau2_16[i]  = v4sigmalapltau2[18 * i + 15];
    v4sigmalapltau2_17[i]  = v4sigmalapltau2[18 * i + 16];
    v4sigmalapltau2_18[i]  = v4sigmalapltau2[18 * i + 17];

// --- v4sigmatau3: 12*np ---
    v4sigmatau3_1[i]  = v4sigmatau3[12 * i + 0];
    v4sigmatau3_2[i]  = v4sigmatau3[12 * i + 1];
    v4sigmatau3_3[i]  = v4sigmatau3[12 * i + 2];
    v4sigmatau3_4[i]  = v4sigmatau3[12 * i + 3];
    v4sigmatau3_5[i]  = v4sigmatau3[12 * i + 4];
    v4sigmatau3_6[i]  = v4sigmatau3[12 * i + 5];
    v4sigmatau3_7[i]  = v4sigmatau3[12 * i + 6];
    v4sigmatau3_8[i]  = v4sigmatau3[12 * i + 7];
    v4sigmatau3_9[i]  = v4sigmatau3[12 * i + 8];
    v4sigmatau3_10[i]  = v4sigmatau3[12 * i + 9];
    v4sigmatau3_11[i]  = v4sigmatau3[12 * i + 10];
    v4sigmatau3_12[i]  = v4sigmatau3[12 * i + 11];

// --- v4lapl4: 5*np ---
    v4lapl4_1[i]  = v4lapl4[5 * i + 0];
    v4lapl4_2[i]  = v4lapl4[5 * i + 1];
    v4lapl4_3[i]  = v4lapl4[5 * i + 2];
    v4lapl4_4[i]  = v4lapl4[5 * i + 3];
    v4lapl4_5[i]  = v4lapl4[5 * i + 4];

// --- v4lapl3tau: 8*np ---
    v4lapl3tau_1[i]  = v4lapl3tau[8 * i + 0];
    v4lapl3tau_2[i]  = v4lapl3tau[8 * i + 1];
    v4lapl3tau_3[i]  = v4lapl3tau[8 * i + 2];
    v4lapl3tau_4[i]  = v4lapl3tau[8 * i + 3];
    v4lapl3tau_5[i]  = v4lapl3tau[8 * i + 4];
    v4lapl3tau_6[i]  = v4lapl3tau[8 * i + 5];
    v4lapl3tau_7[i]  = v4lapl3tau[8 * i + 6];
    v4lapl3tau_8[i]  = v4lapl3tau[8 * i + 7];

// --- v4lapl2tau2: 9*np ---
    v4lapl2tau2_1[i]  = v4lapl2tau2[9 * i + 0];
    v4lapl2tau2_2[i]  = v4lapl2tau2[9 * i + 1];
    v4lapl2tau2_3[i]  = v4lapl2tau2[9 * i + 2];
    v4lapl2tau2_4[i]  = v4lapl2tau2[9 * i + 3];
    v4lapl2tau2_5[i]  = v4lapl2tau2[9 * i + 4];
    v4lapl2tau2_6[i]  = v4lapl2tau2[9 * i + 5];
    v4lapl2tau2_7[i]  = v4lapl2tau2[9 * i + 6];
    v4lapl2tau2_8[i]  = v4lapl2tau2[9 * i + 7];
    v4lapl2tau2_9[i]  = v4lapl2tau2[9 * i + 8];

// --- v4lapltau3: 8*np ---
    v4lapltau3_1[i]  = v4lapltau3[8 * i + 0];
    v4lapltau3_2[i]  = v4lapltau3[8 * i + 1];
    v4lapltau3_3[i]  = v4lapltau3[8 * i + 2];
    v4lapltau3_4[i]  = v4lapltau3[8 * i + 3];
    v4lapltau3_5[i]  = v4lapltau3[8 * i + 4];
    v4lapltau3_6[i]  = v4lapltau3[8 * i + 5];
    v4lapltau3_7[i]  = v4lapltau3[8 * i + 6];
    v4lapltau3_8[i]  = v4lapltau3[8 * i + 7];

// --- v4tau4: 5*np ---
    v4tau4_1[i]  = v4tau4[5 * i + 0];
    v4tau4_2[i]  = v4tau4[5 * i + 1];
    v4tau4_3[i]  = v4tau4[5 * i + 2];
    v4tau4_4[i]  = v4tau4[5 * i + 3];
    v4tau4_5[i]  = v4tau4[5 * i + 4];


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

    int np = rho_up.size();

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

    // Compute third derivatives of the exchange-correlation potential
    std::vector<double> v3rho3_1(np), v3rho3_2(np), v3rho3_3(np), v3rho3_4(np);
    std::vector<double> v3rho2sigma_1(np), v3rho2sigma_2(np), v3rho2sigma_3(np), v3rho2sigma_4(np), v3rho2sigma_5(np), v3rho2sigma_6(np), v3rho2sigma_7(np), v3rho2sigma_8(np), v3rho2sigma_9(np);
    std::vector<double> v3rho2lapl_1(np), v3rho2lapl_2(np), v3rho2lapl_3(np), v3rho2lapl_4(np), v3rho2lapl_5(np), v3rho2lapl_6(np);
    std::vector<double> v3rho2tau_1(np), v3rho2tau_2(np), v3rho2tau_3(np), v3rho2tau_4(np), v3rho2tau_5(np), v3rho2tau_6(np);
    std::vector<double> v3rhosigma2_1(np), v3rhosigma2_2(np), v3rhosigma2_3(np), v3rhosigma2_4(np), v3rhosigma2_5(np), v3rhosigma2_6(np), v3rhosigma2_7(np), v3rhosigma2_8(np), v3rhosigma2_9(np), v3rhosigma2_10(np), v3rhosigma2_11(np), v3rhosigma2_12(np);
    std::vector<double> v3rhosigmalapl_1(np), v3rhosigmalapl_2(np), v3rhosigmalapl_3(np), v3rhosigmalapl_4(np), v3rhosigmalapl_5(np), v3rhosigmalapl_6(np), v3rhosigmalapl_7(np), v3rhosigmalapl_8(np), v3rhosigmalapl_9(np), v3rhosigmalapl_10(np), v3rhosigmalapl_11(np), v3rhosigmalapl_12(np);
    std::vector<double> v3rhosigmatau_1(np), v3rhosigmatau_2(np), v3rhosigmatau_3(np), v3rhosigmatau_4(np), v3rhosigmatau_5(np), v3rhosigmatau_6(np), v3rhosigmatau_7(np), v3rhosigmatau_8(np), v3rhosigmatau_9(np), v3rhosigmatau_10(np), v3rhosigmatau_11(np), v3rhosigmatau_12(np);
    std::vector<double> v3rholapl2_1(np), v3rholapl2_2(np), v3rholapl2_3(np), v3rholapl2_4(np), v3rholapl2_5(np), v3rholapl2_6(np);
    std::vector<double> v3rholapltau_1(np), v3rholapltau_2(np), v3rholapltau_3(np), v3rholapltau_4(np), v3rholapltau_5(np), v3rholapltau_6(np), v3rholapltau_7(np), v3rholapltau_8(np);
    std::vector<double> v3rhotau2_1(np), v3rhotau2_2(np), v3rhotau2_3(np), v3rhotau2_4(np), v3rhotau2_5(np), v3rhotau2_6(np);
    std::vector<double> v3sigma3_1(np), v3sigma3_2(np), v3sigma3_3(np), v3sigma3_4(np), v3sigma3_5(np), v3sigma3_6(np), v3sigma3_7(np), v3sigma3_8(np), v3sigma3_9(np), v3sigma3_10(np);
    std::vector<double> v3sigma2lapl_1(np), v3sigma2lapl_2(np), v3sigma2lapl_3(np), v3sigma2lapl_4(np), v3sigma2lapl_5(np), v3sigma2lapl_6(np), v3sigma2lapl_7(np), v3sigma2lapl_8(np), v3sigma2lapl_9(np), v3sigma2lapl_10(np), v3sigma2lapl_11(np), v3sigma2lapl_12(np);
    std::vector<double> v3sigma2tau_1(np), v3sigma2tau_2(np), v3sigma2tau_3(np), v3sigma2tau_4(np), v3sigma2tau_5(np), v3sigma2tau_6(np), v3sigma2tau_7(np), v3sigma2tau_8(np), v3sigma2tau_9(np), v3sigma2tau_10(np), v3sigma2tau_11(np), v3sigma2tau_12(np);
    std::vector<double> v3sigmalapl2_1(np), v3sigmalapl2_2(np), v3sigmalapl2_3(np), v3sigmalapl2_4(np), v3sigmalapl2_5(np), v3sigmalapl2_6(np), v3sigmalapl2_7(np), v3sigmalapl2_8(np), v3sigmalapl2_9(np);
    std::vector<double> v3sigmalapltau_1(np), v3sigmalapltau_2(np), v3sigmalapltau_3(np), v3sigmalapltau_4(np), v3sigmalapltau_5(np), v3sigmalapltau_6(np), v3sigmalapltau_7(np), v3sigmalapltau_8(np), v3sigmalapltau_9(np), v3sigmalapltau_10(np), v3sigmalapltau_11(np), v3sigmalapltau_12(np);
    std::vector<double> v3sigmatau2_1(np), v3sigmatau2_2(np), v3sigmatau2_3(np), v3sigmatau2_4(np), v3sigmatau2_5(np), v3sigmatau2_6(np), v3sigmatau2_7(np), v3sigmatau2_8(np), v3sigmatau2_9(np);
    std::vector<double> v3lapl3_1(np), v3lapl3_2(np), v3lapl3_3(np), v3lapl3_4(np);
    std::vector<double> v3lapl2tau_1(np), v3lapl2tau_2(np), v3lapl2tau_3(np), v3lapl2tau_4(np), v3lapl2tau_5(np), v3lapl2tau_6(np);
    std::vector<double> v3lapltau2_1(np), v3lapltau2_2(np), v3lapltau2_3(np), v3lapltau2_4(np), v3lapltau2_5(np), v3lapltau2_6(np);
    std::vector<double> v3tau3_1(np), v3tau3_2(np), v3tau3_3(np), v3tau3_4(np);

    mgga_kxc(
        rho_up, rho_down,
        sigma_1, sigma_2, sigma_3,
        lapl_up, lapl_down,
        tau_up, tau_down,
        v3rho3_1, v3rho3_2, v3rho3_3, v3rho3_4,
        v3rho2sigma_1, v3rho2sigma_2, v3rho2sigma_3,
        v3rho2sigma_4, v3rho2sigma_5, v3rho2sigma_6,
        v3rho2sigma_7, v3rho2sigma_8, v3rho2sigma_9,
        v3rho2lapl_1, v3rho2lapl_2, v3rho2lapl_3,
        v3rho2lapl_4, v3rho2lapl_5, v3rho2lapl_6,
        v3rho2tau_1, v3rho2tau_2, v3rho2tau_3,
        v3rho2tau_4, v3rho2tau_5, v3rho2tau_6,
        v3rhosigma2_1, v3rhosigma2_2, v3rhosigma2_3,
        v3rhosigma2_4, v3rhosigma2_5, v3rhosigma2_6,
        v3rhosigma2_7, v3rhosigma2_8, v3rhosigma2_9,
        v3rhosigma2_10, v3rhosigma2_11, v3rhosigma2_12,
        v3rhosigmalapl_1, v3rhosigmalapl_2, v3rhosigmalapl_3,
        v3rhosigmalapl_4, v3rhosigmalapl_5, v3rhosigmalapl_6,
        v3rhosigmalapl_7, v3rhosigmalapl_8, v3rhosigmalapl_9,
        v3rhosigmalapl_10, v3rhosigmalapl_11, v3rhosigmalapl_12,
        v3rhosigmatau_1, v3rhosigmatau_2, v3rhosigmatau_3,
        v3rhosigmatau_4, v3rhosigmatau_5, v3rhosigmatau_6,
        v3rhosigmatau_7, v3rhosigmatau_8, v3rhosigmatau_9,
        v3rhosigmatau_10, v3rhosigmatau_11, v3rhosigmatau_12,
        v3rholapl2_1, v3rholapl2_2, v3rholapl2_3,
        v3rholapl2_4, v3rholapl2_5, v3rholapl2_6,
        v3rholapltau_1, v3rholapltau_2, v3rholapltau_3,
        v3rholapltau_4, v3rholapltau_5, v3rholapltau_6,
        v3rholapltau_7, v3rholapltau_8,
        v3rhotau2_1, v3rhotau2_2, v3rhotau2_3,
        v3rhotau2_4, v3rhotau2_5, v3rhotau2_6,
        v3sigma3_1, v3sigma3_2, v3sigma3_3,
        v3sigma3_4, v3sigma3_5, v3sigma3_6,
        v3sigma3_7, v3sigma3_8, v3sigma3_9,
        v3sigma3_10,
        v3sigma2lapl_1, v3sigma2lapl_2, v3sigma2lapl_3,
        v3sigma2lapl_4, v3sigma2lapl_5, v3sigma2lapl_6,
        v3sigma2lapl_7, v3sigma2lapl_8, v3sigma2lapl_9,
        v3sigma2lapl_10, v3sigma2lapl_11, v3sigma2lapl_12,
        v3sigma2tau_1, v3sigma2tau_2, v3sigma2tau_3,
        v3sigma2tau_4, v3sigma2tau_5, v3sigma2tau_6,
        v3sigma2tau_7, v3sigma2tau_8, v3sigma2tau_9,
        v3sigma2tau_10, v3sigma2tau_11, v3sigma2tau_12,
        v3sigmalapl2_1, v3sigmalapl2_2, v3sigmalapl2_3,
        v3sigmalapl2_4, v3sigmalapl2_5, v3sigmalapl2_6,
        v3sigmalapl2_7, v3sigmalapl2_8, v3sigmalapl2_9,
        v3sigmalapltau_1, v3sigmalapltau_2, v3sigmalapltau_3,
        v3sigmalapltau_4, v3sigmalapltau_5, v3sigmalapltau_6,
        v3sigmalapltau_7, v3sigmalapltau_8, v3sigmalapltau_9,
        v3sigmalapltau_10, v3sigmalapltau_11, v3sigmalapltau_12,
        v3sigmatau2_1, v3sigmatau2_2, v3sigmatau2_3,
        v3sigmatau2_4, v3sigmatau2_5, v3sigmatau2_6,
        v3sigmatau2_7, v3sigmatau2_8, v3sigmatau2_9,
        v3lapl3_1, v3lapl3_2, v3lapl3_3, v3lapl3_4,
        v3lapl2tau_1, v3lapl2tau_2, v3lapl2tau_3,
        v3lapl2tau_4, v3lapl2tau_5, v3lapl2tau_6,
        v3lapltau2_1, v3lapltau2_2, v3lapltau2_3,
        v3lapltau2_4, v3lapltau2_5, v3lapltau2_6,
        v3tau3_1, v3tau3_2, v3tau3_3, v3tau3_4
    );

    // compute the fourth derivative
    std::vector<double> v4rho4_1(np), v4rho4_2(np), v4rho4_3(np), v4rho4_4(np), v4rho4_5(np);
std::vector<double> v4rho3sigma_1(np), v4rho3sigma_2(np), v4rho3sigma_3(np), v4rho3sigma_4(np), v4rho3sigma_5(np), v4rho3sigma_6(np), v4rho3sigma_7(np), v4rho3sigma_8(np), v4rho3sigma_9(np), v4rho3sigma_10(np), v4rho3sigma_11(np), v4rho3sigma_12(np);
std::vector<double> v4rho3lapl_1(np), v4rho3lapl_2(np), v4rho3lapl_3(np), v4rho3lapl_4(np), v4rho3lapl_5(np), v4rho3lapl_6(np), v4rho3lapl_7(np), v4rho3lapl_8(np);
std::vector<double> v4rho3tau_1(np), v4rho3tau_2(np), v4rho3tau_3(np), v4rho3tau_4(np), v4rho3tau_5(np), v4rho3tau_6(np), v4rho3tau_7(np), v4rho3tau_8(np);
std::vector<double> v4rho2sigma2_1(np), v4rho2sigma2_2(np), v4rho2sigma2_3(np), v4rho2sigma2_4(np), v4rho2sigma2_5(np), v4rho2sigma2_6(np), v4rho2sigma2_7(np), v4rho2sigma2_8(np), v4rho2sigma2_9(np), v4rho2sigma2_10(np), v4rho2sigma2_11(np), v4rho2sigma2_12(np), v4rho2sigma2_13(np), v4rho2sigma2_14(np), v4rho2sigma2_15(np), v4rho2sigma2_16(np), v4rho2sigma2_17(np), v4rho2sigma2_18(np);  
std::vector<double> v4rho2sigmalapl_1(np), v4rho2sigmalapl_2(np), v4rho2sigmalapl_3(np), v4rho2sigmalapl_4(np), v4rho2sigmalapl_5(np), v4rho2sigmalapl_6(np), v4rho2sigmalapl_7(np), v4rho2sigmalapl_8(np), v4rho2sigmalapl_9(np), v4rho2sigmalapl_10(np), v4rho2sigmalapl_11(np), v4rho2sigmalapl_12(np), v4rho2sigmalapl_13(np), v4rho2sigmalapl_14(np), v4rho2sigmalapl_15(np), v4rho2sigmalapl_16(np), v4rho2sigmalapl_17(np), v4rho2sigmalapl_18(np);
std::vector<double> v4rho2sigmatau_1(np), v4rho2sigmatau_2(np), v4rho2sigmatau_3(np), v4rho2sigmatau_4(np), v4rho2sigmatau_5(np), v4rho2sigmatau_6(np), v4rho2sigmatau_7(np), v4rho2sigmatau_8(np), v4rho2sigmatau_9(np), v4rho2sigmatau_10(np), v4rho2sigmatau_11(np), v4rho2sigmatau_12(np), v4rho2sigmatau_13(np), v4rho2sigmatau_14(np), v4rho2sigmatau_15(np), v4rho2sigmatau_16(np), v4rho2sigmatau_17(np), v4rho2sigmatau_18(np);
std::vector<double> v4rho2lapl2_1(np), v4rho2lapl2_2(np), v4rho2lapl2_3(np), v4rho2lapl2_4(np), v4rho2lapl2_5(np), v4rho2lapl2_6(np), v4rho2lapl2_7(np), v4rho2lapl2_8(np), v4rho2lapl2_9(np);     
std::vector<double> v4rho2lapltau_1(np), v4rho2lapltau_2(np), v4rho2lapltau_3(np), v4rho2lapltau_4(np), v4rho2lapltau_5(np), v4rho2lapltau_6(np), v4rho2lapltau_7(np), v4rho2lapltau_8(np), v4rho2lapltau_9(np), v4rho2lapltau_10(np), v4rho2lapltau_11(np), v4rho2lapltau_12(np);
std::vector<double> v4rho2tau2_1(np), v4rho2tau2_2(np), v4rho2tau2_3(np), v4rho2tau2_4(np), v4rho2tau2_5(np), v4rho2tau2_6(np), v4rho2tau2_7(np), v4rho2tau2_8(np), v4rho2tau2_9(np);
std::vector<double> v4rhosigma3_1(np), v4rhosigma3_2(np), v4rhosigma3_3(np), v4rhosigma3_4(np), v4rhosigma3_5(np), v4rhosigma3_6(np), v4rhosigma3_7(np), v4rhosigma3_8(np), v4rhosigma3_9(np), v4rhosigma3_10(np), v4rhosigma3_11(np), v4rhosigma3_12(np), v4rhosigma3_13(np), v4rhosigma3_14(np), v4rhosigma3_15(np), v4rhosigma3_16(np), v4rhosigma3_17(np), v4rhosigma3_18(np), v4rhosigma3_19(np), v4rhosigma3_20(np);
std::vector<double> v4rhosigma2lapl_1(np), v4rhosigma2lapl_2(np), v4rhosigma2lapl_3(np), v4rhosigma2lapl_4(np), v4rhosigma2lapl_5(np), v4rhosigma2lapl_6(np), v4rhosigma2lapl_7(np), v4rhosigma2lapl_8(np), v4rhosigma2lapl_9(np), v4rhosigma2lapl_10(np), v4rhosigma2lapl_11(np), v4rhosigma2lapl_12(np), v4rhosigma2lapl_13(np), v4rhosigma2lapl_14(np), v4rhosigma2lapl_15(np), v4rhosigma2lapl_16(np), v4rhosigma2lapl_17(np), v4rhosigma2lapl_18(np), v4rhosigma2lapl_19(np), v4rhosigma2lapl_20(np), v4rhosigma2lapl_21(np), v4rhosigma2lapl_22(np), v4rhosigma2lapl_23(np), v4rhosigma2lapl_24(np);
std::vector<double> v4rhosigma2tau_1(np), v4rhosigma2tau_2(np), v4rhosigma2tau_3(np), v4rhosigma2tau_4(np), v4rhosigma2tau_5(np), v4rhosigma2tau_6(np), v4rhosigma2tau_7(np), v4rhosigma2tau_8(np), v4rhosigma2tau_9(np), v4rhosigma2tau_10(np), v4rhosigma2tau_11(np), v4rhosigma2tau_12(np), v4rhosigma2tau_13(np), v4rhosigma2tau_14(np), v4rhosigma2tau_15(np), v4rhosigma2tau_16(np), v4rhosigma2tau_17(np), v4rhosigma2tau_18(np), v4rhosigma2tau_19(np), v4rhosigma2tau_20(np), v4rhosigma2tau_21(np), v4rhosigma2tau_22(np), v4rhosigma2tau_23(np), v4rhosigma2tau_24(np);
std::vector<double> v4rhosigmalapl2_1(np), v4rhosigmalapl2_2(np), v4rhosigmalapl2_3(np), v4rhosigmalapl2_4(np), v4rhosigmalapl2_5(np), v4rhosigmalapl2_6(np), v4rhosigmalapl2_7(np), v4rhosigmalapl2_8(np), v4rhosigmalapl2_9(np), v4rhosigmalapl2_10(np), v4rhosigmalapl2_11(np), v4rhosigmalapl2_12(np), v4rhosigmalapl2_13(np), v4rhosigmalapl2_14(np), v4rhosigmalapl2_15(np), v4rhosigmalapl2_16(np), v4rhosigmalapl2_17(np), v4rhosigmalapl2_18(np);
std::vector<double> v4rhosigmalapltau_1(np), v4rhosigmalapltau_2(np), v4rhosigmalapltau_3(np), v4rhosigmalapltau_4(np), v4rhosigmalapltau_5(np), v4rhosigmalapltau_6(np), v4rhosigmalapltau_7(np), v4rhosigmalapltau_8(np), v4rhosigmalapltau_9(np), v4rhosigmalapltau_10(np), v4rhosigmalapltau_11(np), v4rhosigmalapltau_12(np), v4rhosigmalapltau_13(np), v4rhosigmalapltau_14(np), v4rhosigmalapltau_15(np), v4rhosigmalapltau_16(np), v4rhosigmalapltau_17(np), v4rhosigmalapltau_18(np), v4rhosigmalapltau_19(np), v4rhosigmalapltau_20(np), v4rhosigmalapltau_21(np), v4rhosigmalapltau_22(np), v4rhosigmalapltau_23(np), v4rhosigmalapltau_24(np);
std::vector<double> v4rhosigmatau2_1(np), v4rhosigmatau2_2(np), v4rhosigmatau2_3(np), v4rhosigmatau2_4(np), v4rhosigmatau2_5(np), v4rhosigmatau2_6(np), v4rhosigmatau2_7(np), v4rhosigmatau2_8(np), v4rhosigmatau2_9(np), v4rhosigmatau2_10(np), v4rhosigmatau2_11(np), v4rhosigmatau2_12(np), v4rhosigmatau2_13(np), v4rhosigmatau2_14(np), v4rhosigmatau2_15(np), v4rhosigmatau2_16(np), v4rhosigmatau2_17(np), v4rhosigmatau2_18(np);
std::vector<double> v4rholapl3_1(np), v4rholapl3_2(np), v4rholapl3_3(np), v4rholapl3_4(np), v4rholapl3_5(np), v4rholapl3_6(np), v4rholapl3_7(np), v4rholapl3_8(np);
std::vector<double> v4rholapl2tau_1(np), v4rholapl2tau_2(np), v4rholapl2tau_3(np), v4rholapl2tau_4(np), v4rholapl2tau_5(np), v4rholapl2tau_6(np), v4rholapl2tau_7(np), v4rholapl2tau_8(np), v4rholapl2tau_9(np), v4rholapl2tau_10(np), v4rholapl2tau_11(np), v4rholapl2tau_12(np);
std::vector<double> v4rholapltau2_1(np), v4rholapltau2_2(np), v4rholapltau2_3(np), v4rholapltau2_4(np), v4rholapltau2_5(np), v4rholapltau2_6(np), v4rholapltau2_7(np), v4rholapltau2_8(np), v4rholapltau2_9(np), v4rholapltau2_10(np), v4rholapltau2_11(np), v4rholapltau2_12(np);
std::vector<double> v4rhotau3_1(np), v4rhotau3_2(np), v4rhotau3_3(np), v4rhotau3_4(np), v4rhotau3_5(np), v4rhotau3_6(np), v4rhotau3_7(np), v4rhotau3_8(np);
std::vector<double> v4sigma4_1(np), v4sigma4_2(np), v4sigma4_3(np), v4sigma4_4(np), v4sigma4_5(np), v4sigma4_6(np), v4sigma4_7(np), v4sigma4_8(np), v4sigma4_9(np), v4sigma4_10(np), v4sigma4_11(np), v4sigma4_12(np), v4sigma4_13(np), v4sigma4_14(np), v4sigma4_15(np);
std::vector<double> v4sigma3lapl_1(np), v4sigma3lapl_2(np), v4sigma3lapl_3(np), v4sigma3lapl_4(np), v4sigma3lapl_5(np), v4sigma3lapl_6(np), v4sigma3lapl_7(np), v4sigma3lapl_8(np), v4sigma3lapl_9(np), v4sigma3lapl_10(np), v4sigma3lapl_11(np), v4sigma3lapl_12(np), v4sigma3lapl_13(np), v4sigma3lapl_14(np), v4sigma3lapl_15(np), v4sigma3lapl_16(np), v4sigma3lapl_17(np), v4sigma3lapl_18(np), v4sigma3lapl_19(np), v4sigma3lapl_20(np);
std::vector<double> v4sigma3tau_1(np), v4sigma3tau_2(np), v4sigma3tau_3(np), v4sigma3tau_4(np), v4sigma3tau_5(np), v4sigma3tau_6(np), v4sigma3tau_7(np), v4sigma3tau_8(np), v4sigma3tau_9(np), v4sigma3tau_10(np), v4sigma3tau_11(np), v4sigma3tau_12(np), v4sigma3tau_13(np), v4sigma3tau_14(np), v4sigma3tau_15(np), v4sigma3tau_16(np), v4sigma3tau_17(np), v4sigma3tau_18(np), v4sigma3tau_19(np), v4sigma3tau_20(np);
std::vector<double> v4sigma2lapl2_1(np), v4sigma2lapl2_2(np), v4sigma2lapl2_3(np), v4sigma2lapl2_4(np), v4sigma2lapl2_5(np), v4sigma2lapl2_6(np), v4sigma2lapl2_7(np), v4sigma2lapl2_8(np), v4sigma2lapl2_9(np), v4sigma2lapl2_10(np), v4sigma2lapl2_11(np), v4sigma2lapl2_12(np), v4sigma2lapl2_13(np), v4sigma2lapl2_14(np), v4sigma2lapl2_15(np), v4sigma2lapl2_16(np), v4sigma2lapl2_17(np), v4sigma2lapl2_18(np);
std::vector<double> v4sigma2lapltau_1(np), v4sigma2lapltau_2(np), v4sigma2lapltau_3(np), v4sigma2lapltau_4(np), v4sigma2lapltau_5(np), v4sigma2lapltau_6(np), v4sigma2lapltau_7(np), v4sigma2lapltau_8(np), v4sigma2lapltau_9(np), v4sigma2lapltau_10(np), v4sigma2lapltau_11(np), v4sigma2lapltau_12(np), v4sigma2lapltau_13(np), v4sigma2lapltau_14(np), v4sigma2lapltau_15(np), v4sigma2lapltau_16(np), v4sigma2lapltau_17(np), v4sigma2lapltau_18(np), v4sigma2lapltau_19(np), v4sigma2lapltau_20(np), v4sigma2lapltau_21(np), v4sigma2lapltau_22(np), v4sigma2lapltau_23(np), v4sigma2lapltau_24(np);
std::vector<double> v4sigma2tau2_1(np), v4sigma2tau2_2(np), v4sigma2tau2_3(np), v4sigma2tau2_4(np), v4sigma2tau2_5(np), v4sigma2tau2_6(np), v4sigma2tau2_7(np), v4sigma2tau2_8(np), v4sigma2tau2_9(np), v4sigma2tau2_10(np), v4sigma2tau2_11(np), v4sigma2tau2_12(np), v4sigma2tau2_13(np), v4sigma2tau2_14(np), v4sigma2tau2_15(np), v4sigma2tau2_16(np), v4sigma2tau2_17(np), v4sigma2tau2_18(np);  
std::vector<double> v4sigmalapl3_1(np), v4sigmalapl3_2(np), v4sigmalapl3_3(np), v4sigmalapl3_4(np), v4sigmalapl3_5(np), v4sigmalapl3_6(np), v4sigmalapl3_7(np), v4sigmalapl3_8(np), v4sigmalapl3_9(np), v4sigmalapl3_10(np), v4sigmalapl3_11(np), v4sigmalapl3_12(np);
std::vector<double> v4sigmalapl2tau_1(np), v4sigmalapl2tau_2(np), v4sigmalapl2tau_3(np), v4sigmalapl2tau_4(np), v4sigmalapl2tau_5(np), v4sigmalapl2tau_6(np), v4sigmalapl2tau_7(np), v4sigmalapl2tau_8(np), v4sigmalapl2tau_9(np), v4sigmalapl2tau_10(np), v4sigmalapl2tau_11(np), v4sigmalapl2tau_12(np), v4sigmalapl2tau_13(np), v4sigmalapl2tau_14(np), v4sigmalapl2tau_15(np), v4sigmalapl2tau_16(np), v4sigmalapl2tau_17(np), v4sigmalapl2tau_18(np);
std::vector<double> v4sigmalapltau2_1(np), v4sigmalapltau2_2(np), v4sigmalapltau2_3(np), v4sigmalapltau2_4(np), v4sigmalapltau2_5(np), v4sigmalapltau2_6(np), v4sigmalapltau2_7(np), v4sigmalapltau2_8(np), v4sigmalapltau2_9(np), v4sigmalapltau2_10(np), v4sigmalapltau2_11(np), v4sigmalapltau2_12(np), v4sigmalapltau2_13(np), v4sigmalapltau2_14(np), v4sigmalapltau2_15(np), v4sigmalapltau2_16(np), v4sigmalapltau2_17(np), v4sigmalapltau2_18(np);
std::vector<double> v4sigmatau3_1(np), v4sigmatau3_2(np), v4sigmatau3_3(np), v4sigmatau3_4(np), v4sigmatau3_5(np), v4sigmatau3_6(np), v4sigmatau3_7(np), v4sigmatau3_8(np), v4sigmatau3_9(np), v4sigmatau3_10(np), v4sigmatau3_11(np), v4sigmatau3_12(np);
std::vector<double> v4lapl4_1(np), v4lapl4_2(np), v4lapl4_3(np), v4lapl4_4(np), v4lapl4_5(np);
std::vector<double> v4lapl3tau_1(np), v4lapl3tau_2(np), v4lapl3tau_3(np), v4lapl3tau_4(np), v4lapl3tau_5(np), v4lapl3tau_6(np), v4lapl3tau_7(np), v4lapl3tau_8(np);
std::vector<double> v4lapl2tau2_1(np), v4lapl2tau2_2(np), v4lapl2tau2_3(np), v4lapl2tau2_4(np), v4lapl2tau2_5(np), v4lapl2tau2_6(np), v4lapl2tau2_7(np), v4lapl2tau2_8(np), v4lapl2tau2_9(np);     
std::vector<double> v4lapltau3_1(np), v4lapltau3_2(np), v4lapltau3_3(np), v4lapltau3_4(np), v4lapltau3_5(np), v4lapltau3_6(np), v4lapltau3_7(np), v4lapltau3_8(np);
std::vector<double> v4tau4_1(np), v4tau4_2(np), v4tau4_3(np), v4tau4_4(np), v4tau4_5(np);

mgga_lxc(rho_up,rho_down,sigma_1,sigma_2,sigma_3,lapl_up,lapl_down,tau_up,tau_down,v4rho4_1,v4rho4_2,v4rho4_3,v4rho4_4,v4rho4_5,v4rho3sigma_1,v4rho3sigma_2,v4rho3sigma_3,v4rho3sigma_4,v4rho3sigma_5,v4rho3sigma_6,v4rho3sigma_7,v4rho3sigma_8,v4rho3sigma_9,v4rho3sigma_10,v4rho3sigma_11,v4rho3sigma_12,v4rho3lapl_1,v4rho3lapl_2,v4rho3lapl_3,v4rho3lapl_4,v4rho3lapl_5,v4rho3lapl_6,v4rho3lapl_7,v4rho3lapl_8,v4rho3tau_1,v4rho3tau_2,v4rho3tau_3,v4rho3tau_4,v4rho3tau_5,v4rho3tau_6,v4rho3tau_7,v4rho3tau_8,v4rho2sigma2_1,v4rho2sigma2_2,v4rho2sigma2_3,v4rho2sigma2_4,v4rho2sigma2_5,v4rho2sigma2_6,v4rho2sigma2_7,v4rho2sigma2_8,v4rho2sigma2_9,v4rho2sigma2_10,v4rho2sigma2_11,v4rho2sigma2_12,v4rho2sigma2_13,v4rho2sigma2_14,v4rho2sigma2_15,v4rho2sigma2_16,v4rho2sigma2_17,v4rho2sigma2_18,v4rho2sigmalapl_1,v4rho2sigmalapl_2,v4rho2sigmalapl_3,v4rho2sigmalapl_4,v4rho2sigmalapl_5,v4rho2sigmalapl_6,v4rho2sigmalapl_7,v4rho2sigmalapl_8,v4rho2sigmalapl_9,v4rho2sigmalapl_10,v4rho2sigmalapl_11,v4rho2sigmalapl_12,v4rho2sigmalapl_13,v4rho2sigmalapl_14,v4rho2sigmalapl_15,v4rho2sigmalapl_16,v4rho2sigmalapl_17,v4rho2sigmalapl_18,v4rho2sigmatau_1,v4rho2sigmatau_2,v4rho2sigmatau_3,v4rho2sigmatau_4,v4rho2sigmatau_5,v4rho2sigmatau_6,v4rho2sigmatau_7,v4rho2sigmatau_8,v4rho2sigmatau_9,v4rho2sigmatau_10,v4rho2sigmatau_11,v4rho2sigmatau_12,v4rho2sigmatau_13,v4rho2sigmatau_14,v4rho2sigmatau_15,v4rho2sigmatau_16,v4rho2sigmatau_17,v4rho2sigmatau_18,v4rho2lapl2_1,v4rho2lapl2_2,v4rho2lapl2_3,v4rho2lapl2_4,v4rho2lapl2_5,v4rho2lapl2_6,v4rho2lapl2_7,v4rho2lapl2_8,v4rho2lapl2_9,v4rho2lapltau_1,v4rho2lapltau_2,v4rho2lapltau_3,v4rho2lapltau_4,v4rho2lapltau_5,v4rho2lapltau_6,v4rho2lapltau_7,v4rho2lapltau_8,v4rho2lapltau_9,v4rho2lapltau_10,v4rho2lapltau_11,v4rho2lapltau_12,v4rho2tau2_1,v4rho2tau2_2,v4rho2tau2_3,v4rho2tau2_4,v4rho2tau2_5,v4rho2tau2_6,v4rho2tau2_7,v4rho2tau2_8,v4rho2tau2_9,v4rhosigma3_1,v4rhosigma3_2,v4rhosigma3_3,v4rhosigma3_4,v4rhosigma3_5,v4rhosigma3_6,v4rhosigma3_7,v4rhosigma3_8,v4rhosigma3_9,v4rhosigma3_10,v4rhosigma3_11,v4rhosigma3_12,v4rhosigma3_13,v4rhosigma3_14,v4rhosigma3_15,v4rhosigma3_16,v4rhosigma3_17,v4rhosigma3_18,v4rhosigma3_19,v4rhosigma3_20,v4rhosigma2lapl_1,v4rhosigma2lapl_2,v4rhosigma2lapl_3,v4rhosigma2lapl_4,v4rhosigma2lapl_5,v4rhosigma2lapl_6,v4rhosigma2lapl_7,v4rhosigma2lapl_8,v4rhosigma2lapl_9,v4rhosigma2lapl_10,v4rhosigma2lapl_11,v4rhosigma2lapl_12,v4rhosigma2lapl_13,v4rhosigma2lapl_14,v4rhosigma2lapl_15,v4rhosigma2lapl_16,v4rhosigma2lapl_17,v4rhosigma2lapl_18,v4rhosigma2lapl_19,v4rhosigma2lapl_20,v4rhosigma2lapl_21,v4rhosigma2lapl_22,v4rhosigma2lapl_23,v4rhosigma2lapl_24,v4rhosigma2tau_1,v4rhosigma2tau_2,v4rhosigma2tau_3,v4rhosigma2tau_4,v4rhosigma2tau_5,v4rhosigma2tau_6,v4rhosigma2tau_7,v4rhosigma2tau_8,v4rhosigma2tau_9,v4rhosigma2tau_10,v4rhosigma2tau_11,v4rhosigma2tau_12,v4rhosigma2tau_13,v4rhosigma2tau_14,v4rhosigma2tau_15,v4rhosigma2tau_16,v4rhosigma2tau_17,v4rhosigma2tau_18,v4rhosigma2tau_19,v4rhosigma2tau_20,v4rhosigma2tau_21,v4rhosigma2tau_22,v4rhosigma2tau_23,v4rhosigma2tau_24,v4rhosigmalapl2_1,v4rhosigmalapl2_2,v4rhosigmalapl2_3,v4rhosigmalapl2_4,v4rhosigmalapl2_5,v4rhosigmalapl2_6,v4rhosigmalapl2_7,v4rhosigmalapl2_8,v4rhosigmalapl2_9,v4rhosigmalapl2_10,v4rhosigmalapl2_11,v4rhosigmalapl2_12,v4rhosigmalapl2_13,v4rhosigmalapl2_14,v4rhosigmalapl2_15,v4rhosigmalapl2_16,v4rhosigmalapl2_17,v4rhosigmalapl2_18,v4rhosigmalapltau_1,v4rhosigmalapltau_2,v4rhosigmalapltau_3,v4rhosigmalapltau_4,v4rhosigmalapltau_5,v4rhosigmalapltau_6,v4rhosigmalapltau_7,v4rhosigmalapltau_8,v4rhosigmalapltau_9,v4rhosigmalapltau_10,v4rhosigmalapltau_11,v4rhosigmalapltau_12,v4rhosigmalapltau_13,v4rhosigmalapltau_14,v4rhosigmalapltau_15,v4rhosigmalapltau_16,v4rhosigmalapltau_17,v4rhosigmalapltau_18,v4rhosigmalapltau_19,v4rhosigmalapltau_20,v4rhosigmalapltau_21,v4rhosigmalapltau_22,v4rhosigmalapltau_23,v4rhosigmalapltau_24,v4rhosigmatau2_1,v4rhosigmatau2_2,v4rhosigmatau2_3,v4rhosigmatau2_4,v4rhosigmatau2_5,v4rhosigmatau2_6,v4rhosigmatau2_7,v4rhosigmatau2_8,v4rhosigmatau2_9,v4rhosigmatau2_10,v4rhosigmatau2_11,v4rhosigmatau2_12,v4rhosigmatau2_13,v4rhosigmatau2_14,v4rhosigmatau2_15,v4rhosigmatau2_16,v4rhosigmatau2_17,v4rhosigmatau2_18,v4rholapl3_1,v4rholapl3_2,v4rholapl3_3,v4rholapl3_4,v4rholapl3_5,v4rholapl3_6,v4rholapl3_7,v4rholapl3_8,v4rholapl2tau_1,v4rholapl2tau_2,v4rholapl2tau_3,v4rholapl2tau_4,v4rholapl2tau_5,v4rholapl2tau_6,v4rholapl2tau_7,v4rholapl2tau_8,v4rholapl2tau_9,v4rholapl2tau_10,v4rholapl2tau_11,v4rholapl2tau_12,v4rholapltau2_1,v4rholapltau2_2,v4rholapltau2_3,v4rholapltau2_4,v4rholapltau2_5,v4rholapltau2_6,v4rholapltau2_7,v4rholapltau2_8,v4rholapltau2_9,v4rholapltau2_10,v4rholapltau2_11,v4rholapltau2_12,v4rhotau3_1,v4rhotau3_2,v4rhotau3_3,v4rhotau3_4,v4rhotau3_5,v4rhotau3_6,v4rhotau3_7,v4rhotau3_8,v4sigma4_1,v4sigma4_2,v4sigma4_3,v4sigma4_4,v4sigma4_5,v4sigma4_6,v4sigma4_7,v4sigma4_8,v4sigma4_9,v4sigma4_10,v4sigma4_11,v4sigma4_12,v4sigma4_13,v4sigma4_14,v4sigma4_15,v4sigma3lapl_1,v4sigma3lapl_2,v4sigma3lapl_3,v4sigma3lapl_4,v4sigma3lapl_5,v4sigma3lapl_6,v4sigma3lapl_7,v4sigma3lapl_8,v4sigma3lapl_9,v4sigma3lapl_10,v4sigma3lapl_11,v4sigma3lapl_12,v4sigma3lapl_13,v4sigma3lapl_14,v4sigma3lapl_15,v4sigma3lapl_16,v4sigma3lapl_17,v4sigma3lapl_18,v4sigma3lapl_19,v4sigma3lapl_20,v4sigma3tau_1,v4sigma3tau_2,v4sigma3tau_3,v4sigma3tau_4,v4sigma3tau_5,v4sigma3tau_6,v4sigma3tau_7,v4sigma3tau_8,v4sigma3tau_9,v4sigma3tau_10,v4sigma3tau_11,v4sigma3tau_12,v4sigma3tau_13,v4sigma3tau_14,v4sigma3tau_15,v4sigma3tau_16,v4sigma3tau_17,v4sigma3tau_18,v4sigma3tau_19,v4sigma3tau_20,v4sigma2lapl2_1,v4sigma2lapl2_2,v4sigma2lapl2_3,v4sigma2lapl2_4,v4sigma2lapl2_5,v4sigma2lapl2_6,v4sigma2lapl2_7,v4sigma2lapl2_8,v4sigma2lapl2_9,v4sigma2lapl2_10,v4sigma2lapl2_11,v4sigma2lapl2_12,v4sigma2lapl2_13,v4sigma2lapl2_14,v4sigma2lapl2_15,v4sigma2lapl2_16,v4sigma2lapl2_17,v4sigma2lapl2_18,v4sigma2lapltau_1,v4sigma2lapltau_2,v4sigma2lapltau_3,v4sigma2lapltau_4,v4sigma2lapltau_5,v4sigma2lapltau_6,v4sigma2lapltau_7,v4sigma2lapltau_8,v4sigma2lapltau_9,v4sigma2lapltau_10,v4sigma2lapltau_11,v4sigma2lapltau_12,v4sigma2lapltau_13,v4sigma2lapltau_14,v4sigma2lapltau_15,v4sigma2lapltau_16,v4sigma2lapltau_17,v4sigma2lapltau_18,v4sigma2lapltau_19,v4sigma2lapltau_20,v4sigma2lapltau_21,v4sigma2lapltau_22,v4sigma2lapltau_23,v4sigma2lapltau_24,v4sigma2tau2_1,v4sigma2tau2_2,v4sigma2tau2_3,v4sigma2tau2_4,v4sigma2tau2_5,v4sigma2tau2_6,v4sigma2tau2_7,v4sigma2tau2_8,v4sigma2tau2_9,v4sigma2tau2_10,v4sigma2tau2_11,v4sigma2tau2_12,v4sigma2tau2_13,v4sigma2tau2_14,v4sigma2tau2_15,v4sigma2tau2_16,v4sigma2tau2_17,v4sigma2tau2_18,v4sigmalapl3_1,v4sigmalapl3_2,v4sigmalapl3_3,v4sigmalapl3_4,v4sigmalapl3_5,v4sigmalapl3_6,v4sigmalapl3_7,v4sigmalapl3_8,v4sigmalapl3_9,v4sigmalapl3_10,v4sigmalapl3_11,v4sigmalapl3_12,v4sigmalapl2tau_1,v4sigmalapl2tau_2,v4sigmalapl2tau_3,v4sigmalapl2tau_4,v4sigmalapl2tau_5,v4sigmalapl2tau_6,v4sigmalapl2tau_7,v4sigmalapl2tau_8,v4sigmalapl2tau_9,v4sigmalapl2tau_10,v4sigmalapl2tau_11,v4sigmalapl2tau_12,v4sigmalapl2tau_13,v4sigmalapl2tau_14,v4sigmalapl2tau_15,v4sigmalapl2tau_16,v4sigmalapl2tau_17,v4sigmalapl2tau_18,v4sigmalapltau2_1,v4sigmalapltau2_2,v4sigmalapltau2_3,v4sigmalapltau2_4,v4sigmalapltau2_5,v4sigmalapltau2_6,v4sigmalapltau2_7,v4sigmalapltau2_8,v4sigmalapltau2_9,v4sigmalapltau2_10,v4sigmalapltau2_11,v4sigmalapltau2_12,v4sigmalapltau2_13,v4sigmalapltau2_14,v4sigmalapltau2_15,v4sigmalapltau2_16,v4sigmalapltau2_17,v4sigmalapltau2_18,v4sigmatau3_1,v4sigmatau3_2,v4sigmatau3_3,v4sigmatau3_4,v4sigmatau3_5,v4sigmatau3_6,v4sigmatau3_7,v4sigmatau3_8,v4sigmatau3_9,v4sigmatau3_10,v4sigmatau3_11,v4sigmatau3_12,v4lapl4_1,v4lapl4_2,v4lapl4_3,v4lapl4_4,v4lapl4_5,v4lapl3tau_1,v4lapl3tau_2,v4lapl3tau_3,v4lapl3tau_4,v4lapl3tau_5,v4lapl3tau_6,v4lapl3tau_7,v4lapl3tau_8,v4lapl2tau2_1,v4lapl2tau2_2,v4lapl2tau2_3,v4lapl2tau2_4,v4lapl2tau2_5,v4lapl2tau2_6,v4lapl2tau2_7,v4lapl2tau2_8,v4lapl2tau2_9,v4lapltau3_1,v4lapltau3_2,v4lapltau3_3,v4lapltau3_4,v4lapltau3_5,v4lapltau3_6,v4lapltau3_7,v4lapltau3_8,v4tau4_1,v4tau4_2,v4tau4_3,v4tau4_4,v4tau4_5);

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

    std::cout << "MGGA V3Tau3_4: ";
    for (auto k : v3tau3_4)
        std::cout << k << " ";

    std::cout << "MGGA V4tau4_5: ";
    for (auto l : v4tau4_5)
        std::cout << l << " ";
}

  
}//end of namespace NCXC

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

