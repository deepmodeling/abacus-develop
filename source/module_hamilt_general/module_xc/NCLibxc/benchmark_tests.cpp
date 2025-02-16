#include "NCLibxc.h"
#include <iostream>
#include <iomanip>

///////////////////////////////////////////////////////////////////////////////////
// collinear limit test for GGA
void NCLibxc::gga_collinear_test()
{
    try {
        NCLibxc::print_NCLibxc();
        
        // Sample 0-th order density
        std::vector<double> n = {1.0, 1.0, 2.0};
        std::vector<double> mx = {0.0, 0.0, 0.0};
        std::vector<double> my = {0.0, 0.0, 0.0};
        std::vector<double> mz = {0.1, -0.1, 0.1};

        // Sample gradients 
        std::vector<double> gradx_n = {0.1, 0.1, 0.1};
        std::vector<double> grady_n = {0.2, 0.2, 0.5};
        std::vector<double> gradz_n = {0.3, 0.3, 0.6};

        std::vector<double> gradx_mx = {0.0, 0.0, 0.0};
        std::vector<double> grady_mx = {0.0, 0.0, 0.0};
        std::vector<double> gradz_mx = {0.0, 0.0, 0.0};

        std::vector<double> gradx_my = {0.0, 0.0, 0.0};
        std::vector<double> grady_my = {0.0, 0.0, 0.0};
        std::vector<double> gradz_my = {0.0, 0.0, 0.0};

        std::vector<double> gradx_mz = {0.2, -0.2, 0.1};
        std::vector<double> grady_mz = {0.3, -0.3, 0.1};
        std::vector<double> gradz_mz = {0.1, -0.1, 0.1};

        // Second derivatives 
        std::vector<double> grad2xx_n = {0.5, 0.5, 0.0};
        std::vector<double> grad2yy_n = {0.2, 0.2, 0.0};
        std::vector<double> grad2zz_n = {0.3, 0.3, 0.0};
        std::vector<double> grad2xy_n = {0.7, 0.7, 0.0};
        std::vector<double> grad2yz_n = {0.8, 0.8, 0.0};
        std::vector<double> grad2xz_n = {0.9, 0.9, 0.0};
        std::vector<double> grad2xx_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2yy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2zz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2xy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2yz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2xz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad2xx_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2yy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2zz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2xy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2yz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2xz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad2xx_mz = {0.01, -0.01, 0.0};
        std::vector<double> grad2yy_mz = {0.02, -0.02, 0.0};
        std::vector<double> grad2zz_mz = {0.03, -0.03, 0.0};
        std::vector<double> grad2xy_mz = {0.06, -0.06, 0.0};
        std::vector<double> grad2yz_mz = {0.07, -0.07, 0.0};
        std::vector<double> grad2xz_mz = {0.08, -0.08, 0.0};

        // Third-order derivatives for n
        std::vector<double> grad3xxx_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3xxy_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3xxz_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3xyy_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3xyz_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3xzz_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3yyy_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3yyz_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3yzz_n = {0.1, 0.1, 0.2};
        std::vector<double> grad3zzz_n = {0.1, 0.1, 0.2};

        // Third-order derivatives for mx (all zeros)
        std::vector<double> grad3xxx_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3xxy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3xxz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3xyy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3xyz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3xzz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3yyy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3yyz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3yzz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad3zzz_mx = {0.0, 0.0, 0.0};

        // Third-order derivatives for my (all zeros)
        std::vector<double> grad3xxx_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3xxy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3xxz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3xyy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3xyz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3xzz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3yyy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3yyz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3yzz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad3zzz_my = {0.0, 0.0, 0.0};

        // Third-order derivatives for mz
        std::vector<double> grad3xxx_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3xxy_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3xxz_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3xyy_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3xyz_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3xzz_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3yyy_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3yyz_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3yzz_mz = {0.05, -0.05, 0.05};
        std::vector<double> grad3zzz_mz = {0.05, -0.05, 0.05};

        // Fourth-order derivatives for n
        std::vector<double> grad4xxxx_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4xxxy_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4xxxz_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4xxyy_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4xxyz_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4xxzz_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4xyyy_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4xyyz_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4xyzz_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4xzzz_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4yyyy_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4yyyz_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4yyzz_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4yzzz_n = {0.05, 0.05, 0.10};
        std::vector<double> grad4zzzz_n = {0.05, 0.05, 0.10};

        // Fourth-order derivatives for mx (all zeros)
        std::vector<double> grad4xxxx_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4xxxy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4xxxz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4xxyy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4xxyz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4xxzz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4xyyy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4xyyz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4xyzz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4xzzz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4yyyy_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4yyyz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4yyzz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4yzzz_mx = {0.0, 0.0, 0.0};
        std::vector<double> grad4zzzz_mx = {0.0, 0.0, 0.0};

        // Fourth-order derivatives for my (all zeros)
        std::vector<double> grad4xxxx_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4xxxy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4xxxz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4xxyy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4xxyz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4xxzz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4xyyy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4xyyz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4xyzz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4xzzz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4yyyy_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4yyyz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4yyzz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4yzzz_my = {0.0, 0.0, 0.0};
        std::vector<double> grad4zzzz_my = {0.0, 0.0, 0.0};

        // Fourth-order derivatives for mz
        std::vector<double> grad4xxxx_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4xxxy_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4xxxz_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4xxyy_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4xxyz_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4xxzz_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4xyyy_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4xyyz_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4xyzz_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4xzzz_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4yyyy_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4yyyz_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4yyzz_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4yzzz_mz = {0.02, -0.02, 0.02};
        std::vector<double> grad4zzzz_mz = {0.02, -0.02, 0.02};


        int xc_id = 106; // Example XC functional ID

        // Call gga_mc
        auto [E_GGA_MC, V_GGA_MC] = NCLibxc::gga_mc(
    xc_id, n, mx, my, mz,
    gradx_n, grady_n, gradz_n,
    gradx_mx, grady_mx, gradz_mx,
    gradx_my, grady_my, gradz_my,
    gradx_mz, grady_mz, gradz_mz,
    grad2xx_n, grad2yy_n, grad2zz_n, grad2xy_n, grad2yz_n, grad2xz_n,
    grad2xx_mx, grad2yy_mx, grad2zz_mx, grad2xy_mx, grad2yz_mx, grad2xz_mx,
    grad2xx_my, grad2yy_my, grad2zz_my, grad2xy_my, grad2yz_my, grad2xz_my,
    grad2xx_mz, grad2yy_mz, grad2zz_mz, grad2xy_mz, grad2yz_mz, grad2xz_mz,
    grad3xxx_n, grad3xxy_n, grad3xxz_n, grad3xyy_n, grad3xyz_n, grad3xzz_n,
    grad3yyy_n, grad3yyz_n, grad3yzz_n, grad3zzz_n,
    grad3xxx_mx, grad3xxy_mx, grad3xxz_mx, grad3xyy_mx, grad3xyz_mx, grad3xzz_mx,
    grad3yyy_mx, grad3yyz_mx, grad3yzz_mx, grad3zzz_mx,
    grad3xxx_my, grad3xxy_my, grad3xxz_my, grad3xyy_my, grad3xyz_my, grad3xzz_my,
    grad3yyy_my, grad3yyz_my, grad3yzz_my, grad3zzz_my,
    grad3xxx_mz, grad3xxy_mz, grad3xxz_mz, grad3xyy_mz, grad3xyz_mz, grad3xzz_mz,
    grad3yyy_mz, grad3yyz_mz, grad3yzz_mz, grad3zzz_mz,
    grad4xxxx_n, grad4xxxy_n, grad4xxxz_n, grad4xxyy_n, grad4xxyz_n, grad4xxzz_n,
    grad4xyyy_n, grad4xyyz_n, grad4xyzz_n, grad4xzzz_n, grad4yyyy_n, grad4yyyz_n,
    grad4yyzz_n, grad4yzzz_n, grad4zzzz_n,
    grad4xxxx_mx, grad4xxxy_mx, grad4xxxz_mx, grad4xxyy_mx, grad4xxyz_mx, grad4xxzz_mx,
    grad4xyyy_mx, grad4xyyz_mx, grad4xyzz_mx, grad4xzzz_mx, grad4yyyy_mx, grad4yyyz_mx,
    grad4yyzz_mx, grad4yzzz_mx, grad4zzzz_mx,
    grad4xxxx_my, grad4xxxy_my, grad4xxxz_my, grad4xxyy_my, grad4xxyz_my, grad4xxzz_my,
    grad4xyyy_my, grad4xyyz_my, grad4xyzz_my, grad4xzzz_my, grad4yyyy_my, grad4yyyz_my,
    grad4yyzz_my, grad4yzzz_my, grad4zzzz_my,
    grad4xxxx_mz, grad4xxxy_mz, grad4xxxz_mz, grad4xxyy_mz, grad4xxyz_mz, grad4xxzz_mz,
    grad4xyyy_mz, grad4xyyz_mz, grad4xyzz_mz, grad4xzzz_mz, grad4yyyy_mz, grad4yyyz_mz,
    grad4yyzz_mz, grad4yzzz_mz, grad4zzzz_mz
);

        // Print results from gga_mc
        std::cout << std::fixed << std::setprecision(8);
        std::cout << "Total E for each real-space grid point (GGA MC):" << std::endl;
        for (const auto &e : E_GGA_MC)
            std::cout << e << " ";
        std::cout << std::endl;
        std::cout << "Total V for each real-space grid point (GGA MC):" << std::endl;
        std::cout << "[ ";
        for (const auto &mat : V_GGA_MC) {
            std::cout << "[ ";
            std::cout << mat[0][0] << " " << mat[0][1] << " "
                      << mat[1][0] << " " << mat[1][1];
            std::cout << "] ";
        }
        std::cout << "]" << std::endl;
        // Prepare inputs for postlibxc_gga
        size_t size = n.size();
        std::vector<double> rho0(size), rho1(size);
        for (size_t i = 0; i < size; ++i) {
            rho0[i] = (n[i] + mz[i]) / 2.0;
            rho1[i] = (n[i] - mz[i]) / 2.0;
        }

        // Gradient calculations
        // First derivatives
        std::vector<double> gradx_rho0(size), grady_rho0(size), gradz_rho0(size);
        std::vector<double> gradx_rho1(size), grady_rho1(size), gradz_rho1(size);

        for (size_t i = 0; i < size; ++i) {
            gradx_rho0[i] = (gradx_n[i] + gradx_mz[i]) / 2.0;
            grady_rho0[i] = (grady_n[i] + grady_mz[i]) / 2.0;
            gradz_rho0[i] = (gradz_n[i] + gradz_mz[i]) / 2.0;

            gradx_rho1[i] = (gradx_n[i] - gradx_mz[i]) / 2.0;
            grady_rho1[i] = (grady_n[i] - grady_mz[i]) / 2.0;
            gradz_rho1[i] = (gradz_n[i] - gradz_mz[i]) / 2.0;
        }

        // Second derivatives
        std::vector<double> grad2xx_rho0(size), grad2yy_rho0(size), grad2zz_rho0(size);
        std::vector<double> grad2xy_rho0(size), grad2yz_rho0(size), grad2xz_rho0(size);
        std::vector<double> grad2xx_rho1(size), grad2yy_rho1(size), grad2zz_rho1(size);
        std::vector<double> grad2xy_rho1(size), grad2yz_rho1(size), grad2xz_rho1(size);

        for (size_t i = 0; i < size; ++i) {
            grad2xx_rho0[i] = (grad2xx_n[i] + grad2xx_mz[i]) / 2.0;
            grad2yy_rho0[i] = (grad2yy_n[i] + grad2yy_mz[i]) / 2.0;
            grad2zz_rho0[i] = (grad2zz_n[i] + grad2zz_mz[i]) / 2.0;
            grad2xy_rho0[i] = (grad2xy_n[i] + grad2xy_mz[i]) / 2.0;
            grad2yz_rho0[i] = (grad2yz_n[i] + grad2yz_mz[i]) / 2.0;
            grad2xz_rho0[i] = (grad2xz_n[i] + grad2xz_mz[i]) / 2.0;

            grad2xx_rho1[i] = (grad2xx_n[i] - grad2xx_mz[i]) / 2.0;
            grad2yy_rho1[i] = (grad2yy_n[i] - grad2yy_mz[i]) / 2.0;
            grad2zz_rho1[i] = (grad2zz_n[i] - grad2zz_mz[i]) / 2.0;
            grad2xy_rho1[i] = (grad2xy_n[i] - grad2xy_mz[i]) / 2.0;
            grad2yz_rho1[i] = (grad2yz_n[i] - grad2yz_mz[i]) / 2.0;
            grad2xz_rho1[i] = (grad2xz_n[i] - grad2xz_mz[i]) / 2.0;
        }

        // Third-order gradients
        std::vector<double> grad3xxx_rho0(size), grad3xxy_rho0(size), grad3xxz_rho0(size);
        std::vector<double> grad3xyy_rho0(size), grad3xyz_rho0(size), grad3xzz_rho0(size);
        std::vector<double> grad3yyy_rho0(size), grad3yyz_rho0(size), grad3yzz_rho0(size);
        std::vector<double> grad3zzz_rho0(size);

        std::vector<double> grad3xxx_rho1(size), grad3xxy_rho1(size), grad3xxz_rho1(size);
        std::vector<double> grad3xyy_rho1(size), grad3xyz_rho1(size), grad3xzz_rho1(size);
        std::vector<double> grad3yyy_rho1(size), grad3yyz_rho1(size), grad3yzz_rho1(size);
        std::vector<double> grad3zzz_rho1(size);

        for (size_t i = 0; i < size; ++i) {
            // Calculate third-order gradients for rho0
            grad3xxx_rho0[i] = (grad3xxx_n[i] + grad3xxx_mz[i]) / 2.0;
            grad3xxy_rho0[i] = (grad3xxy_n[i] + grad3xxy_mz[i]) / 2.0;
            grad3xxz_rho0[i] = (grad3xxz_n[i] + grad3xxz_mz[i]) / 2.0;
            grad3xyy_rho0[i] = (grad3xyy_n[i] + grad3xyy_mz[i]) / 2.0;
            grad3xyz_rho0[i] = (grad3xyz_n[i] + grad3xyz_mz[i]) / 2.0;
            grad3xzz_rho0[i] = (grad3xzz_n[i] + grad3xzz_mz[i]) / 2.0;
            grad3yyy_rho0[i] = (grad3yyy_n[i] + grad3yyy_mz[i]) / 2.0;
            grad3yyz_rho0[i] = (grad3yyz_n[i] + grad3yyz_mz[i]) / 2.0;
            grad3yzz_rho0[i] = (grad3yzz_n[i] + grad3yzz_mz[i]) / 2.0;
            grad3zzz_rho0[i] = (grad3zzz_n[i] + grad3zzz_mz[i]) / 2.0;

            // Calculate third-order gradients for rho1
            grad3xxx_rho1[i] = (grad3xxx_n[i] - grad3xxx_mz[i]) / 2.0;
            grad3xxy_rho1[i] = (grad3xxy_n[i] - grad3xxy_mz[i]) / 2.0;
            grad3xxz_rho1[i] = (grad3xxz_n[i] - grad3xxz_mz[i]) / 2.0;
            grad3xyy_rho1[i] = (grad3xyy_n[i] - grad3xyy_mz[i]) / 2.0;
            grad3xyz_rho1[i] = (grad3xyz_n[i] - grad3xyz_mz[i]) / 2.0;
            grad3xzz_rho1[i] = (grad3xzz_n[i] - grad3xzz_mz[i]) / 2.0;
            grad3yyy_rho1[i] = (grad3yyy_n[i] - grad3yyy_mz[i]) / 2.0;
            grad3yyz_rho1[i] = (grad3yyz_n[i] - grad3yyz_mz[i]) / 2.0;
            grad3yzz_rho1[i] = (grad3yzz_n[i] - grad3yzz_mz[i]) / 2.0;
            grad3zzz_rho1[i] = (grad3zzz_n[i] - grad3zzz_mz[i]) / 2.0;
        }

        // Fourth-order gradients
        std::vector<double> grad4xxxx_rho0(size), grad4xxxy_rho0(size), grad4xxxz_rho0(size);
        std::vector<double> grad4xxyy_rho0(size), grad4xxyz_rho0(size), grad4xxzz_rho0(size);
        std::vector<double> grad4xyyy_rho0(size), grad4xyyz_rho0(size), grad4xyzz_rho0(size);
        std::vector<double> grad4xzzz_rho0(size), grad4yyyy_rho0(size), grad4yyyz_rho0(size);
        std::vector<double> grad4yyzz_rho0(size), grad4yzzz_rho0(size), grad4zzzz_rho0(size);

        std::vector<double> grad4xxxx_rho1(size), grad4xxxy_rho1(size), grad4xxxz_rho1(size);
        std::vector<double> grad4xxyy_rho1(size), grad4xxyz_rho1(size), grad4xxzz_rho1(size);
        std::vector<double> grad4xyyy_rho1(size), grad4xyyz_rho1(size), grad4xyzz_rho1(size);
        std::vector<double> grad4xzzz_rho1(size), grad4yyyy_rho1(size), grad4yyyz_rho1(size);
        std::vector<double> grad4yyzz_rho1(size), grad4yzzz_rho1(size), grad4zzzz_rho1(size);

        for (size_t i = 0; i < size; ++i) {
            // Calculate fourth-order gradients for rho0
            grad4xxxx_rho0[i] = (grad4xxxx_n[i] + grad4xxxx_mz[i]) / 2.0;
            grad4xxxy_rho0[i] = (grad4xxxy_n[i] + grad4xxxy_mz[i]) / 2.0;
            grad4xxxz_rho0[i] = (grad4xxxz_n[i] + grad4xxxz_mz[i]) / 2.0;
            grad4xxyy_rho0[i] = (grad4xxyy_n[i] + grad4xxyy_mz[i]) / 2.0;
            grad4xxyz_rho0[i] = (grad4xxyz_n[i] + grad4xxyz_mz[i]) / 2.0;
            grad4xxzz_rho0[i] = (grad4xxzz_n[i] + grad4xxzz_mz[i]) / 2.0;
            grad4xyyy_rho0[i] = (grad4xyyy_n[i] + grad4xyyy_mz[i]) / 2.0;
            grad4xyyz_rho0[i] = (grad4xyyz_n[i] + grad4xyyz_mz[i]) / 2.0;
            grad4xyzz_rho0[i] = (grad4xyzz_n[i] + grad4xyzz_mz[i]) / 2.0;
            grad4xzzz_rho0[i] = (grad4xzzz_n[i] + grad4xzzz_mz[i]) / 2.0;
            grad4yyyy_rho0[i] = (grad4yyyy_n[i] + grad4yyyy_mz[i]) / 2.0;
            grad4yyyz_rho0[i] = (grad4yyyz_n[i] + grad4yyyz_mz[i]) / 2.0;
            grad4yyzz_rho0[i] = (grad4yyzz_n[i] + grad4yyzz_mz[i]) / 2.0;
            grad4yzzz_rho0[i] = (grad4yzzz_n[i] + grad4yzzz_mz[i]) / 2.0;
            grad4zzzz_rho0[i] = (grad4zzzz_n[i] + grad4zzzz_mz[i]) / 2.0;

            // Calculate fourth-order gradients for rho1
            grad4xxxx_rho1[i] = (grad4xxxx_n[i] - grad4xxxx_mz[i]) / 2.0;
            grad4xxxy_rho1[i] = (grad4xxxy_n[i] - grad4xxxy_mz[i]) / 2.0;
            grad4xxxz_rho1[i] = (grad4xxxz_n[i] - grad4xxxz_mz[i]) / 2.0;
            grad4xxyy_rho1[i] = (grad4xxyy_n[i] - grad4xxyy_mz[i]) / 2.0;
            grad4xxyz_rho1[i] = (grad4xxyz_n[i] - grad4xxyz_mz[i]) / 2.0;
            grad4xxzz_rho1[i] = (grad4xxzz_n[i] - grad4xxzz_mz[i]) / 2.0;
            grad4xyyy_rho1[i] = (grad4xyyy_n[i] - grad4xyyy_mz[i]) / 2.0;
            grad4xyyz_rho1[i] = (grad4xyyz_n[i] - grad4xyyz_mz[i]) / 2.0;
            grad4xyzz_rho1[i] = (grad4xyzz_n[i] - grad4xyzz_mz[i]) / 2.0;
            grad4xzzz_rho1[i] = (grad4xzzz_n[i] - grad4xzzz_mz[i]) / 2.0;
            grad4yyyy_rho1[i] = (grad4yyyy_n[i] - grad4yyyy_mz[i]) / 2.0;
            grad4yyyz_rho1[i] = (grad4yyyz_n[i] - grad4yyyz_mz[i]) / 2.0;
            grad4yyzz_rho1[i] = (grad4yyzz_n[i] - grad4yyzz_mz[i]) / 2.0;
            grad4yzzz_rho1[i] = (grad4yzzz_n[i] - grad4yzzz_mz[i]) / 2.0;
            grad4zzzz_rho1[i] = (grad4zzzz_n[i] - grad4zzzz_mz[i]) / 2.0;
        }

        // Call postlibxc_gga
        std::vector<double> e, v1, v2, f1, f2, f3;
        NCLibxc::postlibxc_gga(
    xc_id,
    rho0, rho1,
    gradx_rho0, grady_rho0, gradz_rho0,
    gradx_rho1, grady_rho1, gradz_rho1,
    grad2xx_rho0, grad2yy_rho0, grad2zz_rho0,
    grad2xy_rho0, grad2yz_rho0, grad2xz_rho0,
    grad2xx_rho1, grad2yy_rho1, grad2zz_rho1,
    grad2xy_rho1, grad2yz_rho1, grad2xz_rho1,
    grad3xxx_rho0, grad3xxy_rho0, grad3xxz_rho0,
    grad3xyy_rho0, grad3xyz_rho0, grad3xzz_rho0,
    grad3yyy_rho0, grad3yyz_rho0, grad3yzz_rho0,
    grad3zzz_rho0,
    grad3xxx_rho1, grad3xxy_rho1, grad3xxz_rho1,
    grad3xyy_rho1, grad3xyz_rho1, grad3xzz_rho1,
    grad3yyy_rho1, grad3yyz_rho1, grad3yzz_rho1,
    grad3zzz_rho1,
    grad4xxxx_rho0, grad4xxxy_rho0, grad4xxxz_rho0,
    grad4xxyy_rho0, grad4xxyz_rho0, grad4xxzz_rho0,
    grad4xyyy_rho0, grad4xyyz_rho0, grad4xyzz_rho0,
    grad4xzzz_rho0, grad4yyyy_rho0, grad4yyyz_rho0,
    grad4yyzz_rho0, grad4yzzz_rho0, grad4zzzz_rho0,
    grad4xxxx_rho1, grad4xxxy_rho1, grad4xxxz_rho1,
    grad4xxyy_rho1, grad4xxyz_rho1, grad4xxzz_rho1,
    grad4xyyy_rho1, grad4xyyz_rho1, grad4xyzz_rho1,
    grad4xzzz_rho1, grad4yyyy_rho1, grad4yyyz_rho1,
    grad4yyzz_rho1, grad4yzzz_rho1, grad4zzzz_rho1,
    e, v1, v2, f1, f2, f3
);

        // Print results from postlibxc_gga
       std::cout << "Total E for each real-space grid point (postlibxc_gga):" << std::endl;
        for (const auto &e : e)
            std::cout << e << " ";
        std::cout << std::endl;

        // Similarly, print v1, v2, f1, f2, f3 if needed
        std::cout << "Total V1 for each real-space grid point (postlibxc_gga):" << std::endl;
        for (const auto &v : v1)
            std::cout << v << " ";
        std::cout << std::endl;

        std::cout << "Total V2 for each real-space grid point (postlibxc_gga):" << std::endl;
        for (const auto &v : v2)
            std::cout << v << " ";
        std::cout << std::endl;

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}


///////////////////////////////////////////////////////////////////////////////////
//consistency between lda_mc and lda_lc
void NCLibxc::lda_mc_lc_test()
{
try {
        NCLibxc::print_NCLibxc();
        // example input data
        std::vector<double> n = {1.0, 1.0, 1.0};
        std::vector<double> mx = {0.1, 0.1, 0.0};
        std::vector<double> my = {0.0, 0.1, 0.1414};
        std::vector<double> mz = {0.1, 0.0, 0.0};
        int xc_id = 1; // for example, set xc_id to 1(lda)

        auto [E_MC, V_MC] = NCLibxc::lda_mc(xc_id, n, mx, my, mz);

        std::cout << "Total E for each real-space grid point(LDA_MC):" << std::endl;
        for (const auto &e : E_MC)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "Total V for each real-space grid point(LDA_MC):" << std::endl;
        for (const auto &matrix : V_MC) {
            for (const auto &row : matrix) {
                for (const auto &elem : row) {
                    std::cout << elem << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        auto [E_LC, V_LC] = NCLibxc::lda_lc(xc_id, n, mx, my, mz);

        std::cout << "Total E for each real-space grid point(LDA_LC):" << std::endl;
        for (const auto &e : E_LC)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "Total V for each real-space grid point(LDA_LC):" << std::endl;
        for (const auto &matrix : V_LC)
        {
            for (const auto &row : matrix)
            {
                for (const auto &elem : row)
                {
                    std::cout << elem << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}

//consistency between lda_mc and postlibxc_lda(collinear)
void NCLibxc::lda_mc_collinear_test()
{
try {
        NCLibxc::print_NCLibxc();
        // example input data
        std::vector<double> n = {1.0, 1.0, 1.0};
        std::vector<double> mx = {0.0, 0.0, 0.0};
        std::vector<double> my = {0.0, 0.0, 0.0};
        std::vector<double> mz = {0.1, 0.2, 0.3};
        int xc_id = 1; // for example, set xc_id to 1(lda)

        auto [E_MC, V_MC] = NCLibxc::lda_mc(xc_id, n, mx, my, mz);
        
        std::cout << "Total E for each real-space grid point(LDA_MC):" << std::endl;
        for (const auto &e : E_MC)
            std::cout << e << " ";
        std::cout << std::endl;

        std::cout << "Total V for each real-space grid point(LDA_MC):" << std::endl;
        for (const auto &matrix : V_MC) {
            for (const auto &row : matrix) {
                for (const auto &elem : row) {
                    std::cout << elem << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

    
    std::vector<double> e;
    std::vector<double> v1, v2, f1, f2, f3;
    
    size_t np = n.size();

    std::vector<double> rho_up(np,0.0) ;     
    std::vector<double> rho_down(np,0.0);

    for(size_t i=0;i<np;i++){
        rho_up[i] = (n[i] + mz[i]) / 2.0;
        rho_down[i] = (n[i] - mz[i]) / 2.0;
    }   

    NCLibxc::postlibxc_lda(xc_id, rho_up, rho_down, e, v1, v2, f1, f2, f3);

    std::cout << "Results from postlibxc_lda:" << std::endl;

    std::cout << "e:" << std::endl;
    for (const auto &val : e)
        std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "v1:" << std::endl;
    for (const auto &val : v1)
        std::cout << val << " ";
    std::cout << std::endl;

    std::cout << "v2:" << std::endl;
    for (const auto &val : v2)
        std::cout << val << " ";
    std::cout << std::endl;

        
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
    }
}

// xc local torque from gga, prove non-vanishing torque
void NCLibxc::gga_local_torque_test()
{
    
    std::vector<double> n  = {1.0, 1.0, 1.0};
    std::vector<double> mx = {0.1, 0.2, 0.3};
    std::vector<double> my = {0.2, 0.3, 0.4};
    std::vector<double> mz = {0.3, 0.4, 0.5};

    // 各阶梯度均采用3个元素的零初始数组（实际计算中请设置合适值）
    std::vector<double> gradx_n(3, 0.1), grady_n(3, 0.0), gradz_n(3, 0.0);
    std::vector<double> gradx_mx(3, 0.0), grady_mx(3, 0.3), gradz_mx(3, 0.0);
    std::vector<double> gradx_my(3, 0.1), grady_my(3, 0.0), gradz_my(3, 0.0);
    std::vector<double> gradx_mz(3, 0.0), grady_mz(3, 0.0), gradz_mz(3, 0.1);

    std::vector<double> grad2xx_n(3, 0.0), grad2yy_n(3, 0.0), grad2zz_n(3, 0.0);
    std::vector<double> grad2xy_n(3, 0.0), grad2yz_n(3, 0.0), grad2xz_n(3, 0.0);
    std::vector<double> grad2xx_mx(3, 0.0), grad2yy_mx(3, 0.0), grad2zz_mx(3, 0.0);
    std::vector<double> grad2xy_mx(3, 0.0), grad2yz_mx(3, 0.0), grad2xz_mx(3, 0.0);
    std::vector<double> grad2xx_my(3, 0.0), grad2yy_my(3, 0.0), grad2zz_my(3, 0.0);
    std::vector<double> grad2xy_my(3, 0.0), grad2yz_my(3, 0.0), grad2xz_my(3, 0.0);
    std::vector<double> grad2xx_mz(3, 0.0), grad2yy_mz(3, 0.0), grad2zz_mz(3, 0.0);
    std::vector<double> grad2xy_mz(3, 0.0), grad2yz_mz(3, 0.0), grad2xz_mz(3, 0.0);

    std::vector<double> grad3xxx_n(3, 0.0), grad3xxy_n(3, 0.0), grad3xxz_n(3, 0.0), 
                        grad3xyy_n(3, 0.0), grad3xyz_n(3, 0.0), grad3xzz_n(3, 0.0);
    std::vector<double> grad3yyy_n(3, 0.0), grad3yyz_n(3, 0.0), grad3yzz_n(3, 0.0), 
                        grad3zzz_n(3, 0.0);
    std::vector<double> grad3xxx_mx(3, 0.0), grad3xxy_mx(3, 0.0), grad3xxz_mx(3, 0.0),
                        grad3xyy_mx(3, 0.0), grad3xyz_mx(3, 0.0), grad3xzz_mx(3, 0.0);
    std::vector<double> grad3yyy_mx(3, 0.0), grad3yyz_mx(3, 0.0), grad3yzz_mx(3, 0.0),
                        grad3zzz_mx(3, 0.0);
    std::vector<double> grad3xxx_my(3, 0.0), grad3xxy_my(3, 0.0), grad3xxz_my(3, 0.0),
                        grad3xyy_my(3, 0.0), grad3xyz_my(3, 0.0), grad3xzz_my(3, 0.0);
    std::vector<double> grad3yyy_my(3, 0.0), grad3yyz_my(3, 0.0), grad3yzz_my(3, 0.0),
                        grad3zzz_my(3, 0.0);
    std::vector<double> grad3xxx_mz(3, 0.0), grad3xxy_mz(3, 0.0), grad3xxz_mz(3, 0.0),
                        grad3xyy_mz(3, 0.0), grad3xyz_mz(3, 0.0), grad3xzz_mz(3, 0.0);
    std::vector<double> grad3yyy_mz(3, 0.0), grad3yyz_mz(3, 0.0), grad3yzz_mz(3, 0.0),
                        grad3zzz_mz(3, 0.0);

    std::vector<double> grad4xxxx_n(3, 0.0), grad4xxxy_n(3, 0.0), grad4xxxz_n(3, 0.0),
                        grad4xxyy_n(3, 0.0), grad4xxyz_n(3, 0.0), grad4xxzz_n(3, 0.0);
    std::vector<double> grad4xyyy_n(3, 0.0), grad4xyyz_n(3, 0.0), grad4xyzz_n(3, 0.0),
                        grad4xzzz_n(3, 0.0), grad4yyyy_n(3, 0.0), grad4yyyz_n(3, 0.0);
    std::vector<double> grad4yyzz_n(3, 0.0), grad4yzzz_n(3, 0.0), grad4zzzz_n(3, 0.0);
    std::vector<double> grad4xxxx_mx(3, 0.0), grad4xxxy_mx(3, 0.0), grad4xxxz_mx(3, 0.0),
                        grad4xxyy_mx(3, 0.0), grad4xxyz_mx(3, 0.0), grad4xxzz_mx(3, 0.0);
    std::vector<double> grad4xyyy_mx(3, 0.0), grad4xyyz_mx(3, 0.0), grad4xyzz_mx(3, 0.0),
                        grad4xzzz_mx(3, 0.0), grad4yyyy_mx(3, 0.0), grad4yyyz_mx(3, 0.0);
    std::vector<double> grad4yyzz_mx(3, 0.0), grad4yzzz_mx(3, 0.0), grad4zzzz_mx(3, 0.0);
    std::vector<double> grad4xxxx_my(3, 0.0), grad4xxxy_my(3, 0.0), grad4xxxz_my(3, 0.0),
                        grad4xxyy_my(3, 0.0), grad4xxyz_my(3, 0.0), grad4xxzz_my(3, 0.0);
    std::vector<double> grad4xyyy_my(3, 0.0), grad4xyyz_my(3, 0.0), grad4xyzz_my(3, 0.0),
                        grad4xzzz_my(3, 0.0), grad4yyyy_my(3, 0.0), grad4yyyz_my(3, 0.0);
    std::vector<double> grad4yyzz_my(3, 0.0), grad4yzzz_my(3, 0.0), grad4zzzz_my(3, 0.0);
    std::vector<double> grad4xxxx_mz(3, 0.0), grad4xxxy_mz(3, 0.0), grad4xxxz_mz(3, 0.0),
                        grad4xxyy_mz(3, 0.0), grad4xxyz_mz(3, 0.0), grad4xxzz_mz(3, 0.0);
    std::vector<double> grad4xyyy_mz(3, 0.0), grad4xyyz_mz(3, 0.0), grad4xyzz_mz(3, 0.0),
                        grad4xzzz_mz(3, 0.0), grad4yyyy_mz(3, 0.0), grad4yyyz_mz(3, 0.0);
    std::vector<double> grad4yyzz_mz(3, 0.0), grad4yzzz_mz(3, 0.0), grad4zzzz_mz(3, 0.0);

    int xc_id = 106;  

    auto torque = gga_torque(xc_id, n, mx, my, mz,
        gradx_n, grady_n, gradz_n,
        gradx_mx, grady_mx, gradz_mx,
        gradx_my, grady_my, gradz_my,
        gradx_mz, grady_mz, gradz_mz,
        grad2xx_n, grad2yy_n, grad2zz_n, grad2xy_n, grad2yz_n, grad2xz_n,
        grad2xx_mx, grad2yy_mx, grad2zz_mx, grad2xy_mx, grad2yz_mx, grad2xz_mx,
        grad2xx_my, grad2yy_my, grad2zz_my, grad2xy_my, grad2yz_my, grad2xz_my,
        grad2xx_mz, grad2yy_mz, grad2zz_mz, grad2xy_mz, grad2yz_mz, grad2xz_mz,
        grad3xxx_n, grad3xxy_n, grad3xxz_n, grad3xyy_n, grad3xyz_n, grad3xzz_n,
        grad3yyy_n, grad3yyz_n, grad3yzz_n, grad3zzz_n,
        grad3xxx_mx, grad3xxy_mx, grad3xxz_mx, grad3xyy_mx, grad3xyz_mx, grad3xzz_mx,
        grad3yyy_mx, grad3yyz_mx, grad3yzz_mx, grad3zzz_mx,
        grad3xxx_my, grad3xxy_my, grad3xxz_my, grad3xyy_my, grad3xyz_my, grad3xzz_my,
        grad3yyy_my, grad3yyz_my, grad3yzz_my, grad3zzz_my,
        grad3xxx_mz, grad3xxy_mz, grad3xxz_mz, grad3xyy_mz, grad3xyz_mz, grad3xzz_mz,
        grad3yyy_mz, grad3yyz_mz, grad3yzz_mz, grad3zzz_mz,
        grad4xxxx_n, grad4xxxy_n, grad4xxxz_n, grad4xxyy_n, grad4xxyz_n, grad4xxzz_n,
        grad4xyyy_n, grad4xyyz_n, grad4xyzz_n, grad4xzzz_n, grad4yyyy_n, grad4yyyz_n,
        grad4yyzz_n, grad4yzzz_n, grad4zzzz_n,
        grad4xxxx_mx, grad4xxxy_mx, grad4xxxz_mx, grad4xxyy_mx, grad4xxyz_mx, grad4xxzz_mx,
        grad4xyyy_mx, grad4xyyz_mx, grad4xyzz_mx, grad4xzzz_mx, grad4yyyy_mx, grad4yyyz_mx,
        grad4yyzz_mx, grad4yzzz_mx, grad4zzzz_mx,
        grad4xxxx_my, grad4xxxy_my, grad4xxxz_my, grad4xxyy_my, grad4xxyz_my, grad4xxzz_my,
        grad4xyyy_my, grad4xyyz_my, grad4xyzz_my, grad4xzzz_my, grad4yyyy_my, grad4yyyz_my,
        grad4yyzz_my, grad4yzzz_my, grad4zzzz_my,
        grad4xxxx_mz, grad4xxxy_mz, grad4xxxz_mz, grad4xxyy_mz, grad4xxyz_mz, grad4xxzz_mz,
        grad4xyyy_mz, grad4xyyz_mz, grad4xyzz_mz, grad4xzzz_mz, grad4yyyy_mz, grad4yyyz_mz,
        grad4yyzz_mz, grad4yzzz_mz, grad4zzzz_mz);

    // print
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "GGA local torque:" << std::endl;
    for (size_t i = 0; i < torque.size()/3; ++i)
    {
        std::cout << "Grid point " << i << ": ("
                  << torque[3*i] << ", " 
                  << torque[3*i+1] << ", " 
                  << torque[3*i+2] << ")" << std::endl;
    }
}