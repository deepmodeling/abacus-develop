// programmed by Xiaoyu Zhang, Peking University, Beijing, China 2025/02/08
// calculate xc local torque = $\mathbf{m}\times\mathbf{B}_{xc}$ and Bxc $\mathbf{B}_{xc}=\frac{\delta E_{xc}[n,\mathbf{m}]}{\delta \mathbf{m}}$
// based on the multi-collinear approach

#include "NCLibxc.h"
#include <vector>
#include <complex>

namespace NCXC{

// gga，input: xc_id, n, mx, my, mz，grad, grad2, grad3, grad4, return Bxc for each real space grid point
 std::vector<double> NCLibxc::gga_Bxc(int xc_id, const std::vector<double>& n, 
    const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz,const std::vector<double> gradx_n,const std::vector<double> grady_n,const std::vector<double> gradz_n, const std::vector<double>& gradx_mx, const std::vector<double>& grady_mx, const std::vector<double>& gradz_mx, const std::vector<double>& gradx_my, const std::vector<double>& grady_my, const std::vector<double>& gradz_my, const std::vector<double>& gradx_mz, const std::vector<double>& grady_mz, const std::vector<double>& gradz_mz, const std::vector<double> grad2xx_n, const std::vector<double> grad2yy_n, const std::vector<double> grad2zz_n, const std::vector<double> grad2xy_n, const std::vector<double> grad2yz_n, const std::vector<double> grad2xz_n, const std::vector<double> grad2xx_mx, const std::vector<double> grad2yy_mx, const std::vector<double> grad2zz_mx, const std::vector<double> grad2xy_mx, const std::vector<double> grad2yz_mx, const std::vector<double> grad2xz_mx, const std::vector<double> grad2xx_my, const std::vector<double> grad2yy_my, const std::vector<double> grad2zz_my, const std::vector<double> grad2xy_my, const std::vector<double> grad2yz_my, const std::vector<double> grad2xz_my, const std::vector<double> grad2xx_mz, const std::vector<double> grad2yy_mz, const std::vector<double> grad2zz_mz, const std::vector<double> grad2xy_mz, const std::vector<double> grad2yz_mz, const std::vector<double> grad2xz_mz, const std::vector<double> grad3xxx_n, const std::vector<double> grad3xxy_n, const std::vector<double> grad3xxz_n, const std::vector<double> grad3xyy_n, const std::vector<double> grad3xyz_n, const std::vector<double> grad3xzz_n, const std::vector<double> grad3yyy_n, const std::vector<double> grad3yyz_n, const std::vector<double> grad3yzz_n, const std::vector<double> grad3zzz_n, const std::vector<double> grad3xxx_mx, const std::vector<double> grad3xxy_mx, const std::vector<double> grad3xxz_mx, const std::vector<double> grad3xyy_mx, const std::vector<double> grad3xyz_mx, const std::vector<double> grad3xzz_mx, const std::vector<double> grad3yyy_mx, const std::vector<double> grad3yyz_mx, const std::vector<double> grad3yzz_mx, const std::vector<double> grad3zzz_mx, const std::vector<double> grad3xxx_my, const std::vector<double> grad3xxy_my, const std::vector<double> grad3xxz_my, const std::vector<double> grad3xyy_my, const std::vector<double> grad3xyz_my, const std::vector<double> grad3xzz_my, const std::vector<double> grad3yyy_my, const std::vector<double> grad3yyz_my, const std::vector<double> grad3yzz_my, const std::vector<double> grad3zzz_my, const std::vector<double> grad3xxx_mz, const std::vector<double> grad3xxy_mz, const std::vector<double> grad3xxz_mz, const std::vector<double> grad3xyy_mz, const std::vector<double> grad3xyz_mz, const std::vector<double> grad3xzz_mz, const std::vector<double> grad3yyy_mz, const std::vector<double> grad3yyz_mz, const std::vector<double> grad3yzz_mz, const std::vector<double> grad3zzz_mz, const std::vector<double> grad4xxxx_n, const std::vector<double> grad4xxxy_n, const std::vector<double> grad4xxxz_n, const std::vector<double> grad4xxyy_n, const std::vector<double> grad4xxyz_n, const std::vector<double> grad4xxzz_n, const std::vector<double> grad4xyyy_n, const std::vector<double> grad4xyyz_n, const std::vector<double> grad4xyzz_n, const std::vector<double> grad4xzzz_n, const std::vector<double> grad4yyyy_n, const std::vector<double> grad4yyyz_n, const std::vector<double> grad4yyzz_n, const std::vector<double> grad4yzzz_n, const std::vector<double> grad4zzzz_n, const std::vector<double> grad4xxxx_mx, const std::vector<double> grad4xxxy_mx, const std::vector<double> grad4xxxz_mx, const std::vector<double> grad4xxyy_mx, const std::vector<double> grad4xxyz_mx, const std::vector<double> grad4xxzz_mx, const std::vector<double> grad4xyyy_mx, const std::vector<double> grad4xyyz_mx, const std::vector<double> grad4xyzz_mx, const std::vector<double> grad4xzzz_mx, const std::vector<double> grad4yyyy_mx, const std::vector<double> grad4yyyz_mx, const std::vector<double> grad4yyzz_mx, const std::vector<double> grad4yzzz_mx, const std::vector<double> grad4zzzz_mx, const std::vector<double> grad4xxxx_my, const std::vector<double> grad4xxxy_my, const std::vector<double> grad4xxxz_my, const std::vector<double> grad4xxyy_my, const std::vector<double> grad4xxyz_my, const std::vector<double> grad4xxzz_my, const std::vector<double> grad4xyyy_my, const std::vector<double> grad4xyyz_my, const std::vector<double> grad4xyzz_my, const std::vector<double> grad4xzzz_my, const std::vector<double> grad4yyyy_my, const std::vector<double> grad4yyyz_my, const std::vector<double> grad4yyzz_my, const std::vector<double> grad4yzzz_my, const std::vector<double> grad4zzzz_my, const std::vector<double> grad4xxxx_mz, const std::vector<double> grad4xxxy_mz, const std::vector<double> grad4xxxz_mz, const std::vector<double> grad4xxyy_mz, const std::vector<double> grad4xxyz_mz, const std::vector<double> grad4xxzz_mz, const std::vector<double> grad4xyyy_mz, const std::vector<double> grad4xyyz_mz, const std::vector<double> grad4xyzz_mz, const std::vector<double> grad4xzzz_mz, const std::vector<double> grad4yyyy_mz, const std::vector<double> grad4yyyz_mz, const std::vector<double> grad4yyzz_mz, const std::vector<double> grad4yzzz_mz, const std::vector<double> grad4zzzz_mz){

        std::complex<double> twoi(0.0, 2.0);
        std::complex<double> two(2.0, 0.0);

        auto result = gga_mc(xc_id, n, mx, my, mz, 
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
    

        auto& E = result.first;
        auto& V = result.second;
    
        size_t num_points = n.size();
        // {Bxc_x, Bxc_y, Bxc_z}
        std::vector<double> Bxc(3 * num_points, 0.0);
        for (size_t i = 0; i < num_points; ++i) {
            Bxc[3*i]   = std::real((V[i][0][1]+V[i][1][0])/two);         // Bxc_x
            Bxc[3*i+1] = std::real((V[i][1][0]-V[i][0][1])/twoi);          // Bxc_y
            Bxc[3*i+2] = std::real((V[i][0][0] - V[i][1][1]) / two); // Bxc_z
        }
    
        return Bxc;
    }

// gga，input: xc_id, n, mx, my, mz，grad, grad2, grad3, grad4, return xc local torque for each real space grid point
// torque = - m × Bxc， note that there is a minus 
std::vector<double> NCLibxc::gga_torque(int xc_id, const std::vector<double>& n, 
    const std::vector<double>& mx, const std::vector<double>& my, const std::vector<double>& mz,const std::vector<double> gradx_n,const std::vector<double> grady_n,const std::vector<double> gradz_n, const std::vector<double>& gradx_mx, const std::vector<double>& grady_mx, const std::vector<double>& gradz_mx, const std::vector<double>& gradx_my, const std::vector<double>& grady_my, const std::vector<double>& gradz_my, const std::vector<double>& gradx_mz, const std::vector<double>& grady_mz, const std::vector<double>& gradz_mz, const std::vector<double> grad2xx_n, const std::vector<double> grad2yy_n, const std::vector<double> grad2zz_n, const std::vector<double> grad2xy_n, const std::vector<double> grad2yz_n, const std::vector<double> grad2xz_n, const std::vector<double> grad2xx_mx, const std::vector<double> grad2yy_mx, const std::vector<double> grad2zz_mx, const std::vector<double> grad2xy_mx, const std::vector<double> grad2yz_mx, const std::vector<double> grad2xz_mx, const std::vector<double> grad2xx_my, const std::vector<double> grad2yy_my, const std::vector<double> grad2zz_my, const std::vector<double> grad2xy_my, const std::vector<double> grad2yz_my, const std::vector<double> grad2xz_my, const std::vector<double> grad2xx_mz, const std::vector<double> grad2yy_mz, const std::vector<double> grad2zz_mz, const std::vector<double> grad2xy_mz, const std::vector<double> grad2yz_mz, const std::vector<double> grad2xz_mz, const std::vector<double> grad3xxx_n, const std::vector<double> grad3xxy_n, const std::vector<double> grad3xxz_n, const std::vector<double> grad3xyy_n, const std::vector<double> grad3xyz_n, const std::vector<double> grad3xzz_n, const std::vector<double> grad3yyy_n, const std::vector<double> grad3yyz_n, const std::vector<double> grad3yzz_n, const std::vector<double> grad3zzz_n, const std::vector<double> grad3xxx_mx, const std::vector<double> grad3xxy_mx, const std::vector<double> grad3xxz_mx, const std::vector<double> grad3xyy_mx, const std::vector<double> grad3xyz_mx, const std::vector<double> grad3xzz_mx, const std::vector<double> grad3yyy_mx, const std::vector<double> grad3yyz_mx, const std::vector<double> grad3yzz_mx, const std::vector<double> grad3zzz_mx, const std::vector<double> grad3xxx_my, const std::vector<double> grad3xxy_my, const std::vector<double> grad3xxz_my, const std::vector<double> grad3xyy_my, const std::vector<double> grad3xyz_my, const std::vector<double> grad3xzz_my, const std::vector<double> grad3yyy_my, const std::vector<double> grad3yyz_my, const std::vector<double> grad3yzz_my, const std::vector<double> grad3zzz_my, const std::vector<double> grad3xxx_mz, const std::vector<double> grad3xxy_mz, const std::vector<double> grad3xxz_mz, const std::vector<double> grad3xyy_mz, const std::vector<double> grad3xyz_mz, const std::vector<double> grad3xzz_mz, const std::vector<double> grad3yyy_mz, const std::vector<double> grad3yyz_mz, const std::vector<double> grad3yzz_mz, const std::vector<double> grad3zzz_mz, const std::vector<double> grad4xxxx_n, const std::vector<double> grad4xxxy_n, const std::vector<double> grad4xxxz_n, const std::vector<double> grad4xxyy_n, const std::vector<double> grad4xxyz_n, const std::vector<double> grad4xxzz_n, const std::vector<double> grad4xyyy_n, const std::vector<double> grad4xyyz_n, const std::vector<double> grad4xyzz_n, const std::vector<double> grad4xzzz_n, const std::vector<double> grad4yyyy_n, const std::vector<double> grad4yyyz_n, const std::vector<double> grad4yyzz_n, const std::vector<double> grad4yzzz_n, const std::vector<double> grad4zzzz_n, const std::vector<double> grad4xxxx_mx, const std::vector<double> grad4xxxy_mx, const std::vector<double> grad4xxxz_mx, const std::vector<double> grad4xxyy_mx, const std::vector<double> grad4xxyz_mx, const std::vector<double> grad4xxzz_mx, const std::vector<double> grad4xyyy_mx, const std::vector<double> grad4xyyz_mx, const std::vector<double> grad4xyzz_mx, const std::vector<double> grad4xzzz_mx, const std::vector<double> grad4yyyy_mx, const std::vector<double> grad4yyyz_mx, const std::vector<double> grad4yyzz_mx, const std::vector<double> grad4yzzz_mx, const std::vector<double> grad4zzzz_mx, const std::vector<double> grad4xxxx_my, const std::vector<double> grad4xxxy_my, const std::vector<double> grad4xxxz_my, const std::vector<double> grad4xxyy_my, const std::vector<double> grad4xxyz_my, const std::vector<double> grad4xxzz_my, const std::vector<double> grad4xyyy_my, const std::vector<double> grad4xyyz_my, const std::vector<double> grad4xyzz_my, const std::vector<double> grad4xzzz_my, const std::vector<double> grad4yyyy_my, const std::vector<double> grad4yyyz_my, const std::vector<double> grad4yyzz_my, const std::vector<double> grad4yzzz_my, const std::vector<double> grad4zzzz_my, const std::vector<double> grad4xxxx_mz, const std::vector<double> grad4xxxy_mz, const std::vector<double> grad4xxxz_mz, const std::vector<double> grad4xxyy_mz, const std::vector<double> grad4xxyz_mz, const std::vector<double> grad4xxzz_mz, const std::vector<double> grad4xyyy_mz, const std::vector<double> grad4xyyz_mz, const std::vector<double> grad4xyzz_mz, const std::vector<double> grad4xzzz_mz, const std::vector<double> grad4yyyy_mz, const std::vector<double> grad4yyyz_mz, const std::vector<double> grad4yyzz_mz, const std::vector<double> grad4yzzz_mz, const std::vector<double> grad4zzzz_mz){
        auto Bxc = gga_Bxc(xc_id, n, mx, my, mz,
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
    
        // torque = m × Bxc
        size_t num_points = n.size();
        std::vector<double> torque(3 * num_points, 0.0);
        for (size_t i = 0; i < num_points; ++i) {
            double Bx = Bxc[3*i];
            double By = Bxc[3*i+1];
            double Bz = Bxc[3*i+2];
      
            // (m×Bxc)_x = my_i*Bz - mz_i*By
            torque[3*i]   = -my[i]*Bz - mz[i]*By;
            // (m×Bxc)_y = mz_i*Bx - mx_i*Bz
            torque[3*i+1] = -mz[i]*Bx - mx[i]*Bz;
            // (m×Bxc)_z = mx_i*By - my_i*Bx
            torque[3*i+2] = -mx[i]*By - my[i]*Bx;
        }
    
        return torque;
    }


}// end of namespace NCXC