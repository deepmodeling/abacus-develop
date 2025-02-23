// programmed by Xiaoyu Zhang, Peking University, Beijing, China 2025/02/08
// calculate xc local torque = $\mathbf{m}\times\mathbf{B}_{xc}$ and Bxc $\mathbf{B}_{xc}=\frac{\delta E_{xc}[n,\mathbf{m}]}{\delta \mathbf{m}}$
// based on the multi-collinear approach

#include "NCLibxc.h"
#include <vector>
#include <complex>

namespace NCXC{

// gga，input: V, return Bxc for each real space grid point
 std::vector<double> NCLibxc::gga_Bxc(const std::vector<Matrix2x2> &V){

        std::complex<double> twoi(0.0, 2.0);
        std::complex<double> two(2.0, 0.0);

       
        size_t num_points = V.size();
        // {Bxc_x, Bxc_y, Bxc_z}
        std::vector<double> Bxc(3 * num_points, 0.0);
        for (size_t i = 0; i < num_points; ++i) {
            Bxc[3*i]   = std::real((V[i][0][1]+V[i][1][0])/two);         // Bxc_x
            Bxc[3*i+1] = std::real((V[i][1][0]-V[i][0][1])/twoi);          // Bxc_y
            Bxc[3*i+2] = std::real((V[i][0][0] - V[i][1][1]) / two); // Bxc_z
        }
    
        return Bxc;
}

// gga，input: m and Bxc, return xc local torque for each real space grid point
// torque = - m × Bxc， note that there is a minus 
std::vector<double> NCLibxc::gga_torque( 
    const std::vector<double>& mx, 
    const std::vector<double>& my, 
    const std::vector<double>& mz, 
    const std::vector<Matrix2x2>& V)
{

    std::vector<double> Bxc = gga_Bxc(V);
    size_t num_points = Bxc.size() / 3;
    std::vector<double> torque(3 * num_points, 0.0);
    for (size_t i = 0; i < num_points; ++i) {
        double Bx = Bxc[3*i];
        double By = Bxc[3*i+1];
        double Bz = Bxc[3*i+2];
  
        // (m×Bxc)_x = my * Bz - mz * By
        // (m×Bxc)_y = mz * Bx - mx * Bz
        // (m×Bxc)_z = mx * By - my * Bx
        torque[3*i]   = - ( my[i] * Bz - mz[i] * By );
        torque[3*i+1] = - ( mz[i] * Bx - mx[i] * Bz );
        torque[3*i+2] = - ( mx[i] * By - my[i] * Bx );
    }
    return torque;
}


}// end of namespace NCXC