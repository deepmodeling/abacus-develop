#ifndef RHOG_IO_H
#define RHOG_IO_H

#include <string>
#include <cassert>
#include "module_basis/module_pw/pw_basis.h"
/**
 * I/O free function of rho(G) in binary format
 * Author: YuLiu98, Kirk0830
 * 
 * The designed foramt of the binary file is kept similar to the one in QuantumESPRESSO
 * Modules/io_base.f90, function write_rhog (non-HDF5):
 * 
 * /3/ /gammaonly/ /ngm_g/ /nspin/ /3/ ! bool, int, int. ngm_g is the global, total number of G-vectors
 * /9/ /b11/ /b12/ /b13/ /b21/ /b22/ /b23/ /b31/ /b32/ /b33/ /9/ ! 9 real numbers, the three lattice vectors
 * /ngm_g*3/
 * /miller(1,1)/ /miller(1,2)/ /miller(1,3)/
 * /miller(2,1)/ /miller(2,2)/ /miller(2,3)/
 * ...
 * /miller(ngm_g,1)/ /miller(ngm_g,2)/ /miller(ngm_g,3)/
 * /ngm_g*3/ ! ngm_g*3 integers, the G-vectors in Miller indices
 * /ngm_g/
 * /rhog(1)/ /rhog(2)/ /rhog(3)/ ...
 * /ngm_g/ ! ngm_g complex numbers, the rho(G) values
 * !/ngm_g/ ! ngm_g complex numbers, the rho(G) values for the second spin component
 * !/rhog(1)/ /rhog(2)/ /rhog(3)/ ...
 * !/ngm_g/ ! ngm_g complex numbers, the rho(G) values for the second spin component
 * 
 * There are some aspects needed to ensure:
 * 1. the correspondence between ABACUS PW and QuantumESPRESSO
 *                      ABACUS PW      QuantumESPRESSO
 *    bij              UnitCell::GT      b1, b2, b3
 *    miller index        ibox[*]           mill
 *    rho              Charge::rhog          rho
 * 
 * 2. the unit of each quantity
 *                      ABACUS PW      QuantumESPRESSO
 *    bij                                   a.u.
 *    miller index
 *    rho   
 */

namespace ModuleIO
{

bool read_rhog(const std::string& filename, const ModulePW::PW_Basis* pw_rhod, std::complex<double>** rhog);

bool write_rhog(const std::string& fchg,
                const bool gamma_only,            // from INPUT
                const ModulePW::PW_Basis* pw_rho, // pw_rho in runtime
                const int nspin,                  // GlobalV
                const ModuleBase::Matrix3& GT,    // from UnitCell, useful for calculating the miller
                std::complex<double>** rhog,
                const int ipool,
                const int irank,
                const int nrank);

} // namespace ModuleIO

#endif
