#include "grid_bigcell.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

Grid_BigCell::Grid_BigCell() {
  this->flag_tib = false;
  this->index_atom = nullptr;
  this->orbital_rmax = 0.0;
  this->nxe = this->nye = this->nze = 0;
  this->bigcell_dx = 0.0;
  this->bigcell_dy = 0.0;
  this->bigcell_dz = 0.0;
  this->dxe = 0;
  this->dye = 0;
  this->dze = 0;
  this->nxe = 0;
  this->nye = 0;
  this->nze = 0;
  this->nxyze = 0;
}

Grid_BigCell::~Grid_BigCell() {
  // delete tau positions.
  if (this->flag_tib) {
    for (int i = 0; i < GlobalC::ucell.nat; i++) {
      delete[] tau_in_bigcell[i];
    }
    delete[] tau_in_bigcell;
  }
  delete[] index_atom;
}

void Grid_BigCell::init_big_latvec(void) {
  ModuleBase::TITLE("Grid_BigCell", "init_big_latvec");
  // initialize the mesh cell vectors.
  assert(nbx > 0);
  assert(nby > 0);
  assert(nbz >= 0);

  // size of each big room (same shape with unitcell)
  this->bigcell_vec1[0] =
      GlobalC::ucell.a1.x / (double)nbx * GlobalC::ucell.lat0;
  this->bigcell_vec1[1] =
      GlobalC::ucell.a1.y / (double)nbx * GlobalC::ucell.lat0;
  this->bigcell_vec1[2] =
      GlobalC::ucell.a1.z / (double)nbx * GlobalC::ucell.lat0;

  this->bigcell_vec2[0] =
      GlobalC::ucell.a2.x / (double)nby * GlobalC::ucell.lat0;
  this->bigcell_vec2[1] =
      GlobalC::ucell.a2.y / (double)nby * GlobalC::ucell.lat0;
  this->bigcell_vec2[2] =
      GlobalC::ucell.a2.z / (double)nby * GlobalC::ucell.lat0;

  this->bigcell_vec3[0] =
      GlobalC::ucell.a3.x / (double)nbz * GlobalC::ucell.lat0;
  this->bigcell_vec3[1] =
      GlobalC::ucell.a3.y / (double)nbz * GlobalC::ucell.lat0;
  this->bigcell_vec3[2] =
      GlobalC::ucell.a3.z / (double)nbz * GlobalC::ucell.lat0;

  this->bigcell_latvec0.e11 = this->bigcell_vec1[0];
  this->bigcell_latvec0.e12 = this->bigcell_vec1[1];
  this->bigcell_latvec0.e13 = this->bigcell_vec1[2];

  this->bigcell_latvec0.e21 = this->bigcell_vec2[0];
  this->bigcell_latvec0.e22 = this->bigcell_vec2[1];
  this->bigcell_latvec0.e23 = this->bigcell_vec2[2];

  this->bigcell_latvec0.e31 = this->bigcell_vec3[0];
  this->bigcell_latvec0.e32 = this->bigcell_vec3[1];
  this->bigcell_latvec0.e33 = this->bigcell_vec3[2];

  // why we need GT = bigcell_latvec0^(-1)?
  // note that (i,j,k) is a grid point.
  // (x,y,z) is the cartesian coordinates.
  // because
  // (x,y,z) = (i,j,k) * bigcell_latvec0
  // once we know (x,y,z) and bigcell_latvec0
  // we need to transform the formula to
  // (x,y,z) * bigcell_latvec0^(-1) = (i,j,k)
  this->bigcell_GT = this->bigcell_latvec0.Inverse();

  if (GlobalV::test_gridt) {
    GlobalV::ofs_running << " the VECTORS of BIGCELL are (Bohr): " << std::endl;
    GlobalV::ofs_running << " vec1( " << std::setw(15) << bigcell_vec1[0]
                         << std::setw(15) << bigcell_vec1[1] << std::setw(15)
                         << bigcell_vec1[2] << ")" << std::endl;

    GlobalV::ofs_running << " vec2( " << std::setw(15) << bigcell_vec2[0]
                         << std::setw(15) << bigcell_vec2[1] << std::setw(15)
                         << bigcell_vec2[2] << ")" << std::endl;

    GlobalV::ofs_running << " vec3( " << std::setw(15) << bigcell_vec3[0]
                         << std::setw(15) << bigcell_vec3[1] << std::setw(15)
                         << bigcell_vec3[2] << ")" << std::endl;
  }
  return;
}

void Grid_BigCell::init_grid_expansion(void) {
  ModuleBase::TITLE("Grid_BigCell", "init_grid_expansion");

  // calculate the max cutoff radius among all orbitals.
  // then we will use this parameter to generate grid expansion.
  for (int T = 0; T < GlobalC::ucell.ntype; T++) {
    this->orbital_rmax =
        std::max(GlobalC::ORB.Phi[T].getRcut(), this->orbital_rmax);
  }
  if (GlobalV::test_gridt)
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                "rmax of periodic grid (bohr)", orbital_rmax);

  // mohan fixed serious bug 2010-03-06
  // G = GT^T
  // g1 = the norm of first std::vector of G
  // g2 = the norm of second std::vector of G
  // g3 = the norm of third std::vector of G
  double g1 =
      sqrt(bigcell_GT.e11 * bigcell_GT.e11 + bigcell_GT.e21 * bigcell_GT.e21 +
           bigcell_GT.e31 * bigcell_GT.e31);

  double g2 =
      sqrt(bigcell_GT.e12 * bigcell_GT.e12 + bigcell_GT.e22 * bigcell_GT.e22 +
           bigcell_GT.e32 * bigcell_GT.e32);

  double g3 =
      sqrt(bigcell_GT.e13 * bigcell_GT.e13 + bigcell_GT.e23 * bigcell_GT.e23 +
           bigcell_GT.e33 * bigcell_GT.e33);

  // we assume the added bigcell can present even the atom
  // is at the edge of the origin grid.
  // mohan add +1, 2011-04-23
  this->dxe = static_cast<int>(this->orbital_rmax * g1) + 1;
  this->dye = static_cast<int>(this->orbital_rmax * g2) + 1;
  this->dze = static_cast<int>(this->orbital_rmax * g3) + 1;

  // xiaohui add 'GlobalV::OUT_LEVEL' line, 2015-09-16
  if (GlobalV::OUT_LEVEL != "m")
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "extended fft grid", dxe,
                                dye, dze);

  // calculate the dimension of expanded grid.
  // +1 in order to cover the spillage atom on the right side.
  assert(nbx > 0);
  assert(nby > 0);
  assert(nbz >= 0);

  this->nxe = nbx + 2 * dxe + 1;
  this->nye = nby + 2 * dye + 1;
  this->nze = nbz + 2 * dze + 1;
  this->nxyze = this->nxe * this->nye * this->nze;

  if (GlobalV::OUT_LEVEL != "m")
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                "dimension of extened grid", nxe, nye, nze);
  return;
}

void Grid_BigCell::init_tau_in_bigcell(void) {
  ModuleBase::TITLE("Grid_BigCell", "init_tau_in_bigcell");

  // allcoate space for atom positions relative
  // to meshcell.

  if (!flag_tib) {
    this->tau_in_bigcell = new double *[GlobalC::ucell.nat];
    for (int i = 0; i < GlobalC::ucell.nat; i++) {
      this->tau_in_bigcell[i] = new double[3];
    }
    this->flag_tib = true;

    // allocate space, these arrays record which meshcell
    // the atom is in.
    delete[] index_atom;
    this->index_atom = new int[GlobalC::ucell.nat];

    ModuleBase::Memory::record("tau_in_bigcell",
                               sizeof(double) * GlobalC::ucell.nat * 3);
  }

  // get the fraction number of (i,j,k)
  ModuleBase::Vector3<double> fraction;
  int iat = 0;
  int ii, jj, kk;
  double delta[3];
  for (int it = 0; it < GlobalC::ucell.ntype; it++) {
    for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++) {
      // direct positions of atoms calculated from cartesian coordinates.
      // not used because the factrion may be <0 (although very small, such as
      // -1.0e-15) mohan note 2012-07-03
      // fraction = ( GlobalC::ucell.atoms[it].tau[ia] * GlobalC::ucell.lat0 )*
      // this->bigcell_GT;

      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // mohan add 2012-07-03,
      // this can make sure faction are always larger than 0.
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      fraction.x = GlobalC::ucell.atoms[it].taud[ia].x / (1.0 / (double)nbx);
      fraction.y = GlobalC::ucell.atoms[it].taud[ia].y / (1.0 / (double)nby);
      fraction.z = GlobalC::ucell.atoms[it].taud[ia].z / (1.0 / (double)nbz);

      // never use the following, especially for k-algorithm,
      // it may move the atom to a cell that it doesn't belong
      // to
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // mohan add 2012-06-07
      // fraction may be very very small, about -1.0e-15,
      // and the fraction must > 0, so I use periodic boundary condition
      //			if( fraction.x < 0.0 ) fraction.x += nxe;
      //			if( fraction.y < 0.0 ) fraction.y += nye;
      //			if( fraction.z < 0.0 ) fraction.z += nze;

      if (fraction.x < 0 || fraction.y < 0 || fraction.z < 0) {
        std::cout << " Atom positions " << std::endl;
        std::cout << GlobalC::ucell.atoms[it].tau[ia].x << " ";
        std::cout << GlobalC::ucell.atoms[it].tau[ia].y << " ";
        std::cout << GlobalC::ucell.atoms[it].tau[ia].z << " ";
        std::cout << " fraction " << std::endl;
        std::cout << fraction.x << " ";
        std::cout << fraction.y << " ";
        std::cout << fraction.z << " ";
        std::cout << std::endl;
        ModuleBase::WARNING_QUIT(
            "Grid_BigCell::init_tau_in_bigcell",
            "fraction.x<0 || fraction.y<0 || fraction.z<0");
      }

      assert(fraction.x >= 0.0);
      assert(fraction.y >= 0.0);
      assert(fraction.z >= 0.0);

      // make clean which meshcell the atom is in.
      ii = static_cast<int>(fraction.x + 1.0e-8);
      jj = static_cast<int>(fraction.y + 1.0e-8);
      kk = static_cast<int>(fraction.z + 1.0e-8);

      // calculate the index of each corresponding meshcell.
      // Notice ! In fact, we need to minus ii,jj,kk by 1.
      // to label the atom belong to which meshcell
      // in a usual way: left, down corner.
      // if we dont' do this, means the start position
      // of atom is another tyep: right,up corner.
      // which cause minus atom position in grid integration.

      // index_atom: atom 'iat' index in extended grid.
      this->index_atom[iat] = (kk + dze) + (jj + dye) * this->nze +
                              (ii + dxe) * this->nye * this->nze;

      /*
      if(index_atom[iat]==3483935)
      {
              std::cout << "\n i=" << kk+dze << " j=" << jj+dye << " k=" <<
      ii+dxe; BLOCK_HERE("check index atom");
      }
      */

      // get the relative position in direct coordinate.
      delta[0] = fraction.x - (double)ii;
      delta[1] = fraction.y - (double)jj;
      delta[2] = fraction.z - (double)kk;

      if (std::abs(delta[0]) < 1.0e-8)
        delta[0] = 0.0;
      if (std::abs(delta[1]) < 1.0e-8)
        delta[1] = 0.0;
      if (std::abs(delta[2]) < 1.0e-8)
        delta[2] = 0.0;

      //			std::cout << " fraction=" << fraction.x << " "
      //<< fraction.y << " " << fraction.z << std::endl;
      //std::cout << " delta="
      //<< delta[0] << " " << delta[1] << " " << delta[2] << std::endl;

      // get the true relative cartesian coordinate of each atom to the
      // coresponding meshcell.
      for (int ic = 0; ic < 3; ic++) {
        this->tau_in_bigcell[iat][ic] = delta[0] * this->bigcell_vec1[ic] +
                                        delta[1] * this->bigcell_vec2[ic] +
                                        delta[2] * this->bigcell_vec3[ic];
      }

      ++iat;
    }
  }

  return;
}

// (3)
// if f2normal == true, calculate the index2normal.
// if f2normal == false, calculate the index2cell.
void Grid_BigCell::grid_expansion_index(bool f2normal, int *target) const {
  ModuleBase::TITLE("Grid_BigCell", "grid_expansion_index");
  ModuleBase::timer::tick("Grid_BigCell", "grid_expansion_index");
  //	std::cout << " ncx=" << ncx << " ncy=" << ncy << " ncz=" << ncz <<
  // std::endl; 	std::stringstream ss; 	ss << GlobalV::global_out_dir <<
  //"expand_grid.dat"; 	std::ofstream ofs(ss.str().c_str());

  int ii, jj, kk, in_ext, in_normal;
  for (int i = 0; i < this->nxe; i++) {
    for (int j = 0; j < this->nye; j++) {
      for (int k = 0; k < this->nze; k++) {
        in_ext = k + j * this->nze + i * this->nye * this->nze;

        // range from [-dxe,ncx+dxe]
        ii = i - this->dxe;
        jj = j - this->dye;
        kk = k - this->dze;

        //---------------------------------------------------
        // mohan add 2010-10-28
        // be careful of the box.
        // it's useful only when k points are used in LCAO.
        // for example, we construct a 2D supercell
        // and using 32 * 32 FFT grid (bigcell ) to do
        // grid integration,
        // then the first cell (0,0) along x is [0,31)
        // others are:
        // cell index: (-2,0)   , (-1,0)  , (0,0),  (0,1)
        // fft index: [-64,-33], [-32,-1], [0,31], [32,63].
        // look at the formulas below,
        // at first, we take grid_index2ucell1=(ii/nbx)
        // but then we found it is wrong if ii < 0.
        // for example, if ii is -31, the box is -1,
        // so we add -1, the formula turns to ii/nbx-1,
        // but if ii is -32, the box is -1-1 = -2, not correct.
        // so we add 1 to ii, the box will be -31/32-1=-1, correct!
        // the formula is (ii+1)/nbx-1,
        // if ii is -1, the box is still -1, correct!
        // if ii is -33, the box is -2, correct!
        //---------------------------------------------------

        int cel1, cel2, cel3;

        if (ii < 0)
          cel1 = (ii + 1) / nbx - 1;
        else
          cel1 = ii / nbx;
        if (jj < 0)
          cel2 = (jj + 1) / nby - 1;
        else
          cel2 = jj / nby;
        if (kk < 0)
          cel3 = (kk + 1) / nbz - 1;
        else
          cel3 = kk / nbz;

        if (!f2normal) {
          // target: index2ucell
          target[in_ext] = this->cal_Rindex(cel1, cel2, cel3);
        } else {
          // std::cout << "i=" << ii << " j=" << jj << " k=" << kk << std::endl;
          // std::cout << "box1=" << box1 << " box2=" << box2 << " box3=" <<
          // box3 << std::endl;

          // if ii < 0, we need to make ii > 0.
          // so we add 10000 layers. It should be enough.
          // ii, jj, kk shoudl -- ?????????????
          ii = (ii + 10000 * nbx) % nbx;
          jj = (jj + 10000 * nby) % nby;
          kk = (kk + 10000 * nbz) % nbz;

          // std::cout << "ii=" << ii << " jj=" << jj << " kk=" << kk <<
          // std::endl; int ok; cin >> ok;

          assert(ii >= 0);
          assert(jj >= 0);
          assert(kk >= 0);

          assert(in_ext < nxyze);

          if (ii < nbx && jj < nby && kk < nbz) {
            in_normal = kk + jj * nbz + ii * nby * nbz;

            // target: index2normal
            target[in_ext] = in_normal;

            /*
               if(in_ext == 5816805)
               {
               std::cout << "\n index2normal[5816805]=" << in_normal <<
               std::endl; BLOCK_HERE("check index2normal");
               }
             */
          } else {
            ModuleBase::WARNING_QUIT("Grid_BigCell::init_grid_expansion_index",
                                     "check ii,jj,kk!");
          }
        } // f2 normal
      }   // k
    }     // j
  }       // i
  ModuleBase::timer::tick("Grid_BigCell", "grid_expansion_index");
  return;
}
