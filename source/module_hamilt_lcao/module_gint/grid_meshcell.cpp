#include "grid_meshcell.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

Grid_MeshCell::Grid_MeshCell() {
  allocate_pos = false;
  bx = by = bz = bxyz = 1;
}

Grid_MeshCell::~Grid_MeshCell() {
  if (allocate_pos) {
    for (int ib = 0; ib < bxyz; ib++) {
      delete[] meshcell_pos[ib];
    }
    delete[] meshcell_pos;
  }
}

void Grid_MeshCell::set_grid_dim(const int &ncx_in, const int &ncy_in,
                                 const int &ncz_in, const int &bx_in,
                                 const int &by_in, const int &bz_in,
                                 const int &nbx_in, const int &nby_in,
                                 const int &nbz_in, const int &nbxx_in,
                                 const int &nbzp_start_in, const int &nbzp_in) {
  this->ncx = ncx_in;
  this->ncy = ncy_in;
  this->ncz = ncz_in;
  this->ncxyz = ncx * ncy * ncz;
  this->bx = bx_in;
  this->by = by_in;
  this->bz = bz_in;
  this->bxyz = bx * by * bz;
  this->nbx = nbx_in;
  this->nby = nby_in;
  this->nbz = nbz_in;
  this->nbxyz = nbx * nby * nbz;
  this->nbxx = nbxx_in;
  this->nbzp_start = nbzp_start_in;
  this->nbzp = nbzp_in;

  // xiaohui add 'GlobalV::OUT_LEVEL' line, 2015-09-16
  if (GlobalV::OUT_LEVEL != "m")
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "real space grid", ncx,
                                ncy, ncz); // real space uniform grid
  if (GlobalV::OUT_LEVEL != "m")
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                "big cell numbers in grid", nbx, nby,
                                nbz); // reduced by BIG_CELL
  if (GlobalV::OUT_LEVEL != "m")
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,
                                "meshcell numbers in big cell", bx, by,
                                bz); // is small integer, typical number 2*2*2

  // std::cout << " bx=" << bx << " by=" << by << " bz=" << bz << std::endl;
  // std::cout << " nbx=" << nbx << " nby=" << nby << " nbz=" << nbz <<
  // std::endl; std::cout << " ncx=" << ncx << " ncy=" << ncy << " ncz=" << ncz
  // << std::endl;

  return;
}

// (1)
void Grid_MeshCell::init_latvec(void) {
  ModuleBase::TITLE("Grid_MeshCell", "init_latvec");
  // initialize the mesh cell vectors.
  assert(ncx > 0);
  assert(ncy > 0);
  assert(ncz > 0);

  // size of each room (same shape with unitcell)
  this->meshcell_vec1[0] =
      GlobalC::ucell.a1.x / (double)ncx * GlobalC::ucell.lat0;
  this->meshcell_vec1[1] =
      GlobalC::ucell.a1.y / (double)ncx * GlobalC::ucell.lat0;
  this->meshcell_vec1[2] =
      GlobalC::ucell.a1.z / (double)ncx * GlobalC::ucell.lat0;

  this->meshcell_vec2[0] =
      GlobalC::ucell.a2.x / (double)ncy * GlobalC::ucell.lat0;
  this->meshcell_vec2[1] =
      GlobalC::ucell.a2.y / (double)ncy * GlobalC::ucell.lat0;
  this->meshcell_vec2[2] =
      GlobalC::ucell.a2.z / (double)ncy * GlobalC::ucell.lat0;

  this->meshcell_vec3[0] =
      GlobalC::ucell.a3.x / (double)ncz * GlobalC::ucell.lat0;
  this->meshcell_vec3[1] =
      GlobalC::ucell.a3.y / (double)ncz * GlobalC::ucell.lat0;
  this->meshcell_vec3[2] =
      GlobalC::ucell.a3.z / (double)ncz * GlobalC::ucell.lat0;

  this->meshcell_latvec0.e11 = this->meshcell_vec1[0];
  this->meshcell_latvec0.e12 = this->meshcell_vec1[1];
  this->meshcell_latvec0.e13 = this->meshcell_vec1[2];

  this->meshcell_latvec0.e21 = this->meshcell_vec2[0];
  this->meshcell_latvec0.e22 = this->meshcell_vec2[1];
  this->meshcell_latvec0.e23 = this->meshcell_vec2[2];

  this->meshcell_latvec0.e31 = this->meshcell_vec3[0];
  this->meshcell_latvec0.e32 = this->meshcell_vec3[1];
  this->meshcell_latvec0.e33 = this->meshcell_vec3[2];

  // why we need GT = meshcell_latvec0^(-1)?
  // note that (i,j,k) is a grid point.
  // (x,y,z) is the cartesian coordinates.
  // because
  // (x,y,z) = (i,j,k) * meshcell_latvec0
  // once we know (x,y,z) and meshcell_latvec0
  // we need to transform the formula to
  // (x,y,z) * meshcell_latvec0^(-1) = (i,j,k)
  this->meshcell_GT = this->meshcell_latvec0.Inverse();

  if (GlobalV::test_gridt) {
    GlobalV::ofs_running << " the VECTORS of MESHCELL are (Bohr): "
                         << std::endl;
    GlobalV::ofs_running << " vec1( " << std::setw(15) << meshcell_vec1[0]
                         << std::setw(15) << meshcell_vec1[1] << std::setw(15)
                         << meshcell_vec1[2] << ")" << std::endl;

    GlobalV::ofs_running << " vec2( " << std::setw(15) << meshcell_vec2[0]
                         << std::setw(15) << meshcell_vec2[1] << std::setw(15)
                         << meshcell_vec2[2] << ")" << std::endl;

    GlobalV::ofs_running << " vec3( " << std::setw(15) << meshcell_vec3[0]
                         << std::setw(15) << meshcell_vec3[1] << std::setw(15)
                         << meshcell_vec3[2] << ")" << std::endl;
  }

  return;
}

void Grid_MeshCell::init_meshcell_pos(void) {
  assert(bx > 0);
  assert(by > 0);
  assert(bz > 0);
  assert(bxyz > 0);

  if (!allocate_pos) {
    meshcell_pos = new double *[bxyz];
    for (int ib = 0; ib < bxyz; ib++) {
      meshcell_pos[ib] = new double[3];
    }
  }

  int index = 0;
  for (int i = 0; i < bx; i++) {
    for (int j = 0; j < by; j++) {
      for (int k = 0; k < bz; k++) {
        //				std::cout << std::setw(5) << i <<
        //std::setw(5)
        //<< j << std::setw(5) << k;
        for (int p = 0; p < 3; p++) {
          meshcell_pos[index][p] = i * meshcell_vec1[p] + j * meshcell_vec2[p] +
                                   k * meshcell_vec3[p];
          // std::cout << std::setw(15) << meshcell_pos[index][p];
        }
        //				std::cout << std::endl;
        ++index;
      }
    }
  }
  allocate_pos = true;

  return;
}

void Grid_MeshCell::cal_extended_cell(const int &dxe, const int &dye,
                                      const int &dze) {
  ModuleBase::TITLE("Grid_MeshK", "cal_extended_cell");

  //--------------------------------------
  // max and min unitcell in expaned grid.
  //--------------------------------------
  this->maxu1 = dxe / this->nbx + 1;
  this->maxu2 = dye / this->nby + 1;
  this->maxu3 = dze / this->nbz + 1;

  this->minu1 = (-dxe + 1) / this->nbx - 1;
  this->minu2 = (-dye + 1) / this->nby - 1;
  this->minu3 = (-dze + 1) / this->nbz - 1;

  if (GlobalV::test_gridt)
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "MaxUnitcell", maxu1,
                                maxu2, maxu3);
  if (GlobalV::test_gridt)
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "MinUnitcell", minu1,
                                minu2, minu3);

  //--------------------------------------
  // number of unitcell in each direction.
  //--------------------------------------
  this->nu1 = maxu1 - minu1 + 1;
  this->nu2 = maxu2 - minu2 + 1;
  this->nu3 = maxu3 - minu3 + 1;
  this->nutot = nu1 * nu2 * nu3;

  if (GlobalV::test_gridt)
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "UnitCellNumber", nu1,
                                nu2, nu3);
  // xiaohui add 'GlobalV::OUT_LEVEL' line, 2015-09-16
  if (GlobalV::OUT_LEVEL != "m")
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "UnitCellTotal", nutot);

  //	std::cout << " nu1 = " << nu1 << " nu2 = " << nu2 << " nu3 = " << nu3 <<
  // std::endl; 	std::cout << " nutot = " << nutot << std::endl;

  delete[] ucell_index2x;
  delete[] ucell_index2y;
  delete[] ucell_index2z;
  this->ucell_index2x = new int[nutot];
  this->ucell_index2y = new int[nutot];
  this->ucell_index2z = new int[nutot];
  ModuleBase::GlobalFunc::ZEROS(ucell_index2x, nutot);
  ModuleBase::GlobalFunc::ZEROS(ucell_index2y, nutot);
  ModuleBase::GlobalFunc::ZEROS(ucell_index2z, nutot);

  this->nutot = nu1 * nu2 * nu3;

  for (int i = minu1; i <= maxu1; i++) {
    for (int j = minu2; j <= maxu2; j++) {
      for (int k = minu3; k <= maxu3; k++) {
        const int cell = cal_Rindex(i, j, k);
        assert(cell < nutot);

        this->ucell_index2x[cell] = i - minu1;
        this->ucell_index2y[cell] = j - minu2;
        this->ucell_index2z[cell] = k - minu3;
      }
    }
  }

  return;
}