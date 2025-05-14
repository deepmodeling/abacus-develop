#include "grid_meshk.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"

Grid_MeshK::Grid_MeshK()
{
}

Grid_MeshK::~Grid_MeshK()
{
}

int Grid_MeshK::cal_Rindex(const int &u1, const int &u2, const int &u3)const
{
	const int x1 = u1 - this->minu1;
	const int x2 = u2 - this->minu2;
	const int x3 = u3 - this->minu3;
	
	if(x1<0 || x2<0 || x3<0)
	{
		std::cout << " u1=" << u1 << " minu1=" << minu1 << std::endl;
		std::cout << " u2=" << u2 << " minu2=" << minu2 << std::endl;
		std::cout << " u3=" << u3 << " minu3=" << minu3 << std::endl;
		ModuleBase::WARNING_QUIT("Grid_MeshK::cal_Rindex","x1<0 || x2<0 || x3<0 !");
	}

	assert(x1>=0);
	assert(x2>=0);
	assert(x3>=0);

	return (x3 + x2 * this->nu3 + x1 * this->nu2 * this->nu3);
}

ModuleBase::Vector3<int> Grid_MeshK::get_ucell_coords(const int &Rindex)const
{
	const int x = ucell_index2x[Rindex];
	const int y = ucell_index2y[Rindex];
	const int z = ucell_index2z[Rindex];

	return ModuleBase::Vector3<int>(x, y, z);
}

void Grid_MeshK::cal_extended_cell(const int &dxe, const int &dye, const int &dze,const int& nbx, const int& nby, const int& nbz)
{
	ModuleBase::TITLE("Grid_MeshK","cal_extended_cell");

	//--------------------------------------
	// max and min unitcell in expaned grid.
	//--------------------------------------
	this->maxu1 = dxe / nbx + 1;
	this->maxu2 = dye / nby + 1;
	this->maxu3 = dze / nbz + 1;

	this->minu1 = (-dxe+1) / nbx - 1; 
	this->minu2 = (-dye+1) / nby - 1; 
	this->minu3 = (-dze+1) / nbz - 1; 

	if(PARAM.inp.test_gridt) {ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MaxUnitcell",maxu1,maxu2,maxu3);
}
	if(PARAM.inp.test_gridt) {ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"MinUnitcell",minu1,minu2,minu3);
}

	//--------------------------------------
	// number of unitcell in each direction.
	//--------------------------------------
	this->nu1 = maxu1 - minu1 + 1;
	this->nu2 = maxu2 - minu2 + 1;
	this->nu3 = maxu3 - minu3 + 1;
	this->nutot = nu1 * nu2 * nu3;

	if(PARAM.inp.test_gridt) {ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"UnitCellNumber",nu1,nu2,nu3);
}
	if(PARAM.inp.out_level != "m") { ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"UnitCellTotal",nutot);
}


    this->ucell_index2x = std::vector<int>(nutot, 0);
    this->ucell_index2y = std::vector<int>(nutot, 0);
    this->ucell_index2z = std::vector<int>(nutot, 0);

	this->nutot = nu1 * nu2 * nu3;

	for(int i=minu1; i<=maxu1; i++)
	{
		for(int j=minu2; j<=maxu2; j++)
		{
			for(int k=minu3; k<=maxu3; k++)
			{
				const int cell = cal_Rindex(i,j,k);	
				assert(cell<nutot);

				this->ucell_index2x[cell] = i;
				this->ucell_index2y[cell] = j;
				this->ucell_index2z[cell] = k;

			}
		}
	}

	return;
}