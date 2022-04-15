#include "./pw_basis.h"
#include "../module_base/constants.h"
#include "../module_base/timer.h"

namespace ModulePW
{
//
//Init the grids for FFT
//Input: lattice vectors of the cell, Energy cut off for G^2/2
//Output: nx, ny, nz, nxyz, latvec, G, GT, GGT
//
void PW_Basis:: initgrids(
        double lat0_in, //unit length (unit in bohr)
        ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
        double gridecut
)
{
    // ModuleBase::timer::start();
    //init latice
    this->lat0 = lat0_in;
    this->latvec = latvec_in;
    this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
    

    //------------------------------------------------------------
    //-------------------------init grids-------------------------
    //------------------------------------------------------------
    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / this->lat0 / this->lat0;
    gridecut = gridecut / tpiba2;
    ModuleBase::Vector3<double> lat;
    int *ibox = new int[3];// ibox[i] are the minimal FFT dimensions,
    
    lat.x = latvec.e11;
    lat.y = latvec.e12;
    lat.z = latvec.e13;
    ibox[0] = 2 * int(sqrt(gridecut) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e21;
    lat.y = latvec.e22;
    lat.z = latvec.e23;
    ibox[1] = 2 * int(sqrt(gridecut) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e31;
    lat.y = latvec.e32;
    lat.z = latvec.e33;
    ibox[2] = 2 * int(sqrt(gridecut) * sqrt(lat * lat)) + 1;
    //lat*lat=lat.x*lat.x+lat.y*lat.y+lat.z+lat.z

    // Find the minimal FFT box size the factors into the primes (2,3,5,7).
    for (int i = 0; i < 3; i++)
    {
    	int b = 0;
        int n2 = 0;
        int n3 = 0;
        int n5 = 0;
        //int n7 = 0;
        bool done_factoring = false;
	
		int ns = 0;
        // increase ibox[i] by 1 until it is totally factorizable by (2,3,5,7) 
        do
        {
            ibox[i] += 1;
			b = ibox[i];          

			//n2 = n3 = n5 = n7 = 0;
			n2 = n3 = n5 = ns = 0;
			done_factoring = false;
			while (!done_factoring)
			{
				if (b % 2 == 0) 
				{
					n2++;
					b /= 2;
					continue;
				}
				if (b % 3 == 0) 
				{
					n3++;
					b /= 3;
					continue;
				}
				if (b % 5 == 0) 
				{
					n5++;
					b /= 5;
					continue;
				}
				//if (b%7==0) { n7++; b /= 7; continue; }
				done_factoring = true;
			}
        }
        while (b != 1);
        //  b==1 means fftbox[i] is (2,3,5,7) factorizable 
    }
    this->nx = ibox[0];
    this->bigny = ibox[1];
    this->nz = ibox[2];
    this->bignxy = this->nx * this->bigny;
    this->bignxyz = this->bignxy * this->nz;

    delete[] ibox;    
    return;
}

//
//Init the grids for FFT
//Input: lattice vectors of the cell, nx, ny, nz
//Output: nx, ny, nz, nxyz, latvec, G, GT, GGT
//
void PW_Basis:: initgrids(
    double lat0_in,
    ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
    int nx_in, int bigny_in, int nz_in
)
{
    this->lat0 = lat0_in;
    this->latvec = latvec_in;
    this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
    this->nx = nx_in;
    this->bigny = bigny_in;
    this->nz = nz_in;
    this->bignxy = this->nx * this->bigny;
    this->bignxyz = this->bignxy * this->nz;

    return;
}


//Init some parameters
void PW_Basis:: initparameters(
    bool gamma_only_in,
    double pwecut_in,
    int poolnproc_in,
    int poolrank_in,
    int distribution_type_in
)
{
    this->gamma_only = gamma_only_in;
    // if use gamma point only, when convert real function f(r) to F(k) = FFT(f),
    // we have F(-k) = F(k)*, so that only half of planewaves are needed.
    if (this->gamma_only)   this->ny = int(this->bigny / 2) + 1;
    else                    this->ny = bigny;
    this->nxy = this->nx * this->ny;
    this->nxyz = this->nxy * this->nz;

    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / this->lat0 / this->lat0;
    this->ggecut = pwecut_in / tpiba2;
    this->poolnproc = poolnproc_in;
    this->poolrank = poolrank_in;
    this->distribution_type = distribution_type_in;
}

}