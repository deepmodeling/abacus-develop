#include "pw_basis.h"
#include "module_base/constants.h"

namespace ModulePW
{
#ifdef __MPI
void PW_Basis:: initmpi(
    const int poolnproc_in,
        const int poolrank_in,
        MPI_Comm pool_world_in
)
{
    this->poolnproc = poolnproc_in;
    this->poolrank = poolrank_in;
    this->pool_world = pool_world_in;
}
#endif

/// 
/// Init the grids for FFT
/// Input: lattice vectors of the cell, Energy cut off for G^2/2
/// Output: fftnx, fftny, fftnz, fftnxyz, latvec, G, GT, GGT
/// 
void PW_Basis:: initgrids(
        const double lat0_in, //unit length (unit in bohr)
        const ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
        const double gridecut
)
{
    //init lattice
    this->lat0 = lat0_in;
    this->tpiba = ModuleBase::TWO_PI / this->lat0;
    this->tpiba2 = this->tpiba*this->tpiba;
    this->latvec = latvec_in;
    this->omega = std::abs(latvec.Det()) * lat0 * lat0 * lat0;
    this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;

    //------------------------------------------------------------
    //-------------------------init grids-------------------------
    //-----------------------------------------------------------
    this->gridecut_lat = gridecut / this->tpiba2;
    ModuleBase::Vector3<double> lat;
    int *ibox = new int[3];// ibox[i] are the minimal FFT dimensions,
    
    lat.x = latvec.e11;
    lat.y = latvec.e12;
    lat.z = latvec.e13;
    ibox[0] = int(sqrt(this->gridecut_lat) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e21;
    lat.y = latvec.e22;
    lat.z = latvec.e23;
    ibox[1] = int(sqrt(this->gridecut_lat) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e31;
    lat.y = latvec.e32;
    lat.z = latvec.e33;
    ibox[2] = int(sqrt(this->gridecut_lat) * sqrt(lat * lat)) + 1;
    
    // We should check if ibox is the minimum number to cover the planewave ball. 
    // Find the minimum number of ibox by traveling all possible ibox
    int n1,n2,n3; 
    n1 = n2 = n3 = 0;
    for(int igz = -ibox[2]+this->poolrank; igz <= ibox[2]; igz += this->poolnproc)
    {
        for(int igy = -ibox[1]; igy <= ibox[1]; ++igy)
        {
            for(int igx = -ibox[0]; igx <= ibox[0]; ++igx)
            {
                ModuleBase::Vector3<double> f;
                f.x = igx;
                f.y = igy;
                f.z = igz;
                double modulus = f * (this->GGT * f);
                if(modulus <= this->gridecut_lat)
                {
                    if(n1 < abs(igx)) n1 = abs(igx);
                    if(n2 < abs(igy)) n2 = abs(igy);
                    if(n3 < abs(igz)) n3 = abs(igz);
                }
            }
        }
    }
    ibox[0] = 2*n1+1;
    ibox[1] = 2*n2+1;
    ibox[2] = 2*n3+1;
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, ibox, 3, MPI_INT, MPI_MAX , this->pool_world);
#endif

    // Find the minimal FFT box size the factors into the primes (2,3,5,7).
    for (int i = 0; i < 3; i++)
    {
    	int b = 0;
        int n2 = 0;
        int n3 = 0;
        int n5 = 0;
        //int n7 = 0;
        bool done_factoring = false;
	
        // increase ibox[i] by 1 until it is totally factorizable by (2,3,5,7) 
        do
        {
			b = ibox[i];          

			//n2 = n3 = n5 = n7 = 0;
			n2 = n3 = n5 = 0;
			done_factoring = false;
            if ((this->full_pw && this->full_pw_dim == 2) && b % 2 != 0) done_factoring = true; // full_pw_dim = 2 means FFT dimensions should be even.
			while (!done_factoring)
			{
				if (b % 2 == 0 && (!this->full_pw || this->full_pw_dim != 1)) // full_pw_dim = 1 means FFT dimension should be odd.
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
            ibox[i] += 1;
        }
        while (b != 1);
        ibox[i] -= 1;
        //  b==1 means fftbox[i] is (2,3,5,7) factorizable 
    }
    this->nx = ibox[0];
    this->ny = ibox[1];
    this->nz = ibox[2];
    this->nxy =this->nx * this->ny;
    this->nxyz = this->nxy * this->nz;

    delete[] ibox;    
    return;
}

/// 
/// Init the grids for FFT
/// Input: lattice vectors of the cell, nx, ny, nz
/// Output: nx, ny, nz, nxyz, latvec, G, GT, GGT
/// 
void PW_Basis:: initgrids(
    const double lat0_in,
    const ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
    const int nx_in, int ny_in, int nz_in
)
{
    this->lat0 = lat0_in;
    this->tpiba = ModuleBase::TWO_PI / this->lat0;
    this->tpiba2 = this->tpiba*this->tpiba;
    this->latvec = latvec_in;
    this->omega = std::abs(latvec.Det()) * lat0 * lat0 * lat0;
    this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
    this->nxy = this->nx * this->ny;
    this->nxyz = this->nxy * this->nz;

    int *ibox = new int[3];
    ibox[0] = int((this->nx-1)/2)+1;
    ibox[1] = int((this->ny-1)/2)+1;
    ibox[2] = int((this->nz-1)/2)+1;
    this->gridecut_lat = 1e20;
    int count = 0;
    for(int igz = -ibox[2]; igz <= ibox[2]; ++igz)
    {
        for(int igy = -ibox[1]; igy <= ibox[1]; ++igy)
        {
            for(int igx = -ibox[0]; igx <= ibox[0]; ++igx)
            {
                ++count;
                if(count%this->poolnproc != this->poolrank) continue;
                if(abs(igx)<=ibox[0]-1 && abs(igy)<=ibox[1]-1 && abs(igz)<=ibox[2]-1 ) continue;
                ModuleBase::Vector3<double> f;
                f.x = igx;
                f.y = igy;
                f.z = igz;
                double modulus = f * (this->GGT * f);
                if(modulus < this->gridecut_lat)
                {
                    this->gridecut_lat = modulus;
                }
            }
        }
    }
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &this->gridecut_lat, 1, MPI_DOUBLE, MPI_MIN , this->pool_world);
#endif
    this->gridecut_lat -= 1e-6;

    delete[] ibox;
    return;
}


//Init some parameters
void PW_Basis:: initparameters(
    const bool gamma_only_in,
    const double pwecut_in,
    const int distribution_type_in,
    const bool xprime_in
)
{
    this->xprime = xprime_in;
    this->gamma_only = gamma_only_in;
    // if use gamma point only, when convert real function f(r) to F(k) = FFT(f),
    // we have F(-k) = F(k)*, so that only half of planewaves are needed.
    this->fftny = this->ny;
    this->fftnx = this->nx;
    if (this->gamma_only)   
    {
        if(this->xprime) this->fftnx = int(this->nx / 2) + 1;
        else            this->fftny = int(this->ny / 2) + 1;
    }      
    
    this->fftnz = this->nz;
    this->fftnxy = this->fftnx * this->fftny;
    this->fftnxyz = this->fftnxy * this->fftnz;

    this->ggecut = pwecut_in / this->tpiba2;
    //ggecut should be no larger than gridecut
    if(this->ggecut > this->gridecut_lat) 
    {
        this->ggecut = this->gridecut_lat;
    }
    this->distribution_type = distribution_type_in;
}

// Set parameters about full planewave, used only in OFDFT for now. sunliang added 2022-08-30
void PW_Basis::setfullpw(
    const bool inpt_full_pw,
    const int inpt_full_pw_dim
)
{
    this->full_pw = inpt_full_pw;
    this->full_pw_dim = inpt_full_pw_dim;
    if (!this->full_pw) this->full_pw_dim = 0;
}
}