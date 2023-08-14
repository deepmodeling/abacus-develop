#ifndef PW_BASIS_BIG_H
#define PW_BASIS_BIG_H
#include "module_base/constants.h"
#include "module_base/global_function.h"
#ifdef __MPI
#include "mpi.h"
#endif

// temporary class, because previous ABACUS consider big grid for fft grids 
// which are used for grid integration in LCAO.
// In fact, it is unnecessary. It will be moved after grid integration is refactored.
namespace ModulePW
{

class PW_Basis_Big: public PW_Basis
{
public:
    
    // combine [bx,by,bz] FFT grids into a big one
	// typical values are bx=2, by=2, bz=2
	// nbx=nx/bx, nby=ny/by, nbz=nz/bz, 
    PW_Basis_Big(){
        bx = 1;
        by = 1;
        bz = 1;
    }
    PW_Basis_Big(std::string device_, std::string precision_) : PW_Basis(device_, precision_) {}

    ~PW_Basis_Big(){};
    void setbxyz(const int bx_in, const int by_in, const int bz_in)
    {
        bx = bx_in;
        by = by_in;
        bz = bz_in;
        bxyz = bx * by * bz;
    }
    int bx,by,bz,bxyz;
    int nbx, nby, nbz;
    int nbzp;
    int nbxx;
    int nbzp_start;

    void autoset_big_cell_size(int& b_size, const int& nc_size, const int nproc = 0)
    {
        //original default setting is 4
        b_size = 4;
        //only for bz
        if(nproc > 0)
        {
            int candidate_lists[4] = {4, 3, 5, 2};
            int max_bz[4];
            for(int i=0;i<4;i++)
            {
                int tmp = candidate_lists[i];
                max_bz[i] = nc_size / tmp;
                if(nc_size % tmp!=0)
                {//ignore candidates which can't be factored by nc_size
                    max_bz[i]=0;
                    continue;
                } 
                if(max_bz[i] % nproc == 0)
                {
                    b_size = tmp;
                    return;
                }
            }

            //choose maximum residual
            double res = 0.0;
            double res_temp = 0.0;
            for(int i=0;i<4;i++)
            {
                if(max_bz[i]==0) continue;
                res_temp = double(max_bz[i] % nproc) / nproc;
                if(res < res_temp)
                {
                    res = res_temp;
                    b_size = candidate_lists[i];
                }
            }
            return;
        }
        //for bx and by, choose maximum residual of (5,4,3)
        else
        {
            int res = 0;
            int res_temp = 0;
            for(int i=5;i>2;i--)
            {
                res_temp = nc_size % i;
                if(res_temp == 0)
                {
                    b_size = i;
                    return;
                }
                else if(res < res_temp)
                {
                    res = res_temp;
                    b_size = i;
                } 
            }
            return;
        }
    }


    virtual void initgrids(const double lat0_in,const ModuleBase::Matrix3 latvec_in,
        const double gridecut){
        //init lattice
    this->lat0 = lat0_in;
    this->latvec = latvec_in;
    this->omega = std::abs(latvec.Det()) * lat0 * lat0 * lat0;
    this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;

    //------------------------------------------------------------
    //-------------------------init grids-------------------------
    //------------------------------------------------------------
    this->tpiba = ModuleBase::TWO_PI / this->lat0;
    this->tpiba2 = this->tpiba * this->tpiba;
    this->gridecut_lat = gridecut / tpiba2;
    ModuleBase::Vector3<double> lat;
    int *ibox = new int[3];
    
    lat.x = latvec.e11;
    lat.y = latvec.e12;
    lat.z = latvec.e13;
    ibox[0] = 2 * int(sqrt(gridecut_lat) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e21;
    lat.y = latvec.e22;
    lat.z = latvec.e23;
    ibox[1] = 2 * int(sqrt(gridecut_lat) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e31;
    lat.y = latvec.e32;
    lat.z = latvec.e33;
    ibox[2] = 2 * int(sqrt(gridecut_lat) * sqrt(lat * lat)) + 1;

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
                    if(n1 < std::abs(igx)) n1 = std::abs(igx);
                    if(n2 < std::abs(igy)) n2 = std::abs(igy);
                    if(n3 < std::abs(igz)) n3 = std::abs(igz);
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
    //autoset bx/by/bz if not set in INPUT
    if(!this->bz)
    {
        this->autoset_big_cell_size(this->bz, ibox[2], this->poolnproc);
    }
    if(!this->bx)
    {
        //if cz == cx, autoset bx==bz for keeping same symmetry 
        if(ibox[0] == ibox[2])
        {
            this->bx = this->bz;
        }
        else
        {
            this->autoset_big_cell_size(this->bx, ibox[0]);
        }
    }
    if(!this->by)
    {
        //if cz == cy, autoset by==bz for keeping same symmetry 
        if(ibox[1] == ibox[2])
        {
            this->by = this->bz;
        }
        else
        {
            this->autoset_big_cell_size(this->by, ibox[1]);
        }
    }
    this->bxyz = this->bx * this->by * this->bz;
    if(ibox[0]%this->bx != 0) ibox[0] += (this->bx - ibox[0] % this->bx);
    if(ibox[1]%this->by != 0) ibox[1] += (this->by - ibox[1] % this->by);
    if(ibox[2]%this->bz != 0) ibox[2] += (this->bz - ibox[2] % this->bz);

    this->nx = ibox[0];
    this->ny = ibox[1];
    this->nz = ibox[2];
    this->nxy =this->nx * this->ny;
    this->nxyz = this->nxy * this->nz;
    this->nbx = this->nx / bx;
    this->nby = this->ny / by;
    this->nbz = this->nz / bz;

    delete[] ibox;    
    return;

    }

    virtual void initgrids(
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
        // autoset bx/by/bz if not set in INPUT
        if (!this->bz)
        {
        this->autoset_big_cell_size(this->bz, nz, this->poolnproc);
        }
        if (!this->bx)
        {
        // if cz == cx, autoset bx==bz for keeping same symmetry
        if (nx == nz)
        {
            this->bx = this->bz;
        }
        else
        {
            this->autoset_big_cell_size(this->bx, nx);
        }
        }
        if (!this->by)
        {
        // if cz == cy, autoset by==bz for keeping same symmetry
        if (ny == nz)
        {
            this->by = this->bz;
        }
        else
        {
            this->autoset_big_cell_size(this->by, ny);
        }
        }
        this->bxyz = this->bx * this->by * this->bz;
        if(this->nx%this->bx != 0) this->nx += (this->bx - this->nx % this->bx);
        if(this->ny%this->by != 0) this->ny += (this->by - this->ny % this->by);
        if(this->nz%this->bz != 0) this->nz += (this->bz - this->nz % this->bz);
        this->nbx = this->nx / bx;
        this->nby = this->ny / by;
        this->nbz = this->nz / bz;
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
                    if(std::abs(igx)<=ibox[0]-1 && std::abs(igy)<=ibox[1]-1 && std::abs(igz)<=ibox[2]-1 ) continue;
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

    virtual void distribute_r()
    {   
        delete[] this->numz; this->numz = new int[this->poolnproc];
        delete[] this->startz; this->startz = new int[this->poolnproc];
        ModuleBase::GlobalFunc::ZEROS(this->numz, this->poolnproc);
        ModuleBase::GlobalFunc::ZEROS(this->startz, this->poolnproc);

        int npbz = this->nbz / this->poolnproc;
        int modbz = this->nbz % this->poolnproc;
        this->startz[0] = 0;
        for(int ip = 0 ; ip < this->poolnproc ; ++ip)
        {
            this->numz[ip] = npbz*this->bz;
            if(ip < modbz)   this->numz[ip]+=this->bz;
            if(ip < this->poolnproc - 1)   this->startz[ip+1] = this->startz[ip] + numz[ip];
            if(ip == this->poolrank) 
            {
                this->nplane = numz[ip];
                this->startz_current = startz[ip];
            }
        }
        this->nbzp = this->nplane / this->bz;
        this->nrxx = this->numz[this->poolrank] * this->nxy;
        this->nbxx = this->nbzp * this->nbx * this->nby;
        this->nbzp_start = this->startz[this->poolrank] / this->bz;
        return;
    }

};
}
#endif