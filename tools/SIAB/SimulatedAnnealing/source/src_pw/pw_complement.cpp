#include "pw_complement.h"
#include "heapsort.h"

int PW_complement::get_total_pw_number(
    const double& ggcut_start,
    const double& ggcut_end,
    const int& nx,
    const int& ny,
    const int& nz,
    const Matrix3& GGT)
{
    TITLE("PW_complement","get_total_pw_number");
    if (ggcut_end<=0.0)
    {
        WARNING_QUIT("PW_complement::get_total_pw_number","ggcut <= 0.0");
    }
    // First, figure out the number of G-vectors within the cutoff G2max.
    int ibox[3]={0,0,0};

    // mohan modify 2008-3-25 nx -> ncx ,ny-> ncy, nz->ncz
    ibox[0] = int(nx / 2) + 1;
    ibox[1] = int(ny / 2) + 1;
    ibox[2] = int(nz / 2) + 1;

    // first create the auxiliary arrays for the 1D G vectors
    Vector3<double> f;
    int ngm = 0 ;
    for (int i = -ibox[0]; i <= ibox[0]; i++)
    {
        for (int j = -ibox[1]; j <= ibox[1]; j++)
        {
            for (int k = -ibox[2]; k <= ibox[2]; k++)
            {
                f.x = i;
                f.y = j;
                f.z = k;
                double g2 = f * (GGT * f);  //G2= |f|^2 in the unit of (2Pi/lat0)^2
                
				if (g2 < ggcut_end && g2 >= ggcut_start)
                {
                    ngm++;
                }
            }
        }
    }
    return ngm;
}

void PW_complement::get_total_pw(
    double* gg,
    Vector3<double> *ig,
    const double& ggcut_start,
    const double& ggcut_end,
    const int& nx,
    const int& ny,
    const int& nz,
    const Matrix3& GGT, // GGT = G*GT.
    int& ngm// number of total plane waves.
)
{
    TITLE("PW_complement","get_total_pw");
    timer::tick("PW_complement","get_total_pw");
    if (ggcut_end<=0.0)
    {
        WARNING_QUIT("PW_complement::get_total_pw","ggcut <= 0.0");
    }
    // First, figure out the number of G-vectors within the cutoff G2max.
    int ibox[3]={0,0,0};

    // mohan modify 2008-3-25 nx -> ncx ,ny-> ncy, nz->ncz
    ibox[0] = int(nx / 2) + 1;
    ibox[1] = int(ny / 2) + 1;
    ibox[2] = int(nz / 2) + 1;

    // first create the auxiliary arrays for the 1D G vectors
    Vector3<double> f;
    int ng = 0;
    for (int i = -ibox[0]; i <= ibox[0]; i++)
    {
        for (int j = -ibox[1]; j <= ibox[1]; j++)
        {
            for (int k = -ibox[2]; k <= ibox[2]; k++)
            {
                f.x = i;
                f.y = j;
                f.z = k;
                double g2 = f * (GGT * f);  //G2= |f|^2 in the unit of (2Pi/lat0)^2
                if (g2 < ggcut_end && g2 >= ggcut_start)
                {
                    /*** g vector indices f=(i,j,k) ***/
                    ig[ng] = f ;
                    gg[ng] = g2;
                    //cout<<setw(12)<<f.x<<setw(12)<<f.y<<setw(12)<<f.z<<setw(12)<<g2<<endl;
                    ++ng;
                }
            }
        }
    }

    //cout << "\n ng = " << ng;
    timer::tick("PW_complement","get_total_pw");
    return;
}

void PW_complement::get_FFT_dimension(
    const Matrix3 &latvec,
    const double &ggcut,
    int &nx_tmp,
    int &ny_tmp,
    int &nz_tmp)
{
    TITLE("PW_complement","get_FFT_dimension");
    // read in the FFT dimension (Nx, Ny, Nz) from paratab,
    // if Nx*Ny*Nz >0,use the input
    // FFT grid, otherwise generate the FFT grid in the code.

    int i = 0;
    Vector3<double> lat;
    int ibox[3]={0,0,0};

    // ibox[i] are the minimal FFT dimensions,
    lat.x = latvec.e11;
    lat.y = latvec.e12;
    lat.z = latvec.e13;
    ibox[0] = 2 * int(sqrt(ggcut) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e21;
    lat.y = latvec.e22;
    lat.z = latvec.e23;
    ibox[1] = 2 * int(sqrt(ggcut) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e31;
    lat.y = latvec.e32;
    lat.z = latvec.e33;
    ibox[2] = 2 * int(sqrt(ggcut) * sqrt(lat * lat)) + 1;
    //lat*lat=lat.x*lat.x+lat.y*lat.y+lat.z+lat.z

    /* Find the minimal FFT box size the factors into the primes (2,3,5,7). */

//  log << "\n ibox[0]= "<< ibox[0] << " ibox[1] = "<<ibox[1]
//      << " ibox[2] = "<< ibox[2];

    for (i = 0; i < 3; i++)
    {
        int b = 0;
        int n2 = 0;
        int n3 = 0;
        int n5 = 0;
        int n7 = 0;
        bool done_factoring = false;

        /* increase ibox[i] by 1 until it is totally factorizable
           by (2,3,5,7) */

        do
        {
            ibox[i] += 1;
            b = ibox[i];
            //n2 = n3 = n5 = n7 = 0;
            n2 = n3 = n5 = 0;
            done_factoring = false;

            while (!done_factoring)
            {
                if (b % 2 == 0) {
                    n2++;
                    b /= 2;
                    continue;
                }
                if (b % 3 == 0) {
                    n3++;
                    b /= 3;
                    continue;
                }
                if (b % 5 == 0) {
                    n5++;
                    b /= 5;
                    continue;
                }
                //if (b%7==0) { n7++; b /= 7; continue; }
                done_factoring = true;
            }
        }
        while (b != 1);
        /*  b==1 means fftbox[i] is (2,3,5,7) factorizable */
    }

    // Nx, Ny, Nz are the FFT dimensions

    nx_tmp = ibox[0];
    ny_tmp = ibox[1];
    nz_tmp = ibox[2];

    int nxyz_tmp = 0;
    nxyz_tmp = nx_tmp * ny_tmp * nz_tmp;

    return;
}


//==========================================================
// MEMBER FUNCTION :
// NAME : PW_Basis::setup_GVectors
// Second part of the initialization : find out all G
// vectors that |G|^2<=G2max and map it into a one
// dimentional array G1d in the increase order of |G|^2.
// Next generate the indices between the 1D array and
// the 3D G-grid and the FFT grid.
// generate : gg_global, g_global, ig_global
//==========================================================
void PW_complement::setup_GVectors(
    const Matrix3& G,
    const int &ngmc_g_in,
    double* gg,
    Vector3<double>* ig,
    Vector3<double>* g)
{
    TITLE("PW_complement","setup_GVectors");
    timer::tick("PW_complement","setup_GVectors");

    int *ind = new int[ngmc_g_in];// auxiliary array for the 1d G vector index
    ZEROS( ind, ngmc_g_in );
    ind[0] = 0;//ind[0]=0, meaning ind[] is not initialized.

    //----------------------------------------------------------
    // CALL GLOBAL FUNCTION :
    // sort the norm of the G vectors of G1d[] in the ascending
    // order
    //----------------------------------------------------------
    heapsort(ngmc_g_in, gg, ind);

    /***************************************************************/
    // adjust G1d[NG] according to the new order.
    // establish the link between the 1d array and 3d grid,
    // i.e. give the position "i" in 1d array,
    // G1d[i] is the vector in 3d grid, G1d2[i] is its norm.
    /***************************************************************/

    Vector3<double> *igsort = new Vector3<double>[ngmc_g_in];
    for (int i=0;i<ngmc_g_in;i++)
    {
        igsort[i] = ig[ind[i]];
//		cout << i << " " << ind[i] << " " << ig[ind[i]].x << " " << ig[ind[i]].y << " " << ig[ind[i]].z << endl;
    }

    for (int i=0; i<ngmc_g_in; i++)
    {
        ig[i] = igsort[i];
    }

    delete[] igsort;
    delete[] ind;

    //----------------------------------------------------------
    // EXPLAIN : store the G vectors in 1d array (Cartian
    // coordinate)
    //----------------------------------------------------------
    for (int i=0; i<ngmc_g_in; i++)
    {
        g[i] = ig[i] * G;
//		cout << i << " " << ig[i].x << " " << ig[i].y << " " << ig[i].z << endl;
        // g vector (in 2*pi/lat0), E_k=g*g*(2*pi/lat0)^2
    }

    timer::tick("PW_complement","setup_GVectors");
    return;
}
