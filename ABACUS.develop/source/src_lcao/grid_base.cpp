#include "src_pw/global.h"
#include "src_global/ylm.h"
#include "ORB_read.h"
#include "grid_base.h"
#include "global_fp.h" // mohan add 2021-01-30

Grid_Base::Grid_Base()
{ 
	this->Rcut_max = new double[1];
	this->Rcut_max_direct = new Vector3<double>[1];
	this->ijk_index = new int[1];
	this->norm1 = new double[1];
	this->norm2 = new double[1];
	this->dR1 = new Vector3<double>[1];
	this->dR2 = new Vector3<double>[1];
	this->grid_number_last = 0; // must initialized!!
	this->cartesian = new Vector3<double>[1];
	this->yy1 = new double*[1];
	this->yy2 = new double*[1];
	this->yy1[0] = new double[1];
	this->yy2[0] = new double[1];
	// must initialized!!
	this->n1_last = 1; 
	this->n2_last = 1;
	this->test = 1;
}

Grid_Base::~Grid_Base()
{
	delete[] Rcut_max;
	delete[] Rcut_max_direct;
	delete[] ijk_index;
	delete[] norm1;
	delete[] norm2;
	delete[] dR1;
	delete[] dR2;
	delete[] cartesian;
	for(int i=0; i<n1_last; i++) delete[] yy1[i];
	for(int i=0; i<n2_last; i++) delete[] yy2[i]; // so can delete yy2[0]
	delete[] yy1;
	delete[] yy2;
}

void Grid_Base::init(
    const Matrix3 &latvec_in,
    const double &lat0_in,
    const int &nx_in,
    const int &ny_in,
    const int &nz_in,
    const int &nxyz_in)
{
	this->latvec = latvec_in;
	this->lat0 = lat0_in;

	this->latvec0 = this->latvec;
	this->latvec0 *= this->lat0;

	this->a1.x = latvec0.e11;
	this->a1.y = latvec0.e12;
	this->a1.z = latvec0.e13;

	this->a2.x = latvec0.e21;
	this->a2.y = latvec0.e22;
	this->a2.z = latvec0.e23;

	this->a3.x = latvec0.e31;
	this->a3.y = latvec0.e32;
	this->a3.z = latvec0.e33;

	this->a1_len = a1.norm();
	this->a2_len = a2.norm();
	this->a3_len = a3.norm();

	this->nx = nx_in;
	this->ny = ny_in;
	this->nz = nz_in;
	this->nxyz = nxyz_in;

	this->da1 = a1_len / nx;
	this->da2 = a2_len / ny;
	this->da3 = a3_len / nz;

	this->get_rcut_max();

	OUT(ofs_running,"lat0 = ", lat0);
	OUT(ofs_running,"|a1| = ", a1_len);
	OUT(ofs_running,"|a2| = ", a2_len);
	OUT(ofs_running,"|a3| = ", a3_len);
	OUT(ofs_running,"da_d.x = ", da_d.x);
	OUT(ofs_running,"da_d.y = ", da_d.y);
	OUT(ofs_running,"da_d.z = ", da_d.z);
	OUT(ofs_running,"nx = ", nx);
	OUT(ofs_running,"ny = ", ny);
	OUT(ofs_running,"nz = ", nz);
	OUT(ofs_running,"nxyz = ", nxyz);
	OUT(ofs_running,"da1 = ", da1);
	OUT(ofs_running,"da2 = ", da2);
	OUT(ofs_running,"da3 = ", da3);

	delete[] this->cartesian;
	this->cartesian = new Vector3<double>[ this->nxyz ];
	Memory::record("Grid_Base","cartesian",nxyz*3,"double");

	Vector3<double> direct;
	for (int i = 0; i < this->nx; i++)
	{
		direct.x = static_cast<double>(i) / this->nx;
		for (int j = 0; j < this->ny; j++)
		{
			direct.y = static_cast<double>(j) / this->ny;
			for (int k = 0; k < this->nz; k++)
			{
				// unit : a.u.
				direct.z = static_cast<double>(k) / this->nz;

				const int index = j * ny * nz + i  * nz + k;

				Mathzone::Direct_to_Cartesian(
				    direct.x, direct.y, direct.z,
				    latvec0.e11, latvec0.e12, latvec0.e13,
				    latvec0.e21, latvec0.e22, latvec0.e23,
				    latvec0.e31, latvec0.e32, latvec0.e33,
				    cartesian[index].x, cartesian[index].y, cartesian[index].z);
			}
		}
	}

	int tot_pairs = 0;
	for(int T1=0; T1<ucell.ntype; T1++)
	{
		for(int I1=0; I1<ucell.atoms[T1].na; I1++)
		{
			//GridD.Find_atom(ucell.atoms[T1].tau[I1]);
			GridD.Find_atom(ucell.atoms[T1].tau[I1], T1, I1);
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
			{
				++tot_pairs;
			}
		}
	}
	OUT(ofs_running,"total atom pairs = ",tot_pairs);	

	return;
}

void Grid_Base::get_rcut_max(void)
{
	assert( ORB.get_ntype() > 0 );

	this->da_d.x = 1.0 / this->a1_len;
	this->da_d.y = 1.0 / this->a2_len;
	this->da_d.z = 1.0 / this->a3_len;

	delete[] Rcut_max;
	delete[] Rcut_max_direct;
	this->Rcut_max = new double[ ORB.get_ntype() ];
	this->Rcut_max_direct = new Vector3<double>[ ORB.get_ntype() ];

	for(int it=0; it<ORB.get_ntype(); it++)
	{
		this->Rcut_max[it] = 0.0;
		for (int L = 0; L < ucell.atoms[it].nwl + 1; L++)
		{
			for (int N = 0; N < ucell.atoms[it].l_nchi[L]; N++)
			{
				Rcut_max[it] = max(ORB.Phi[it].getRcut(),
			               Rcut_max[it]);
			}
		}
		this->Rcut_max_direct[it].x = this->da_d.x * Rcut_max[it];
		this->Rcut_max_direct[it].y = this->da_d.y * Rcut_max[it];
		this->Rcut_max_direct[it].z = this->da_d.z * Rcut_max[it];
		
		ofs_running << "\n" << "Type=" << it << " Rcut_max(a.u.)=" << Rcut_max[it] 
			<< " Direct = " 
			<< Rcut_max_direct[it].x << " "
			<< Rcut_max_direct[it].y << " "
			<< Rcut_max_direct[it].z;
	}

	return;
}

// be called by cal_region.
void Grid_Base::get_small_box( 
		const Vector3<double> &tau, 
		const int &T,
		Vector3<double> &tau_max_direct,
		Vector3<double> &tau_min_direct)
{
	Vector3<double> tau_dir;

	// NOTE: tau_dir can be calculated and stored firstly.
	Mathzone::Cartesian_to_Direct
	(
	    tau.x, tau.y, tau.z,
	    latvec.e11, latvec.e12, latvec.e13,
	    latvec.e21, latvec.e22, latvec.e23,
	    latvec.e31, latvec.e32, latvec.e33,
	    tau_dir.x, tau_dir.y, tau_dir.z
	);

	if (tau_dir.x < 0.0 || tau_dir.x > 1.0 || tau_dir.y < 0.0 || tau_dir.y > 1.0 || tau_dir.z < 0.0 || tau_dir.z > 1.0)
	{
		cout << "\n tau.x = " << tau.x;
		cout << "\n tau.y = " << tau.y;
		cout << "\n tau.z = " << tau.z;
		
		cout << "\n tau_dir.x = " << tau_dir.x;
		cout << "\n tau_dir.y = " << tau_dir.y;
		cout << "\n tau_dir.z = " << tau_dir.z;
		WARNING_QUIT("Grid_Base::get_small_box","Positions(x,y,z) Of tau and R2 in Direct Coordinates should be between 0 and 1!");
	}

	tau_max_direct = tau_dir + this->Rcut_max_direct[T];
	tau_min_direct = tau_dir - this->Rcut_max_direct[T];
	return;
}


// be called by cal_region
void Grid_Base::edge_grid_points(
	const Vector3<double> &R10,
	const Vector3<double> &R20,
	const Vector3<double> &max_direct_coodinate,
	const Vector3<double> &min_direct_coodinate)
{
	//index of the 8 edge points in FFT box
	timer::tick("Grid_Base","edge_grid_points");
	edge_min.x = static_cast<int>(min_direct_coodinate.x * this->nx) ;
	edge_max.x = static_cast<int>(max_direct_coodinate.x * this->nx) + 1;
	edge_min.y = static_cast<int>(min_direct_coodinate.y * this->ny) ;
	edge_max.y = static_cast<int>(max_direct_coodinate.y * this->ny) + 1;
	edge_min.z = static_cast<int>(min_direct_coodinate.z * this->nz) ;
	edge_max.z = static_cast<int>(max_direct_coodinate.z * this->nz) + 1;

	assert( edge_min.x >= 0);
	assert( edge_max.x < nx );
	assert( edge_min.y >= 0);
	assert( edge_max.y < ny );
	assert( edge_min.z >= 0);
	assert( edge_max.z < nz );

	this->grid_number = (edge_max.x-edge_min.x) * (edge_max.y-edge_min.y) * (edge_max.z-edge_min.z);

	if(grid_number < 0)
	{
		cout << "\n edge_x = " << edge_min.x << " " << edge_max.x;
		cout << "\n edge_y = " << edge_min.y << " " << edge_max.y;
		cout << "\n edge_z = " << edge_min.z << " " << edge_max.z << endl;
	}


	if( grid_number > grid_number_last )
	{
		delete[] norm1;
		delete[] norm2;
		delete[] dR1;
		delete[] dR2;
		delete[] ijk_index;
		this->norm1 = new double[grid_number];
		this->norm2 = new double[grid_number];
		this->dR1 = new Vector3<double>[grid_number];
		this->dR2 = new Vector3<double>[grid_number];
		this->ijk_index = new int[grid_number];
	}

	if( (grid_number > grid_number_last) || (n1 > n1_last) || (n2 > n2_last)  )
	{
		ofs_running << "\n Renew!" << " n1=" << n1 << " n1_last=" << n1_last 
			<< " n2=" << n2 << " n2_last=" << n2
			<< " grid_number=" << grid_number << " grid_number_last=" << grid_number_last;
		for(int i=0; i<n1_last; i++) delete[] yy1[i];
		for(int i=0; i<n2_last; i++) delete[] yy2[i];
		delete[] yy1;
		delete[] yy2;

		this->yy1 = new double*[ this->n1 ];
		this->yy2 = new double*[ this->n2 ];

		const int bigger = std::max( grid_number, grid_number_last );
		for(int i=0; i<n1; i++) yy1[i] = new double[bigger];
		for(int i=0; i<n2; i++) yy2[i] = new double[bigger];

		this->n1_last = n1;
		this->n2_last = n2;
	}

	if( grid_number > grid_number_last )
	{
		this->grid_number_last = grid_number;
	}

	int count=0;
	timer::tick("mohan_test","1");
	for (int i = edge_min.x; i < edge_max.x; i++)
	{
		for (int j = edge_min.y; j < edge_max.y; j++)
		{
			for (int k = edge_min.z; k < edge_max.z; k++)
			{
				// unit : a.u.
				int index = j * ny * nz + i  * nz + k;

				// be careful!! the order between R1 and cartesian
				// make the two atom has the same center : cartesian
				this->dR1[count] = R10 - cartesian[index];
				this->dR2[count] = R20 - cartesian[index];
				this->norm1[count] = dR1[count].norm();
				this->norm2[count] = dR2[count].norm();
				this->ijk_index[count] = index;
				++count;
			}
		}
	}
	timer::tick("mohan_test","1");

	timer::tick("mohan_test","2");
	static double yy_tmp[400];
	count = 0;
	for (int i = edge_min.x; i < edge_max.x; i++)
	{
		for (int j = edge_min.y; j < edge_max.y; j++)
		{
			for (int k = edge_min.z; k < edge_max.z; k++)
			{
				Ylm::get_ylm_real(this->lmax1, dR1[count], yy_tmp );
				for(int in=0; in<n1; in++) yy1[in][count] = yy_tmp[in];
				++count;
			}
		}
	}

	count = 0;
	for (int i = edge_min.x; i < edge_max.x; i++)
	{
		for (int j = edge_min.y; j < edge_max.y; j++)
		{
			for (int k = edge_min.z; k < edge_max.z; k++)
			{
				Ylm::get_ylm_real(this->lmax2, dR2[count], yy_tmp );
				for(int in=0; in<n2; in++) yy2[in][count] = yy_tmp[in];
				++count;
			}
		}
	}
	timer::tick("mohan_test","2");

	if(count!=grid_number)
	{
		ofs_warning << "\n count = " << count;
		ofs_warning << "\n grid_number = " << grid_number;
		WARNING_QUIT("Grid_Base::edge_grid_points","count!=grid_number");
	}
	timer::tick("Grid_Base","edge_grid_points");
	return;
}
