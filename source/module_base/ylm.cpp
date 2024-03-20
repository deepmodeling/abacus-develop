#include "ylm.h"

#include <cassert>
#include <iomanip>

#include "constants.h"
#include "timer.h"
#include "tool_quit.h"

namespace ModuleBase
{

int Ylm::nlm = 0;
std::vector<double> Ylm::ylmcoef(100);

// here Lmax == max angular momentum + 1
void Ylm::get_ylm_real( const int &Lmax, const ModuleBase::Vector3<double> &vec, double ylmr[] )
{
	ModuleBase::timer::tick ("Ylm","get_ylm_real");
	//1e-9 is too large
	const double cut0 = 1e-12;
	// allocate space.
	Ylm::nlm = Lmax * Lmax;
	if(Lmax==1)
	{
		for(int i=0; i<Ylm::nlm; i++)
		{
			ylmr[i] = ModuleBase::SQRT_INVERSE_FOUR_PI;
		}
	}

	double cost = 0.0; // must initialized.
	double sint = cut0;
	double phi = 0.0;

	double vnorm = vec.norm();

	if(vnorm < cut0)
	{
		vnorm += cut0;
	}

	cost = vec.z / vnorm;

	if(fabs(cost) > 1.0 - cut0)
	{
		cost = sgn(cost) * (1.0 - cut0);
		//		std::cout << "\n" << "cost = " << cost << std::endl;
	}
	sint = sqrt(1.0 - cost*cost);

	if(vec.x > cut0)
	{
		phi = std::atan( vec.y / vec.x );
	}
	else if( vec.x < -cut0 )
	{
		phi = std::atan( vec.y / vec.x ) + ModuleBase::PI;
	}
	else
	{
		phi = ModuleBase::PI_HALF * ((vec.y >= 0.0) ? 1.0 : -1.0);
	}

	//===============================
	// NAME : p(Legendre Polynomials)
	//===============================
	double p[20][20];
	assert(Lmax <= 20);

	int m=0;
	double x1=0.0;	// x2;
	int lm = -1; // must initialized!

	for (int l=0; l<Lmax; l++)
	{
		const double c = sqrt((2*l+1) / ModuleBase::FOUR_PI);
		if (l == 0)
		{
			p[0][0] = 1.0;
		}
		else if (l == 1)
		{
			p[0][1] = cost;
			p[1][1] = -sint;
		}
		else
		{
			const int l1 = l-1;
			const int l2 = l-2;
			const int l3 = 2*l-1;
			//  recursion on l for P(:,l,m)
			for (m=0; m<=l2; m++)  // do m = 0, l - 2//mohan modify 2007-10-13
			{
				p[m][l] = (cost * l3 * p[m][l1] -
					      (l1 + m ) * p[m][l2]) / (l - m);
			} // end do
			p[l1][l] = cost * l3 * p[l1][l1];
			x1 = Ylm::Semi_Fact(l3) * pow(sint, static_cast<double>(l)) ;//liaochen modify 2009-09-06
			if (l%2 == 1)
			{
				x1 = -x1;
			}
			p[l][l] = x1;
		} // end if

		// Y_lm, m = 0
		++lm;
		ylmr[lm] = c*p[0][l];

		for(m=1;m<=l;m++)
		{
			// Y_lm, m > 0
			const double same = c * sqrt
					(
					static_cast<double>( Ylm::Fact(l - m)) /
					static_cast<double>( Ylm::Fact(l + m))
					)
					*ModuleBase::SQRT2;

			++lm;
			ylmr[lm] = same * p[m][l] * cos(m * phi);

			// Y_lm, m < 0
			++lm;
			ylmr[lm] = same * p[m][l] * sin(m * phi);
		}
	}// end do

	ModuleBase::timer::tick ("Ylm", "get_ylm_real");
	return;
}

void Ylm::get_ylm_real( const int &Lmax, const ModuleBase::Vector3<double> &vec, double ylmr[], double dylmdr[][3] )
{
	//1e-9 is too large
	const double cut0 = 1e-12;
	// allocate space.
	Ylm::nlm = Lmax * Lmax;
	if(Lmax==1)
	{
		for(int i=0; i<Ylm::nlm; i++)
		{
			ylmr[i] = ModuleBase::SQRT_INVERSE_FOUR_PI;
			for(int j = 0; j < 3; j++)
			{
				dylmdr[i][j] = 0.0;
			}
		}
	}

	double cost = 0.0; // must initialized.
	double sint = 0.0;
	double phi = 0.0;

	double vnorm = vec.norm();

	if(vnorm < cut0)
	{
		vnorm += cut0;
	}

	cost = vec.z / vnorm;

	if(fabs(cost) > 1.0-cut0)
	{
		cost = sgn(cost) * (1.0 - cut0);
	}
	sint = sqrt(1.0 - cost*cost);

	if(vec.x > cut0)
	{
		phi = std::atan( vec.y / vec.x );
	}
	else if( vec.x < -cut0 )
	{
		phi = std::atan( vec.y / vec.x ) + ModuleBase::PI;
	}
	else
	{
		phi = ModuleBase::PI_HALF * ((vec.y >= 0.0) ? 1.0 : -1.0);
	}

	//===============================
	// NAME : p(Legendre Polynomials)
	//===============================
	double p[20][20];
	double dp[20][20];
	assert(Lmax <= 20);

	int m = 0;
	int lm = -1; // must initialized!
	for (int l=0; l<Lmax; l++)
	{
		const double c = sqrt((2*l+1) / ModuleBase::FOUR_PI);
		if (l == 0)
		{
			p[0][0] = 1.0;
			dp[0][0] = 0.0;
		}
		else if (l == 1)
		{
			p[0][1] = cost;
			dp[0][1] = -sint;

			p[1][1] = -sint;
			dp[1][1] = -cost;
		}
		else
		{
			const int l1 = l-1;
			const int l2 = l-2;
			const int l3 = 2*l-1;
			//  recursion on l for P(:,l,m)
			for (m=0; m<=l2; m++)  // do m = 0, l - 2//mohan modify 2007-10-13
			{
				p[m][l] = (cost * l3 * p[m][l1] -
					      (l1 + m ) * p[m][l2]) / (l - m);
			}

			p[l1][l] = cost * l3 * p[l1][l1];

			p[l][l] = Ylm::Semi_Fact(l3) * pow(sint, static_cast<double>(l)) ;
			if (l%2 == 1)
			{
				p[l][l] = -p[l][l];
			}
		}

		for(m=0; m <= l; m++)
		{
			if( m == l )
			{
				dp[l][l] = l * cost * p[l][l] / sint;
			}
			else
			{
				dp[m][l] = (l * cost * p[m][l] - (l+m) * p[m][l-1]) / sint;
			}
		}

		// Y_lm, m = 0
		++lm;
		ylmr[lm] = c*p[0][l];

		dylmdr[lm][0] = c * dp[0][l] * cost * cos(phi) / vnorm;
		dylmdr[lm][1] = c * dp[0][l] * cost * sin(phi) / vnorm;
		dylmdr[lm][2] = -c * dp[0][l] * sint / vnorm;

		for(m=1; m <= l; m++)
		{
			// Y_lm, m > 0
			const double same = c * sqrt
					(
					static_cast<double>( Ylm::Fact(l - m)) /
					static_cast<double>( Ylm::Fact(l + m))
					)
					*ModuleBase::SQRT2;

			++lm;
			ylmr[lm] = same * p[m][l] * cos(m * phi);

			dylmdr[lm][0] = same * dp[m][l] * cos(m * phi) * cost * cos(phi) / vnorm
							+ same * p[m][l] * m * sin(m * phi) * sin(phi) / sint / vnorm;
			dylmdr[lm][1] = same * dp[m][l] * cos(m * phi) * cost * sin(phi) / vnorm
							- same * p[m][l] * m * sin(m * phi) * cos(phi) / sint / vnorm;
			dylmdr[lm][2] = -same * dp[m][l] * cos(m * phi) * sint / vnorm;

			// Y_lm, m < 0
			++lm;
			ylmr[lm] = same * p[m][l] * sin(m * phi);

			dylmdr[lm][0] = same * dp[m][l] * sin(m * phi) * cost * cos(phi) / vnorm
						 - same * p[m][l] * m * cos(m * phi) * sin(phi) / sint / vnorm;
			dylmdr[lm][1] = same * dp[m][l] * sin(m * phi) * cost * sin(phi) / vnorm
						 + same * p[m][l] * m * cos(m * phi) * cos(phi) / sint / vnorm;
			dylmdr[lm][2] = -same * dp[m][l] * sin(m * phi) * sint / vnorm;
		}
	}// end do

	return;
}

/***************************
 * Solid Spherical Harmonic
 * *************************/
void Ylm::rlylm
(
 	const int& Lmax, //max momentum of l + 1
 	const double& x,
	const double& y,
	const double& z,
	double rly[]
)
{
//	ModuleBase::TITLE("Ylm","rlylm");
//	ModuleBase::timer::tick("Ylm","rlylm");

	int MaxL = Lmax - 1;

	assert(MaxL >= 0);

	//get xy_dependence
	assert(MaxL <= 19);

	double Am[20];
	double Bm[20];

//	ZEROS(Am, 20);
//	ZEROS(Bm, 20);

	double x2, x3, x4, x5;
	double y2, y3, y4, y5;

	x2 = x * x;
	x3 = x2 * x;
	x4 = x3 * x;
	x5 = x4 * x;

	y2 = y * y;
	y3 = y2 * y;
	y4 = y3 * y;
	y5 = y4 * y;

	//x-y dependence
	//Am
	//Bm
	for(int im = 0; im < MaxL+1; im++)
	{
		if(im == 0)
		{
			Am[0] = 1.0;
			Bm[0] = 0.0;
		}
		else if(im == 1)
		{
			Am[1] = x;
			Bm[1] = y;
		}
		else if(im == 2)
		{
			Am[2] = x2- y2;
			Bm[2] = 2.0 * x * y;
		}
		else if(im == 3)
		{
			Am[3] = x3 - 3.0 * x * y2;
			Bm[3] = 3.0 * x2 * y - y3;
		}
		else if(im == 4)
		{
			Am[4] = x4 - 6.0 * x2 * y2 + y4;
			Bm[4] = 4.0 * (x3 * y - x * y3);
		}
		else if(im == 5)
		{
			Am[5] = x5 - 10.0 * x3 * y2 + 5.0 * x * y4;
			Bm[5] = 5.0 * x4 * y - 10.0 * x2 * y3 + y5;
		}
		else
		{
			for(int ip = 0; ip <= im; ip++)
			{
				double aux = Fact(im) / Fact(ip) / Fact(im - ip);
				Am[im] += aux * pow(x, ip) * pow(y, im-ip) * cos( (im-ip) * ModuleBase::PI / 2.0 );
				Bm[im] += aux * pow(x, ip) * pow(y, im-ip) * sin( (im-ip) * ModuleBase::PI / 2.0 );
			}
		}
	}

	//z dependence
	double zdep[20][20];

//	for(int il = 0; il < 20; il++)
//	{
//		ZEROS(zdep[il], 20);
//	}

	double z2 = z * z;
	double z3 = z2 * z;
	double z4 = z3 * z;
	//double z5 = z4 * z;

	double r = sqrt(x*x + y*y + z*z);
	double r2 = r * r;
	double r3 = r2 * r;
	double r4 = r3 * r;

	for(int il = 0; il < MaxL+1; il++)
	{
		if(il == 0)
		{
			zdep[0][0] = 1.0;
		}
		else if(il == 1)
		{
			zdep[1][0] = z;
			zdep[1][1] = 1.0;
		}
		else if(il == 2)
		{
			zdep[2][0] = 0.5 * (3.0 * z2 - r2);
			zdep[2][1] = sqrt(3.0) * z;
			zdep[2][2] = sqrt(3.0) * 0.5;
		}
		else if(il == 3)
		{
			zdep[3][0] = 2.5 * z3 - 1.5 * z * r2;
			zdep[3][1] = 0.25 * sqrt(6.0) * (5.0 * z2 - r2);
			zdep[3][2] = 0.5 * sqrt(15.0) * z;
			zdep[3][3] = 0.25 * sqrt(10.0);
		}
		else if(il == 4)
		{
			zdep[4][0] = 0.125 * (35.0 * z4 - 30.0 * r2 * z2 + 3.0 * r4);
			zdep[4][1] = sqrt(10.0) * 0.25 * z * (7.0 * z2 - 3.0 * r2);
			zdep[4][2] = sqrt(5.0) * 0.25 * (7.0 * z2 - r2);
			zdep[4][3] = sqrt(70.0) * 0.25 * z;
			zdep[4][4] = sqrt(35.0) * 0.125;
		}
		else if(il == 5)
		{
			zdep[5][0] = 0.125 * z *( 63.0 * z4 - 70.0 * z2 * r2 + 15.0 * r4);
			zdep[5][1] = 0.125 * sqrt(15.0) * (21.0 * z4 - 14.0 * z2 * r2 + r4);
			zdep[5][2] = 0.25 * sqrt(105.0) * z * (3.0 * z2 - r2);
			zdep[5][3] = 0.0625 * sqrt(70.0) * (9.0 * z2 - r2);
			zdep[5][4] = 0.375 * sqrt(35.0) * z;
			zdep[5][5] = 0.1875 * sqrt(14.0);
		}
		else
		{
			for(int im = 0; im <= il; im++)
			{
				int kmax = static_cast<int>( (il - im) / 2 );
				for(int ik = 0; ik <= kmax; ik++)
				{
					int twok = 2 * ik;

					double gamma;
					double aux0, aux1, aux2, aux3;

					aux0 = pow(-1.0, ik) * pow(2.0, -il);
					aux1 = Fact(il) / Fact(ik) / Fact(il-ik);
					aux2 = Fact(2*il - twok) / Fact(il) / Fact(il - twok);
					aux3 = Fact(il - twok) / Fact(il - twok - im);

					gamma = aux0 * aux1 * aux2 * aux3;

					assert(il - twok - im >= 0);
					zdep[il][im] += pow(r, twok) * pow(z, il-twok-im) * gamma;
				}

				if(im >= 1)
				{
					zdep[il][im] *= sqrt(2 * Fact(il - im) / Fact(il + im));
				}
			}
		}
	}

	//calc
	int ic = 0;
	for(int il = 0; il <= MaxL; il++)
	{
		double fac = sqrt( (2.0 * il + 1.0) / ModuleBase::FOUR_PI );

		//m=0
		rly[ic] = Am[0] * zdep[il][0] * fac;

		ic++;

		//m ! = 0
		for(int im = 1; im <= il; im++)
		{
			//m>0
			rly[ic] = Am[im] * zdep[il][im] * pow(-1.0, im) * fac;

			ic++;

			//m<0
			rly[ic] = Bm[im] * zdep[il][im] * pow(-1.0, im) * fac;

			ic++;
		}
	}

//	ModuleBase::timer::tick("Ylm", "rlylm");
	return;
}

//return ylm, not rlylm
void Ylm::sph_harm
(
 	const int& Lmax, //max momentum of l
 	const double& xdr,
	const double& ydr,
	const double& zdr,
	std::vector<double> &rly
)
{
	rly.reserve( (Lmax+1)*(Lmax+1) );

	//begin calculation
	/***************************
			 L = 0
	***************************/
	rly[0] = Ylm::ylmcoef[0]; //l=0, m=0
	if (Lmax == 0) return;

	/***************************
			 L = 1
	***************************/
	rly[1] = Ylm::ylmcoef[1]*zdr; //l=1, m=0
	rly[2] = -Ylm::ylmcoef[1]*xdr; //l=1, m=1
	rly[3] = -Ylm::ylmcoef[1]*ydr; //l=1, m=-1
	if (Lmax == 1) return;

	/***************************
			 L = 2
	***************************/
	rly[4] = Ylm::ylmcoef[2]*zdr*rly[1]-Ylm::ylmcoef[3]*rly[0];//l=2, m=0

	double tmp0 = Ylm::ylmcoef[4]*zdr;
	rly[5] = tmp0*rly[2];//l=2,m=1
	rly[6] = tmp0*rly[3];//l=2,m=-1

	double tmp2 = Ylm::ylmcoef[4]*xdr;
	rly[7]= Ylm::ylmcoef[5]*rly[4]-Ylm::ylmcoef[6]*rly[0] - tmp2*rly[2];//l=2,m=2
	rly[8] = -tmp2*rly[3];
//	rly[8] = tmp1+tmp2*rly[3];//l=2,m=-2
	if (Lmax == 2) return;

	/***************************
			 L = 3
	***************************/
	rly[9] = Ylm::ylmcoef[7]*zdr*rly[4]-Ylm::ylmcoef[8]*rly[1]; //l=3, m=0

	double tmp3 = Ylm::ylmcoef[9]*zdr;
	rly[10] = tmp3*rly[5]-Ylm::ylmcoef[10]*rly[2];//l=3,m=1
	rly[11] = tmp3*rly[6]-Ylm::ylmcoef[10]*rly[3];//l=3,m=-1

	double tmp4 = Ylm::ylmcoef[11]*zdr;
	rly[12] = tmp4*rly[7];//l=3,m=2
	rly[13] = tmp4*rly[8];//l=3,m=-2

	double tmp5 = Ylm::ylmcoef[14]*xdr;
	rly[14] = Ylm::ylmcoef[12]*rly[10]-Ylm::ylmcoef[13]*rly[2]-tmp5*rly[7];//l=3,m=3
	rly[15] = Ylm::ylmcoef[12]*rly[11]-Ylm::ylmcoef[13]*rly[3]-tmp5*rly[8];//l=3,m=-3
	if (Lmax == 3) return;

	/***************************
			 L = 4
	***************************/
	rly[16] = Ylm::ylmcoef[15]*zdr*rly[9]-Ylm::ylmcoef[16]*rly[4];//l=4,m=0

	double tmp6 = Ylm::ylmcoef[17]*zdr;
	rly[17] = tmp6*rly[10]-Ylm::ylmcoef[18]*rly[5];//l=4,m=1
	rly[18] = tmp6*rly[11]-Ylm::ylmcoef[18]*rly[6];//l=4,m=-1

	double tmp7 = Ylm::ylmcoef[19]*zdr;
	rly[19] = tmp7*rly[12]-Ylm::ylmcoef[20]*rly[7];//l=4,m=2
	rly[20] = tmp7*rly[13]-Ylm::ylmcoef[20]*rly[8];//l=4,m=-2

	double tmp8 = 3.0*zdr;
	rly[21] = tmp8*rly[14];//l=4,m=3
	rly[22] = tmp8*rly[15];//l=4,m=-3

	double tmp9 = Ylm::ylmcoef[23]*xdr;
	rly[23] = Ylm::ylmcoef[21]*rly[19]-Ylm::ylmcoef[22]*rly[7]-tmp9*rly[14];//l=4,m=4
	rly[24] = Ylm::ylmcoef[21]*rly[20]-Ylm::ylmcoef[22]*rly[8]-tmp9*rly[15];//l=4,m=-4
	if (Lmax == 4) return;

	/***************************
			 L = 5
	***************************/
	rly[25] = Ylm::ylmcoef[24]*zdr*rly[16]-Ylm::ylmcoef[25]*rly[9];//l=5,m=0

	double tmp10 = Ylm::ylmcoef[26]*zdr;
	rly[26] = tmp10*rly[17]-Ylm::ylmcoef[27]*rly[10];//l=5,m=1
	rly[27] = tmp10*rly[18]-Ylm::ylmcoef[27]*rly[11];//l=5,m=-1

	double tmp11 = Ylm::ylmcoef[28]*zdr;
	rly[28] = tmp11*rly[19]-Ylm::ylmcoef[29]*rly[12];//l=5,m=2
	rly[29] = tmp11*rly[20]-Ylm::ylmcoef[29]*rly[13];//l=5,m=-2

	double tmp12 = Ylm::ylmcoef[30]*zdr;
	rly[30] = tmp12*rly[21]-Ylm::ylmcoef[31]*rly[14];//l=5,m=3
	rly[31] = tmp12*rly[22]-Ylm::ylmcoef[31]*rly[15];//l=5,m=-3

	double tmp13 = Ylm::ylmcoef[32]*zdr;
	rly[32] = tmp13*rly[23];//l=5,m=4
	rly[33] = tmp13*rly[24];//l=5,m=-4

	double tmp14 = Ylm::ylmcoef[35]*xdr;
	rly[34] = Ylm::ylmcoef[33]*rly[30]-Ylm::ylmcoef[34]*rly[14]-tmp14*rly[23];//l=5,m=5
	rly[35] = Ylm::ylmcoef[33]*rly[31]-Ylm::ylmcoef[34]*rly[15]-tmp14*rly[24];//l=5,m=-5
	if (Lmax == 5) return;

	//if Lmax > 5
	for (int il = 6; il <= Lmax; il++)
	{
		int istart = il*il;
		int istart1 = (il-1)*(il-1);
		int istart2 = (il-2)*(il-2);

		double fac2 = sqrt(4.0*istart-1.0);
		double fac4 = sqrt(4.0*istart1-1.0);

		for (int im = 0; im < 2*il-1; im++)
		{
			int imm = (im+1)/2;
//			if (im % 2 == 0) imm *= -1;

			rly[istart+im] = fac2/sqrt((double)istart-imm*imm)*
								(zdr*rly[istart1+im] - sqrt((double)istart1-imm*imm)/fac4*rly[istart2+im]);
		}

		double bl1 = sqrt(2.0*il/(2.0*il+1.0));
		double bl2 = sqrt((2.0*il-2.0)/(2.0*il-1.0));
		double bl3 = sqrt(2.0)/fac2;

		rly[istart+2*il-1] = (bl3*rly[istart+2*il-5]-bl2*rly[istart2+2*il-5]-2.0*xdr*rly[istart1+2*il-3]) / bl1;
		rly[istart+2*il] = (bl3*rly[istart+2*il-4]-bl2*rly[istart2+2*il-4]-2.0*xdr*rly[istart1+2*il-2]) / bl1;
	}


	return;
}

// Peize Lin change rly 2016-08-26
void Ylm::rl_sph_harm
(
 	const int& Lmax, //max momentum of L
 	const double& x,
	const double& y,
	const double& z,
	std::vector<double>& rly
)
{
	rly.resize( (Lmax+1)*(Lmax+1) );

	double radius2 = x*x+y*y+z*z;

	//begin calculation
	/***************************
			 L = 0
	***************************/
	rly[0] = Ylm::ylmcoef[0]; //l=0, m=0
	if (Lmax == 0) return;

	/***************************
			 L = 1
	***************************/
	rly[1] = Ylm::ylmcoef[1]*z; //l=1, m=0
	rly[2] = -Ylm::ylmcoef[1]*x; //l=1, m=1
	rly[3] = -Ylm::ylmcoef[1]*y; //l=1, m=-1
	if (Lmax == 1) return;

	/***************************
			 L = 2
	***************************/
	rly[4] = Ylm::ylmcoef[2]*z*rly[1]-Ylm::ylmcoef[3]*rly[0]*radius2;//l=2, m=0

	double tmp0 = Ylm::ylmcoef[4]*z;
	rly[5] = tmp0*rly[2];//l=2,m=1
	rly[6] = tmp0*rly[3];//l=2,m=-1

	double tmp2 = Ylm::ylmcoef[4]*x;
	rly[7]= Ylm::ylmcoef[5]*rly[4]-Ylm::ylmcoef[6]*rly[0]*radius2 - tmp2*rly[2];//l=2,m=2
	rly[8] = -tmp2*rly[3];
//	rly[8] = tmp1+tmp2*rly[3];//l=2,m=-2
	if (Lmax == 2) return;

	/***************************
			 L = 3
	***************************/
	rly[9] = Ylm::ylmcoef[7]*z*rly[4]-Ylm::ylmcoef[8]*rly[1]*radius2; //l=3, m=0

	double tmp3 = Ylm::ylmcoef[9]*z;
	rly[10] = tmp3*rly[5]-Ylm::ylmcoef[10]*rly[2]*radius2;//l=3,m=1
	rly[11] = tmp3*rly[6]-Ylm::ylmcoef[10]*rly[3]*radius2;//l=3,m=-1

	double tmp4 = Ylm::ylmcoef[11]*z;
	rly[12] = tmp4*rly[7];//l=3,m=2
	rly[13] = tmp4*rly[8];//l=3,m=-2

	double tmp5 = Ylm::ylmcoef[14]*x;
	rly[14] = Ylm::ylmcoef[12]*rly[10]-Ylm::ylmcoef[13]*rly[2]*radius2-tmp5*rly[7];//l=3,m=3
	rly[15] = Ylm::ylmcoef[12]*rly[11]-Ylm::ylmcoef[13]*rly[3]*radius2-tmp5*rly[8];//l=3,m=-3
	if (Lmax == 3) return;

	/***************************
			 L = 4
	***************************/
	rly[16] = Ylm::ylmcoef[15]*z*rly[9]-Ylm::ylmcoef[16]*rly[4]*radius2;//l=4,m=0

	double tmp6 = Ylm::ylmcoef[17]*z;
	rly[17] = tmp6*rly[10]-Ylm::ylmcoef[18]*rly[5]*radius2;//l=4,m=1
	rly[18] = tmp6*rly[11]-Ylm::ylmcoef[18]*rly[6]*radius2;//l=4,m=-1

	double tmp7 = Ylm::ylmcoef[19]*z;
	rly[19] = tmp7*rly[12]-Ylm::ylmcoef[20]*rly[7]*radius2;//l=4,m=2
	rly[20] = tmp7*rly[13]-Ylm::ylmcoef[20]*rly[8]*radius2;//l=4,m=-2

	double tmp8 = 3.0*z;
	rly[21] = tmp8*rly[14];//l=4,m=3
	rly[22] = tmp8*rly[15];//l=4,m=-3

	double tmp9 = Ylm::ylmcoef[23]*x;
	rly[23] = Ylm::ylmcoef[21]*rly[19]-Ylm::ylmcoef[22]*rly[7]*radius2-tmp9*rly[14];//l=4,m=4
	rly[24] = Ylm::ylmcoef[21]*rly[20]-Ylm::ylmcoef[22]*rly[8]*radius2-tmp9*rly[15];//l=4,m=-4
	if (Lmax == 4) return;

	/***************************
			 L = 5
	***************************/
	rly[25] = Ylm::ylmcoef[24]*z*rly[16]-Ylm::ylmcoef[25]*rly[9]*radius2;//l=5,m=0

	double tmp10 = Ylm::ylmcoef[26]*z;
	rly[26] = tmp10*rly[17]-Ylm::ylmcoef[27]*rly[10]*radius2;//l=5,m=1
	rly[27] = tmp10*rly[18]-Ylm::ylmcoef[27]*rly[11]*radius2;//l=5,m=-1

	double tmp11 = Ylm::ylmcoef[28]*z;
	rly[28] = tmp11*rly[19]-Ylm::ylmcoef[29]*rly[12]*radius2;//l=5,m=2
	rly[29] = tmp11*rly[20]-Ylm::ylmcoef[29]*rly[13]*radius2;//l=5,m=-2

	double tmp12 = Ylm::ylmcoef[30]*z;
	rly[30] = tmp12*rly[21]-Ylm::ylmcoef[31]*rly[14]*radius2;//l=5,m=3
	rly[31] = tmp12*rly[22]-Ylm::ylmcoef[31]*rly[15]*radius2;//l=5,m=-3

	double tmp13 = Ylm::ylmcoef[32]*z;
	rly[32] = tmp13*rly[23];//l=5,m=4
	rly[33] = tmp13*rly[24];//l=5,m=-4

	double tmp14 = Ylm::ylmcoef[35]*x;
	rly[34] = Ylm::ylmcoef[33]*rly[30]-Ylm::ylmcoef[34]*rly[14]*radius2-tmp14*rly[23];//l=5,m=5
	rly[35] = Ylm::ylmcoef[33]*rly[31]-Ylm::ylmcoef[34]*rly[15]*radius2-tmp14*rly[24];//l=5,m=-5
	if (Lmax == 5) return;

	//if Lmax > 5
	for (int il = 6; il <= Lmax; il++)
	{
		int istart = il*il;
		int istart1 = (il-1)*(il-1);
		int istart2 = (il-2)*(il-2);

		double fac2 = sqrt(4.0*istart-1);
		double fac4 = sqrt(4.0*istart1-1);

		for (int im = 0; im < 2*il-1; im++)
		{
			int imm = (im+1)/2;
//			if (im % 2 == 0) imm *= -1;

			rly[istart+im] = fac2/sqrt((double)istart-imm*imm)*
								(z*rly[istart1+im] - sqrt((double)istart1-imm*imm)/fac4*rly[istart2+im]*radius2);
		}

		double bl1 = sqrt(2.0*il/(2.0*il+1.0));
		double bl2 = sqrt((2.0*il-2.0)/(2.0*il-1.0));
		double bl3 = sqrt(2.0)/fac2;

		rly[istart+2*il-1] = (bl3*rly[istart+2*il-5]-bl2*rly[istart2+2*il-5]*radius2-2.0*x*rly[istart1+2*il-3]) / bl1;
		rly[istart+2*il] = (bl3*rly[istart+2*il-4]-bl2*rly[istart2+2*il-4]*radius2-2.0*x*rly[istart1+2*il-2]) / bl1;
	}

	return;
}

void Ylm::grad_rl_sph_harm
(
 	const int& Lmax, //max momentum of L
 	const double& x,
	const double& y,
	const double& z,
	std::vector<double>& rly,
	std::vector<std::vector<double>>& grly
)
{
	rly.resize( (Lmax+1)*(Lmax+1) );
	grly.resize( (Lmax+1)*(Lmax+1), std::vector<double>(3) );

	double radius2 = x*x+y*y+z*z;
	double tx = 2.0*x;
	double ty = 2.0*y;
	double tz = 2.0*z;

	//begin calculation
	/***************************
			 L = 0
	***************************/
	rly[0] = Ylm::ylmcoef[0]; //l=0, m=0
	grly[0][0] = grly[0][1] = grly[0][2] = 0.0;
	if (Lmax == 0) return;

	/***************************
			 L = 1
	***************************/
	rly[1] = Ylm::ylmcoef[1]*z; //l=1, m=0
	grly[1][0] = grly[1][1] = 0.0;
	grly[1][2] = Ylm::ylmcoef[1];

	rly[2] = -Ylm::ylmcoef[1]*x; //l=1, m=1
	grly[2][1] = grly[2][2] = 0.0;
	grly[2][0] = -Ylm::ylmcoef[1];

	rly[3] = -Ylm::ylmcoef[1]*y; //l=1, m=-1
	grly[3][0] = grly[3][2] = 0.0;
	grly[3][1] = -Ylm::ylmcoef[1];

	if (Lmax == 1) return;

	/***************************
			 L = 2
	***************************/
	rly[4] = Ylm::ylmcoef[2]*z*rly[1]-Ylm::ylmcoef[3]*rly[0]*radius2;//l=2, m=0
	grly[4][0] = Ylm::ylmcoef[2]*z*grly[1][0]-Ylm::ylmcoef[3]*(grly[0][0]*radius2+rly[0]*tx);//l=2, m=0
	grly[4][1] = Ylm::ylmcoef[2]*z*grly[1][1]-Ylm::ylmcoef[3]*(grly[0][1]*radius2+rly[0]*ty);//l=2, m=0
	grly[4][2] = Ylm::ylmcoef[2]*(z*grly[1][2]+rly[1])-Ylm::ylmcoef[3]*(grly[0][2]*radius2+rly[0]*tz);//l=2, m=0


	double tmp0 = Ylm::ylmcoef[4]*z;
	rly[5] = tmp0*rly[2];//l=2,m=1
	grly[5][0] = tmp0*grly[2][0];
	grly[5][1] = tmp0*grly[2][1];
	grly[5][2] = Ylm::ylmcoef[4]*(rly[2]+z*grly[2][2]);

	rly[6] = tmp0*rly[3];//l=2,m=-1
	grly[6][0] = tmp0*grly[3][0];
	grly[6][1] = tmp0*grly[3][1];
	grly[6][2] = Ylm::ylmcoef[4]*(rly[3]+z*grly[3][2]);

	double tmp2 = Ylm::ylmcoef[4]*x;
	rly[7]= Ylm::ylmcoef[5]*rly[4]-Ylm::ylmcoef[6]*rly[0]*radius2 - tmp2*rly[2];//l=2,m=2
	grly[7][0] = Ylm::ylmcoef[5]*grly[4][0]-Ylm::ylmcoef[6]*(rly[0]*tx+grly[0][0]*radius2)-Ylm::ylmcoef[4]*(x*grly[2][0]+rly[2]);

//	std::cout << "\np1 = "<< Ylm::ylmcoef[5]*grly[4][0] << " p2 = " << -Ylm::ylmcoef[6]*rly[0]*tx
//						<< " p3 = " << -Ylm::ylmcoef[4]*x*grly[2][0] << " p4 = " << -Ylm::ylmcoef[4]*rly[2] << std::endl;

	grly[7][1] = Ylm::ylmcoef[5]*grly[4][1]-Ylm::ylmcoef[6]*(rly[0]*ty+grly[0][1]*radius2)-tmp2*grly[2][1];
	grly[7][2] = Ylm::ylmcoef[5]*grly[4][2]-Ylm::ylmcoef[6]*(rly[0]*tz+grly[0][2]*radius2)-tmp2*grly[2][2];

	rly[8] = -tmp2*rly[3];
	grly[8][0] = -Ylm::ylmcoef[4]*(rly[3]+x*grly[3][0]);
	grly[8][1] = -tmp2*grly[3][1];
	grly[8][2] = -tmp2*grly[3][2];
//	rly[8] = tmp1+tmp2*rly[3];//l=2,m=-2
	if (Lmax == 2) return;

	/***************************
			 L = 3
	***************************/
	rly[9] = Ylm::ylmcoef[7]*z*rly[4]-Ylm::ylmcoef[8]*rly[1]*radius2; //l=3, m=0
	grly[9][0] = Ylm::ylmcoef[7]*z*grly[4][0]-Ylm::ylmcoef[8]*(rly[1]*tx+grly[1][0]*radius2);
	grly[9][1] = Ylm::ylmcoef[7]*z*grly[4][1]-Ylm::ylmcoef[8]*(rly[1]*ty+grly[1][1]*radius2);
	grly[9][2] = Ylm::ylmcoef[7]*(rly[4]+z*grly[4][2])-Ylm::ylmcoef[8]*(rly[1]*tz+grly[1][2]*radius2);

	double tmp3 = Ylm::ylmcoef[9]*z;
	rly[10] = tmp3*rly[5]-Ylm::ylmcoef[10]*rly[2]*radius2;//l=3,m=1
	grly[10][0] = tmp3*grly[5][0]-Ylm::ylmcoef[10]*(grly[2][0]*radius2+rly[2]*tx);
	grly[10][1] = tmp3*grly[5][1]-Ylm::ylmcoef[10]*(grly[2][1]*radius2+rly[2]*ty);
	grly[10][2] = Ylm::ylmcoef[9]*(z*grly[5][2]+rly[5])-Ylm::ylmcoef[10]*(grly[2][2]*radius2+rly[2]*tz);

	rly[11] = tmp3*rly[6]-Ylm::ylmcoef[10]*rly[3]*radius2;//l=3,m=-1
	grly[11][0] = tmp3*grly[6][0]-Ylm::ylmcoef[10]*(grly[3][0]*radius2+rly[3]*tx);
	grly[11][1] = tmp3*grly[6][1]-Ylm::ylmcoef[10]*(grly[3][1]*radius2+rly[3]*ty);
	grly[11][2] = Ylm::ylmcoef[9]*(z*grly[6][2]+rly[6])-Ylm::ylmcoef[10]*(grly[3][2]*radius2+rly[3]*tz);

	double tmp4 = Ylm::ylmcoef[11]*z;
	rly[12] = tmp4*rly[7];//l=3,m=2
	grly[12][0] = tmp4*grly[7][0];
	grly[12][1] = tmp4*grly[7][1];
	grly[12][2] = Ylm::ylmcoef[11]*(z*grly[7][2]+rly[7]);

	rly[13] = tmp4*rly[8];//l=3,m=-2
	grly[13][0] = tmp4*grly[8][0];
	grly[13][1] = tmp4*grly[8][1];
	grly[13][2] = Ylm::ylmcoef[11]*(z*grly[8][2]+rly[8]);

	double tmp5 = Ylm::ylmcoef[14]*x;
	rly[14] = Ylm::ylmcoef[12]*rly[10]-Ylm::ylmcoef[13]*rly[2]*radius2-tmp5*rly[7];//l=3,m=3
	grly[14][0] = Ylm::ylmcoef[12]*grly[10][0]-Ylm::ylmcoef[13]*(rly[2]*tx+grly[2][0]*radius2)-Ylm::ylmcoef[14]*(rly[7]+x*grly[7][0]);
	grly[14][1] = Ylm::ylmcoef[12]*grly[10][1]-Ylm::ylmcoef[13]*(rly[2]*ty+grly[2][1]*radius2)-tmp5*grly[7][1];
	grly[14][2] = Ylm::ylmcoef[12]*grly[10][2]-Ylm::ylmcoef[13]*(rly[2]*tz+grly[2][2]*radius2)-tmp5*grly[7][2];

	rly[15] = Ylm::ylmcoef[12]*rly[11]-Ylm::ylmcoef[13]*rly[3]*radius2-tmp5*rly[8];//l=3,m=-3
	grly[15][0] = Ylm::ylmcoef[12]*grly[11][0]-Ylm::ylmcoef[13]*(rly[3]*tx+grly[3][0]*radius2)-Ylm::ylmcoef[14]*(rly[8]+x*grly[8][0]);
	grly[15][1] = Ylm::ylmcoef[12]*grly[11][1]-Ylm::ylmcoef[13]*(rly[3]*ty+grly[3][1]*radius2)-tmp5*grly[8][1];
	grly[15][2] = Ylm::ylmcoef[12]*grly[11][2]-Ylm::ylmcoef[13]*(rly[3]*tz+grly[3][2]*radius2)-tmp5*grly[8][2];
	if (Lmax == 3) return;

	/***************************
			 L = 4
	***************************/
	rly[16] = Ylm::ylmcoef[15]*z*rly[9]-Ylm::ylmcoef[16]*rly[4]*radius2;//l=4,m=0
	grly[16][0] = Ylm::ylmcoef[15]*z*grly[9][0]-Ylm::ylmcoef[16]*(rly[4]*tx+grly[4][0]*radius2);
	grly[16][1] = Ylm::ylmcoef[15]*z*grly[9][1]-Ylm::ylmcoef[16]*(rly[4]*ty+grly[4][1]*radius2);
	grly[16][2] = Ylm::ylmcoef[15]*(z*grly[9][2]+rly[9])-Ylm::ylmcoef[16]*(rly[4]*tz+grly[4][2]*radius2);

	double tmp6 = Ylm::ylmcoef[17]*z;
	rly[17] = tmp6*rly[10]-Ylm::ylmcoef[18]*rly[5]*radius2;//l=4,m=1
	grly[17][0] = tmp6*grly[10][0]-Ylm::ylmcoef[18]*(rly[5]*tx+grly[5][0]*radius2);
	grly[17][1] = tmp6*grly[10][1]-Ylm::ylmcoef[18]*(rly[5]*ty+grly[5][1]*radius2);
	grly[17][2] = Ylm::ylmcoef[17]*(z*grly[10][2]+rly[10])-Ylm::ylmcoef[18]*(rly[5]*tz+grly[5][2]*radius2);

	rly[18] = tmp6*rly[11]-Ylm::ylmcoef[18]*rly[6]*radius2;//l=4,m=-1
	grly[18][0] = tmp6*grly[11][0]-Ylm::ylmcoef[18]*(rly[6]*tx+grly[6][0]*radius2);
	grly[18][1] = tmp6*grly[11][1]-Ylm::ylmcoef[18]*(rly[6]*ty+grly[6][1]*radius2);
	grly[18][2] = Ylm::ylmcoef[17]*(z*grly[11][2]+rly[11])-Ylm::ylmcoef[18]*(rly[6]*tz+grly[6][2]*radius2);

	double tmp7 = Ylm::ylmcoef[19]*z;
	rly[19] = tmp7*rly[12]-Ylm::ylmcoef[20]*rly[7]*radius2;//l=4,m=2
	grly[19][0] = tmp7*grly[12][0]-Ylm::ylmcoef[20]*(rly[7]*tx+grly[7][0]*radius2);
	grly[19][1] = tmp7*grly[12][1]-Ylm::ylmcoef[20]*(rly[7]*ty+grly[7][1]*radius2);
	grly[19][2] = Ylm::ylmcoef[19]*(z*grly[12][2]+rly[12])-Ylm::ylmcoef[20]*(rly[7]*tz+grly[7][2]*radius2);

	rly[20] = tmp7*rly[13]-Ylm::ylmcoef[20]*rly[8]*radius2;//l=4,m=-2
	grly[20][0] = tmp7*grly[13][0]-Ylm::ylmcoef[20]*(rly[8]*tx+grly[8][0]*radius2);
	grly[20][1] = tmp7*grly[13][1]-Ylm::ylmcoef[20]*(rly[8]*ty+grly[8][1]*radius2);
	grly[20][2] = Ylm::ylmcoef[19]*(z*grly[13][2]+rly[13])-Ylm::ylmcoef[20]*(rly[8]*tz+grly[8][2]*radius2);

	double tmp8 = 3.0*z;
	rly[21] = tmp8*rly[14];//l=4,m=3
	grly[21][0] = tmp8*grly[14][0];
	grly[21][1] = tmp8*grly[14][1];
	grly[21][2] = 3.0*(z*grly[14][2]+rly[14]);


	rly[22] = tmp8*rly[15];//l=4,m=-3
	grly[22][0] = tmp8*grly[15][0];
	grly[22][1] = tmp8*grly[15][1];
	grly[22][2] = 3.0*(z*grly[15][2]+rly[15]);

	double tmp9 = Ylm::ylmcoef[23]*x;
	rly[23] = Ylm::ylmcoef[21]*rly[19]-Ylm::ylmcoef[22]*rly[7]*radius2-tmp9*rly[14];//l=4,m=4
	grly[23][0] = Ylm::ylmcoef[21]*grly[19][0]-Ylm::ylmcoef[22]*(rly[7]*tx+grly[7][0]*radius2)-Ylm::ylmcoef[23]*(x*grly[14][0]+rly[14]);
	grly[23][1] = Ylm::ylmcoef[21]*grly[19][1]-Ylm::ylmcoef[22]*(rly[7]*ty+grly[7][1]*radius2)-tmp9*grly[14][1];
	grly[23][2] = Ylm::ylmcoef[21]*grly[19][2]-Ylm::ylmcoef[22]*(rly[7]*tz+grly[7][2]*radius2)-tmp9*grly[14][2];

	rly[24] = Ylm::ylmcoef[21]*rly[20]-Ylm::ylmcoef[22]*rly[8]*radius2-tmp9*rly[15];//l=4,m=-4
	grly[24][0] = Ylm::ylmcoef[21]*grly[20][0]-Ylm::ylmcoef[22]*(rly[8]*tx+grly[8][0]*radius2)-Ylm::ylmcoef[23]*(x*grly[15][0]+rly[15]);
	grly[24][1] = Ylm::ylmcoef[21]*grly[20][1]-Ylm::ylmcoef[22]*(rly[8]*ty+grly[8][1]*radius2)-tmp9*grly[15][1];
	grly[24][2] = Ylm::ylmcoef[21]*grly[20][2]-Ylm::ylmcoef[22]*(rly[8]*tz+grly[8][2]*radius2)-tmp9*grly[15][2];

	if (Lmax == 4) return;

	/***************************
			 L = 5
	***************************/
	rly[25] = Ylm::ylmcoef[24]*z*rly[16]-Ylm::ylmcoef[25]*rly[9]*radius2;//l=5,m=0
	grly[25][0] = Ylm::ylmcoef[24]*z*grly[16][0]-Ylm::ylmcoef[25]*(rly[9]*tx+grly[9][0]*radius2);
	grly[25][1] = Ylm::ylmcoef[24]*z*grly[16][1]-Ylm::ylmcoef[25]*(rly[9]*ty+grly[9][1]*radius2);
	grly[25][2] = Ylm::ylmcoef[24]*(z*grly[16][2]+rly[16])-Ylm::ylmcoef[25]*(rly[9]*tz+grly[9][2]*radius2);

	double tmp10 = Ylm::ylmcoef[26]*z;
	rly[26] = tmp10*rly[17]-Ylm::ylmcoef[27]*rly[10]*radius2;//l=5,m=1
	grly[26][0] = tmp10*grly[17][0]-Ylm::ylmcoef[27]*(rly[10]*tx+grly[10][0]*radius2);
	grly[26][1] = tmp10*grly[17][1]-Ylm::ylmcoef[27]*(rly[10]*ty+grly[10][1]*radius2);
	grly[26][2] = Ylm::ylmcoef[26]*(z*grly[17][2]+rly[17])-Ylm::ylmcoef[27]*(rly[10]*tz+grly[10][2]*radius2);

	rly[27] = tmp10*rly[18]-Ylm::ylmcoef[27]*rly[11]*radius2;//l=5,m=-1
	grly[27][0] = tmp10*grly[18][0]-Ylm::ylmcoef[27]*(rly[11]*tx+grly[11][0]*radius2);
	grly[27][1] = tmp10*grly[18][1]-Ylm::ylmcoef[27]*(rly[11]*ty+grly[11][1]*radius2);
	grly[27][2] = Ylm::ylmcoef[26]*(z*grly[18][2]+rly[18])-Ylm::ylmcoef[27]*(rly[11]*tz+grly[11][2]*radius2);

	double tmp11 = Ylm::ylmcoef[28]*z;
	rly[28] = tmp11*rly[19]-Ylm::ylmcoef[29]*rly[12]*radius2;//l=5,m=2
	grly[28][0] = tmp11*grly[19][0]-Ylm::ylmcoef[29]*(rly[12]*tx+grly[12][0]*radius2);
	grly[28][1] = tmp11*grly[19][1]-Ylm::ylmcoef[29]*(rly[12]*ty+grly[12][1]*radius2);
	grly[28][2] = Ylm::ylmcoef[28]*(z*grly[19][2]+rly[19])-Ylm::ylmcoef[29]*(rly[12]*tz+grly[12][2]*radius2);

	rly[29] = tmp11*rly[20]-Ylm::ylmcoef[29]*rly[13]*radius2;//l=5,m=-2
	grly[29][0] = tmp11*grly[20][0]-Ylm::ylmcoef[29]*(rly[13]*tx+grly[13][0]*radius2);
	grly[29][1] = tmp11*grly[20][1]-Ylm::ylmcoef[29]*(rly[13]*ty+grly[13][1]*radius2);
	grly[29][2] = Ylm::ylmcoef[28]*(z*grly[20][2]+rly[20])-Ylm::ylmcoef[29]*(rly[13]*tz+grly[13][2]*radius2);

	double tmp12 = Ylm::ylmcoef[30]*z;
	rly[30] = tmp12*rly[21]-Ylm::ylmcoef[31]*rly[14]*radius2;//l=5,m=3
	grly[30][0] = tmp12*grly[21][0]-Ylm::ylmcoef[31]*(grly[14][0]*radius2+rly[14]*tx);
	grly[30][1] = tmp12*grly[21][1]-Ylm::ylmcoef[31]*(grly[14][1]*radius2+rly[14]*ty);
	grly[30][2] = Ylm::ylmcoef[30]*(z*grly[21][2]+rly[21])-Ylm::ylmcoef[31]*(grly[14][2]*radius2+rly[14]*tz);

	rly[31] = tmp12*rly[22]-Ylm::ylmcoef[31]*rly[15]*radius2;//l=5,m=-3
	grly[31][0] = tmp12*grly[22][0]-Ylm::ylmcoef[31]*(grly[15][0]*radius2+rly[15]*tx);
	grly[31][1] = tmp12*grly[22][1]-Ylm::ylmcoef[31]*(grly[15][1]*radius2+rly[15]*ty);
	grly[31][2] = Ylm::ylmcoef[30]*(z*grly[22][2]+rly[22])-Ylm::ylmcoef[31]*(grly[15][2]*radius2+rly[15]*tz);

	double tmp13 = Ylm::ylmcoef[32]*z;
	rly[32] = tmp13*rly[23];//l=5,m=4
	grly[32][0] = tmp13*grly[23][0];
	grly[32][1] = tmp13*grly[23][1];
	grly[32][2] = Ylm::ylmcoef[32]*(rly[23]+z*grly[23][2]);

	rly[33] = tmp13*rly[24];//l=5,m=-4
	grly[33][0] = tmp13*grly[24][0];
	grly[33][1] = tmp13*grly[24][1];
	grly[33][2] = Ylm::ylmcoef[32]*(rly[24]+z*grly[24][2]);

	double tmp14 = Ylm::ylmcoef[35]*x;
	rly[34] = Ylm::ylmcoef[33]*rly[30]-Ylm::ylmcoef[34]*rly[14]*radius2-tmp14*rly[23];//l=5,m=5
	grly[34][0] = Ylm::ylmcoef[33]*grly[30][0]-Ylm::ylmcoef[34]*(rly[14]*tx+grly[14][0]*radius2)-Ylm::ylmcoef[35]*(x*grly[23][0]+rly[23]);
	grly[34][1] = Ylm::ylmcoef[33]*grly[30][1]-Ylm::ylmcoef[34]*(rly[14]*ty+grly[14][1]*radius2)-tmp14*grly[23][1];
	grly[34][2] = Ylm::ylmcoef[33]*grly[30][2]-Ylm::ylmcoef[34]*(rly[14]*tz+grly[14][2]*radius2)-tmp14*grly[23][2];

	rly[35] = Ylm::ylmcoef[33]*rly[31]-Ylm::ylmcoef[34]*rly[15]*radius2-tmp14*rly[24];//l=5,m=-5
	grly[35][0] = Ylm::ylmcoef[33]*grly[31][0]-Ylm::ylmcoef[34]*(rly[15]*tx+grly[15][0]*radius2)-Ylm::ylmcoef[35]*(x*grly[24][0]+rly[24]);
	grly[35][1] = Ylm::ylmcoef[33]*grly[31][1]-Ylm::ylmcoef[34]*(rly[15]*ty+grly[15][1]*radius2)-tmp14*grly[24][1];
	grly[35][2] = Ylm::ylmcoef[33]*grly[31][2]-Ylm::ylmcoef[34]*(rly[15]*tz+grly[15][2]*radius2)-tmp14*grly[24][2];

	if (Lmax == 5) return;

	//if Lmax > 5
	for (int il = 6; il <= Lmax; il++)
	{
		int istart = il*il;
		int istart1 = (il-1)*(il-1);
		int istart2 = (il-2)*(il-2);

		double fac2 = sqrt(4.0*istart-1.0);
		double fac4 = sqrt(4.0*istart1-1.0);

		for (int im = 0; im < 2*il-1; im++)
		{
			int imm = (im+1)/2;
//			if (im % 2 == 0) imm *= -1;

			double var1 = fac2/sqrt((double)istart-imm*imm);
			double var2 = sqrt((double)istart1-imm*imm)/fac4;

			rly[istart+im] = var1*(z*rly[istart1+im] - var2*rly[istart2+im]*radius2);

			grly[istart+im][0]=var1*(z*grly[istart1+im][0]-var2*(rly[istart2+im]*tx+grly[istart2+im][0]*radius2));
			grly[istart+im][1]=var1*(z*grly[istart1+im][1]-var2*(rly[istart2+im]*ty+grly[istart2+im][1]*radius2));
			grly[istart+im][2]=var1*(z*grly[istart1+im][2]+rly[istart1+im]-var2*(rly[istart2+im]*tz+grly[istart2+im][2]*radius2));

		}

		double bl1 = sqrt(2.0*il/(2.0*il+1.0));
		double bl2 = sqrt((2.0*il-2.0)/(2.0*il-1.0));
		double bl3 = sqrt(2.0)/fac2;

		int id1 = istart+2*il-1;
		int id2 = istart+2*il-5;
		int id3 = istart2+2*il-5;
		int id4 = istart1+2*il-3;

		rly[id1] = (bl3*rly[id2]-bl2*rly[id3]*radius2-2.0*x*rly[id4]) / bl1;
		grly[id1][0] = (bl3*grly[id2][0]-bl2*(grly[id3][0]*radius2+rly[id3]*tx)-2.0*(rly[id4]+x*grly[id4][0]))/bl1;
		grly[id1][1] = (bl3*grly[id2][1]-bl2*(grly[id3][1]*radius2+rly[id3]*ty)-2.0*x*grly[id4][1])/bl1;
		grly[id1][2] = (bl3*grly[id2][2]-bl2*(grly[id3][2]*radius2+rly[id3]*tz)-2.0*x*grly[id4][2])/bl1;


		rly[id1+1] = (bl3*rly[id2+1]-bl2*rly[id3+1]*radius2-2.0*x*rly[id4+1]) / bl1;
		grly[id1+1][0] = (bl3*grly[id2+1][0]-bl2*(grly[id3+1][0]*radius2+rly[id3+1]*tx)-2.0*(rly[id4+1]+x*grly[id4+1][0]))/bl1;
		grly[id1+1][1] = (bl3*grly[id2+1][1]-bl2*(grly[id3+1][1]*radius2+rly[id3+1]*ty)-2.0*x*grly[id4+1][1])/bl1;
		grly[id1+1][2] = (bl3*grly[id2+1][2]-bl2*(grly[id3+1][2]*radius2+rly[id3+1]*tz)-2.0*x*grly[id4+1][2])/bl1;
	}


	return;
}

void Ylm::hes_rl_sph_harm
(
 	const int& Lmax, //max momentum of L
 	const double& x,
	const double& y,
	const double& z,
	std::vector<std::vector<double>>& hrly
)
{
	hrly.resize( (Lmax+1)*(Lmax+1), std::vector<double>(6) );

	double radius2 = x*x+y*y+z*z;
	double coeff;

	//begin calculation
	/***************************
			 L = 0
	***************************/
	hrly[0][0] = hrly[0][1] = hrly[0][2] = 0.0;
	hrly[0][3] = hrly[0][4] = hrly[0][5] = 0.0;
	if (Lmax == 0) return;

	/***************************
			 L = 1
	***************************/
	hrly[1][0] = hrly[1][1] = hrly[1][2] = 0.0;
	hrly[1][3] = hrly[1][4] = hrly[1][5] = 0.0;

	hrly[2][0] = hrly[2][1] = hrly[2][2] = 0.0;
	hrly[2][3] = hrly[2][4] = hrly[2][5] = 0.0;

	hrly[3][0] = hrly[3][1] = hrly[3][2] = 0.0;
	hrly[3][3] = hrly[3][4] = hrly[3][5] = 0.0;

	if (Lmax == 1) return;

	/***************************
			 L = 2
	***************************/
	//m=0 : 3z^2-r^2
	coeff = sqrt(5.0 / ModuleBase::PI) / 4.0;
	hrly[4][0] = hrly[4][3] = -2.0 * coeff;
	hrly[4][5] = 4.0 * coeff;
	hrly[4][1] = hrly[4][2] = hrly[4][4] = 0.0;

	//m=1 : xz
	coeff = sqrt(15.0 / ModuleBase::PI) / 2.0;
	hrly[5][2] = coeff;
	hrly[5][0] = hrly[5][1] = 0.0;
	hrly[5][3] = hrly[5][4] = hrly[5][5] = 0.0;

	//m=-1 : yz
	hrly[6][4] = coeff;
	hrly[6][0] = hrly[6][1] = 0.0;
	hrly[6][2] = hrly[6][3] = hrly[6][5] = 0.0;

	//m=-2 : xy
	hrly[8][1] = coeff;
	hrly[8][0] = hrly[8][2] = 0.0;
	hrly[8][3] = hrly[8][4] = hrly[8][5] = 0.0;

	//m=2 : (x^2-y^2)
	coeff = sqrt(15.0 / ModuleBase::PI) / 4.0;
	hrly[7][0] =  2.0 * coeff;
	hrly[7][3] = -2.0 * coeff;
	hrly[7][1] = hrly[7][2] = 0.0;
	hrly[7][4] = hrly[7][5] = 0.0;

	if (Lmax == 2) return;

	/***************************
			 L = 3
	***************************/
	//m=0 : (5z^3-3zr^2)
	coeff = sqrt(7.0 / ModuleBase::PI) / 4.0;
	hrly[9][0] = hrly[9][3] = -6.0 * z * coeff;
	hrly[9][1] =  0.0;
	hrly[9][2] = -6.0 * x * coeff;
	hrly[9][4] = -6.0 * y * coeff;
	hrly[9][5] = 12.0 * z * coeff;

	//m=1 : x(5z^2-r^2)
	coeff = sqrt(21.0 / 2.0 / ModuleBase::PI) / 4.0;
	hrly[10][0] = -6.0 * x * coeff;
	hrly[10][1] = -2.0 * y * coeff;
	hrly[10][2] =  8.0 * z * coeff;
	hrly[10][3] = -2.0 * x * coeff;
	hrly[10][4] =  0.0;
	hrly[10][5] =  8.0 * x * coeff;

	//m=-1 : y(5z^2-r^2)
	hrly[11][0] = -2.0 * y * coeff;
	hrly[11][1] = -2.0 * x * coeff;
	hrly[11][2] =  0.0;
	hrly[11][3] = -6.0 * y * coeff;
	hrly[11][4] =  8.0 * z * coeff;
	hrly[11][5] =  8.0 * y * coeff;

	//m=2 : (x^2-y^2)z
	coeff = sqrt(105.0 / ModuleBase::PI) / 4.0;
	hrly[12][0] =  2.0 * z * coeff;
	hrly[12][1] =  0.0;
	hrly[12][2] =  2.0 * x * coeff;
	hrly[12][3] = -hrly[12][0];
	hrly[12][4] = -2.0 * y * coeff;
	hrly[12][5] =  0.0;

	//m=-2 : xyz
	coeff = sqrt(105.0 / ModuleBase::PI) / 2.0;
	hrly[13][0] = 0.0;
	hrly[13][1] = z * coeff;
	hrly[13][2] = y * coeff;
	hrly[13][3] = 0.0;
	hrly[13][4] = x * coeff;
	hrly[13][5] = 0.0;

	//m=3 : x(x^2-3y^2)
	coeff = sqrt(35.0 / 2.0 / ModuleBase::PI) / 4.0;
	hrly[14][0] =  6.0 * x * coeff;
	hrly[14][1] = -6.0 * y * coeff;
	hrly[14][2] = 0.0;
	hrly[14][3] = -hrly[14][0];
	hrly[14][4] = 0.0;
	hrly[14][5] = 0.0;

	//m=-3 : y(3x^2-y^2)
	hrly[15][0] =  6.0 * y * coeff;
	hrly[15][1] =  6.0 * x * coeff;
	hrly[15][2] =  0.0;
	hrly[15][3] = -hrly[15][0];
	hrly[15][4] =  0.0;
	hrly[15][5] =  0.0;

	if (Lmax == 3) return;

	/***************************
			 L = 4
	***************************/
	//m=0 : (35z^4 - 30z^2r^2 + 3r^4)
	coeff = sqrt(1.0 / ModuleBase::PI) * 3.0 / 16.0;
	hrly[16][0] =  12.0 * (3.0 * x*x + y*y - 4.0 * z*z) * coeff;
	hrly[16][1] =  24.0 * x * y * coeff;
	hrly[16][2] = -96.0 * x * z * coeff;
	hrly[16][3] =  12.0 * (x*x + 3.0 * y*y - 4.0 * z*z) * coeff;
	hrly[16][4] = -96.0 * y * z * coeff;
	hrly[16][5] = -48.0 * (x*x + y*y -2.0 * z*z) * coeff;

	//m=1 : x(7z^3 - 3zr^2)
	coeff = 3.0 / 4.0 * sqrt(5.0 / 2.0 / ModuleBase::PI);
	hrly[17][0] = -18.0 * x * z * coeff;
	hrly[17][1] =  -6.0 * y * z * coeff;
	hrly[17][2] =  -3.0 * (3.0 * x*x + y*y  - 4.0 * z*z) * coeff;
	hrly[17][3] =  -6.0 * x * z * coeff;
	hrly[17][4] =  -6.0 * x * y * coeff;
	hrly[17][5] =  24.0 * x * z * coeff;

	//m=-1 : y(7z^3 - 3zr^2)
	hrly[18][0] =  -6.0 * y * z * coeff;
	hrly[18][1] =  -6.0 * x * z * coeff;
	hrly[18][2] =  -6.0 * x * y * coeff;
	hrly[18][3] = -18.0 * y * z * coeff;
	hrly[18][4] =  -3.0 * (x*x + 3.0 * y*y  - 4.0 * z*z) * coeff;
	hrly[18][5] =  24.0 * y * z * coeff;

	//m=2 : (x^2 - y^2)(7z^2 - r^2)
	coeff = 3.0 / 8.0 * sqrt(5.0 / ModuleBase::PI);
	hrly[19][0] = -12.0 * (x*x - z*z) * coeff;
	hrly[19][1] =  0.0;
	hrly[19][2] =  24.0 * x * z * coeff;
	hrly[19][3] =  12.0 * (y*y - z*z) * coeff;
	hrly[19][4] = -24.0 * y * z * coeff;
	hrly[19][5] =  12.0 * (x*x - y*y) * coeff;

	//m=-2 : xy(7z^2 - r^2)
	coeff = 3.0 / 4.0 * sqrt(5.0 / ModuleBase::PI);
	hrly[20][0] = -6.0 * x * y * coeff;
	hrly[20][1] = -3.0 * (x*x + y*y - 2.0 * z*z) * coeff;
	hrly[20][2] =  2.0 * y * z * coeff;
	hrly[20][3] =  hrly[20][0];
	hrly[20][4] = 12.0 * x * z * coeff;
	hrly[20][5] = 12.0 * x * y * coeff;

	//m=3 : x(x^2-3y^2)z
	coeff = 3.0 / 4.0 * sqrt(35.0 / 2.0 / ModuleBase::PI);
	hrly[21][0] =  6.0 * x * z * coeff;
	hrly[21][1] = -6.0 * y * z * coeff;
	hrly[21][2] =  3.0 * (x*x - y*y) * coeff;
	hrly[21][3] = -6.0 * x * z * coeff;
	hrly[21][4] = -6.0 * x * y * coeff;
	hrly[21][5] =  0.0;

	//m=-3 : y(3x^2-y^2)z
	hrly[22][0] =  6.0 * y * z * coeff;
	hrly[22][1] =  6.0 * x * z * coeff;
	hrly[22][2] =  6.0 * x * y * coeff;
	hrly[22][3] = -6.0 * y * z * coeff;
	hrly[22][4] =  3.0 * (x*x - y*y) * coeff;
	hrly[22][5] =  0.0;

	//m=4 : x^4 + y^4 - 6 x^2y^2
	coeff = 3.0 / 16.0 * sqrt(35.0 / ModuleBase::PI);
	hrly[23][0] =  12.0 * (x*x - y*y) * coeff;
	hrly[23][1] = -24.0 * x * y * coeff;
	hrly[23][2] =   0.0;
	hrly[23][3] =  -hrly[23][0];
	hrly[23][4] =   0.0;
	hrly[23][5] =   0.0;

	//m=-4 : xy(x^2 - y^2)
	coeff = 3.0 / 4.0 * sqrt(35.0 / ModuleBase::PI);
	hrly[24][0] =  6.0 * x * y * coeff;
	hrly[24][1] =  3.0 * (x*x - y*y) * coeff;
	hrly[24][2] =  0.0;
	hrly[24][3] = -hrly[24][0];
	hrly[24][4] =  0.0;
	hrly[24][5] =  0.0;

	if (Lmax == 4) return;

	/***************************
			 L > 4
	***************************/
	ModuleBase::WARNING_QUIT("hes_rl_sph_harm","l>4 not implemented!");


	return;
}

void Ylm::set_coefficients(void)
{
	Ylm::ylmcoef[0] = 1.0 / sqrt(ModuleBase::FOUR_PI);
	Ylm::ylmcoef[1] = sqrt (3.0 / ModuleBase::FOUR_PI);
	Ylm::ylmcoef[2] = sqrt (15.0) / 2.0;
	Ylm::ylmcoef[3] = sqrt (5.0) / 2.0;
	Ylm::ylmcoef[4] = sqrt (5.0);
	Ylm::ylmcoef[5] = 1.0 / sqrt(3.0);
	Ylm::ylmcoef[6] = sqrt (5.0 / 3.0);
	Ylm::ylmcoef[7] = sqrt (35.0 / 9.0);
	Ylm::ylmcoef[8] = sqrt (7.0/3.0)/1.5;
	Ylm::ylmcoef[9] = sqrt (35.0 / 8.0);
	Ylm::ylmcoef[10] = sqrt (7.0 / 8.0);
	Ylm::ylmcoef[11] = sqrt (7.0);
	Ylm::ylmcoef[12] = 1.0 / sqrt (15.0);
	Ylm::ylmcoef[13] = sqrt (14.0 / 15.0);
	Ylm::ylmcoef[14] = sqrt (14.0 / 3.0);
	Ylm::ylmcoef[15] = sqrt(7.0)*3.0/4.0;
	Ylm::ylmcoef[16] = 9.0/4.0/sqrt(5.0);
	Ylm::ylmcoef[17] = sqrt(21.0/5.0);
	Ylm::ylmcoef[18] = sqrt(24.0/25.0);
	Ylm::ylmcoef[19] = sqrt(21.0)/2.0;
	Ylm::ylmcoef[20] = sqrt(3.0)/2.0;
	Ylm::ylmcoef[21] = 0.5/sqrt(7.0);
	Ylm::ylmcoef[22] = 1.5*sqrt(3.0/7.0);
	Ylm::ylmcoef[23] = 3.0/sqrt(2.0);
	Ylm::ylmcoef[24] = 0.6*sqrt(11.0);
	Ylm::ylmcoef[25] = 0.8*sqrt(11.0/7.0);
	Ylm::ylmcoef[26] = sqrt (33.0/8.0);
	Ylm::ylmcoef[27] = sqrt (55.0/56.0);
	Ylm::ylmcoef[28] = sqrt (33.0/7.0);
	Ylm::ylmcoef[29] = sqrt (11.0)*2.0/7.0;
	Ylm::ylmcoef[30] = sqrt (11.0)*0.75;
	Ylm::ylmcoef[31] = sqrt (11.0)*0.25;
	Ylm::ylmcoef[32] = sqrt (11.0);
	Ylm::ylmcoef[33] = 1.0/3.0/sqrt(5.0);
	Ylm::ylmcoef[34] = 2.0/3.0*sqrt(11.0/5.0);
	Ylm::ylmcoef[35] = sqrt(22.0/5.0);
	return;
}

/*
void Ylm::test1 (void)
{
	ModuleBase::Vector3<double> R (20.0, 0.0, 0.0);
	double xdr = R.x/R.norm();
	double ydr = R.y/R.norm();
	double zdr = R.z/R.norm();
	const int L = 9;
	const double rl = std::pow( R.norm(), L);
	std::cout << " rl=" << rl << std::endl;
	Ylm::set_coefficients();

	int nu = 100;

	// Peize Lin change rlya 2016-08-26
	std::vector<double> rlya;
	double rlyb[400];
	ZEROS( rlyb, 400);

//	Ylm::sph_harm (9, xdr, ydr, zdr, rlya);
	Ylm::rl_sph_harm (L, xdr, ydr, zdr, rlya);
//	Ylm::rlylm (10, R.x, R.y, R.z, rlyb);
	Ylm::get_ylm_real (L+1, R, rlyb);

	for (int i=0; i < nu; i++)
	{
	//	std::cout << "\ni= " << i << " rlya = " << rlya[i] << " rlyb = " << rlyb[i] << std::endl;
		double diff = fabs(rlya[i]-rlyb[i]);
		if (diff > 1e-8)
		{
			std::cout << "Ylm::test1, error is too large!" << std::endl;
			//WARNING_QUIT ("Ylm::test1","error is too large!");
			exit(0);
		}
	}
	return;
}
*/
/*
void Ylm::test2 (void)
{
	ModuleBase::Vector3<double> R (0.1,-0.2,0.5);
	Ylm::set_coefficients();

	//int nu = 100;

	std::vector<double> rlya;
	double rlyb[400];

	std::vector<std::vector<double>> grlya;
	double grlyb[400][3];

	Ylm::grad_rl_sph_harm (9, R.x, R.y, R.z, rlya, grlya);
	Ylm::rlylm (10, R.x, R.y, R.z, rlyb, grlyb);

	for (int i = 0; i < 100; i++)
	{
		double diffx = fabs(grlya[i][2]-grlyb[i][2]);
		if (diffx > 1e-8)
		{
			std::cout << "Ylm::test2, Large error in Direv X!" << std::endl;
			//WARNING_QUIT ("Ylm::test2","Large error in Direv X!");
			exit(0);
		}
	}
	return;
}
*/

void Ylm::rlylm
(
 	const int& Lmax, //max momentum of l + 1
 	const double& x,
	const double& y,
	const double& z,
	double rly[],
	double grly[][3]
)
{
	int MaxL = Lmax - 1;

	assert(MaxL >= 0);

	//get xy_dependence
	assert(MaxL <= 19);

	double Am[20];
	double Bm[20];
	double Gx_Am[20];
	double Gx_Bm[20];
	double Gy_Am[20];
	double Gy_Bm[20];

	ZEROS(Am, 20);
	ZEROS(Bm, 20);
	ZEROS(Gx_Am, 20);
	ZEROS(Gy_Am, 20);

	double x2, x3, x4, x5;
	double y2, y3, y4, y5;

	x2 = x * x;
	x3 = x2 * x;
	x4 = x3 * x;
	x5 = x4 * x;

	y2 = y * y;
	y3 = y2 * y;
	y4 = y3 * y;
	y5 = y4 * y;

	//x-y dependence
	//Am
	//Bm
	for(int im = 0; im < MaxL+1; im++)
	{
		if(im == 0)
		{
			Am[0] = 1.0;
			Bm[0] = 0.0;

			Gx_Am[0] = 0.0;
			Gy_Am[0] = 0.0;

			Gx_Bm[0] = 0.0;
			Gy_Bm[0] = 0.0;
		}
		else if(im == 1)
		{
			Am[1] = x;
			Bm[1] = y;

			Gx_Am[1] = 1.0;
			Gy_Am[1] = 0.0;

			Gx_Bm[1] = 0.0;
			Gy_Bm[1] = 1.0;
		}
		else if(im == 2)
		{
			Am[2] = x2- y2;
			Bm[2] = 2.0 * x * y;

			Gx_Am[2] = 2.0 * x;
			Gy_Am[2] = -2.0 * y;

			Gx_Bm[2] = 2.0 * y;
			Gy_Bm[2] = 2.0 * x;
		}
		else if(im == 3)
		{
			Am[3] = x3 - 3.0 * x * y2;
			Bm[3] = 3.0 * x2 * y - y3;

			Gx_Am[3] = 3.0 * (x2 - y2);
			Gy_Am[3] = -6.0 * x * y;

			Gx_Bm[3] = 6.0 * x * y;
			Gy_Bm[3] = 3.0 * (x2 - y2);
		}
		else if(im == 4)
		{
			Am[4] = x4 - 6.0 * x2 * y2 + y4;
			Bm[4] = 4.0 * (x3 * y - x * y3);

			Gx_Am[4] = 4.0 * x3 - 12.0 * x * y2;
			Gy_Am[4] = -12.0 * x2 * y + 4.0 * y3;

			Gx_Bm[4] = 12.0 * x2 * y - 4.0 * y3;
			Gy_Bm[4] = 4.0 * x3 - 12.0 * x * y2;
		}
		else if(im == 5)
		{
			Am[5] = x5 - 10.0 * x3 * y2 + 5.0 * x * y4;
			Bm[5] = 5.0 * x4 * y - 10.0 * x2 * y3 + y5;

			Gx_Am[5] = 5.0 * x4 - 30.0 * x2 * y2 + 5.0 * y4;
			Gy_Am[5] = 20.0 * (x * y3 - x3 * y);

			Gx_Bm[5] = 20.0 * (x3 * y - x * y3);
			Gy_Bm[5] = 5.0 * x4 - 30.0 * x2 * y2 + 5.0 * y4;
		}
		else
		{
			for(int ip = 0; ip <= im; ip++)
			{
				double aux = Fact(im) / Fact(ip) / Fact(im - ip);
				Am[im] += aux * pow(x, ip) * pow(y, im-ip) * cos( (im-ip) * ModuleBase::PI / 2.0 );
				Bm[im] += aux * pow(x, ip) * pow(y, im-ip) * sin( (im-ip) * ModuleBase::PI / 2.0 );

				if(ip > 0)
				{
					Gx_Am[im] += aux * ip * pow(x, ip-1) * pow(y, im-ip) * cos( (im-ip) * ModuleBase::PI / 2.0 );
					Gx_Bm[im] += aux * ip * pow(x, ip-1) * pow(y, im-ip) * sin( (im-ip) * ModuleBase::PI / 2.0 );
				}

				if(ip < im)
				{
					Gy_Am[im] += aux * pow(x, ip) * (im - ip) * pow(y, im-ip-1) * cos( (im-ip) * ModuleBase::PI / 2.0 );
					Gy_Bm[im] += aux * pow(x, ip) * (im - ip) * pow(y, im-ip-1) * sin( (im-ip) * ModuleBase::PI / 2.0 );
				}
			}
		}
	}

	//z dependence
	double zdep[20][20];
	double Gx_dep[20][20];
	double Gy_dep[20][20];
	double Gz_dep[20][20];

	for(int il = 0; il < 20; il++)
	{
		ZEROS(zdep[il], 20);
		ZEROS(Gx_dep[il], 20);
		ZEROS(Gy_dep[il], 20);
		ZEROS(Gz_dep[il], 20);
	}

	double z2 = z * z;
	double z3 = z2 * z;
	double z4 = z3 * z;
	//double z5 = z4 * z;

	double r = sqrt(x*x + y*y + z*z);
	double r2 = r * r;
	double r3 = r2 * r;
	double r4 = r3 * r;

	for(int il = 0; il < MaxL+1; il++)
	{
		if(il == 0)
		{
			zdep[0][0] = 1.0;
		}
		else if(il == 1)
		{
			zdep[1][0] = z;
			zdep[1][1] = 1.0;

			Gz_dep[1][0] = 1.0;
		}
		else if(il == 2)
		{
			zdep[2][0] = 0.5 * (3.0 * z2 - r2);
			Gx_dep[2][0] = -x;
			Gy_dep[2][0] = -y;
			Gz_dep[2][0] = 2.0 * z;

			zdep[2][1] = sqrt(3.0) * z;
			Gz_dep[2][1] = sqrt(3.0);

			zdep[2][2] = sqrt(3.0) * 0.5;
		}
		else if(il == 3)
		{
			zdep[3][0] = 2.5 * z3 - 1.5 * z * r2;
			Gx_dep[3][0] = -3.0 * x * z;
			Gy_dep[3][0] = -3.0 * y * z;
			Gz_dep[3][0] = 1.5 * (3.0 * z2 - r2);

			zdep[3][1] = 0.25 * sqrt(6.0) * (5.0 * z2 - r2);
			Gx_dep[3][1] = -0.5 * sqrt(6.0) * x;
			Gy_dep[3][1] = -0.5 * sqrt(6.0) * y;
			Gz_dep[3][1] = sqrt(6.0) * 2.0 * z;

			zdep[3][2] = 0.5 * sqrt(15.0) * z;
			Gz_dep[3][2] = 0.5 * sqrt(15.0);

			zdep[3][3] = 0.25 * sqrt(10.0);
		}
		else if(il == 4)
		{
			zdep[4][0] = 0.125 * (35.0 * z4 - 30.0 * r2 * z2 + 3.0 * r4);
			Gx_dep[4][0] = -7.5 * x * z2 + 1.5 * x * r2;
			Gy_dep[4][0] = -7.5 * y * z2 + 1.5 * y * r2;
			Gz_dep[4][0] = 10.0 * z3 - 6.0 * r2 * z;

			zdep[4][1] = sqrt(10.0) * 0.25 * z * (7.0 * z2 - 3.0 * r2);
			Gx_dep[4][1] = -1.5 * sqrt(10.0) * x * z;
			Gy_dep[4][1] = -1.5 * sqrt(10.0) * y * z;
			Gz_dep[4][1] = 0.75 * sqrt(10.0) * (5.0 * z2 - r2);

			zdep[4][2] = sqrt(5.0) * 0.25 * (7.0 * z2 - r2);
			Gx_dep[4][2] = -0.5 * sqrt(5.0) * x;
			Gy_dep[4][2] = -0.5 * sqrt(5.0) * y;
			Gz_dep[4][2] = 3.0 * sqrt(5.0) * z;

			zdep[4][3] = sqrt(70.0) * 0.25 * z;
			Gz_dep[4][3] = 0.25 * sqrt(70.0);

			zdep[4][4] = sqrt(35.0) * 0.125;
		}
		else if(il == 5)
		{
			zdep[5][0] = 0.125 * z *( 63.0 * z4 - 70.0 * z2 * r2 + 15.0 * r4);
			Gx_dep[5][0] = -17.5 * x * z3 + 7.5 * x * z * r2;
			Gy_dep[5][0] = -17.5 * y * z3 + 7.5 * y * z * r2;
			Gz_dep[5][0] = 175.0 * 0.125 * z4 + 15.0 * 0.125 * r4 - 150.0 * 0.125 * r2 * z2;

			zdep[5][1] = 0.125 * sqrt(15.0) * (21.0 * z4 - 14.0 * z2 * r2 + r4);
			Gx_dep[5][1] = -3.5 * sqrt(15.0) * x * z2 + 0.5 * sqrt(15.0) * x * r2;
			Gy_dep[5][1] = -3.5 * sqrt(15.0) * y * z2 + 0.5 * sqrt(15.0) * y * r2;
			Gz_dep[5][1] = 7.0 * sqrt(15.0) * z3 - 3.0 * sqrt(15.0) * r2 * z;

			zdep[5][2] = 0.25 * sqrt(105.0) * z * (3.0 * z2 - r2);
			Gx_dep[5][2] = -0.5 * sqrt(105.0) * x * z;
			Gy_dep[5][2] = -0.5 * sqrt(105.0) * y * z;
			Gz_dep[5][2] = 0.25 * sqrt(105.0) * (7.0 * z2 - r2);

			zdep[5][3] = 0.0625 * sqrt(70.0) * (9.0 * z2 - r2);
			Gx_dep[5][3] = -0.125 * sqrt(70.0) * x;
			Gy_dep[5][3] = -0.125 * sqrt(70.0) * y;
			Gz_dep[5][3] = sqrt(70.0) * z;

			zdep[5][4] = 0.375 * sqrt(35.0) * z;
			Gz_dep[5][4] = 0.375 * sqrt(35.0);

			zdep[5][5] = 0.1875 * sqrt(14.0);
		}
		else
		{
			for(int im = 0; im <= il; im++)
			{
				int kmax = static_cast<int>( (il - im) / 2 );
				for(int ik = 0; ik <= kmax; ik++)
				{
					int twok = 2 * ik;

					double gamma;
					double aux0, aux1, aux2, aux3;

					aux0 = pow(-1.0, ik) * pow(2.0, -il);
					aux1 = Fact(il) / Fact(ik) / Fact(il-ik);
					aux2 = Fact(2*il - twok) / Fact(il) / Fact(il - twok);
					aux3 = Fact(il - twok) / Fact(il - twok - im);

					gamma = aux0 * aux1 * aux2 * aux3;

					assert(il - twok - im >= 0);
					zdep[il][im] += pow(r, twok) * pow(z, il-twok-im) * gamma;

					if(ik > 0)
					{
						Gx_dep[il][im] += (ik * pow(r2, ik-1) * 2.0 * x) * pow(z, il-twok-im) * gamma;
						Gy_dep[il][im] += (ik * pow(r2, ik-1) * 2.0 * y) * pow(z, il-twok-im) * gamma;
					}

					if(ik == 0)
					{
						if(il > im)
						{
							Gz_dep[il][im] += (il-im) * pow(z, il-im-1) * gamma;
						}
					}
					else
					{
						if(il - twok - im == 0)
						{
							Gz_dep[il][im] += gamma * ik * pow(r2, ik-1) * 2.0 * z;
						}
						else
						{
							Gz_dep[il][im] += gamma * (ik * pow(r2, ik-1) * 2.0 * z * pow(z, il-twok-im)
													+ pow(r, twok) * (il-twok-im) * pow(z, il-twok-im-1));
						}
					}
				}

				if(im >= 1)
				{
					zdep[il][im] *= sqrt(2 * Fact(il - im) / Fact(il + im));
					Gx_dep[il][im] *= sqrt(2 * Fact(il - im) / Fact(il + im));
					Gy_dep[il][im] *= sqrt(2 * Fact(il - im) / Fact(il + im));
					Gz_dep[il][im] *= sqrt(2 * Fact(il - im) / Fact(il + im));

				}
			}
		}
	}

	//calc
	int ic = 0;
	for(int il = 0; il <= MaxL; il++)
	{
		double fac = sqrt( (2.0 * il + 1.0) / ModuleBase::FOUR_PI );

		//m=0
		rly[ic] = Am[0] * zdep[il][0] * fac;
		grly[ic][0] = (Gx_dep[il][0] * Am[0] + zdep[il][0] * Gx_Am[0]) * fac;
		grly[ic][1] = (Gy_dep[il][0] * Am[0] + zdep[il][0] * Gy_Am[0]) * fac;
		grly[ic][2] = Gz_dep[il][0] * Am[0] * fac;

		ic++;

		//m ! = 0
		for(int im = 1; im <= il; im++)
		{
			//m>0
			rly[ic] = Am[im] * zdep[il][im] * pow(-1.0, im) * fac;
			grly[ic][0] = (Gx_dep[il][im] * Am[im] + zdep[il][im] * Gx_Am[im]) * pow(-1.0, im) * fac;
			grly[ic][1] = (Gy_dep[il][im] * Am[im] + zdep[il][im] * Gy_Am[im]) * pow(-1.0, im) * fac;
			grly[ic][2] = Gz_dep[il][im] * Am[im] * pow(-1.0, im) * fac;

			ic++;

			//m<0
			rly[ic] = Bm[im] * zdep[il][im] * pow(-1.0, im) * fac;
			grly[ic][0] = (Gx_dep[il][im] * Bm[im] + zdep[il][im] * Gx_Bm[im]) * pow(-1.0, im) * fac;
			grly[ic][1] = (Gy_dep[il][im] * Bm[im] + zdep[il][im] * Gy_Bm[im]) * pow(-1.0, im) * fac;
			grly[ic][2] = Gz_dep[il][im] * Bm[im] * pow(-1.0, im) * fac;

			ic++;
		}
	}

	return;
}

/*
void Ylm::test(void)
{
	ModuleBase::Vector3<double> R(0.0, 0.0, 1.0);

	double r,r2,r3,r4,r5,r6,r7;
	r = R.norm();
	r2 = r * r;
	r3 = r2 * r;
	r4 = r3 * r;
	r5 = r4 * r;
	r6 = r5 * r;
	r7 = r6 * r;

	//Max L = 7;
	double ylm[64];
	double dylmdr[64][3];

	double rly[64];
	double grly[64][3];
//	std::cout << R.x << " " << R.y << " " << R.z << std::endl;
	get_ylm_real(8, R, ylm, dylmdr);
	rlylm(8, R.x, R.y, R.z, rly, grly);

//	std::cout << R.x << " " << R.y << " " << R.z << std::endl;
	for(int i = 0; i < 64; i++)
	{
		if(i >= 1 && i <= 3)
		{
			dylmdr[i][0] = dylmdr[i][0] * r + ylm[i] * R.x / r;
			dylmdr[i][1] = dylmdr[i][1] * r + ylm[i] * R.y / r;
			dylmdr[i][2] = dylmdr[i][2] * r + ylm[i] * R.z / r;

			ylm[i] *= r;
		}
		if(i >= 4 && i <= 8)
		{
			dylmdr[i][0] = dylmdr[i][0] * r2 + ylm[i] * R.x * 2.0;
			dylmdr[i][1] = dylmdr[i][1] * r2 + ylm[i] * R.y * 2.0;
			dylmdr[i][2] = dylmdr[i][2] * r2 + ylm[i] * R.z * 2.0;

			ylm[i] *= r2;
		}
		if(i >= 9 && i <= 15)
		{
			dylmdr[i][0] = dylmdr[i][0] * r3 + ylm[i] * R.x * 3.0 * r;
			dylmdr[i][1] = dylmdr[i][1] * r3 + ylm[i] * R.y * 3.0 * r;
			dylmdr[i][2] = dylmdr[i][2] * r3 + ylm[i] * R.z * 3.0 * r;

			ylm[i] *= pow(R.norm(),3);
		}
		if(i >= 16 && i <=24)
		{
			dylmdr[i][0] = dylmdr[i][0] * r4 + ylm[i] * R.x * 4.0 * r2;
			dylmdr[i][1] = dylmdr[i][1] * r4 + ylm[i] * R.y * 4.0 * r2;
			dylmdr[i][2] = dylmdr[i][2] * r4 + ylm[i] * R.z * 4.0 * r2;

			ylm[i] *= pow(R.norm(), 4);
		}
		if(i >= 25 &&  i <= 35)
		{
			dylmdr[i][0] = dylmdr[i][0] * r5 + ylm[i] * R.x * 5.0 * r3;
			dylmdr[i][1] = dylmdr[i][1] * r5 + ylm[i] * R.y * 5.0 * r3;
			dylmdr[i][2] = dylmdr[i][2] * r5 + ylm[i] * R.z * 5.0 * r3;

			ylm[i] *= pow(R.norm(), 5);
		}
		if(i >= 36 && i <= 48)
		{
			dylmdr[i][0] = dylmdr[i][0] * r6 + ylm[i] * R.x * 6.0 * r4;
			dylmdr[i][1] = dylmdr[i][1] * r6 + ylm[i] * R.y * 6.0 * r4;
			dylmdr[i][2] = dylmdr[i][2] * r6 + ylm[i] * R.z * 6.0 * r4;
			ylm[i] *= pow(R.norm(), 6);
		}
		if(i >= 49 && i <= 63)
		{
			dylmdr[i][0] = dylmdr[i][0] * r7 + ylm[i] * R.x * 7.0 * r5;
			dylmdr[i][1] = dylmdr[i][1] * r7 + ylm[i] * R.y * 7.0 * r5;
			dylmdr[i][2] = dylmdr[i][2] * r7 + ylm[i] * R.z * 7.0 * r5;
			ylm[i] *= pow(R.norm(), 7);
		}

		std::cout << grly[i][0] << std::setw(20) << grly[i][1] << std::setw(20) << grly[i][2] << std::endl;
	}

	return;
}
*/

void Ylm::ZEROS(double u[], const int& n)
{
	for(int i = 0; i < n; i++)
	{
		u[i] = 0.0;
	}
	return;
}


//==========================================================
// MEMBER FUNCTION :
// NAME : Fact ( n! )
// NAME : Semi_Fact ( n!! )
//==========================================================
long double Ylm::Fact(const int n)
{
	long double f = 1;
	for(int i=n; i>1; i--)
	{
		f *= i;
	}
	return f;
}


int Ylm::Semi_Fact(const int n)
{
	int semif = 1;
	for(int i=n; i>2; i -= 2)
	{
		semif *= i;
	}
	return semif;
}


double Ylm::sgn(const double x)
{
	if(x < 0.0) return -1.0;
	if(x > 0.0) return 1.0;
	return 0.0;
}

}
