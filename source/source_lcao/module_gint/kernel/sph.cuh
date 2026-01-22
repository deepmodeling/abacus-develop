#pragma once

namespace ModuleGint
{

static __device__ void sph_harm(
    const int nwl,
    const double* __restrict__ ylmcoef,
    const double x,
	const double y,
	const double z,
    double* __restrict__ ylma
)
{
   /***************************
   L = 0
   ***************************/
   ylma[0] = ylmcoef[0]; // l=0, m=0
   double tmp0;
   if (nwl == 0)
       return;

   /***************************
   L = 1
   ***************************/
   ylma[1] = ylmcoef[1] * z;  // l=1, m=0
   ylma[2] = -ylmcoef[1] * x; // l=1, m=1
   ylma[3] = -ylmcoef[1] * y; // l=1, m=-1
   if (nwl == 1)
       return;

   /***************************
   L = 2
   ***************************/
   tmp0=ylmcoef[3] * ylma[0];
   ylma[4] = ylmcoef[2] * z * ylma[1] - tmp0 ; // l=2, m=0
   tmp0 = ylmcoef[4] * z;
   ylma[5] = tmp0 * ylma[2]; // l=2,m=1
   ylma[6] = tmp0 * ylma[3]; // l=2,m=-1

   tmp0 = ylmcoef[4] * x;
   ylma[7] = ylmcoef[5] * ylma[4] - ylmcoef[6] * ylma[0]
             - tmp0 * ylma[2]; // l=2,m=2
   ylma[8] = -tmp0 * ylma[3];
   if (nwl == 2)
       return;

   /***************************
   L = 3
   ***************************/
   tmp0=ylmcoef[8] * ylma[1];
   ylma[9] = ylmcoef[7] * z * ylma[4] - tmp0; // l=3, m=0

   tmp0 = ylmcoef[9] * z;
   ylma[10] = tmp0 * ylma[5] - ylmcoef[10] * ylma[2]; // l=3,m=1
   ylma[11] = tmp0 * ylma[6] - ylmcoef[10] * ylma[3]; // l=3,m=-1

   tmp0 = ylmcoef[11] * z;
   ylma[12] = tmp0 * ylma[7]; // l=3,m=2
   ylma[13] = tmp0 * ylma[8]; // l=3,m=-2

   tmp0 = ylmcoef[14] * x;
   ylma[14] = ylmcoef[12] * ylma[10] - ylmcoef[13] * ylma[2]
              - tmp0 * ylma[7]; // l=3,m=3
   ylma[15] = ylmcoef[12] * ylma[11] - ylmcoef[13] * ylma[3]
              - tmp0 * ylma[8]; // l=3,m=-3
   if (nwl == 3)
       return;

   /***************************
   L = 4
   ***************************/
   tmp0=ylmcoef[16] * ylma[4];
   ylma[16] = ylmcoef[15] * z * ylma[9] - tmp0; // l=4,m=0

   tmp0 = ylmcoef[17] * z;
   ylma[17] = tmp0 * ylma[10] - ylmcoef[18] * ylma[5]; // l=4,m=1
   ylma[18] = tmp0 * ylma[11] - ylmcoef[18] * ylma[6]; // l=4,m=-1

   tmp0 = ylmcoef[19] * z;
   ylma[19] = tmp0 * ylma[12] - ylmcoef[20] * ylma[7]; // l=4,m=2
   ylma[20] = tmp0 * ylma[13] - ylmcoef[20] * ylma[8]; // l=4,m=-2

   tmp0 = 3.0 * z;
   ylma[21] = tmp0 * ylma[14]; // l=4,m=3
   ylma[22] = tmp0 * ylma[15]; // l=4,m=-3

   tmp0 = ylmcoef[23] * x;
   ylma[23] = ylmcoef[21] * ylma[19] - ylmcoef[22] * ylma[7]
              - tmp0 * ylma[14]; // l=4,m=4
   ylma[24] = ylmcoef[21] * ylma[20] - ylmcoef[22] * ylma[8]
              - tmp0 * ylma[15]; // l=4,m=-4
   if (nwl == 4)
       return;

   /***************************
   L = 5
   ***************************/
   tmp0=ylmcoef[25] * ylma[9];
   ylma[25]
       = ylmcoef[24] * z * ylma[16] - tmp0; // l=5,m=0

   tmp0 = ylmcoef[26] * z;
   ylma[26] = tmp0 * ylma[17] - ylmcoef[27] * ylma[10]; // l=5,m=1
   ylma[27] = tmp0 * ylma[18] - ylmcoef[27] * ylma[11]; // l=5,m=-1

   tmp0 = ylmcoef[28] * z;
   ylma[28] = tmp0 * ylma[19] - ylmcoef[29] * ylma[12]; // l=5,m=2
   ylma[29] = tmp0 * ylma[20] - ylmcoef[29] * ylma[13]; // l=5,m=-2

   tmp0 = ylmcoef[30] * z;
   ylma[30] = tmp0 * ylma[21] - ylmcoef[31] * ylma[14]; // l=5,m=3
   ylma[31] = tmp0 * ylma[22] - ylmcoef[31] * ylma[15]; // l=5,m=-3

   tmp0 = ylmcoef[32] * z;
   ylma[32] = tmp0 * ylma[23]; // l=5,m=4
   ylma[33] = tmp0 * ylma[24]; // l=5,m=-4

   tmp0 = ylmcoef[35] * x;
   ylma[34] = ylmcoef[33] * ylma[30] - ylmcoef[34] * ylma[14]
              - tmp0 * ylma[23]; // l=5,m=5
   ylma[35] = ylmcoef[33] * ylma[31] - ylmcoef[34] * ylma[15]
              - tmp0 * ylma[24]; // l=5,m=-5
   if (nwl == 5)
       return;
   /*
   // if nwl > 5
   for (int il = 6; il <= nwl; il++)
   {
       int istart = il * il;
       int istart1 = (il - 1) * (il - 1);
       int istart2 = (il - 2) * (il - 2);

       double fac2 = sqrt(4.0 * istart - 1.0);
       double fac4 = sqrt(4.0 * istart1 - 1.0);

       for (int im = 0; im < 2 * il - 1; im++)
       {
           int imm = (im + 1) / 2;
           ylma[istart + im] = fac2 / sqrt((double)istart - imm * imm) * (z
   * ylma[istart1 + im] - sqrt((double)istart1 - imm * imm) / fac4 *
   ylma[istart2 + im]);
       }

       double bl1 = sqrt(2.0 * il / (2.0 * il + 1.0));
       double bl2 = sqrt((2.0 * il - 2.0) / (2.0 * il - 1.0));
       double bl3 = sqrt(2.0) / fac2;

       ylma[istart + 2 * il - 1] = (bl3 * ylma[istart + 2 * il - 5] - bl2 *
   ylma[istart2 + 2 * il - 5] - 2.0 * x * ylma[istart1 + 2 * il - 3]) /
   bl1; ylma[istart + 2 * il] = (bl3 * ylma[istart + 2 * il - 4] - bl2 *
   ylma[istart2 + 2 * il - 4] - 2.0 * x * ylma[istart1 + 2 * il - 2]) /
   bl1;
   }*/
}

static __device__ void grad_rl_sph_harm(
    const int nwl,
    const double* __restrict__ ylmcoef,
    const double x,
	const double y,
	const double z,
    double* __restrict__ rly,
    double* __restrict__ grly
)
{
    double r2 = x * x + y * y + z * z;
    double tx = x * 2;
    double ty = y * 2;
    double tz = z * 2;

    //begin calculation
	/***************************
			 L = 0
	***************************/
	rly[0] = ylmcoef[0]; //l=0, m=0
	grly[0] = grly[1] = grly[2] = 0.0;
	if (nwl == 0) return;

	/***************************
			 L = 1
	***************************/
	rly[1] = ylmcoef[1]*z; //l=1, m=0
	grly[3] = grly[4] = 0.0;
	grly[5] = ylmcoef[1];

	rly[2] = -ylmcoef[1]*x; //l=1, m=1
	grly[7] = grly[8] = 0.0;
	grly[6] = -ylmcoef[1];

	rly[3] = -ylmcoef[1]*y; //l=1, m=-1
	grly[9] = grly[11] = 0.0;
	grly[10] = -ylmcoef[1];

	if (nwl == 1) return;

	/***************************
			 L = 2
	***************************/
	rly[4] = ylmcoef[2]*z*rly[1]-ylmcoef[3]*rly[0]*r2;//l=2, m=0
	grly[12] = ylmcoef[2]*z*grly[3]-ylmcoef[3]*(grly[0]*r2+rly[0]*tx);//l=2, m=0
	grly[13] = ylmcoef[2]*z*grly[4]-ylmcoef[3]*(grly[1]*r2+rly[0]*ty);//l=2, m=0
	grly[14] = ylmcoef[2]*(z*grly[5]+rly[1])-ylmcoef[3]*(grly[2]*r2+rly[0]*tz);//l=2, m=0


	double tmp0 = ylmcoef[4]*z;
	rly[5] = tmp0*rly[2];//l=2,m=1
	grly[15] = tmp0*grly[6];
	grly[16] = tmp0*grly[7];
	grly[17] = ylmcoef[4]*(rly[2]+z*grly[8]);

	rly[6] = tmp0*rly[3];//l=2,m=-1
	grly[18] = tmp0*grly[9];
	grly[19] = tmp0*grly[10];
	grly[20] = ylmcoef[4]*(rly[3]+z*grly[11]);

	double tmp2 = ylmcoef[4]*x;
	rly[7]= ylmcoef[5]*rly[4]-ylmcoef[6]*rly[0]*r2 - tmp2*rly[2];//l=2,m=2
	grly[21] = ylmcoef[5]*grly[12]-ylmcoef[6]*(rly[0]*tx+grly[0]*r2)-ylmcoef[4]*(x*grly[6]+rly[2]);

//	std::cout << "\np1 = "<< ylmcoef[5]*grly[12] << " p2 = " << -ylmcoef[6]*rly[0]*tx
//						<< " p3 = " << -ylmcoef[4]*x*grly[6] << " p4 = " << -ylmcoef[4]*rly[2] << std::endl;

	grly[22] = ylmcoef[5]*grly[13]-ylmcoef[6]*(rly[0]*ty+grly[1]*r2)-tmp2*grly[7];
	grly[23] = ylmcoef[5]*grly[14]-ylmcoef[6]*(rly[0]*tz+grly[2]*r2)-tmp2*grly[8];

	rly[8] = -tmp2*rly[3];
	grly[24] = -ylmcoef[4]*(rly[3]+x*grly[9]);
	grly[25] = -tmp2*grly[10];
	grly[26] = -tmp2*grly[11];
//	rly[8] = tmp1+tmp2*rly[3];//l=2,m=-2
	if (nwl == 2) return;

	/***************************
			 L = 3
	***************************/
	rly[9] = ylmcoef[7]*z*rly[4]-ylmcoef[8]*rly[1]*r2; //l=3, m=0
	grly[27] = ylmcoef[7]*z*grly[12]-ylmcoef[8]*(rly[1]*tx+grly[3]*r2);
	grly[28] = ylmcoef[7]*z*grly[13]-ylmcoef[8]*(rly[1]*ty+grly[4]*r2);
	grly[29] = ylmcoef[7]*(rly[4]+z*grly[14])-ylmcoef[8]*(rly[1]*tz+grly[5]*r2);

	double tmp3 = ylmcoef[9]*z;
	rly[10] = tmp3*rly[5]-ylmcoef[10]*rly[2]*r2;//l=3,m=1
	grly[30] = tmp3*grly[15]-ylmcoef[10]*(grly[6]*r2+rly[2]*tx);
	grly[31] = tmp3*grly[16]-ylmcoef[10]*(grly[7]*r2+rly[2]*ty);
	grly[32] = ylmcoef[9]*(z*grly[17]+rly[5])-ylmcoef[10]*(grly[8]*r2+rly[2]*tz);

	rly[11] = tmp3*rly[6]-ylmcoef[10]*rly[3]*r2;//l=3,m=-1
	grly[33] = tmp3*grly[18]-ylmcoef[10]*(grly[9]*r2+rly[3]*tx);
	grly[34] = tmp3*grly[19]-ylmcoef[10]*(grly[10]*r2+rly[3]*ty);
	grly[35] = ylmcoef[9]*(z*grly[20]+rly[6])-ylmcoef[10]*(grly[11]*r2+rly[3]*tz);

	double tmp4 = ylmcoef[11]*z;
	rly[12] = tmp4*rly[7];//l=3,m=2
	grly[36] = tmp4*grly[21];
	grly[37] = tmp4*grly[22];
	grly[38] = ylmcoef[11]*(z*grly[23]+rly[7]);

	rly[13] = tmp4*rly[8];//l=3,m=-2
	grly[39] = tmp4*grly[24];
	grly[40] = tmp4*grly[25];
	grly[41] = ylmcoef[11]*(z*grly[26]+rly[8]);

	double tmp5 = ylmcoef[14]*x;
	rly[14] = ylmcoef[12]*rly[10]-ylmcoef[13]*rly[2]*r2-tmp5*rly[7];//l=3,m=3
	grly[42] = ylmcoef[12]*grly[30]-ylmcoef[13]*(rly[2]*tx+grly[6]*r2)-ylmcoef[14]*(rly[7]+x*grly[21]);
	grly[43] = ylmcoef[12]*grly[31]-ylmcoef[13]*(rly[2]*ty+grly[7]*r2)-tmp5*grly[22];
	grly[44] = ylmcoef[12]*grly[32]-ylmcoef[13]*(rly[2]*tz+grly[8]*r2)-tmp5*grly[23];

	rly[15] = ylmcoef[12]*rly[11]-ylmcoef[13]*rly[3]*r2-tmp5*rly[8];//l=3,m=-3
	grly[45] = ylmcoef[12]*grly[33]-ylmcoef[13]*(rly[3]*tx+grly[9]*r2)-ylmcoef[14]*(rly[8]+x*grly[24]);
	grly[46] = ylmcoef[12]*grly[34]-ylmcoef[13]*(rly[3]*ty+grly[10]*r2)-tmp5*grly[25];
	grly[47] = ylmcoef[12]*grly[35]-ylmcoef[13]*(rly[3]*tz+grly[11]*r2)-tmp5*grly[26];
	if (nwl == 3) return;

	/***************************
			 L = 4
	***************************/
	rly[16] = ylmcoef[15]*z*rly[9]-ylmcoef[16]*rly[4]*r2;//l=4,m=0
	grly[48] = ylmcoef[15]*z*grly[27]-ylmcoef[16]*(rly[4]*tx+grly[12]*r2);
	grly[49] = ylmcoef[15]*z*grly[28]-ylmcoef[16]*(rly[4]*ty+grly[13]*r2);
	grly[50] = ylmcoef[15]*(z*grly[29]+rly[9])-ylmcoef[16]*(rly[4]*tz+grly[14]*r2);

	double tmp6 = ylmcoef[17]*z;
	rly[17] = tmp6*rly[10]-ylmcoef[18]*rly[5]*r2;//l=4,m=1
	grly[51] = tmp6*grly[30]-ylmcoef[18]*(rly[5]*tx+grly[15]*r2);
	grly[52] = tmp6*grly[31]-ylmcoef[18]*(rly[5]*ty+grly[16]*r2);
	grly[53] = ylmcoef[17]*(z*grly[32]+rly[10])-ylmcoef[18]*(rly[5]*tz+grly[17]*r2);

	rly[18] = tmp6*rly[11]-ylmcoef[18]*rly[6]*r2;//l=4,m=-1
	grly[54] = tmp6*grly[33]-ylmcoef[18]*(rly[6]*tx+grly[18]*r2);
	grly[55] = tmp6*grly[34]-ylmcoef[18]*(rly[6]*ty+grly[19]*r2);
	grly[56] = ylmcoef[17]*(z*grly[35]+rly[11])-ylmcoef[18]*(rly[6]*tz+grly[20]*r2);

	double tmp7 = ylmcoef[19]*z;
	rly[19] = tmp7*rly[12]-ylmcoef[20]*rly[7]*r2;//l=4,m=2
	grly[57] = tmp7*grly[36]-ylmcoef[20]*(rly[7]*tx+grly[21]*r2);
	grly[58] = tmp7*grly[37]-ylmcoef[20]*(rly[7]*ty+grly[22]*r2);
	grly[59] = ylmcoef[19]*(z*grly[38]+rly[12])-ylmcoef[20]*(rly[7]*tz+grly[23]*r2);

	rly[20] = tmp7*rly[13]-ylmcoef[20]*rly[8]*r2;//l=4,m=-2
	grly[60] = tmp7*grly[39]-ylmcoef[20]*(rly[8]*tx+grly[24]*r2);
	grly[61] = tmp7*grly[40]-ylmcoef[20]*(rly[8]*ty+grly[25]*r2);
	grly[62] = ylmcoef[19]*(z*grly[41]+rly[13])-ylmcoef[20]*(rly[8]*tz+grly[26]*r2);

	double tmp8 = 3.0*z;
	rly[21] = tmp8*rly[14];//l=4,m=3
	grly[63] = tmp8*grly[42];
	grly[64] = tmp8*grly[43];
	grly[65] = 3.0*(z*grly[44]+rly[14]);


	rly[22] = tmp8*rly[15];//l=4,m=-3
	grly[66] = tmp8*grly[45];
	grly[67] = tmp8*grly[46];
	grly[68] = 3.0*(z*grly[47]+rly[15]);

	double tmp9 = ylmcoef[23]*x;
	rly[23] = ylmcoef[21]*rly[19]-ylmcoef[22]*rly[7]*r2-tmp9*rly[14];//l=4,m=4
	grly[69] = ylmcoef[21]*grly[57]-ylmcoef[22]*(rly[7]*tx+grly[21]*r2)-ylmcoef[23]*(x*grly[42]+rly[14]);
	grly[70] = ylmcoef[21]*grly[58]-ylmcoef[22]*(rly[7]*ty+grly[22]*r2)-tmp9*grly[43];
	grly[71] = ylmcoef[21]*grly[59]-ylmcoef[22]*(rly[7]*tz+grly[23]*r2)-tmp9*grly[44];

	rly[24] = ylmcoef[21]*rly[20]-ylmcoef[22]*rly[8]*r2-tmp9*rly[15];//l=4,m=-4
	grly[72] = ylmcoef[21]*grly[60]-ylmcoef[22]*(rly[8]*tx+grly[24]*r2)-ylmcoef[23]*(x*grly[45]+rly[15]);
	grly[73] = ylmcoef[21]*grly[61]-ylmcoef[22]*(rly[8]*ty+grly[25]*r2)-tmp9*grly[46];
	grly[74] = ylmcoef[21]*grly[62]-ylmcoef[22]*(rly[8]*tz+grly[26]*r2)-tmp9*grly[47];

	if (nwl == 4) return;

	/***************************
			 L = 5
	***************************/
	rly[25] = ylmcoef[24]*z*rly[16]-ylmcoef[25]*rly[9]*r2;//l=5,m=0
	grly[75] = ylmcoef[24]*z*grly[48]-ylmcoef[25]*(rly[9]*tx+grly[27]*r2);
	grly[76] = ylmcoef[24]*z*grly[49]-ylmcoef[25]*(rly[9]*ty+grly[28]*r2);
	grly[77] = ylmcoef[24]*(z*grly[50]+rly[16])-ylmcoef[25]*(rly[9]*tz+grly[29]*r2);

	double tmp10 = ylmcoef[26]*z;
	rly[26] = tmp10*rly[17]-ylmcoef[27]*rly[10]*r2;//l=5,m=1
	grly[78] = tmp10*grly[51]-ylmcoef[27]*(rly[10]*tx+grly[30]*r2);
	grly[79] = tmp10*grly[52]-ylmcoef[27]*(rly[10]*ty+grly[31]*r2);
	grly[80] = ylmcoef[26]*(z*grly[53]+rly[17])-ylmcoef[27]*(rly[10]*tz+grly[32]*r2);

	rly[27] = tmp10*rly[18]-ylmcoef[27]*rly[11]*r2;//l=5,m=-1
	grly[81] = tmp10*grly[54]-ylmcoef[27]*(rly[11]*tx+grly[33]*r2);
	grly[82] = tmp10*grly[55]-ylmcoef[27]*(rly[11]*ty+grly[34]*r2);
	grly[83] = ylmcoef[26]*(z*grly[56]+rly[18])-ylmcoef[27]*(rly[11]*tz+grly[35]*r2);

	double tmp11 = ylmcoef[28]*z;
	rly[28] = tmp11*rly[19]-ylmcoef[29]*rly[12]*r2;//l=5,m=2
	grly[84] = tmp11*grly[57]-ylmcoef[29]*(rly[12]*tx+grly[36]*r2);
	grly[85] = tmp11*grly[58]-ylmcoef[29]*(rly[12]*ty+grly[37]*r2);
	grly[86] = ylmcoef[28]*(z*grly[59]+rly[19])-ylmcoef[29]*(rly[12]*tz+grly[38]*r2);

	rly[29] = tmp11*rly[20]-ylmcoef[29]*rly[13]*r2;//l=5,m=-2
	grly[87] = tmp11*grly[60]-ylmcoef[29]*(rly[13]*tx+grly[39]*r2);
	grly[88] = tmp11*grly[61]-ylmcoef[29]*(rly[13]*ty+grly[40]*r2);
	grly[89] = ylmcoef[28]*(z*grly[62]+rly[20])-ylmcoef[29]*(rly[13]*tz+grly[41]*r2);

	double tmp12 = ylmcoef[30]*z;
	rly[30] = tmp12*rly[21]-ylmcoef[31]*rly[14]*r2;//l=5,m=3
	grly[90] = tmp12*grly[63]-ylmcoef[31]*(grly[42]*r2+rly[14]*tx);
	grly[91] = tmp12*grly[64]-ylmcoef[31]*(grly[43]*r2+rly[14]*ty);
	grly[92] = ylmcoef[30]*(z*grly[65]+rly[21])-ylmcoef[31]*(grly[44]*r2+rly[14]*tz);

	rly[31] = tmp12*rly[22]-ylmcoef[31]*rly[15]*r2;//l=5,m=-3
	grly[93] = tmp12*grly[66]-ylmcoef[31]*(grly[45]*r2+rly[15]*tx);
	grly[94] = tmp12*grly[67]-ylmcoef[31]*(grly[46]*r2+rly[15]*ty);
	grly[95] = ylmcoef[30]*(z*grly[68]+rly[22])-ylmcoef[31]*(grly[47]*r2+rly[15]*tz);

	double tmp13 = ylmcoef[32]*z;
	rly[32] = tmp13*rly[23];//l=5,m=4
	grly[96] = tmp13*grly[69];
	grly[97] = tmp13*grly[70];
	grly[98] = ylmcoef[32]*(rly[23]+z*grly[71]);

	rly[33] = tmp13*rly[24];//l=5,m=-4
	grly[99] = tmp13*grly[72];
	grly[100] = tmp13*grly[73];
	grly[101] = ylmcoef[32]*(rly[24]+z*grly[74]);

	double tmp14 = ylmcoef[35]*x;
	rly[34] = ylmcoef[33]*rly[30]-ylmcoef[34]*rly[14]*r2-tmp14*rly[23];//l=5,m=5
	grly[102] = ylmcoef[33]*grly[90]-ylmcoef[34]*(rly[14]*tx+grly[42]*r2)-ylmcoef[35]*(x*grly[69]+rly[23]);
	grly[103] = ylmcoef[33]*grly[91]-ylmcoef[34]*(rly[14]*ty+grly[43]*r2)-tmp14*grly[70];
	grly[104] = ylmcoef[33]*grly[92]-ylmcoef[34]*(rly[14]*tz+grly[44]*r2)-tmp14*grly[71];

	rly[35] = ylmcoef[33]*rly[31]-ylmcoef[34]*rly[15]*r2-tmp14*rly[24];//l=5,m=-5
	grly[105] = ylmcoef[33]*grly[93]-ylmcoef[34]*(rly[15]*tx+grly[45]*r2)-ylmcoef[35]*(x*grly[72]+rly[24]);
	grly[106] = ylmcoef[33]*grly[94]-ylmcoef[34]*(rly[15]*ty+grly[46]*r2)-tmp14*grly[73];
	grly[107] = ylmcoef[33]*grly[95]-ylmcoef[34]*(rly[15]*tz+grly[47]*r2)-tmp14*grly[74];

	if (nwl == 5) return;
}
}