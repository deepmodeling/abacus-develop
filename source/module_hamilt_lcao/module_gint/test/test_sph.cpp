#include "test_sph.h"
using namespace std;

void sph_harm(const int& Lmax, // max momentum of l
              const double& xdr,
              const double& ydr,
              const double& zdr,
              std::vector<double>& rly,
              double* ylmcoef)
{

    // begin calculation
    /***************************
             L = 0
    ***************************/
    rly[0] = ylmcoef[0]; // l=0, m=0
    if (Lmax == 0)
        return;

    /***************************
             L = 1
    ***************************/
    rly[1] = ylmcoef[1] * zdr;  // l=1, m=0
    rly[2] = -ylmcoef[1] * xdr; // l=1, m=1
    rly[3] = -ylmcoef[1] * ydr; // l=1, m=-1
    if (Lmax == 1)
        return;

    /***************************
             L = 2
    ***************************/
    double tmp0 = ylmcoef[3] * rly[0];
    rly[4] = ylmcoef[2] * zdr * rly[1] - tmp0; // l=2, m=0

    tmp0 = ylmcoef[4] * zdr;
    rly[5] = tmp0 * rly[2]; // l=2,m=1
    rly[6] = tmp0 * rly[3]; // l=2,m=-1

    double tmp2 = ylmcoef[4] * xdr;
    rly[7]
        = ylmcoef[5] * rly[4] - ylmcoef[6] * rly[0] - tmp2 * rly[2]; // l=2,m=2
    rly[8] = -tmp2 * rly[3];
    //	rly[8] = tmp1+tmp2*rly[3];//l=2,m=-2
    if (Lmax == 2)
        return;

    /***************************
             L = 3
    ***************************/
    tmp0 = ylmcoef[8] * rly[1];
    rly[9] = ylmcoef[7] * zdr * rly[4] - tmp0; // l=3, m=0

    double tmp3 = ylmcoef[9] * zdr;
    rly[10] = tmp3 * rly[5] - ylmcoef[10] * rly[2]; // l=3,m=1
    rly[11] = tmp3 * rly[6] - ylmcoef[10] * rly[3]; // l=3,m=-1

    double tmp4 = ylmcoef[11] * zdr;
    rly[12] = tmp4 * rly[7]; // l=3,m=2
    rly[13] = tmp4 * rly[8]; // l=3,m=-2

    double tmp5 = ylmcoef[14] * xdr;
    rly[14] = ylmcoef[12] * rly[10] - ylmcoef[13] * rly[2]
              - tmp5 * rly[7]; // l=3,m=3
    rly[15] = ylmcoef[12] * rly[11] - ylmcoef[13] * rly[3]
              - tmp5 * rly[8]; // l=3,m=-3
    if (Lmax == 3)
        return;

    /***************************
             L = 4
    ***************************/
    tmp0 = ylmcoef[16] * rly[4];
    rly[16] = ylmcoef[15] * zdr * rly[9] - tmp0; // l=4,m=0

    double tmp6 = ylmcoef[17] * zdr;
    rly[17] = tmp6 * rly[10] - ylmcoef[18] * rly[5]; // l=4,m=1
    rly[18] = tmp6 * rly[11] - ylmcoef[18] * rly[6]; // l=4,m=-1

    double tmp7 = ylmcoef[19] * zdr;
    rly[19] = tmp7 * rly[12] - ylmcoef[20] * rly[7]; // l=4,m=2
    rly[20] = tmp7 * rly[13] - ylmcoef[20] * rly[8]; // l=4,m=-2

    double tmp8 = 3.0 * zdr;
    rly[21] = tmp8 * rly[14]; // l=4,m=3
    rly[22] = tmp8 * rly[15]; // l=4,m=-3

    double tmp9 = ylmcoef[23] * xdr;
    rly[23] = ylmcoef[21] * rly[19] - ylmcoef[22] * rly[7]
              - tmp9 * rly[14]; // l=4,m=4
    rly[24] = ylmcoef[21] * rly[20] - ylmcoef[22] * rly[8]
              - tmp9 * rly[15]; // l=4,m=-4
    if (Lmax == 4)
        return;

    /***************************
             L = 5
    ***************************/
    tmp0 = ylmcoef[25] * rly[9];
    rly[25] = ylmcoef[24] * zdr * rly[16] - tmp0; // l=5,m=0

    double tmp10 = ylmcoef[26] * zdr;
    rly[26] = tmp10 * rly[17] - ylmcoef[27] * rly[10]; // l=5,m=1
    rly[27] = tmp10 * rly[18] - ylmcoef[27] * rly[11]; // l=5,m=-1

    double tmp11 = ylmcoef[28] * zdr;
    rly[28] = tmp11 * rly[19] - ylmcoef[29] * rly[12]; // l=5,m=2
    rly[29] = tmp11 * rly[20] - ylmcoef[29] * rly[13]; // l=5,m=-2

    double tmp12 = ylmcoef[30] * zdr;
    rly[30] = tmp12 * rly[21] - ylmcoef[31] * rly[14]; // l=5,m=3
    rly[31] = tmp12 * rly[22] - ylmcoef[31] * rly[15]; // l=5,m=-3

    double tmp13 = ylmcoef[32] * zdr;
    rly[32] = tmp13 * rly[23]; // l=5,m=4
    rly[33] = tmp13 * rly[24]; // l=5,m=-4

    double tmp14 = ylmcoef[35] * xdr;
    rly[34] = ylmcoef[33] * rly[30] - ylmcoef[34] * rly[14]
              - tmp14 * rly[23]; // l=5,m=5
    rly[35] = ylmcoef[33] * rly[31] - ylmcoef[34] * rly[15]
              - tmp14 * rly[24]; // l=5,m=-5
    if (Lmax == 5)
        return;

    // if Lmax > 5
    for (int il = 6; il <= Lmax; il++)
    {
        int istart = il * il;
        int istart1 = (il - 1) * (il - 1);
        int istart2 = (il - 2) * (il - 2);

        double fac2 = sqrt(4.0 * istart - 1.0);
        double fac4 = sqrt(4.0 * istart1 - 1.0);

        for (int im = 0; im < 2 * il - 1; im++)
        {
            int imm = (im + 1) / 2;
            //			if (im % 2 == 0) imm *= -1;

            rly[istart + im] = fac2 / sqrt((double)istart - imm * imm)
                               * (zdr * rly[istart1 + im]
                                  - sqrt((double)istart1 - imm * imm) / fac4
                                        * rly[istart2 + im]);
        }

        double bl1 = sqrt(2.0 * il / (2.0 * il + 1.0));
        double bl2 = sqrt((2.0 * il - 2.0) / (2.0 * il - 1.0));
        double bl3 = sqrt(2.0) / fac2;

        rly[istart + 2 * il - 1]
            = (bl3 * rly[istart + 2 * il - 5] - bl2 * rly[istart2 + 2 * il - 5]
               - 2.0 * xdr * rly[istart1 + 2 * il - 3])
              / bl1;
        rly[istart + 2 * il]
            = (bl3 * rly[istart + 2 * il - 4] - bl2 * rly[istart2 + 2 * il - 4]
               - 2.0 * xdr * rly[istart1 + 2 * il - 2])
              / bl1;
    }

    return;
}
void grad_rl_sph_harm(const int& Lmax, // max momentum of L
                      const double& x,
                      const double& y,
                      const double& z,
                      double* rly,
                      double** grly,
                      const double* ylmcoef)
{
    double radius2 = x * x + y * y + z * z;
    double tx = 2.0 * x;
    double ty = 2.0 * y;
    double tz = 2.0 * z;

    // begin calculation
    /***************************
             L = 0
    ***************************/
    rly[0] = ylmcoef[0]; // l=0, m=0
    grly[0][0] = grly[0][1] = grly[0][2] = 0.0;
    if (Lmax == 0)
        return;

    /***************************
             L = 1
    ***************************/
    rly[1] = ylmcoef[1] * z; // l=1, m=0
    grly[1][0] = grly[1][1] = 0.0;
    grly[1][2] = ylmcoef[1];

    rly[2] = -ylmcoef[1] * x; // l=1, m=1
    grly[2][1] = grly[2][2] = 0.0;
    grly[2][0] = -ylmcoef[1];

    rly[3] = -ylmcoef[1] * y; // l=1, m=-1
    grly[3][0] = grly[3][2] = 0.0;
    grly[3][1] = -ylmcoef[1];

    if (Lmax == 1)
        return;

    /***************************
             L = 2
    ***************************/
    rly[4]
        = ylmcoef[2] * z * rly[1] - ylmcoef[3] * rly[0] * radius2; // l=2, m=0
    grly[4][0]
        = ylmcoef[2] * z * grly[1][0]
          - ylmcoef[3] * (grly[0][0] * radius2 + rly[0] * tx); // l=2, m=0
    grly[4][1]
        = ylmcoef[2] * z * grly[1][1]
          - ylmcoef[3] * (grly[0][1] * radius2 + rly[0] * ty); // l=2, m=0
    grly[4][2]
        = ylmcoef[2] * (z * grly[1][2] + rly[1])
          - ylmcoef[3] * (grly[0][2] * radius2 + rly[0] * tz); // l=2, m=0

    double tmp0 = ylmcoef[4] * z;
    rly[5] = tmp0 * rly[2]; // l=2,m=1
    grly[5][0] = tmp0 * grly[2][0];
    grly[5][1] = tmp0 * grly[2][1];
    grly[5][2] = ylmcoef[4] * (rly[2] + z * grly[2][2]);

    rly[6] = tmp0 * rly[3]; // l=2,m=-1
    grly[6][0] = tmp0 * grly[3][0];
    grly[6][1] = tmp0 * grly[3][1];
    grly[6][2] = ylmcoef[4] * (rly[3] + z * grly[3][2]);

    double tmp2 = ylmcoef[4] * x;
    rly[7] = ylmcoef[5] * rly[4] - ylmcoef[6] * rly[0] * radius2
             - tmp2 * rly[2]; // l=2,m=2
    grly[7][0] = ylmcoef[5] * grly[4][0]
                 - ylmcoef[6] * (rly[0] * tx + grly[0][0] * radius2)
                 - ylmcoef[4] * (x * grly[2][0] + rly[2]);

    //	std::cout << "\np1 = "<< ylmcoef[5]*grly[4][0] << " p2 = " <<
    //-ylmcoef[6]*rly[0]*tx
    //						<< " p3 = " << -ylmcoef[4]*x*grly[2][0] << " p4 = "
    //<< -ylmcoef[4]*rly[2] << std::endl;

    grly[7][1] = ylmcoef[5] * grly[4][1]
                 - ylmcoef[6] * (rly[0] * ty + grly[0][1] * radius2)
                 - tmp2 * grly[2][1];
    grly[7][2] = ylmcoef[5] * grly[4][2]
                 - ylmcoef[6] * (rly[0] * tz + grly[0][2] * radius2)
                 - tmp2 * grly[2][2];

    rly[8] = -tmp2 * rly[3];
    grly[8][0] = -ylmcoef[4] * (rly[3] + x * grly[3][0]);
    grly[8][1] = -tmp2 * grly[3][1];
    grly[8][2] = -tmp2 * grly[3][2];
    //	rly[8] = tmp1+tmp2*rly[3];//l=2,m=-2
    if (Lmax == 2)
        return;

    /***************************
             L = 3
    ***************************/
    rly[9]
        = ylmcoef[7] * z * rly[4] - ylmcoef[8] * rly[1] * radius2; // l=3, m=0
    grly[9][0] = ylmcoef[7] * z * grly[4][0]
                 - ylmcoef[8] * (rly[1] * tx + grly[1][0] * radius2);
    grly[9][1] = ylmcoef[7] * z * grly[4][1]
                 - ylmcoef[8] * (rly[1] * ty + grly[1][1] * radius2);
    grly[9][2] = ylmcoef[7] * (rly[4] + z * grly[4][2])
                 - ylmcoef[8] * (rly[1] * tz + grly[1][2] * radius2);

    double tmp3 = ylmcoef[9] * z;
    rly[10] = tmp3 * rly[5] - ylmcoef[10] * rly[2] * radius2; // l=3,m=1
    grly[10][0] = tmp3 * grly[5][0]
                  - ylmcoef[10] * (grly[2][0] * radius2 + rly[2] * tx);
    grly[10][1] = tmp3 * grly[5][1]
                  - ylmcoef[10] * (grly[2][1] * radius2 + rly[2] * ty);
    grly[10][2] = ylmcoef[9] * (z * grly[5][2] + rly[5])
                  - ylmcoef[10] * (grly[2][2] * radius2 + rly[2] * tz);

    rly[11] = tmp3 * rly[6] - ylmcoef[10] * rly[3] * radius2; // l=3,m=-1
    grly[11][0] = tmp3 * grly[6][0]
                  - ylmcoef[10] * (grly[3][0] * radius2 + rly[3] * tx);
    grly[11][1] = tmp3 * grly[6][1]
                  - ylmcoef[10] * (grly[3][1] * radius2 + rly[3] * ty);
    grly[11][2] = ylmcoef[9] * (z * grly[6][2] + rly[6])
                  - ylmcoef[10] * (grly[3][2] * radius2 + rly[3] * tz);

    double tmp4 = ylmcoef[11] * z;
    rly[12] = tmp4 * rly[7]; // l=3,m=2
    grly[12][0] = tmp4 * grly[7][0];
    grly[12][1] = tmp4 * grly[7][1];
    grly[12][2] = ylmcoef[11] * (z * grly[7][2] + rly[7]);

    rly[13] = tmp4 * rly[8]; // l=3,m=-2
    grly[13][0] = tmp4 * grly[8][0];
    grly[13][1] = tmp4 * grly[8][1];
    grly[13][2] = ylmcoef[11] * (z * grly[8][2] + rly[8]);

    double tmp5 = ylmcoef[14] * x;
    rly[14] = ylmcoef[12] * rly[10] - ylmcoef[13] * rly[2] * radius2
              - tmp5 * rly[7]; // l=3,m=3
    grly[14][0] = ylmcoef[12] * grly[10][0]
                  - ylmcoef[13] * (rly[2] * tx + grly[2][0] * radius2)
                  - ylmcoef[14] * (rly[7] + x * grly[7][0]);
    grly[14][1] = ylmcoef[12] * grly[10][1]
                  - ylmcoef[13] * (rly[2] * ty + grly[2][1] * radius2)
                  - tmp5 * grly[7][1];
    grly[14][2] = ylmcoef[12] * grly[10][2]
                  - ylmcoef[13] * (rly[2] * tz + grly[2][2] * radius2)
                  - tmp5 * grly[7][2];

    rly[15] = ylmcoef[12] * rly[11] - ylmcoef[13] * rly[3] * radius2
              - tmp5 * rly[8]; // l=3,m=-3
    grly[15][0] = ylmcoef[12] * grly[11][0]
                  - ylmcoef[13] * (rly[3] * tx + grly[3][0] * radius2)
                  - ylmcoef[14] * (rly[8] + x * grly[8][0]);
    grly[15][1] = ylmcoef[12] * grly[11][1]
                  - ylmcoef[13] * (rly[3] * ty + grly[3][1] * radius2)
                  - tmp5 * grly[8][1];
    grly[15][2] = ylmcoef[12] * grly[11][2]
                  - ylmcoef[13] * (rly[3] * tz + grly[3][2] * radius2)
                  - tmp5 * grly[8][2];
    if (Lmax == 3)
        return;

    /***************************
             L = 4
    ***************************/
    rly[16]
        = ylmcoef[15] * z * rly[9] - ylmcoef[16] * rly[4] * radius2; // l=4,m=0
    grly[16][0] = ylmcoef[15] * z * grly[9][0]
                  - ylmcoef[16] * (rly[4] * tx + grly[4][0] * radius2);
    grly[16][1] = ylmcoef[15] * z * grly[9][1]
                  - ylmcoef[16] * (rly[4] * ty + grly[4][1] * radius2);
    grly[16][2] = ylmcoef[15] * (z * grly[9][2] + rly[9])
                  - ylmcoef[16] * (rly[4] * tz + grly[4][2] * radius2);

    double tmp6 = ylmcoef[17] * z;
    rly[17] = tmp6 * rly[10] - ylmcoef[18] * rly[5] * radius2; // l=4,m=1
    grly[17][0] = tmp6 * grly[10][0]
                  - ylmcoef[18] * (rly[5] * tx + grly[5][0] * radius2);
    grly[17][1] = tmp6 * grly[10][1]
                  - ylmcoef[18] * (rly[5] * ty + grly[5][1] * radius2);
    grly[17][2] = ylmcoef[17] * (z * grly[10][2] + rly[10])
                  - ylmcoef[18] * (rly[5] * tz + grly[5][2] * radius2);

    rly[18] = tmp6 * rly[11] - ylmcoef[18] * rly[6] * radius2; // l=4,m=-1
    grly[18][0] = tmp6 * grly[11][0]
                  - ylmcoef[18] * (rly[6] * tx + grly[6][0] * radius2);
    grly[18][1] = tmp6 * grly[11][1]
                  - ylmcoef[18] * (rly[6] * ty + grly[6][1] * radius2);
    grly[18][2] = ylmcoef[17] * (z * grly[11][2] + rly[11])
                  - ylmcoef[18] * (rly[6] * tz + grly[6][2] * radius2);

    double tmp7 = ylmcoef[19] * z;
    rly[19] = tmp7 * rly[12] - ylmcoef[20] * rly[7] * radius2; // l=4,m=2
    grly[19][0] = tmp7 * grly[12][0]
                  - ylmcoef[20] * (rly[7] * tx + grly[7][0] * radius2);
    grly[19][1] = tmp7 * grly[12][1]
                  - ylmcoef[20] * (rly[7] * ty + grly[7][1] * radius2);
    grly[19][2] = ylmcoef[19] * (z * grly[12][2] + rly[12])
                  - ylmcoef[20] * (rly[7] * tz + grly[7][2] * radius2);

    rly[20] = tmp7 * rly[13] - ylmcoef[20] * rly[8] * radius2; // l=4,m=-2
    grly[20][0] = tmp7 * grly[13][0]
                  - ylmcoef[20] * (rly[8] * tx + grly[8][0] * radius2);
    grly[20][1] = tmp7 * grly[13][1]
                  - ylmcoef[20] * (rly[8] * ty + grly[8][1] * radius2);
    grly[20][2] = ylmcoef[19] * (z * grly[13][2] + rly[13])
                  - ylmcoef[20] * (rly[8] * tz + grly[8][2] * radius2);

    double tmp8 = 3.0 * z;
    rly[21] = tmp8 * rly[14]; // l=4,m=3
    grly[21][0] = tmp8 * grly[14][0];
    grly[21][1] = tmp8 * grly[14][1];
    grly[21][2] = 3.0 * (z * grly[14][2] + rly[14]);

    rly[22] = tmp8 * rly[15]; // l=4,m=-3
    grly[22][0] = tmp8 * grly[15][0];
    grly[22][1] = tmp8 * grly[15][1];
    grly[22][2] = 3.0 * (z * grly[15][2] + rly[15]);

    double tmp9 = ylmcoef[23] * x;
    rly[23] = ylmcoef[21] * rly[19] - ylmcoef[22] * rly[7] * radius2
              - tmp9 * rly[14]; // l=4,m=4
    grly[23][0] = ylmcoef[21] * grly[19][0]
                  - ylmcoef[22] * (rly[7] * tx + grly[7][0] * radius2)
                  - ylmcoef[23] * (x * grly[14][0] + rly[14]);
    grly[23][1] = ylmcoef[21] * grly[19][1]
                  - ylmcoef[22] * (rly[7] * ty + grly[7][1] * radius2)
                  - tmp9 * grly[14][1];
    grly[23][2] = ylmcoef[21] * grly[19][2]
                  - ylmcoef[22] * (rly[7] * tz + grly[7][2] * radius2)
                  - tmp9 * grly[14][2];

    rly[24] = ylmcoef[21] * rly[20] - ylmcoef[22] * rly[8] * radius2
              - tmp9 * rly[15]; // l=4,m=-4
    grly[24][0] = ylmcoef[21] * grly[20][0]
                  - ylmcoef[22] * (rly[8] * tx + grly[8][0] * radius2)
                  - ylmcoef[23] * (x * grly[15][0] + rly[15]);
    grly[24][1] = ylmcoef[21] * grly[20][1]
                  - ylmcoef[22] * (rly[8] * ty + grly[8][1] * radius2)
                  - tmp9 * grly[15][1];
    grly[24][2] = ylmcoef[21] * grly[20][2]
                  - ylmcoef[22] * (rly[8] * tz + grly[8][2] * radius2)
                  - tmp9 * grly[15][2];

    if (Lmax == 4)
        return;

    /***************************
             L = 5
    ***************************/
    rly[25]
        = ylmcoef[24] * z * rly[16] - ylmcoef[25] * rly[9] * radius2; // l=5,m=0
    grly[25][0] = ylmcoef[24] * z * grly[16][0]
                  - ylmcoef[25] * (rly[9] * tx + grly[9][0] * radius2);
    grly[25][1] = ylmcoef[24] * z * grly[16][1]
                  - ylmcoef[25] * (rly[9] * ty + grly[9][1] * radius2);
    grly[25][2] = ylmcoef[24] * (z * grly[16][2] + rly[16])
                  - ylmcoef[25] * (rly[9] * tz + grly[9][2] * radius2);

    double tmp10 = ylmcoef[26] * z;
    rly[26] = tmp10 * rly[17] - ylmcoef[27] * rly[10] * radius2; // l=5,m=1
    grly[26][0] = tmp10 * grly[17][0]
                  - ylmcoef[27] * (rly[10] * tx + grly[10][0] * radius2);
    grly[26][1] = tmp10 * grly[17][1]
                  - ylmcoef[27] * (rly[10] * ty + grly[10][1] * radius2);
    grly[26][2] = ylmcoef[26] * (z * grly[17][2] + rly[17])
                  - ylmcoef[27] * (rly[10] * tz + grly[10][2] * radius2);

    rly[27] = tmp10 * rly[18] - ylmcoef[27] * rly[11] * radius2; // l=5,m=-1
    grly[27][0] = tmp10 * grly[18][0]
                  - ylmcoef[27] * (rly[11] * tx + grly[11][0] * radius2);
    grly[27][1] = tmp10 * grly[18][1]
                  - ylmcoef[27] * (rly[11] * ty + grly[11][1] * radius2);
    grly[27][2] = ylmcoef[26] * (z * grly[18][2] + rly[18])
                  - ylmcoef[27] * (rly[11] * tz + grly[11][2] * radius2);

    double tmp11 = ylmcoef[28] * z;
    rly[28] = tmp11 * rly[19] - ylmcoef[29] * rly[12] * radius2; // l=5,m=2
    grly[28][0] = tmp11 * grly[19][0]
                  - ylmcoef[29] * (rly[12] * tx + grly[12][0] * radius2);
    grly[28][1] = tmp11 * grly[19][1]
                  - ylmcoef[29] * (rly[12] * ty + grly[12][1] * radius2);
    grly[28][2] = ylmcoef[28] * (z * grly[19][2] + rly[19])
                  - ylmcoef[29] * (rly[12] * tz + grly[12][2] * radius2);

    rly[29] = tmp11 * rly[20] - ylmcoef[29] * rly[13] * radius2; // l=5,m=-2
    grly[29][0] = tmp11 * grly[20][0]
                  - ylmcoef[29] * (rly[13] * tx + grly[13][0] * radius2);
    grly[29][1] = tmp11 * grly[20][1]
                  - ylmcoef[29] * (rly[13] * ty + grly[13][1] * radius2);
    grly[29][2] = ylmcoef[28] * (z * grly[20][2] + rly[20])
                  - ylmcoef[29] * (rly[13] * tz + grly[13][2] * radius2);

    double tmp12 = ylmcoef[30] * z;
    rly[30] = tmp12 * rly[21] - ylmcoef[31] * rly[14] * radius2; // l=5,m=3
    grly[30][0] = tmp12 * grly[21][0]
                  - ylmcoef[31] * (grly[14][0] * radius2 + rly[14] * tx);
    grly[30][1] = tmp12 * grly[21][1]
                  - ylmcoef[31] * (grly[14][1] * radius2 + rly[14] * ty);
    grly[30][2] = ylmcoef[30] * (z * grly[21][2] + rly[21])
                  - ylmcoef[31] * (grly[14][2] * radius2 + rly[14] * tz);

    rly[31] = tmp12 * rly[22] - ylmcoef[31] * rly[15] * radius2; // l=5,m=-3
    grly[31][0] = tmp12 * grly[22][0]
                  - ylmcoef[31] * (grly[15][0] * radius2 + rly[15] * tx);
    grly[31][1] = tmp12 * grly[22][1]
                  - ylmcoef[31] * (grly[15][1] * radius2 + rly[15] * ty);
    grly[31][2] = ylmcoef[30] * (z * grly[22][2] + rly[22])
                  - ylmcoef[31] * (grly[15][2] * radius2 + rly[15] * tz);

    double tmp13 = ylmcoef[32] * z;
    rly[32] = tmp13 * rly[23]; // l=5,m=4
    grly[32][0] = tmp13 * grly[23][0];
    grly[32][1] = tmp13 * grly[23][1];
    grly[32][2] = ylmcoef[32] * (rly[23] + z * grly[23][2]);

    rly[33] = tmp13 * rly[24]; // l=5,m=-4
    grly[33][0] = tmp13 * grly[24][0];
    grly[33][1] = tmp13 * grly[24][1];
    grly[33][2] = ylmcoef[32] * (rly[24] + z * grly[24][2]);

    double tmp14 = ylmcoef[35] * x;
    rly[34] = ylmcoef[33] * rly[30] - ylmcoef[34] * rly[14] * radius2
              - tmp14 * rly[23]; // l=5,m=5
    grly[34][0] = ylmcoef[33] * grly[30][0]
                  - ylmcoef[34] * (rly[14] * tx + grly[14][0] * radius2)
                  - ylmcoef[35] * (x * grly[23][0] + rly[23]);
    grly[34][1] = ylmcoef[33] * grly[30][1]
                  - ylmcoef[34] * (rly[14] * ty + grly[14][1] * radius2)
                  - tmp14 * grly[23][1];
    grly[34][2] = ylmcoef[33] * grly[30][2]
                  - ylmcoef[34] * (rly[14] * tz + grly[14][2] * radius2)
                  - tmp14 * grly[23][2];

    rly[35] = ylmcoef[33] * rly[31] - ylmcoef[34] * rly[15] * radius2
              - tmp14 * rly[24]; // l=5,m=-5
    grly[35][0] = ylmcoef[33] * grly[31][0]
                  - ylmcoef[34] * (rly[15] * tx + grly[15][0] * radius2)
                  - ylmcoef[35] * (x * grly[24][0] + rly[24]);
    grly[35][1] = ylmcoef[33] * grly[31][1]
                  - ylmcoef[34] * (rly[15] * ty + grly[15][1] * radius2)
                  - tmp14 * grly[24][1];
    grly[35][2] = ylmcoef[33] * grly[31][2]
                  - ylmcoef[34] * (rly[15] * tz + grly[15][2] * radius2)
                  - tmp14 * grly[24][2];

    if (Lmax == 5)
        return;

    // if Lmax > 5
    for (int il = 6; il <= Lmax; il++)
    {
        int istart = il * il;
        int istart1 = (il - 1) * (il - 1);
        int istart2 = (il - 2) * (il - 2);

        double fac2 = sqrt(4.0 * istart - 1.0);
        double fac4 = sqrt(4.0 * istart1 - 1.0);

        for (int im = 0; im < 2 * il - 1; im++)
        {
            int imm = (im + 1) / 2;
            //			if (im % 2 == 0) imm *= -1;

            double var1 = fac2 / sqrt((double)istart - imm * imm);
            double var2 = sqrt((double)istart1 - imm * imm) / fac4;

            rly[istart + im] = var1
                               * (z * rly[istart1 + im]
                                  - var2 * rly[istart2 + im] * radius2);

            grly[istart + im][0]
                = var1
                  * (z * grly[istart1 + im][0]
                     - var2
                           * (rly[istart2 + im] * tx
                              + grly[istart2 + im][0] * radius2));
            grly[istart + im][1]
                = var1
                  * (z * grly[istart1 + im][1]
                     - var2
                           * (rly[istart2 + im] * ty
                              + grly[istart2 + im][1] * radius2));
            grly[istart + im][2]
                = var1
                  * (z * grly[istart1 + im][2] + rly[istart1 + im]
                     - var2
                           * (rly[istart2 + im] * tz
                              + grly[istart2 + im][2] * radius2));
        }

        double bl1 = sqrt(2.0 * il / (2.0 * il + 1.0));
        double bl2 = sqrt((2.0 * il - 2.0) / (2.0 * il - 1.0));
        double bl3 = sqrt(2.0) / fac2;

        int id1 = istart + 2 * il - 1;
        int id2 = istart + 2 * il - 5;
        int id3 = istart2 + 2 * il - 5;
        int id4 = istart1 + 2 * il - 3;

        rly[id1]
            = (bl3 * rly[id2] - bl2 * rly[id3] * radius2 - 2.0 * x * rly[id4])
              / bl1;
        grly[id1][0] = (bl3 * grly[id2][0]
                        - bl2 * (grly[id3][0] * radius2 + rly[id3] * tx)
                        - 2.0 * (rly[id4] + x * grly[id4][0]))
                       / bl1;
        grly[id1][1] = (bl3 * grly[id2][1]
                        - bl2 * (grly[id3][1] * radius2 + rly[id3] * ty)
                        - 2.0 * x * grly[id4][1])
                       / bl1;
        grly[id1][2] = (bl3 * grly[id2][2]
                        - bl2 * (grly[id3][2] * radius2 + rly[id3] * tz)
                        - 2.0 * x * grly[id4][2])
                       / bl1;

        rly[id1 + 1] = (bl3 * rly[id2 + 1] - bl2 * rly[id3 + 1] * radius2
                        - 2.0 * x * rly[id4 + 1])
                       / bl1;
        grly[id1 + 1][0]
            = (bl3 * grly[id2 + 1][0]
               - bl2 * (grly[id3 + 1][0] * radius2 + rly[id3 + 1] * tx)
               - 2.0 * (rly[id4 + 1] + x * grly[id4 + 1][0]))
              / bl1;
        grly[id1 + 1][1]
            = (bl3 * grly[id2 + 1][1]
               - bl2 * (grly[id3 + 1][1] * radius2 + rly[id3 + 1] * ty)
               - 2.0 * x * grly[id4 + 1][1])
              / bl1;
        grly[id1 + 1][2]
            = (bl3 * grly[id2 + 1][2]
               - bl2 * (grly[id3 + 1][2] * radius2 + rly[id3 + 1] * tz)
               - 2.0 * x * grly[id4 + 1][2])
              / bl1;
    }

    return;
}