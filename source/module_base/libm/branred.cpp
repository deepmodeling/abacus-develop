//==========================================================
// AUTHOR : alcanderian@gmail.com
// DATE : 2023-01-06
//==========================================================

#include <math.h>
#include <endian.h>

namespace ModuleBase
{
namespace libm
{

typedef int int4;
typedef union { unsigned int u[2]; int4 i[2]; double x; double d; } mynumber;

#define max(x, y)  (((y) > (x)) ? (y) : (x))
#define min(x, y)  (((y) < (x)) ? (y) : (x))

#if (__BYTE_ORDER == __BIG_ENDIAN)

#define HIGH_HALF 0
#define  LOW_HALF 1

static const mynumber

/**/           t576 = {{0x63f00000, 0x00000000}}, /* 2 ^ 576  */
/**/          tm600 = {{0x1a700000, 0x00000000}}, /* 2 ^- 600 */
/**/           tm24 = {{0x3e700000, 0x00000000}}, /* 2 ^- 24  */
/**/            big = {{0x43380000, 0x00000000}}, /*  6755399441055744      */
/**/           big1 = {{0x43580000, 0x00000000}}, /* 27021597764222976      */
/**/            hp0 = {{0x3FF921FB, 0x54442D18}} ,/* 1.5707963267948966     */
/**/            hp1 = {{0x3C91A626, 0x33145C07}} ,/* 6.123233995736766e-17  */
/**/            mp1 = {{0x3FF921FB, 0x58000000}}, /* 1.5707963407039642     */
/**/            mp2 = {{0xBE4DDE97, 0x40000000}}; /*-1.3909067675399456e-08 */

#endif

#if (__BYTE_ORDER == __LITTLE_ENDIAN)

#define HIGH_HALF 1
#define  LOW_HALF 0

static const mynumber

/**/           t576 = {{0x00000000, 0x63f00000}},  /* 2 ^ 576  */
/**/          tm600 = {{0x00000000, 0x1a700000}},  /* 2 ^- 600 */
/**/           tm24 = {{0x00000000, 0x3e700000}},  /* 2 ^- 24  */
/**/            big = {{0x00000000, 0x43380000}},  /*  6755399441055744      */
/**/           big1 = {{0x00000000, 0x43580000}},  /* 27021597764222976      */
/**/            hp0 = {{0x54442D18, 0x3FF921FB}},  /* 1.5707963267948966     */
/**/            hp1 = {{0x33145C07, 0x3C91A626}},  /* 6.123233995736766e-17  */
/**/            mp1 = {{0x58000000, 0x3FF921FB}},  /* 1.5707963407039642     */
/**/            mp2 = {{0x40000000, 0xBE4DDE97}};  /*-1.3909067675399456e-08 */

#endif

static const double toverp[75] = { /*  2/ PI base 24*/
  10680707.0,  7228996.0,  1387004.0,  2578385.0, 16069853.0,
  12639074.0,  9804092.0,  4427841.0, 16666979.0, 11263675.0,
  12935607.0,  2387514.0,  4345298.0, 14681673.0,  3074569.0,
  13734428.0, 16653803.0,  1880361.0, 10960616.0,  8533493.0,
   3062596.0,  8710556.0,  7349940.0,  6258241.0,  3772886.0,
   3769171.0,  3798172.0,  8675211.0, 12450088.0,  3874808.0,
   9961438.0,   366607.0, 15675153.0,  9132554.0,  7151469.0,
   3571407.0,  2607881.0, 12013382.0,  4155038.0,  6285869.0,
   7677882.0, 13102053.0, 15825725.0,   473591.0,  9065106.0,
  15363067.0,  6271263.0,  9264392.0,  5636912.0,  4652155.0,
   7056368.0, 13614112.0, 10155062.0,  1944035.0,  9527646.0,
  15080200.0,  6658437.0,  6231200.0,  6832269.0, 16767104.0,
   5075751.0,  3212806.0,  1398474.0,  7579849.0,  6349435.0,
  12618859.0,  4703257.0, 12806093.0, 14477321.0,  2786137.0,
  12875403.0,  9837734.0, 14528324.0, 13719321.0,   343717.0 };

/* CN = 1+2**27 = '41a0000002000000' IEEE double format.  Use it to split a
   double for better accuracy.  */
#define  CN   134217729.0
static const double split =  CN;	/* 2^27 + 1 */

/*******************************************************************/
/* Routine  branred() performs range  reduction of a double number */
/* x into Double length number a+aa,such that                      */
/* x=n*pi/2+(a+aa), abs(a+aa)<pi/4, n=0,+-1,+-2,....               */
/* Routine return integer (n mod 4)                                */
/*******************************************************************/
int
__branred(double x, double *a, double *aa)
{
  int i=0,k=0;
  mynumber  u,gor;
  double r[6],s,t,sum,b,bb,sum1,sum2,b1,bb1,b2,bb2,x1,x2,t1,t2;

  x*=tm600.x;
  t=x*split;   /* split x to two numbers */
  x1=t-(t-x);
  x2=x-x1;
  sum=0;
  u.x = x1;
  k = (u.i[HIGH_HALF]>>20)&2047;
  k = (k-450)/24;
  if (k<0)
    k=0;
  gor.x = t576.x;
  gor.i[HIGH_HALF] -= ((k*24)<<20);
  for (i=0;i<6;i++)
    { r[i] = x1*toverp[k+i]*gor.x; gor.x *= tm24.x; }
  for (i=0;i<3;i++) {
    s=(r[i]+big.x)-big.x;
    sum+=s;
    r[i]-=s;
  }
  t=0;
  for (i=0;i<6;i++)
    t+=r[5-i];
  bb=(((((r[0]-t)+r[1])+r[2])+r[3])+r[4])+r[5];
  s=(t+big.x)-big.x;
  sum+=s;
  t-=s;
  b=t+bb;
  bb=(t-b)+bb;
  s=(sum+big1.x)-big1.x;
  sum-=s;
  b1=b;
  bb1=bb;
  sum1=sum;
  sum=0;

  u.x = x2;
  k = (u.i[HIGH_HALF]>>20)&2047;
  k = (k-450)/24;
  if (k<0)
    k=0;
  gor.x = t576.x;
  gor.i[HIGH_HALF] -= ((k*24)<<20);
  for (i=0;i<6;i++)
    { r[i] = x2*toverp[k+i]*gor.x; gor.x *= tm24.x; }
  for (i=0;i<3;i++) {
    s=(r[i]+big.x)-big.x;
    sum+=s;
    r[i]-=s;
  }
  t=0;
  for (i=0;i<6;i++)
    t+=r[5-i];
  bb=(((((r[0]-t)+r[1])+r[2])+r[3])+r[4])+r[5];
  s=(t+big.x)-big.x;
 sum+=s;
 t-=s;
 b=t+bb;
 bb=(t-b)+bb;
 s=(sum+big1.x)-big1.x;
 sum-=s;

 b2=b;
 bb2=bb;
 sum2=sum;

 sum=sum1+sum2;
 b=b1+b2;
 bb = (fabs(b1)>fabs(b2))? (b1-b)+b2 : (b2-b)+b1;
 if (b > 0.5)
   {b-=1.0; sum+=1.0;}
 else if (b < -0.5)
   {b+=1.0; sum-=1.0;}
 s=b+(bb+bb1+bb2);
 t=((b-s)+bb)+(bb1+bb2);
 b=s*split;
 t1=b-(b-s);
 t2=s-t1;
 b=s*hp0.x;
 bb=(((t1*mp1.x-b)+t1*mp2.x)+t2*mp1.x)+(t2*mp2.x+s*hp1.x+t*hp0.x);
 s=b+bb;
 t=(b-s)+bb;
 *a=s;
 *aa=t;
 return ((int) sum)&3; /* return quater of unit circle */
}

};
};
