#ifndef MRRR_H_
#define MRRR_H_
#ifdef __MPI
#include "mpi.h"
#endif
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>

//typedef int bool;
//#define true 1
//#define false 0
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)
#define dlamch dlamch_

extern "C" double dlamch(char *);
/*
 extern int pdstemr(char *jobz, char *range, int *n, double * d__, double *e,
 double *vl, double *vu, int *il, int *iu, int *m, double *w,
 double *z__, int *ldz, int *nzc, int *isuppz, bool *tryrac,
 double *work, int *lwork, int *iwork, int *liwork, int *info);
extern int pdstemr_mpi(MPI_Comm comm1D, char *jobz, char *range, int *n,
		double * d__, double *e, double *vl, double *vu, int *il, int *iu,
		int *m, double *w, double *z__, int *ldz, int *nzc, int *isuppz,
		bool *tryrac, double *work, int *lwork, int *iwork, int *liwork,
		int *info);
*/
extern void merge(double *x, double *xe, double *y, double *ye, double *z,
		int *a, int *b, int *a1);
extern int sort_w(double *w, int n, int il, int iu, int nsplit, int *blk_begins,
		int *blk_sizes, int *a);
extern int picopy(int *n, int *dx, int *incx, int *dy, int *incy);
extern int pdcopy(int *n, double *dx, int *incx, double *dy, int *incy);
extern int pdlarnv(int *idist, int *iseed, int *n, double *x);
extern int pdlar1v(int *n, int *b1, int *bn, double *lambda, double *d__,
		double *l, double *ld, double * lld, double *pivmin, double *gaptol,
		double *z__, bool *wantnc, int *negcnt, double *ztz, double *mingma,
		int *r__, int *isuppz, double *nrminv, double *resid, double *rqcorr,
		double *work);
extern int pdlaneg(int *n, double *d__, double *lld, double * sigma,
		double *pivmin, int *r__);
extern double pdlanst(char *norm, int *n, double *d__, double *e);
extern int pdlaev2(double *a, double *b, double *c__, double *rt1, double *rt2,
		double *cs1, double *sn1);
extern int pdlaebz(int *ijob, int *nitmax, int *n, int *mmax, int *minp,
		int *nbmin, double *abstol, double *reltol, double *pivmin, double *d__,
		double * e, double *e2, int *nval, double *ab, double *c__, int *mout,
		int *nab, double *work, int *iwork, int *info);
extern int pdlae2(double *a, double *b, double *c__, double *rt1, double *rt2);
extern int pdlarra(int *n, double *d__, double *e, double *e2, double *spltol,
		double *tnrm, int *nsplit, int *isplit, int *info);
extern int pdlarrb(int *n, double *d__, double *lld, int *ifirst, int *ilast,
		double *rtol1, double *rtol2, int *offset, double *w, double *wgap,
		double *werr, double *work, int *iwork, double *pivmin, double * spdiam,
		int *twist, int *info);
extern int pdlarrc(char *jobt, int *n, double *vl, double *vu, double *d__,
		double *e, double *pivmin, int *eigcnt, int *lcnt, int *rcnt, int *info);
extern int pdlarrd(char *range, char *order, int *n, double *vl, double *vu,
		int *il, int *iu, double *gers, double *reltol, double *d__, double *e,
		double *e2, double *pivmin, int *nsplit, int *isplit, int *m, double *w,
		double *werr, double *wl, double *wu, int *iblock, int *indexw,
		double *work, int *iwork, int *info);
extern int pdlarre(char *range, int *n, double *vl, double *vu, int *il,
		int *iu, double *d__, double *e, double *e2, double *rtol1,
		double *rtol2, double * spltol, int *nsplit, int *isplit, int *m,
		double *w, double *werr, double *wgap, int *iblock, int *indexw,
		double *gers, double *pivmin, double *work, int * iwork, int *info);
extern int pdlarrf(int *n, double *d__, double *l, double *ld, int *clstrt,
		int *clend, double *w, double *wgap, double *werr, double *spdiam,
		double * clgapl, double *clgapr, double *pivmin, double *sigma,
		double *dplus, double *lplus, double *work, int *info);
extern int pdlarrj(int *n, double *d__, double *e2, int *ifirst, int *ilast,
		double *rtol, int *offset, double *w, double *werr, double *work,
		int *iwork, double *pivmin, double *spdiam, int *info);
extern int pdlarrk(int *n, int *iw, double *gl, double *gu, double *d__,
		double *e2, double *pivmin, double *reltol, double *w, double *werr,
		int *info);
extern int pdlarrr(int *n, double *d__, double *e, int *info);
extern int pdlarrv(int *n, double *vl, double *vu, double *d__, double *l,
		double *pivmin, int *isplit, int *m, int *dol, int *dou, double *minrgp,
		double *rtol1, double *rtol2, double *w, double *werr, double *wgap,
		int *iblock, int *indexw, double *gers, double *z__, int *ldz,
		int *isuppz, double *work, int *iwork, int *info);
int pdlarrv_mpi(int *il, int *iu, int *n, double *vl, double *vu, double *d__,
		double *l, double *pivmin, int *isplit, int *m, int *dol, int *dou,
		double *minrgp, double *rtol1, double *rtol2, double *w, double *werr,
		double *wgap, int *iblock, int *indexw, double *gers, double *z__,
		int *ldz, int *isuppz, double *work, int *iwork, int *info);
extern int pdlaruv(int *iseed, int *n, double *x);
extern int pdlas2(double *f, double *g, double *h__, double *ssmin,
		double *ssmax);
extern int pdlascl(char *type__, int *kl, int *ku, double *cfrom, double *cto,
		int *m, int *n, double *a, int *lda, int *info);
extern int pdlaset(char *uplo, int *m, int *n, double * alpha, double *beta,
		double *a, int *lda);
extern int pdlasq2(int *n, double *z__, int *info);
extern int pdlasq3(int *i0, int *n0, double *z__, int *pp, double *dmin__,
		double *sigma, double *desig, double *qmax, int *nfail, int *iter,
		int *ndiv, bool *ieee, int *ttype, double *dmin1, double *dmin2,
		double *dn, double *dn1, double *dn2, double *g, double *tau);
extern int pdlasq4(int *i0, int *n0, double *z__, int *pp, int *n0in,
		double *dmin__, double *dmin1, double *dmin2, double *dn, double *dn1,
		double *dn2, double *tau, int *ttype, double *g);
extern int pdlasq5(int *i0, int *n0, double *z__, int *pp, double *tau,
		double *dmin__, double *dmin1, double *dmin2, double *dn, double *dnm1,
		double *dnm2, bool *ieee);
extern int pdlasq6(int *i0, int *n0, double *z__, int *pp, double *dmin__,
		double *dmin1, double *dmin2, double *dn, double *dnm1, double *dnm2);
extern int pdlasrt(char *id, int *n, double *d__, int * info);
extern int pdlassq(int *n, double *x, int *incx, double *scale, double *sumsq);
extern int pdscal(int *n, double *da, double *dx, int *incx);
extern int pilaenv(int *ispec, char *name__, char *opts, int *n1, int *n2,
		int *n3, int *n4);
extern int pdswap(int *n, double *dx, int *incx, double *dy, int *incy);
extern int i_nint(double *x);
extern int pieeeck(int *ispec, double *zero, double *one);
extern int piparmq(int *ispec, char *name__, char *opts, int *n, int *ilo,
		int *ihi, int *lwork);
extern bool plsame(char *ca, char *cb);
extern int s_cmp(char *a0, char *b0, int la, int lb);
extern void s_copy(register char *a, register char *b, int la, int lb);
extern int xerbla(char *srname, int *info);
#endif /* MRRR_H_ */
