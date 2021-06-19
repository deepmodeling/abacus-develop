//==========================================================
// AUTHOR : Linxin He, Kaimin Duan, Wei Fang, Mohan Chen
// START: 2006
// LAST UPDATE: 2009-05-18
//==========================================================
#ifndef K_VECTORS_H
#define K_VECTORS_H

#include "tools.h"
#include "symmetry.h"

class K_Vectors
{
public:

    Vector3<double> *kvec_c;		// Cartesian coordinates of k points
    Vector3<double> *kvec_d;		// Direct coordinates of k points
    Vector3<double> *kvec_d_ibz;	// ibz Direct coordinates of k points

    double *wk;						// wk, weight of k points
    double *wk_ibz;					// ibz kpoint wk ,weight of k points

    int *ngk;						// ngk, number of plane waves for each k point
    int *isk;						// distinguish spin up and down k points
    int *ibz2bz;					// mohan added 2009-05-18

    int nks;						// number of k points in this pool(processor, up+dw)
    int nkstot;						// total number of k points
    int nkstot_ibz;

    int nmp[3];						// Number of Monhorst-Pack

    K_Vectors();
    ~K_Vectors();

    void set(
        const Symmetry &symm,
        const string &k_file_name,
        const int& nspin,
        const Matrix3 &reciprocal_vec,
        const Matrix3 &latvec);

    void ibz_kpoint( const Symmetry &symm);
    //LiuXh add 20180515
    void set_after_vc(
            const Symmetry &symm,
            const string &k_file_name,
            const int& nspin,
            const Matrix3 &reciprocal_vec,
            const Matrix3 &latvec);

private:
    int nspin;
    bool kc_done;
    bool kd_done;
    double koffset[3];     			// used only in automatic k-points.
    string k_kword; //LiuXh add 20180619
    int k_nkstot; //LiuXh add 20180619

    void renew( const int &kpoint_number );

    // step 1 : generate kpoints
    bool read_kpoints(const string &fn); // return 0: something wrong.
    void Monkhorst_Pack(const int *nmp_in,const double *koffset_in,const int tipo);
    double Monkhorst_Pack_formula( const int &k_type, const double &offset,
                                   const int& n, const int &dim);

    // step 2 : set both kvec and kved; normalize weight
    void update_use_ibz( void );
    void set_both_kvec(const Matrix3 &G,const Matrix3 &R);
    void normalize_wk( const int &degspin );

    // step 3 : mpi kpoints information.
    void mpi_k();

    // step 4 : *2 or *4 kpoints.
    // *2 for LSDA
    // *4 for non-collinear
    void set_kup_and_kdw();

    // step 5
    // print k lists.
    void print_klists(ofstream &fn);
    //bool read_kpoints_after_vc(const string &fn); //LiuXh add 20180515
    //void Monkhorst_Pack_after_vc(const int *nmp_in,const double *koffset_in,const int tipo); //LiuXh add 20180515
    void mpi_k_after_vc(); //LiuXh add 20180515
    void set_both_kvec_after_vc(const Matrix3 &G,const Matrix3 &R);
    void set_kup_and_kdw_after_vc();
};

#endif // KVECT_H
