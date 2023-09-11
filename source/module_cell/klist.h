#ifndef K_VECTORS_H
#define K_VECTORS_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix3.h"
#include "module_cell/module_symmetry/symmetry.h"
#include <vector>

class K_Vectors
{
public:

    std::vector<ModuleBase::Vector3<double>> kvec_c;		// Cartesian coordinates of k points
    std::vector<ModuleBase::Vector3<double>> kvec_d;		// Direct coordinates of k points
    std::vector<ModuleBase::Vector3<double>> kvec_d_ibz;	// ibz Direct coordinates of k points

    std::vector<double> wk;						// wk, weight of k points
    std::vector<double> wk_ibz;					// ibz kpoint wk ,weight of k points

    std::vector<int> ngk;						// ngk, number of plane waves for each k point
    std::vector<int> isk;						// distinguish spin up and down k points
    std::vector<int> ibz2bz;					// mohan added 2009-05-18

    int nks;						// number of k points in this pool(processor, up+dw)
    int nkstot;						/// total number of k points, equal to nkstot_ibz after reducing k points
    int nkstot_ibz;             /// number of k points in IBZ
    int nkstot_full;    /// number of k points in full k mesh

    int nmp[3];						// Number of Monhorst-Pack

    K_Vectors();
    ~K_Vectors();

    void set(
        const ModuleSymmetry::Symmetry &symm,
        const std::string &k_file_name,
        const int& nspin,
        const ModuleBase::Matrix3 &reciprocal_vec,
        const ModuleBase::Matrix3 &latvec);

    void ibz_kpoint(const ModuleSymmetry::Symmetry &symm, bool use_symm,std::string& skpt, const UnitCell &ucell, bool& match);
    //LiuXh add 20180515
    void set_after_vc(
            const ModuleSymmetry::Symmetry &symm,
            const std::string &k_file_name,
            const int& nspin,
            const ModuleBase::Matrix3 &reciprocal_vec,
            const ModuleBase::Matrix3 &latvec);
    //get global index for ik
    inline int getik_global(const int& ik) const;

private:
    int nspin;
    bool kc_done;
    bool kd_done;
    double koffset[3];     			// used only in automatic k-points.
    std::string k_kword; //LiuXh add 20180619
    int k_nkstot; //LiuXh add 20180619
    bool is_mp = false; //Monkhorst-Pack
    
    void renew( const int &kpoint_number );

    // step 1 : generate kpoints
    bool read_kpoints(const std::string &fn); // return 0: something wrong.
    void Monkhorst_Pack(const int *nmp_in,const double *koffset_in,const int tipo);
    double Monkhorst_Pack_formula( const int &k_type, const double &offset,
                                   const int& n, const int &dim);

    // step 2 : set both kvec and kved; normalize weight
    void update_use_ibz( void );
    void set_both_kvec(const ModuleBase::Matrix3 &G,const ModuleBase::Matrix3 &R, std::string& skpt);
    void normalize_wk( const int &degspin );

    // step 3 : mpi kpoints information.
    void mpi_k();

    // step 4 : *2 or *4 kpoints.
    // *2 for LSDA
    // *4 for non-collinear
    void set_kup_and_kdw();

    // step 5
    // print k lists.
    void print_klists(std::ofstream &fn);
    //bool read_kpoints_after_vc(const std::string &fn); //LiuXh add 20180515
    //void Monkhorst_Pack_after_vc(const int *nmp_in,const double *koffset_in,const int tipo); //LiuXh add 20180515
    void mpi_k_after_vc(); //LiuXh add 20180515
    void set_both_kvec_after_vc(const ModuleBase::Matrix3 &G,const ModuleBase::Matrix3 &R);
    void set_kup_and_kdw_after_vc();
};

inline int K_Vectors:: getik_global(const int& ik) const
{
    int nkp = this->nkstot / GlobalV::KPAR;
    int rem = this->nkstot % GlobalV::KPAR;
    if(GlobalV::MY_POOL < rem)
    {
        return GlobalV::MY_POOL*nkp + GlobalV::MY_POOL + ik;
    }
    else
    {
        return GlobalV::MY_POOL*nkp + rem + ik;       
    }
}

#endif // KVECT_H
