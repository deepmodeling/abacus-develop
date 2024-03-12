///klist : adapted from klist from module_hamilt_pw/hamilt_pwdft
///deals with k point sampling

#include "module_base/vector3.h"
#include "module_base/matrix3.h"
#include "module_base/memory.h"
#include "module_base/global_function.h"
#include <iostream>
#include <fstream>

namespace Test_Deepks
{
class K_Vectors
{
public:

    ModuleBase::Vector3<double> *kvec_c;		// Cartesian coordinates of k points
    std::vector<ModuleBase::Vector3<double>> kvec_d;		// Direct coordinates of k points

    double *wk;						// wk, weight of k points

    int *isk;						// distinguish spin up and down k points

    int nkstot;						// total number of k points

    int nmp[3];						// Number of Monhorst-Pack

    K_Vectors();
    ~K_Vectors();

    void set(
        const std::string &k_file_name,
        const int& nspin,
        const ModuleBase::Matrix3 &reciprocal_vec,
        const ModuleBase::Matrix3 &latvec,
		bool &GAMMA_ONLY_LOCAL,
		std::ofstream &ofs_running,
		std::ofstream &ofs_warning);

private:
    int nspin;
    bool kc_done;
    bool kd_done;
    double koffset[3];     			// used only in automatic k-points.
    std::string k_kword; //LiuXh add 20180619
    int k_nkstot; //LiuXh add 20180619

    // step 1 : generate kpoints
    bool read_kpoints(
		const std::string &fn,
		bool &GAMMA_ONLY_LOCAL,
		std::ofstream &ofs_warning,
		std::ofstream &ofs_running);
    void Monkhorst_Pack(const int *nmp_in,const double *koffset_in,const int tipo);
    double Monkhorst_Pack_formula( const int &k_type, const double &offset,
                                   const int& n, const int &dim);

    // step 2 : set both kvec and kved; normalize weight
    void set_both_kvec(const ModuleBase::Matrix3 &G,const ModuleBase::Matrix3 &Rm, std::ofstream &ofs_running);
	void renew(const int &kpoint_number);
    void normalize_wk( const int &degspin );

    // step 3 : *2 or *4 kpoints.
    // *2 for LSDA
    // *4 for non-collinear
    void set_kup_and_kdw(std::ofstream &ofs_running);

    // step 4
    // print k lists.
    void print_klists(std::ofstream &fn_running);
};
}
