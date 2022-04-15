#ifndef EXX_ABFS_DM_H
#define EXX_ABFS_DM_H

#include "exx_abfs.h"
#include "abfs-vector3_order.h"

#include <vector>
#include <map>
#include <deque>

class Charge_Broyden;

class Exx_Abfs::DM
{
public:
	std::map<size_t,std::map<size_t,std::vector<ModuleBase::ComplexMatrix>>> DMk;								// DMk[iat1][iat2][ik](iw1,iw2)
	std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::vector<ModuleBase::matrix>>>> DMr;		// DMr[iat1][iat2][box][is](iw1,iw2)
	std::deque<std::map<size_t,std::map<size_t,std::vector<ModuleBase::ComplexMatrix>>>> DMk_pulay_seq;				// DMk_pulay_seq[istep][iat1][iat2][ik](iw1,iw2)
//	double DM_threshold = std::numeric_limits<double>::min();
	bool flag_mix;

public:
	void cal_DM(
		const std::set<std::pair<size_t,size_t>> &atom_pairs,
        const std::vector<Abfs::Vector3_Order<int>>& Born_von_Karman_boxes,
        std::complex<double>*** wfc_k_grid);
    std::map<size_t, std::map<size_t, std::vector<ModuleBase::ComplexMatrix>>> cal_DMk_raw(const std::set<std::pair<size_t, size_t>>& atom_pairs, std::complex<double>*** wfc_k_grid) const;
		
private:
	void cal_DMk_mixing(
		const Charge_Broyden &charge,
        const std::set<std::pair<size_t, size_t>>& atom_pairs,
        std::complex<double>*** wfc_k_grid);
    void plain_mixing(
		const Charge_Broyden &charge,
        const std::set<std::pair<size_t, size_t>>& atom_pairs,
        std::complex<double>*** wfc_k_grid);
    void pulay_mixing(
		const Charge_Broyden &charge,
        const std::set<std::pair<size_t, size_t>>& atom_pairs,
        complex<double>*** wfc_k_grid);

//	double cal_DM_delta();

//	void   set_DM_threshold(double DM_threshold_in) { DM_threshold = DM_threshold_in; }
//	double get_DM_threshold() const { return DM_threshold; }
};

#endif	// EXX_ABFS_DM_H