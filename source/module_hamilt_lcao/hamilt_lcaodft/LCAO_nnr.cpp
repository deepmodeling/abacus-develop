#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "record_adj.h" //mohan add 2012-07-06
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif
#include "module_base/libm/libm.h"
#include <algorithm>

// This is for cell R dependent part. 
void Grid_Technique::cal_nnrg(Parallel_Orbitals* pv, const std::vector<double>& orb_cutoff)
{
	ModuleBase::TITLE("LCAO_nnr","cal_nnrg");

	this->cal_max_box_index();

	this->nnrg = 0;

	nlocdimg.clear();
	nlocstartg.clear();
	nad.clear();
	this->nnrg_index.resize(0);
	
	this->nad = std::vector<int>(GlobalC::ucell.nat, 0);
	this->nlocdimg = std::vector<int>(GlobalC::ucell.nat, 0);
	this->nlocstartg =std::vector<int>(GlobalC::ucell.nat, 0);

	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
	for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
	{
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];

			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);

			const int iat = GlobalC::ucell.itia2iat(T1,I1);

			// for grid integration (on FFT box),
			// we only need to consider <phi_i | phi_j>,
			// which is different from non-local term,
			// which we need to consdier <phi_i|beta_k><beta_k|phi_j>

			// whether this atom is in this processor.
			if(this->in_this_processor[iat])
			{
				// starting index of adjacents.
				this->nlocstartg[iat] = this->nnrg;

				// number of adjacent atoms for atom 'iat'
				this->nad[iat] = 0;

				int count = 0;
				for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
				{
					const int T2 = GlobalC::GridD.getType(ad);
					const int I2 = GlobalC::GridD.getNatom(ad);
					const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
					Atom* atom2 = &GlobalC::ucell.atoms[T2]; 

					// if the adjacent atom is in this processor.
					if(this->in_this_processor[iat2])
					{
						tau2 = GlobalC::GridD.getAdjacentTau(ad);
						dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
						double distance = dtau.norm() * GlobalC::ucell.lat0;
						double rcut = orb_cutoff[T1] + orb_cutoff[T2];


						//if(distance < rcut)
						// mohan reset this 2013-07-02 in Princeton
						// we should make absolutely sure that the distance is smaller than orb_cutoff[T1] 
						// this should be consistant with LCAO_nnr::cal_nnrg function 
						// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
						// distance = 7.0000000000000000
						// orb_cutoff[T1] = 7.0000000000000008
						if(distance < rcut - 1.0e-15)
						{
							//storing the indexed for nnrg
							const int mu = pv->global2local_row(iat);
							const int nu = pv->global2local_col(iat2);
							this->nnrg_index.push_back(gridIntegral::gridIndex{this->nnrg, mu, nu, GlobalC::GridD.getBox(ad), atom1->nw, atom2->nw});
							
							const int nelement = atom1->nw * atom2->nw;
							this->nnrg += nelement;
							this->nlocdimg[iat] += nelement; 
							this->nad[iat]++;
							++count;
						}
					}// end iat2
				}// end ad
//				GlobalV::ofs_running << " iat=" << iat << " nlocstartg=" << nlocstartg[iat] << " nad=" << nad[iat] << std::endl;
			}// end iat
		}// end I1
	}// end T1

	if(GlobalV::OUT_LEVEL != "m") { ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nnrg",this->nnrg);
}

	//--------------------------------------------------
	// search again, to order each (iat2, b1, b2, b3)
	// find_R2 is used to target DM_R.
	// because DM_R is allocated with nnrg.
	// So once we had dR = R2 - R1 and iat2,
	// we need to find out the corresponding positions
	// in DM_R
	//--------------------------------------------------
	if(allocate_find_R2)
	{
		find_R2.clear();
		find_R2st.clear();
		find_R2_sorted_index.clear();
		allocate_find_R2 = false;
	}
	this->find_R2.resize(GlobalC::ucell.nat);
	this->find_R2st.resize(GlobalC::ucell.nat);
	this->find_R2_sorted_index.resize(GlobalC::ucell.nat);
	for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
	{
		this->find_R2_sorted_index[iat].resize(nad[iat]);
		std::fill(this->find_R2_sorted_index[iat].begin(), this->find_R2_sorted_index[iat].end(), 0);
		this->find_R2st[iat].resize(nad[iat]);
		std::fill(this->find_R2st[iat].begin(), this->find_R2st[iat].end(), 0);
		this->find_R2[iat].resize(nad[iat]);
		std::fill(this->find_R2[iat].begin(), this->find_R2[iat].end(), 0);
	}

	allocate_find_R2 = true;

	for (int T1 = 0; T1 < GlobalC::ucell.ntype; T1++)
	{
		for (int I1 = 0; I1 < GlobalC::ucell.atoms[T1].na; I1++)
		{
//			std::cout << " T1=" << T1 << " I1=" << I1 << std::endl; 
			tau1 = GlobalC::ucell.atoms[T1].tau[I1];
			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
			const int iat = GlobalC::ucell.itia2iat(T1,I1);

//			std::cout << " Number of adjacent = " << GlobalC::GridD.getAdjacentNum()+1 << std::endl;
			
			int count=0;
			for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ad++)
			{
		//		std::cout << " ad=" << ad << std::endl;
				const int T2 = GlobalC::GridD.getType(ad);
				const int I2 = GlobalC::GridD.getNatom(ad);
				const int iat2 = GlobalC::ucell.itia2iat(T2,I2);

				// if this atom is in this processor.
				if(this->in_this_processor[iat])
				{
					if(this->in_this_processor[iat2])
					{
						dtau = GlobalC::GridD.getAdjacentTau(ad) - tau1;
                        double distance = dtau.norm() * GlobalC::ucell.lat0;
                        double rcut = orb_cutoff[T1] + orb_cutoff[T2];

						const int b1 = GlobalC::GridD.getBox(ad).x;
						const int b2 = GlobalC::GridD.getBox(ad).y;
						const int b3 = GlobalC::GridD.getBox(ad).z;
					
						// mohan fix bug 2011-06-26, should be '<', not '<='	
						//			if(distance < rcut)

						// mohan reset this 2013-07-02 in Princeton
						// we should make absolutely sure that the distance is smaller than orb_cutoff[it] 
						// this should be consistant with LCAO_nnr::cal_nnrg function 
						// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
						// distance = 7.0000000000000000
						// orb_cutoff[it] = 7.0000000000000008
						if(distance < rcut - 1.0e-15)
						{
						//	assert( count < nad[iat] );
							//--------------------------------------------------------------
							// start positions of adjacent atom of 'iat'
							// note: the first is not zero.
							//--------------------------------------------------------------
							find_R2[iat][count] = this->cal_RindexAtom(b1, b2, b3, iat2);
							find_R2_sorted_index[iat][count] = count;


							// find_R2st
							// note: the first must be zero.
							// find_R2st: start position of each adjacen atom.
							if( count + 1 < nad[iat] )
							{
								find_R2st[iat][count+1] = find_R2st[iat][count] 
								+ GlobalC::ucell.atoms[T1].nw * GlobalC::ucell.atoms[T2].nw; //modified by zhengdy-soc
							}
							++count;
						}
					}
				}
			}
 // Include the necessary header file

			std::stable_sort(find_R2_sorted_index[iat].begin(), find_R2_sorted_index[iat].begin() + nad[iat],
				[&](int pos1, int pos2){ return find_R2[iat][pos1] < find_R2[iat][pos2]; });
		}
	}

	return;
}

void Grid_Technique::cal_max_box_index()
{
	ModuleBase::TITLE("LCAO_nnr","cal_max_box_index");
	this->maxB1 = this->maxB2 = this->maxB3 = -10000;
	this->minB1 = this->minB2 = this->minB3 = 10000;
	for (int T1 = 0; T1 < GlobalC::ucell.ntype; T1++)
	{
		for (int I1 = 0; I1 < GlobalC::ucell.atoms[T1].na; I1++)
		{
			ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[T1].tau[I1];
			//GlobalC::GridD.Find_atom(tau1);
			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
			for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ad++)
			{
				this->maxB1 = std::max( GlobalC::GridD.getBox(ad).x, maxB1 ); 
				this->maxB2 = std::max( GlobalC::GridD.getBox(ad).y, maxB2 ); 
				this->maxB3 = std::max( GlobalC::GridD.getBox(ad).z, maxB3 ); 

				this->minB1 = std::min( GlobalC::GridD.getBox(ad).x, minB1 ); 
				this->minB2 = std::min( GlobalC::GridD.getBox(ad).y, minB2 ); 
				this->minB3 = std::min( GlobalC::GridD.getBox(ad).z, minB3 ); 
			}
		}
	}

	nB1 = maxB1-minB1+1;
	nB2 = maxB2-minB2+1;
	nB3 = maxB3-minB3+1;

	nbox = nB1 * nB2 * nB3;
	
	//ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nbox",nbox);

	return;
}

int Grid_Technique::cal_RindexAtom(const int &u1, const int &u2, const int &u3, const int &iat2) const
{
	const int x1 = u1 - this->minB1;
	const int x2 = u2 - this->minB2;
	const int x3 = u3 - this->minB3;
	
	if(x1<0 || x2<0 || x3<0)
	{
		std::cout << " u1=" << u1 << " minB1=" << minB1 << std::endl;
		std::cout << " u2=" << u2 << " minB2=" << minB2 << std::endl;
		std::cout << " u3=" << u3 << " minB3=" << minB3 << std::endl;
		ModuleBase::WARNING_QUIT("LCAO_nnr::cal_Rindex","x1<0 || x2<0 || x3<0 !");
	}

	assert(x1>=0);
	assert(x2>=0);
	assert(x3>=0);

	return (iat2 + (x3 + x2 * this->nB3 + x1 * this->nB2 * this->nB3) * GlobalC::ucell.nat);
}

int Grid_Technique::binary_search_find_R2_offset(int val, int iat) const
{
    auto findR2 = this->find_R2[iat];
	auto findR2_index = this->find_R2_sorted_index[iat];

	int left = 0;
    int right = nad[iat] - 1;
    while(left <= right)
    {
        int mid = left + ((right - left) >> 1);
		int idx = findR2_index[mid];
        if(val == findR2[idx])
        {
            return idx;
        }
        if(val < findR2[idx])
        {
            right = mid - 1;
        }
        else
        {
            left = mid + 1;
        }
    }
    return -1;
}
