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

// This is for cell R dependent part. 
void Grid_Technique::cal_nnrg(Parallel_Orbitals* pv)
{
	ModuleBase::TITLE("LCAO_nnr","cal_nnrg");

	this->cal_max_box_index();

	this->nnrg = 0;

	delete[] nlocdimg;
	delete[] nlocstartg;
	delete[] nad;
	this->nnrg_index.resize(0);
	
	this->nad = new int[GlobalC::ucell.nat];
	this->nlocdimg = new int[GlobalC::ucell.nat];	
	this->nlocstartg = new int[GlobalC::ucell.nat];
	
	ModuleBase::GlobalFunc::ZEROS(nad, GlobalC::ucell.nat);
	ModuleBase::GlobalFunc::ZEROS(nlocdimg, GlobalC::ucell.nat);
	ModuleBase::GlobalFunc::ZEROS(nlocstartg, GlobalC::ucell.nat);


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
						double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();


						//if(distance < rcut)
						// mohan reset this 2013-07-02 in Princeton
						// we should make absolutely sure that the distance is smaller than GlobalC::ORB.Phi[it].getRcut
						// this should be consistant with LCAO_nnr::cal_nnrg function 
						// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
						// distance = 7.0000000000000000
						// GlobalC::ORB.Phi[it].getRcut = 7.0000000000000008
						if(distance < rcut - 1.0e-15)
						{
							//storing the indexed for nnrg
							const int mu = pv->trace_loc_row[iat];
							const int nu = pv->trace_loc_col[iat2];
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

	if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nnrg",this->nnrg);

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
		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			delete[] find_R2_sorted_index[iat];
		}
		delete[] find_R2_sorted_index;

		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			delete[] find_R2[iat];
		}
		delete[] find_R2;

		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			delete[] find_R2st[iat];
		}
		delete[] find_R2st;
		allocate_find_R2 = false;
	}

	this->find_R2_sorted_index = new int*[GlobalC::ucell.nat];
	for(int iat=0; iat<GlobalC::ucell.nat; iat++)
	{
		// at least nad contains itself, so nad[iat] can not be 0.
		this->find_R2_sorted_index[iat] = new int[nad[iat]];
		ModuleBase::GlobalFunc::ZEROS(find_R2_sorted_index[iat], nad[iat]);
	}

	this->find_R2 = new int*[GlobalC::ucell.nat];
	for(int iat=0; iat<GlobalC::ucell.nat; iat++)
	{
		// at least nad contains itself, so nad[iat] can not be 0.
		this->find_R2[iat] = new int[nad[iat]];
		ModuleBase::GlobalFunc::ZEROS(find_R2[iat], nad[iat]);
	}

	this->find_R2st = new int*[GlobalC::ucell.nat];
	for(int iat=0; iat<GlobalC::ucell.nat; iat++)
	{
		this->find_R2st[iat] = new int[nad[iat]];
		ModuleBase::GlobalFunc::ZEROS(find_R2st[iat], nad[iat]);
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
                        double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

						const int b1 = GlobalC::GridD.getBox(ad).x;
						const int b2 = GlobalC::GridD.getBox(ad).y;
						const int b3 = GlobalC::GridD.getBox(ad).z;
					
						// mohan fix bug 2011-06-26, should be '<', not '<='	
						//			if(distance < rcut)

						// mohan reset this 2013-07-02 in Princeton
						// we should make absolutely sure that the distance is smaller than GlobalC::ORB.Phi[it].getRcut
						// this should be consistant with LCAO_nnr::cal_nnrg function 
						// typical example : 7 Bohr cutoff Si orbital in 14 Bohr length of cell.
						// distance = 7.0000000000000000
						// GlobalC::ORB.Phi[it].getRcut = 7.0000000000000008
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
			std::stable_sort(find_R2_sorted_index[iat], find_R2_sorted_index[iat] + nad[iat],
				[&](int pos1, int pos2){ return find_R2[iat][pos1] < find_R2[iat][pos2]; });
		}
	}

	return;
}

void Grid_Technique::cal_max_box_index(void)
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
				this->maxB1 = max( GlobalC::GridD.getBox(ad).x, maxB1 ); 
				this->maxB2 = max( GlobalC::GridD.getBox(ad).y, maxB2 ); 
				this->maxB3 = max( GlobalC::GridD.getBox(ad).z, maxB3 ); 

				this->minB1 = min( GlobalC::GridD.getBox(ad).x, minB1 ); 
				this->minB2 = min( GlobalC::GridD.getBox(ad).y, minB2 ); 
				this->minB3 = min( GlobalC::GridD.getBox(ad).z, minB3 ); 
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

// be called in LCAO_Hamilt::calculate_Hk.
void LCAO_Matrix::folding_fixedH(
						const int &ik, 
						const std::vector<ModuleBase::Vector3<double>>& kvec_d,
						bool cal_syns)
{
	ModuleBase::TITLE("LCAO_nnr","folding_fixedH");
    ModuleBase::timer::tick("LCAO_nnr", "folding_fixedH");
    const Parallel_Orbitals* pv = this->ParaV;

	int tot_index = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:tot_index)
{
	const int num_threads = omp_get_num_threads();
	const int thread_id = omp_get_thread_num();
#else
	const int num_threads = 1;
	const int thread_id = 0;
#endif
	ModuleBase::Vector3<double> dtau;
	ModuleBase::Vector3<double> tau1;
	ModuleBase::Vector3<double> tau2;

	ModuleBase::Vector3<double> dtau1;
	ModuleBase::Vector3<double> dtau2;
	ModuleBase::Vector3<double> tau0;

#ifdef __DEEPKS
	if (GlobalV::deepks_scf)
    {
		int beg, len;
		ModuleBase::BLOCK_TASK_DIST_1D(num_threads, thread_id, (int)pv->nloc, 1024, beg, len);
		ModuleBase::GlobalFunc::ZEROS(GlobalC::ld.H_V_delta_k[ik] + beg, len);
	}
#endif

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for (int iat=0; iat<GlobalC::ucell.nat; ++iat)
	{
		const int T1 = GlobalC::ucell.iat2it[iat];
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		const int I1 = GlobalC::ucell.iat2ia[iat];
		{
			tau1 = atom1->tau[I1];
			//GlobalC::GridD.Find_atom(tau1);
			AdjacentAtomInfo adjs;
			GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1, &adjs);
			Atom* atom1 = &GlobalC::ucell.atoms[T1];
			const int start = GlobalC::ucell.itiaiw2iwt(T1,I1,0);
			int index = pv->nlocstart[iat];

			if (cal_syns)
            {
                for (int k = 0; k < 3; k++)
                {
                    tau1[k] = tau1[k] - atom1->vel[I1][k] * INPUT.mdp.md_dt / GlobalC::ucell.lat0 ;
                }
            }

			// (2) search among all adjacent atoms.
			for (int ad = 0; ad < adjs.adj_num+1; ++ad)
			{
				const int T2 = adjs.ntype[ad];
				const int I2 = adjs.natom[ad];
				Atom* atom2 = &GlobalC::ucell.atoms[T2];

				tau2 = adjs.adjacent_tau[ad];
				dtau = tau2 - tau1;
				double distance = dtau.norm() * GlobalC::ucell.lat0;
				double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

				bool adj = false;

				if(distance < rcut) 
				{
					adj = true;
				}
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < adjs.adj_num+1; ++ad0)
					{
						const int T0 = adjs.ntype[ad0]; 
						const int I0 = adjs.natom[ad0]; 
						//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = adjs.adjacent_tau[ad0];
						dtau1 = tau0 - tau1;
						dtau2 = tau0 - tau2;

						double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
						double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

						double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
						double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							adj = true;
							break;
						}
					}
				}

				if(adj) // mohan fix bug 2011-06-26, should not be '<='
				{
					// (3) calculate the nu of atom (T2, I2)
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);
					//------------------------------------------------
					// exp(k dot dR)
					// dR is the index of box in Crystal coordinates
					//------------------------------------------------
					ModuleBase::Vector3<double> dR(adjs.box[ad].x, adjs.box[ad].y, adjs.box[ad].z); 
					const double arg = ( kvec_d[ik] * dR ) * ModuleBase::TWO_PI;
					//const double arg = ( kvec_d[ik] * GlobalC::GridD.getBox(ad) ) * ModuleBase::TWO_PI;
					double sinp, cosp;
					ModuleBase::libm::sincos(arg, &sinp, &cosp);
					const std::complex<double> kphase = std::complex <double> ( cosp,  sinp );

					//--------------------------------------------------
					// calculate how many matrix elements are in 
					// this processor.
					//--------------------------------------------------
					//########################### EXPLAIN ###############################
					// 1. overlap matrix with k point
					// this->SlocR = < phi_0i | phi_Rj >, where 0, R are the cell index
					// while i,j are the orbital index.

					// 2. H_fixed=T+Vnl matrix element with k point (if Vna is not used).
					// H_fixed=T+Vnl+Vna matrix element with k point (if Vna is used).
					// this->Hloc_fixed = < phi_0i | H_fixed | phi_Rj>

					// 3. H(k) |psi(k)> = S(k) | psi(k)> 
					// Sloc2 is used to diagonalize for a give k point.
					// Hloc_fixed2 is used to diagonalize (eliminate index R).
					//###################################################################
#define NW_IJ_LOOP_BEGIN() \
					for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)				\
					{															\
						/* the index of orbitals in this processor */			\
						const int iw1_all = start + ii;							\
						const int mu = pv->trace_loc_row[iw1_all];				\
						if(mu<0) {continue;}									\
						for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)			\
						{														\
							int iw2_all = start2 + jj;							\
							const int nu = pv->trace_loc_col[iw2_all];			\
							if(nu<0) {continue;}

								// DO STH.

#define NW_IJ_LOOP_END() \
							++index;											\
							++tot_index;										\
						}/*end jj*/												\
					}/*end ii*/

					if(GlobalV::NSPIN!=4)
					{
						auto compute_kernel = [&](int iic, int index) {
							this->Sloc2[iic] += this->SlocR[index] * kphase;
							this->Hloc_fixed2[iic] += this->Hloc_fixedR[index] * kphase;
#ifdef __DEEPKS
							if(GlobalV::deepks_scf)
							{
								GlobalC::ld.H_V_delta_k[ik][iic] += GlobalC::ld.H_V_deltaR[index] * kphase;
							}
#endif
						};
						if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
						{
							NW_IJ_LOOP_BEGIN()
							compute_kernel(mu+nu*pv->nrow, index);
							NW_IJ_LOOP_END()
						} else {
							NW_IJ_LOOP_BEGIN()
							compute_kernel(mu*pv->ncol+nu, index);
							NW_IJ_LOOP_END()
						}
					}
					else
					{
						auto compute_kernel = [&](int iic, int index, int iw1_all, int iw2_all) {
							this->Sloc2[iic] += this->SlocR_soc[index] * kphase;
							this->Hloc_fixed2[iic] += this->Hloc_fixedR_soc[index] * kphase;
#ifdef __DEEPKS
							if(GlobalV::deepks_scf)
							{
								if (iw1_all % 2 == iw2_all % 2)
								{
									GlobalC::ld.H_V_delta_k[ik][iic] += GlobalC::ld.H_V_deltaR[index] * kphase;
								}
							}
#endif
						};
						if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
						{
							NW_IJ_LOOP_BEGIN()
							compute_kernel(mu+nu*pv->nrow, index, iw1_all, iw2_all);
							NW_IJ_LOOP_END()
						} else {
							NW_IJ_LOOP_BEGIN()
							compute_kernel(mu*pv->ncol+nu, index, iw1_all, iw2_all);
							NW_IJ_LOOP_END()
						}
					}
				}
			}// end ad
		}// end I1
	} // end T1
#ifdef _OPENMP
}
#endif
	assert(tot_index==this->ParaV->nnr);

	ModuleBase::timer::tick("LCAO_nnr","folding_fixedH");
	return;
}
