#include "record_adj.h"

#include "module_base/timer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
Record_adj::Record_adj() : iat2ca(nullptr) {}
Record_adj::~Record_adj(){
	if(info_modified)
		this->delete_grid();
}

void Record_adj::delete_grid(void)
{
	for(int i=0; i<na_proc; i++)
	{
		// how many 'numerical orbital' adjacents
		// for each atom in this process.
		for(int j=0; j< na_each[i]; j++)
		{
			delete[] info[i][j];
		}
		delete[] info[i];
	}
	delete[] info;
	delete[] na_each;
	if (iat2ca) delete[] iat2ca;
	info_modified=false;
}



//--------------------------------------------
// This will record the orbitals according to
// HPSEPS's 2D block division.
// If multi-k, calculate nnr at the same time.
// be called only once in an ion-step.
//--------------------------------------------
void Record_adj::for_2d(Parallel_Orbitals &pv, bool gamma_only)
{
	ModuleBase::TITLE("Record_adj","for_2d");
	ModuleBase::timer::tick("Record_adj","for_2d");

    assert(GlobalC::ucell.nat > 0);
    if (!gamma_only)
    {
        delete[] pv.nlocdim;
        delete[] pv.nlocstart;
        pv.nlocdim = new int[GlobalC::ucell.nat];	
        pv.nlocstart = new int[GlobalC::ucell.nat];
        ModuleBase::GlobalFunc::ZEROS(pv.nlocdim, GlobalC::ucell.nat);
        ModuleBase::GlobalFunc::ZEROS(pv.nlocstart, GlobalC::ucell.nat);
        pv.nnr = 0;
    }
{
	// (1) find the adjacent atoms of atom[T1,I1];
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;

	this->na_proc = GlobalC::ucell.nat;

	// number of adjacents for each atom.	
	this->na_each = new int[na_proc];
	ModuleBase::GlobalFunc::ZEROS(na_each, na_proc);
	int iat = 0;
	
//	std::cout << " in for_2d" << std::endl;

	for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
	{
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			//GlobalC::GridD.Find_atom( tau1 );
			GlobalC::GridD.Find_atom(GlobalC::ucell,  tau1 ,T1, I1);
			const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            if(!gamma_only) pv.nlocstart[iat] = pv.nnr;
            
            // (2) search among all adjacent atoms.
			for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GlobalC::GridD.getType(ad);
				const int I2 = GlobalC::GridD.getNatom(ad);
				const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
				tau2 = GlobalC::GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * GlobalC::ucell.lat0;
				double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

				bool is_adj = false;
                if (distance < rcut) is_adj = true;
                // there is another possibility that i and j are adjacent atoms.
				// which is that <i|beta> are adjacents while <beta|j> are also
				// adjacents, these considerations are only considered in k-point
				// algorithm, 
                else if (distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GlobalC::GridD.getType(ad0);
						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GlobalC::GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
						double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * GlobalC::ucell.lat0;
						double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							is_adj = true;
							break;
						} // dis1, dis2
					}
				}

				if(is_adj)
				{
                    ++na_each[iat];
                    if (!gamma_only)
                    {
                        for(int ii=0; ii<atom1->nw * GlobalV::NPOL; ++ii)
                        {
                            // the index of orbitals in this processor
                            const int iw1_all = start1 + ii;
                            const int mu = pv.trace_loc_row[iw1_all];
                            if(mu<0)continue;

                            for(int jj=0; jj<GlobalC::ucell.atoms[T2].nw * GlobalV::NPOL; ++jj)
                            {
                                const int iw2_all = start2 + jj;
                                const int nu = pv.trace_loc_col[iw2_all];
                                if(nu<0)continue;
                                
                                pv.nlocdim[iat]++;
                                ++(pv.nnr);
                            }
                        }
                    }
                }//end is_adj
            }//end ad
            ++iat;
		}//end I1
	}//end T1
}
    //xiaohui add "OUT_LEVEL", 2015-09-16
    if (GlobalV::OUT_LEVEL != "m" && !gamma_only) ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "ParaV.nnr", pv.nnr);

	//------------------------------------------------
	// info will identify each atom in each unitcell.
	//------------------------------------------------
	this->info = new int**[na_proc];
#ifdef _OPENMP
#pragma omp parallel
{
#endif

	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for(int i=0; i<na_proc; i++)
	{
//		GlobalV::ofs_running << " atom" << std::setw(5) << i << std::setw(10) << na_each[i] << std::endl;
		if( na_each[i] > 0)
		{
			info[i] = new int*[ na_each[i] ];
			for(int j=0; j<na_each[i]; j++)
			{
				// (Rx, Ry, Rz, T, I)
				info[i][j] = new int[5];
				ModuleBase::GlobalFunc::ZEROS(info[i][j],5);
			}
		}
	}

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
	{
		const int T1 = GlobalC::ucell.iat2it[iat];
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		const int I1 = GlobalC::ucell.iat2ia[iat];
		{
			tau1 = atom1->tau[I1];
			//GlobalC::GridD.Find_atom( tau1 );
			AdjacentAtomInfo adjs;
			GlobalC::GridD.Find_atom(GlobalC::ucell,  tau1 ,T1, I1, &adjs);

			// (2) search among all adjacent atoms.
			int cb = 0;
			for (int ad = 0; ad < adjs.adj_num+1; ++ad)
			{
				const int T2 = adjs.ntype[ad];
				const int I2 = adjs.natom[ad];
				tau2 = adjs.adjacent_tau[ad];
				dtau = tau2 - tau1;
				double distance = dtau.norm() * GlobalC::ucell.lat0;
				double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();


				bool is_adj = false;
				if(distance < rcut) is_adj = true;
				else if(distance >= rcut) 
				{
					for (int ad0 = 0; ad0 < adjs.adj_num+1; ++ad0)
					{
						const int T0 = adjs.ntype[ad0];
						//const int I0 = GlobalC::GridD.getNatom(ad0);
						//const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
						//const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = adjs.adjacent_tau[ad0];
						dtau1 = tau0 - tau1;
						double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
						double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * GlobalC::ucell.lat0;
						double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							is_adj = true;
							break;
						} // dis1, dis2
					}
				}

				if(is_adj)
				{
					info[iat][cb][0] = adjs.box[ad].x; 
					info[iat][cb][1] = adjs.box[ad].y; 
					info[iat][cb][2] = adjs.box[ad].z; 
					info[iat][cb][3] = T2;
                    info[iat][cb][4] = I2;
                    ++cb;
				}
			}//end ad
//			GlobalV::ofs_running << " nadj = " << cb << std::endl;
		}//end I1
	}//end T1
#ifdef _OPENMP
}
#endif
	ModuleBase::timer::tick("Record_adj","for_2d");
	info_modified=true;
	return;
}


//--------------------------------------------
// This will record the orbitals according to
// grid division (cut along z direction) 
//--------------------------------------------
void Record_adj::for_grid(const Grid_Technique &gt)
{
    ModuleBase::TITLE("Record_adj","for_grid");
	ModuleBase::timer::tick("Record_adj","for_grid");
	
	this->na_proc = 0;
	this->iat2ca = new int[GlobalC::ucell.nat];
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
	{
		{
			if(gt.in_this_processor[iat])
			{
				iat2ca[iat] = na_proc;
				++na_proc;
			} else {
				iat2ca[iat] = -1;
			}
		}
	}

	// number of adjacents for each atom.	
	this->na_each = new int[na_proc];
	ModuleBase::GlobalFunc::ZEROS(na_each, na_proc);
	this->info = new int**[na_proc];
#ifdef _OPENMP
#pragma omp parallel
{
#endif
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> tau0, dtau1, dtau2;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
	{
		const int T1 = GlobalC::ucell.iat2it[iat];
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		const int I1 = GlobalC::ucell.iat2ia[iat];
		{
			const int ca = iat2ca[iat];
			// key in this function
			if(gt.in_this_processor[iat])
			{	
				tau1 = atom1->tau[I1];
				//GlobalC::GridD.Find_atom(tau1);
				AdjacentAtomInfo adjs;
				GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1, &adjs);
				for (int ad = 0; ad < adjs.adj_num+1; ad++)
				{
					const int T2 = adjs.ntype[ad];
					const int I2 = adjs.natom[ad];
					const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
					if(gt.in_this_processor[iat2])
					{
						//Atom* atom2 = &GlobalC::ucell.atoms[T2];
						tau2 = adjs.adjacent_tau[ad];
						dtau = tau2 - tau1;
						double distance = dtau.norm() * GlobalC::ucell.lat0;
						double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();


						bool is_adj = false;
						if(distance < rcut) is_adj = true;
						/*
						else if(distance >= rcut)
						{
                            for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
                            {
                                const int T0 = GlobalC::GridD.getType(ad0);
                                const int I0 = GlobalC::GridD.getNatom(ad0);
                                const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                                const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

                                tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                                dtau1 = tau0 - tau1;
                                dtau2 = tau0 - tau2;

                                double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                                double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                                double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                                double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                                if( distance1 < rcut1 && distance2 < rcut2 )
                                {
                                   	is_adj = true; 
                                    break;
								} // dis1, dis2
							}
						}
						*/

						// check the distance
						if(is_adj)
						{
							++na_each[ca];
						}
					}//end judge 2
				}//end ad
			}//end judge 1
		}//end I1
	}//end T1

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for(int i=0; i<na_proc; i++)
	{
		assert(na_each[i]>0);
		info[i] = new int*[ na_each[i] ];
		for(int j=0; j<na_each[i]; j++)
		{
			// (Rx, Ry, Rz, T, I)
			info[i][j] = new int[5];
			ModuleBase::GlobalFunc::ZEROS(info[i][j],5);
		}
	}

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
	{
		const int T1 = GlobalC::ucell.iat2it[iat];
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		const int I1 = GlobalC::ucell.iat2ia[iat];
		{
			const int ca = iat2ca[iat];

			// key of this function
			if(gt.in_this_processor[iat])
			{
				tau1 = atom1->tau[I1];
				//GlobalC::GridD.Find_atom(tau1);
				AdjacentAtomInfo adjs;
				GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1, &adjs);

				int cb = 0;
				for (int ad = 0; ad < adjs.adj_num+1; ad++)
				{
					const int T2 = adjs.ntype[ad];
					const int I2 = adjs.natom[ad];
					const int iat2 = GlobalC::ucell.itia2iat(T2, I2);

					// key of this function
					if(gt.in_this_processor[iat2])
					{
						//Atom* atom2 = &GlobalC::ucell.atoms[T2];
						tau2 = adjs.adjacent_tau[ad];
						dtau = tau2 - tau1;
						double distance = dtau.norm() * GlobalC::ucell.lat0;
						double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

						// check the distance
						if(distance < rcut)
						{
							info[ca][cb][0] = adjs.box[ad].x; 
							info[ca][cb][1] = adjs.box[ad].y; 
							info[ca][cb][2] = adjs.box[ad].z; 
							info[ca][cb][3] = T2;
							info[ca][cb][4] = I2;
							++cb;
						}
						/*
                        else if(distance >= rcut)
                        {
                            for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum()+1; ++ad0)
                            {
                                const int T0 = GlobalC::GridD.getType(ad0);
                                const int I0 = GlobalC::GridD.getNatom(ad0);
                                const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                                const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

                                tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                                dtau1 = tau0 - tau1;
                                dtau2 = tau0 - tau2;

                                double distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                                double distance2 = dtau2.norm() * GlobalC::ucell.lat0;

                                double rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();
                                double rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                                if( distance1 < rcut1 && distance2 < rcut2 )
                                {
									info[ca][cb][0] = GlobalC::GridD.getBox(ad).x; 
									info[ca][cb][1] = GlobalC::GridD.getBox(ad).y; 
									info[ca][cb][2] = GlobalC::GridD.getBox(ad).z; 
									info[ca][cb][3] = T2;
									info[ca][cb][4] = I2;
									++cb;
									break;
								} // dis1, dis2
							}
                        }
						*/
					}
				}// end ad

				assert(cb == na_each[ca]);
			}
		}
	}
#ifdef _OPENMP
}
#endif
	ModuleBase::timer::tick("Record_adj","for_grid");
	info_modified=true;
//	std::cout << " after for_grid" << std::endl;
	return;
}

