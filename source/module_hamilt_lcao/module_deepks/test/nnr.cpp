#include "LCAO_deepks_test.h"

void test_deepks::cal_nnr(void)
{
	ModuleBase::TITLE("test_deepks","cal_nnr");

	this->nnr = 0;
	int start = 0;
	//int ind1 = 0;
	int iat = 0;

	// (1) find the adjacent atoms of atom[T1,I1];
	ModuleBase::Vector3<double> tau1;
	ModuleBase::Vector3<double> tau2;
	ModuleBase::Vector3<double> dtau;
	ModuleBase::Vector3<double> tau0;
	ModuleBase::Vector3<double> dtau1;
	ModuleBase::Vector3<double> dtau2;

	for (int T1 = 0; T1 < ucell.ntype; T1++)
	{
		for (int I1 = 0; I1 < ucell.atoms[T1].na; I1++)
		{
			tau1 = ucell.atoms[T1].tau[I1];
			Test_Deepks::GridD.Find_atom(ucell,  tau1 ,T1, I1);
			const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
			int nw1 = ucell.atoms[T1].nw;

			// (2) search among all adjacent atoms.
			for (int ad = 0; ad < Test_Deepks::GridD.getAdjacentNum()+1; ad++)
			{
				const int T2 = Test_Deepks::GridD.getType(ad);
				const int I2 = Test_Deepks::GridD.getNatom(ad);
				const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
				int nw2 = ucell.atoms[T2].nw;

				tau2 = Test_Deepks::GridD.getAdjacentTau(ad);

				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
				
				if(distance < rcut)
				{
					//--------------------------------------------------
					// calculate how many matrix elements are in 
					// this processor.
					for(int ii=0; ii<nw1; ii++)
					{
						// the index of orbitals in this processor
						// according to HPSEPS's division method.
						const int iw1_all = start1 + ii;
						const int mu = ParaO.global2local_row(iw1_all);
						if(mu<0)continue;

						for(int jj=0; jj<nw2; jj++)
						{
							const int iw2_all = start2 + jj;
							const int nu = ParaO.global2local_col(iw2_all);
							if(nu<0)continue;
							++nnr;
						}// end jj
					}// end ii
				}//end distance
				// there is another possibility that i and j are adjacent atoms.
				// which is that <i|beta> are adjacents while <beta|j> are also
				// adjacents, these considerations are only considered in k-point
				// algorithm, 
				// mohan fix bug 2012-07-03
				else if(distance >= rcut)
				{
					for (int ad0 = 0; ad0 < Test_Deepks::GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = Test_Deepks::GridD.getType(ad0);
						const int I0 = Test_Deepks::GridD.getNatom(ad0);
					
						tau0 = Test_Deepks::GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1; 
						double distance1 = dtau1.norm() * ucell.lat0;
						double rcut1 = ORB.Phi[T1].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max();

						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * ucell.lat0;
						double rcut2 = ORB.Phi[T2].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							for(int ii=0; ii<nw1; ++ii)
							{
								const int iw1_all = start1 + ii;
								const int mu = ParaO.global2local_row(iw1_all);
								if(mu<0)continue;

								for(int jj=0; jj<nw2; ++jj)
								{
									const int iw2_all = start2 + jj;
									const int nu = ParaO.global2local_col(iw2_all);
									if(nu<0)continue;
									++nnr;
								}
							}
							break;	
						} // dis1, dis2
					}//ad0
				}
			}// end ad
		}// end I1
	} // end T1

    GlobalV::ofs_running << "nnr : " << nnr << std::endl; 

	return;
}

void test_deepks::folding_nnr(const Test_Deepks::K_Vectors &kv)
{

	int iat = 0;
	int index = 0;
	ModuleBase::Vector3<double> dtau;
	ModuleBase::Vector3<double> tau1;
	ModuleBase::Vector3<double> tau2;

	ModuleBase::Vector3<double> dtau1;
	ModuleBase::Vector3<double> dtau2;
	ModuleBase::Vector3<double> tau0;

	for(int ik=0;ik<kv.nkstot;ik++)
	{
		ModuleBase::GlobalFunc::ZEROS(ld.H_V_delta_k[ik].data(), ParaO.nloc);
		index=0;
		for (int T1 = 0; T1 < ucell.ntype; ++T1)
		{
			Atom* atom1 = &ucell.atoms[T1];
			for (int I1 = 0; I1 < atom1->na; ++I1)
			{
				tau1 = atom1->tau[I1];
				Test_Deepks::GridD.Find_atom(ucell, tau1, T1, I1);
				Atom* atom1 = &ucell.atoms[T1];
				const int start = ucell.itiaiw2iwt(T1,I1,0);

				// (2) search among all adjacent atoms.
				for (int ad = 0; ad < Test_Deepks::GridD.getAdjacentNum()+1; ++ad)
				{
					const int T2 = Test_Deepks::GridD.getType(ad);
					const int I2 = Test_Deepks::GridD.getNatom(ad);
					Atom* atom2 = &ucell.atoms[T2];

					tau2 = Test_Deepks::GridD.getAdjacentTau(ad);
					dtau = tau2 - tau1;
					double distance = dtau.norm() * ucell.lat0;
					double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

					bool adj = false;

					if(distance < rcut) 
					{
						adj = true;
					}
					else if(distance >= rcut)
					{
						for (int ad0 = 0; ad0 < Test_Deepks::GridD.getAdjacentNum()+1; ++ad0)
						{
							const int T0 = Test_Deepks::GridD.getType(ad0); 
							const int I0 = Test_Deepks::GridD.getNatom(ad0); 

							tau0 = Test_Deepks::GridD.getAdjacentTau(ad0);
							dtau1 = tau0 - tau1;
							dtau2 = tau0 - tau2;

							double distance1 = dtau1.norm() * ucell.lat0;
							double distance2 = dtau2.norm() * ucell.lat0;

							double rcut1 = ORB.Phi[T1].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max();
							double rcut2 = ORB.Phi[T2].getRcut() + ucell.infoNL.Beta[T0].get_rcut_max();

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
						const int start2 = ucell.itiaiw2iwt(T2,I2,0);
						//------------------------------------------------
						// exp(k dot dR)
						// dR is the index of box in Crystal coordinates
						//------------------------------------------------
						ModuleBase::Vector3<double> dR(Test_Deepks::GridD.getBox(ad).x, Test_Deepks::GridD.getBox(ad).y, Test_Deepks::GridD.getBox(ad).z); 
						const double arg = ( kv.kvec_d[ik] * dR ) * ModuleBase::TWO_PI;
						const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );

						//--------------------------------------------------
						// calculate how many matrix elements are in 
						// this processor.
						//--------------------------------------------------
						for(int ii=0; ii<atom1->nw*GlobalV::NPOL; ii++)
						{
							// the index of orbitals in this processor
							const int iw1_all = start + ii;
							const int mu = ParaO.global2local_row(iw1_all);
							if(mu<0)continue;

							for(int jj=0; jj<atom2->nw*GlobalV::NPOL; jj++)
							{
								int iw2_all = start2 + jj;
								const int nu = ParaO.global2local_col(iw2_all);

								if(nu<0)continue;
								int iic;
                                if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER() )
								{
									iic=mu+nu*ParaO.nrow;
								}
								else
								{
									iic=mu*ParaO.ncol+nu;
								}
								ld.H_V_delta_k[ik][iic] += ld.H_V_deltaR[index] * kphase;
								
								++index;
							}//end jj
						}//end ii
					}
				}// end ad
				++iat;
			}// end I1
		} // end T1
		assert(index==this->nnr);
	}
}
