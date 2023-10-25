//==========================================================
// AUTHOR : Zhengpan , mohan , spshu
// DATE : 2007-9
//==========================================================
#include "symmetry.h"
#include "module_base/mymath.h"
bool ModuleSymmetry::test_brav = 0;

namespace ModuleSymmetry
{
Symmetry_Basic::Symmetry_Basic()
{
	this->epsilon = 1e-6;
}

Symmetry_Basic::~Symmetry_Basic()
{
}


// Find the type of bravais lattice.
std::string Symmetry_Basic::get_brav_name(const int ibrav) const
{
	switch(ibrav)
	{
		case 1: return "01. Cubic P (simple)";
		case 2: return "02. Cubic I (body-centered)";
		case 3: return "03. Cubic F (face-centered)";
		case 4: return "04. Hexagonal cell";
		case 5: return "05. Tetrogonal P (simple)";
		case 6: return "06. Tetrogonal I (body-centered)";
		case 7: return "07. Rhombohedral (Trigonal) cell";
		case 8: return "08. Orthorhombic P(simple)";
		case 9: return "09. Orthorhombic I (body-centered)";
		case 10: return "10. Orthorhombic F (face-centered)";
		case 11: return "11. Orthorhombic C (base-centered)";
		case 12: return "12. Monoclinic P (simple)";
		case 13: return "13. Monoclinic A (base-center)";
		case 14: return "14. Triclinic cell";
		case 15: return "wrong !! ";
	}
	return "Congratulations!You have found a bravais lattice that never existed!";
}

// Control the accuracy
bool Symmetry_Basic::equal(const double &m,const double &n)const
{
	//if( fabs(m-n) < 1.0e-5 )
    if( fabs(m-n) < epsilon ) //LiuXh add 2021-08-12, use accuracy for symmetry
	{
		return true;
	}
	return false;
}

// check the boundary condition of atom positions.
void Symmetry_Basic::check_boundary(double &x)const
{
	if(equal(x,-0.5) || equal(x,0.5)) x=-0.5;
	return;
}

double Symmetry_Basic::get_translation_vector(const double &x1, const double &x2)
{
	double t=0.0; // "t"ranslation
	t = x2 - x1;
	t = fmod(t+100.0, 1.0);
	if( fabs(t-1) < epsilon * 0.5) t = 0.0;
	return t;
}

void Symmetry_Basic::check_translation(double &x, const double &t) const
{
	x += t;
	//impose the periodic boundary condition
	x = fmod(x + 100.5,1) - 0.5;
	return;
}

double Symmetry_Basic::check_diff(const double &x1, const double &x2)
{
	double diff = x1 - x2;
	diff = fmod(diff + 100,1);
	//for reasons of safety
	if(fabs(diff - 1.0) < epsilon * 0.5)
	{
		diff = 0;
	}
	return diff;
}


void Symmetry_Basic::order_atoms(double* pos, const int& nat, const int* index) const
{
	double** tmp = new double*[nat];
	for(int ia=0; ia<nat; ia++)
	{
		tmp[ia] = new double[3];
	}

	for(int ia=0; ia<nat; ia++)
	{
		tmp[ia][0] = pos[index[ia]*3+0];
		tmp[ia][1] = pos[index[ia]*3+1];
		tmp[ia][2] = pos[index[ia]*3+2];
	}

	for(int ia=0; ia<nat; ia++)
	{
		pos[ia*3+0] = tmp[ia][0];
		pos[ia*3+1] = tmp[ia][1];
		pos[ia*3+2] = tmp[ia][2];
	}

	for(int ia=0; ia<nat; ia++)
	{
		delete[] tmp[ia];
	}
	delete[] tmp;

	return;
}
/*
// atom ordering for each atom type.
void Symmetry_Basic::atom_ordering(double *posi, const int natom, int *subindex)
{
	//order the atomic positions inside a supercell by a unique ordering scheme	
	subindex[0] = 0;

	if(natom == 1)
	{
		//if there is only one atom, it is not necessary to order
		return;
	}

	//order the x-dimension of the atomic position, 
	//only get a permutation table, data in momery are not changed
	//this->heapsort_pos(natom, posi, subindex);

	double *tmppos = new double[natom];
	for(int i=0; i<natom; i++)
	{
		tmppos[i] = posi[i*3];
	}
	ModuleBase::heapsort(natom, tmppos, subindex);

//	for(int i=0; i<natom; i++)
//	{
//		std::cout << "\n subindex[" << i << "]=" << subindex[i];
//	}
	
	this->order_atoms(posi, natom, subindex);

//	std::cout << "\n after Order X:" << std::endl;
//	print_pos(posi, natom);
	
	// order y
	// the old x value 
	// and the position
	double oldv = posi[0];
	int oldpos = 0;

	// the new x value
	// and the position
	double vnew;
	int newpos;
	for(int ia=0; ia<natom; ia++)
	{
		vnew = posi[ia*3];
		newpos = ia;
		if( vnew != oldv )
		{
			this->order_y(&posi[oldpos*3], oldpos, newpos);
			oldv = vnew;
			oldpos = ia;
		}
		else if(ia==natom-1)
		{
			// in this case, the new one is the last one,
			// also the same as the old one,
			// so we need to add 1
			this->order_y(&posi[oldpos*3], oldpos, newpos+1);	
		}
	}
	delete[] tmppos;
	return;
}

void Symmetry_Basic::order_y(double *pos1, const int &oldpos1, const int &newpos1)
{
//	ModuleBase::TITLE("Symmetry_Basic","order_y");
	// how many atoms need to be reordered according to same x value.
	const int nat1 = newpos1 - oldpos1;
//	std::cout << "\n nat1=" << nat1 << std::endl; 
	if(nat1 == 1) return;

	double* tmp1 = new double[nat1];
	int* index1 = new int[nat1];
	ModuleBase::GlobalFunc::ZEROS(index1, nat1);
	
	for(int ia=0; ia<nat1; ia++)
	{
		//+1 means y
		tmp1[ia] = pos1[3*ia+1];
//		std::cout << "\n ia=" << ia << " y=" << tmp1[ia];
	}

	ModuleBase::heapsort(nat1, tmp1, index1);

	this->order_atoms(pos1,nat1,index1);
	
//	for(int ia=0; ia<nat1; ia++)
//	{
//		std::cout << "\n index[" << ia << "]=" << index1[ia];
//	}

	double oldv2 = pos1[1];
	int oldpos2 = 0;

	double newv2;
	int newpos2;

	for(int ia=0; ia<nat1; ia++)
	{
		newv2 = pos1[ia*3+1];
		newpos2 = ia;
		
		if(newv2 != oldv2)
		{
			this->order_z(&pos1[oldpos2*3], oldpos2, newpos2);
			oldv2 = newv2;
			oldpos2 = newpos2;
		}
		else if(ia==nat1-1)
		{
			this->order_z(&pos1[oldpos2*3], oldpos2, newpos2+1);
		}
	}

	delete[] index1;
	delete[] tmp1;

	return;
}

void Symmetry_Basic::order_z(double* pos2, const int &oldpos2, const int &newpos2)
{
//	ModuleBase::TITLE("Symmetry_Basic","order_z");
	const int nat2 = newpos2 - oldpos2;
//	std::cout << "\n nat2=" << nat2; 
//	if(nat2==1) return;
	double* tmp2 = new double[nat2];
	int* index2 = new int[nat2];
	ModuleBase::GlobalFunc::ZEROS(index2, nat2);
	for(int ia=0; ia<nat2; ia++)
	{
		//+2 means z
		tmp2[ia] = pos2[3*ia+2];
	}

	ModuleBase::heapsort(nat2, tmp2, index2);

//	for(int ia=0; ia<nat2; ia++)
//	{
//		std::cout << "\n Z1 pos2[" << ia << "]=" << pos2[3*ia] << " " << pos2[3*ia+1] << " " << pos2[3*ia+2];
//	}
	
	this->order_atoms(pos2,nat2,index2);

//	for(int ia=0; ia<nat2; ia++)
//	{
//		std::cout << "\n Z2 pos2[" << ia << "]=" << pos2[3*ia] << " " << pos2[3*ia+1] << " " << pos2[3*ia+2];
//	}

	delete[] index2;
	delete[] tmp2;
}
*/
//convert a set of vectors (va) represented in the basis vectors old1, old2, old3 
//to a set of vectors (vb) represented in the basis vectors new1, new2, new3
void Symmetry_Basic::veccon(
		double *carpos, 
		double *rotpos, 
		const int num, 
		const ModuleBase::Vector3<double> &old1, 
		const ModuleBase::Vector3<double> &old2, 
		const ModuleBase::Vector3<double> &old3, 
		const ModuleBase::Vector3<double> &new1, 
		const ModuleBase::Vector3<double> &new2, 
		const ModuleBase::Vector3<double> &new3
		)
{

	GlobalV::ofs_running << "\n old1:" << old1.x << " " << old1.y << " " << old1.z;
	GlobalV::ofs_running << "\n old2:" << old2.x << " " << old2.y << " " << old2.z;
	GlobalV::ofs_running << "\n old3:" << old3.x << " " << old3.y << " " << old3.z;

	GlobalV::ofs_running << "\n new1:" << new1.x << " " << new1.y << " " << new1.z;
	GlobalV::ofs_running << "\n new2:" << new2.x << " " << new2.y << " " << new2.z;
	GlobalV::ofs_running << "\n new3:" << new3.x << " " << new3.y << " " << new3.z;

	ModuleBase::Matrix3 oldlat;
	oldlat.e11 = old1.x;
	oldlat.e12 = old1.y;
	oldlat.e13 = old1.z;
	oldlat.e21 = old2.x;
	oldlat.e22 = old2.y;
	oldlat.e23 = old2.z;
	oldlat.e31 = old3.x;
	oldlat.e32 = old3.y;
	oldlat.e33 = old3.z;

	ModuleBase::Matrix3 newlat;
	newlat.e11 = new1.x;
	newlat.e12 = new1.y;
	newlat.e13 = new1.z;
	newlat.e21 = new2.x;
	newlat.e22 = new2.y;
	newlat.e23 = new2.z;
	newlat.e31 = new3.x;
	newlat.e32 = new3.y;
	newlat.e33 = new3.z;

	ModuleBase::Matrix3 GT = newlat.Inverse();
	
	ModuleBase::Vector3<double> car;
	ModuleBase::Vector3<double> direct_old;
	ModuleBase::Vector3<double> direct_new;

	//calculate the reciprocal vectors rb1, rb2, rb3 for the vectors new1, new2, new3
	//this->recip(1.0, new1, new2, new3, rb1, rb2, rb3);
 	
	for(int i = 0; i < num; ++i)
	{
		direct_old.x = carpos[i * 3 + 0];
		direct_old.y = carpos[i * 3 + 1];
		direct_old.z = carpos[i * 3 + 2];

		car = direct_old * oldlat;

		
		
		direct_new = car * GT;

		//std::cout << "\n car[" << i << "]=" << car.x << " " << car.y << " " << car.z;
		//std::cout << "\n dir[" << i << "]=" << direct_new.x << " " << direct_new.y << " " << direct_new.z;

		rotpos[i * 3 + 0] = direct_new.x;
		rotpos[i * 3 + 1] = direct_new.y;
		rotpos[i * 3 + 2] = direct_new.z;
	}
	return;
}


//generate all point group symmetry operations from the generation group
void Symmetry_Basic::matrigen(ModuleBase::Matrix3 *symgen, const int ngen, ModuleBase::Matrix3* symop, int &nop) const
{
	int m1, m2;
	int n;
	ModuleBase::Matrix3 iden(1,0,0,0,1,0,0,0,1);
	ModuleBase::Matrix3 sig(1,0,0,0,1,0,0,0,1);
	ModuleBase::Matrix3 temp1(1,0,0,0,1,0,0,0,1);
	ModuleBase::Matrix3 temp2(1,0,0,0,1,0,0,0,1);
	bool flag = 0;
	int order = 0;
	int now = 0;

	symop[0] = iden;	//identity (the trivial element)
	nop = 1;

	//take all generators
	for(int i = 0; i < ngen; ++i)
	{
//		std::cout<<"\n symgen = "<<i<<std::endl;
		sig = symgen[i];
		flag = 1;
		for(int j = 0; j < nop; ++j)
		{
			if(symop[j] ==  sig)
			{
				flag = 0;
				break;
			}
		}
		if(flag == 0)
		{
			continue;
		}

		//Determine the order of the operation
		temp1 = sig;
		for(int j = 1; j < 100; ++j)
		{
			order = j;
			if(temp1 == iden)
			{
				break;
			}
			temp1= sig * temp1;
		//	temp1 = temp2;
		}
//		std::cout<<"\n order = "<<order<<std::endl;
		now = nop;
		for(int j = 0; j < nop; ++j)
		{	
			temp1 = symop[j];
		//	std::cout<<"\n j = "<<j<<std::endl<<"================================================="<<std::endl;
		//	out.printM3("temp1",temp1);
			for(int k = 1; k < order; ++k)
			{
				temp1 = sig * temp1;
			//	std::cout<<"\n k = "<<k<<std::endl<<"================================================="<<std::endl;
			//	out.printM3("temp1",temp1);
			//	out.printM3("temp2",temp2);

				//if(i==0 && j==0 && k==1){
					//std::cout<<"sig:"<<std::endl;
					//for(l=0; l<3; l++)
						//for(m=0; m<3; m++) std::cout<<sig[l][m]<<"\t";
					//std::cout<<std::endl;
				//}

				for(int l = 0; l < nop; ++l)
				{
					temp2 = symop[l] * temp1;
				//	std::cout<<"\n l = "<<l<<std::endl<<"================================================="<<std::endl;
				//	out.printM3("temp1",temp1);
				//	out.printM3("temp2",temp2);
					flag = 1;
					for(int m = 0; m < now; ++m)
					{
						if(symop[m] == temp2)
						{
							flag = 0;
							break;
						}
					}
					if(flag == 0)
					{
						continue;	//the new-found element has already existed.
					}

					//if(i==0 && j==0 && k==1 && l==0){
						//std::cout<<"symop[0]:"<<std::endl;
						//for(m=0; m<3; m++)
							//for(t=0; t<3; t++) std::cout<<temp2[m][t]<<"\t";
						//std::cout<<std::endl;
					//}

					++now;	//the number of elements we found
					if(now > 48)
					{
						std::cout<<"\n a: now= "<<now<<std::endl;
						std::cout<<"\n There are too many symmetrical matrices!"<<std::endl;
						return;
					}
					symop[now - 1] = temp2;
				//	std::cout<<"\n symop[now]:  "<<now<<std::endl;
				//	out.printM3("",temp2);
				}
			}
			if(j == 0)
			{
				n = now;
			}
		}

		m1 = nop;
		m2 = now;
		for(int j = 1; j < 50; ++j)
		{
			for(int k = nop; k < n; ++k)
			{
				for(int m = m1; m < m2; ++m)
				{	
					temp1 = symop[k] * symop[m];
					flag = 1;
					for(int l = 0; l < now; ++l)	
					{
						if(symop[l] == temp1)
						{
							flag = 0;
							break;
						}
					}
					if(flag == 0)	 
					{
						continue;	//the new-found element has already existed
					}

					++now;
					if(now > 48)
					{
						std::cout<<"\n b: now= "<<now<<std::endl;
						std::cout<<"\n There are too many symmetrical matrices!"<<std::endl;
						return;
					}
					symop[now - 1] = temp1;
				}
			}
			if(now == m2)
			{
				break;	//if no more new element could be found, stop the loop
			}
			m1 = m2;
			m2 = now;
		}
		nop = now;
//		for(int i = 0; i < nop; ++i)
//		{
//			std::cout<<"\n i = "<<i + 1<<std::endl;
//			out.printM3("",symop[i]);
//		}
	}
//	std::cout<<"\n in matrigen:"<<std::endl;
//	for(int i = 0; i < nop; ++i)
//	{
//		std::cout<<"\n i = "<<i<<std::endl;
//		out.printM3("",symop[i]);
//	}
}

//--------------------------------------------------------------
// set up all possible space group operators 
// (integer rotation matrices and nontrivial translations
// given in crystal coordinates) 
// of a lattice with some arbitrary basis (atomic arrangement).
//--------------------------------------------------------------
void Symmetry_Basic::setgroup(ModuleBase::Matrix3* symop, int &nop, const int &ibrav) const
{
	if(GlobalV::test_symmetry) ModuleBase::TITLE("Symmetry_Basic","setgroup");

	ModuleBase::Matrix3 symgen[3];

	ModuleBase::Matrix3 inv(-1,0,0,0,-1,0,0,0,-1);
	ModuleBase::Matrix3 r3d(0,1,0,0,0,1,1,0,0);
	ModuleBase::Matrix3 r6z(1,1,0,-1,0,0,0,0,1);
	ModuleBase::Matrix3 r2hex(1,0,0,-1,-1,0,0,0,-1);
	ModuleBase::Matrix3 r2tri(-1,0,0,0,0,-1,0,-1,0);
	ModuleBase::Matrix3 r4zp(0,1,0,-1,0,0,0,0,1);
	ModuleBase::Matrix3 r2yp(-1,0,0,0,1,0,0,0,-1);
	ModuleBase::Matrix3 r4zbc(0,0,-1,1,1,1,0,-1,0);
	ModuleBase::Matrix3 r4zfc(1,0,-1,1,0,0,1,-1,0);
	ModuleBase::Matrix3 r2zp(-1,0,0,0,-1,0,0,0,1);
	ModuleBase::Matrix3 r2ybc(0,0,1,-1,-1,-1,1,0,0);
	ModuleBase::Matrix3 r2zbc(0,1,0,1,0,0,-1,-1,-1);
	ModuleBase::Matrix3 r2ybas(0,-1,0,-1,0,0,0,0,-1);
	ModuleBase::Matrix3 r2yfc(0,-1,1,0,-1,0,1,-1,0);
	ModuleBase::Matrix3 r2zfc(0,1,-1,1,0,-1,0,0,-1);

	//the pure translation lattice (bravais lattice) has some maximum symmetry
	//set first up the point group operations for this symmetry.
	symgen[0] = inv;

	if(ibrav == 14)
	{
		this->matrigen(symgen, 1, symop, nop);
	}
	else if(ibrav == 13)
	{
		symgen[1] = r2ybas;
		this->matrigen(symgen, 2, symop, nop);
	}
	else if(ibrav == 12)
	{
		symgen[1] = r2yp;
		this->matrigen(symgen, 2, symop, nop);
	}
	else if(ibrav == 11)
	{
		symgen[1] = r2zp;
		symgen[2] = r2ybas;
		this->matrigen(symgen, 3, symop, nop);
	}
	else if(ibrav == 10)
	{
		symgen[1] = r2zfc;
		symgen[2] = r2yfc;
		this->matrigen(symgen, 3, symop, nop);
	}
	else if(ibrav == 9)
	{
		symgen[1] = r2zbc;
		symgen[2] = r2ybc;
		this->matrigen(symgen, 3, symop, nop);
	}
	else if(ibrav == 8)
	{
		symgen[1] = r2zp;
		symgen[2] = r2yp;
		this->matrigen(symgen, 3, symop, nop);
	}
	else if(ibrav == 7)
	{
		symgen[1] = r2tri;
		symgen[2] = r3d;
		this->matrigen(symgen, 3, symop, nop);
	}
	else if(ibrav == 6)
	{
		symgen[1] = r4zbc;
		symgen[2] = r2ybc;
		this->matrigen(symgen, 3, symop, nop);
	}
	else if(ibrav == 5)
	{
		symgen[1] = r4zp;
		symgen[2] = r2yp;
		this->matrigen(symgen, 3, symop, nop);
	}
	else if(ibrav == 4)
	{
		symgen[1] = r6z;
		symgen[2] = r2hex;
		this->matrigen(symgen, 3, symop, nop);
	}
	else if(ibrav == 3)
	{
		symgen[1] = r3d;
		symgen[2] = r4zfc;
		this->matrigen(symgen, 3, symop, nop);
	}
	else if(ibrav == 2)
	{
		symgen[1] = r3d;
		symgen[2] = r4zbc;
		this->matrigen(symgen, 3, symop, nop);
	}
	else if(ibrav == 1)
	{
		symgen[1] = r3d;
		symgen[2] = r4zp;
		this->matrigen(symgen, 3, symop, nop);
	}

	if(test_brav)
	{
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"ROTATION MATRICES",nop);
	}

	if(GlobalV::test_symmetry > 1)
	{
		GlobalV::ofs_running<<" THERE ARE " << nop << " ROTATION MATRICES FOR THE PURE BRAVAIS LATTICE"<<std::endl;
		GlobalV::ofs_running<<"    E11 E12 E13 E21 E22 E23 E31 E32 E33"<<std::endl;
		for(int i = 0; i < nop; ++i)
		{
			GlobalV::ofs_running << " " 
			<< std::setw(3) << i + 1
			<< std::setw(4) << symop[i].e11
			<< std::setw(4) << symop[i].e12
			<< std::setw(4) << symop[i].e13
			<< std::setw(4) << symop[i].e21
			<< std::setw(4) << symop[i].e22
			<< std::setw(4) << symop[i].e23
			<< std::setw(4) << symop[i].e31
			<< std::setw(4) << symop[i].e32
			<< std::setw(4) << symop[i].e33 << std::endl;
//			out.printM3("",symop[i]);
//			GlobalV::ofs_running<<std::endl;
		}
	}

	return;
}

int Symmetry_Basic::subgroup(const int& nrot, const int& ninv, const int& nc2, const int& nc3, const int& nc4, const int& nc6,
    const int& ns1, const int& ns3, const int& ns4, const int& ns6)
{
    if (nrot > 24)
    {
        if (ninv)
        {
            // if (nc2 >= 7 && nc3 >= 2 && nc6 >= 2 && ns1 >= 7 && ns3 >= 2 && ns6 >= 2) return 27;//D_6h
            if (nc2 >= 3 && nc3 >= 8 && ns1 >= 3 && ns6 >= 8) return 29; //T_h
        }
        else
        {
            if (nc2 >= 9 && nc3 >= 8 && nc4 >= 6) return 30;    //O
            if (nc2 >= 3 && nc3 >= 8 && ns1 >= 6 && ns4 >= 6) return 31; //T_d
        }
    }
    if (nrot > 16)//not else if: nrot>24 can also fall in this part and below
    {
        if (ninv && nc2 >= 5 && nc4 >= 2 && ns1 >= 5 && ns4 >= 2) return 20; //D_4h
    }
    if (nrot > 12)
    {
        if (ninv)
        {
            if (nc2 >= 1 && nc3 >= 2 && nc6 >= 2 && ns1 >= 1 && ns3 >= 2 && ns6 >= 2) return 23;//C_6h
            if (nc2 >= 3 && nc3 >= 2 && ns1 >= 3 && ns6 >= 2)return 13;//D_3d
        }
        else
        {
            if (nc2 >= 3 && nc3 >= 8)return 28; //T
            if (nc2 >= 3 && nc3 >= 2 && ns1 >= 4 && ns3 >= 2) return 26;//D_3h
            if (nc2 >= 1 && nc3 >= 2 && nc6 >= 2 && ns1 >= 6) return 25;//C_6v
            if (nc2 >= 7 && nc3 >= 2 && nc6 >= 2) return 24;//D_6
        }
    }
    if (nrot > 8)
    {
        if (ninv)
        {
            if (nc2 >= 1 && nc4 >= 2 && ns1 >= 1 && ns4 >= 2) return 16;//C_4h
            if (nc2 >= 3 && ns1 >= 3)return 8;//D_2h
        }
        else
        {
            if (nc2 >= 3 && ns1 >= 2 && ns4 >= 2)return 19;//D_2d
            if (nc2 >= 1 && nc4 >= 2 && ns1 >= 4) return 18;//C_4v
            if (nc2 >= 5 && nc4 >= 2)return 17;//D_4
        }
    }
    if (nrot > 6)
    {
        if (nc3 >= 2 && ns1 >= 1 && ns3 >= 2)return 22;//C_3h
        if (nc2 >= 1 && nc3 >= 2 && nc6 >= 2)return 21;//C_6
        if (nc3 >= 2 && ns1 >= 3)return 12;//C_3v
        if (nc2 >= 3 && nc3 >= 2)return 11;//D_3
        if (ninv && nc3 >= 2 && ns3 >= 2)return 10;//S_6
    }
    if (nrot > 4)
    {
        if (nc2 >= 1 && ns4 >= 2)return 15;//S_4
        if (nc2 >= 1 && nc4 >= 2)return 14;//C_4
        if (nc2 >= 1 && ns1 >= 2)return 7;//C_2v
        if (nc2 >= 3)return 6;//D_2
        if (ninv && nc2 >= 1 && ns1 >= 1)return 5;//C_2h
    }
    if (nrot > 3)
    {
        if (nc3 >= 2)return 9;//C_3
    }
    if (nrot > 2)
    {
        if (ns1 >= 1)return 4;//C_1h
        if (nc2 >= 1)return 3;//C_2
        if (ninv)return 2;//S_2
    }
    return 1;//C_1
}
void Symmetry_Basic::pointgroup(const int& nrot, int& pgnumber, std::string& pgname, const ModuleBase::Matrix3* gmatrix, std::ofstream& ofs_running)
{
	//-------------------------------------------------------------------------
	//return the name of the point group
	//the "name" (Schoenflies mark) of the group defined by following key:
	//       1 --> C_1       9 --> C_3      17 --> D_4      25 --> C_6v     *
	//       2 --> S_2      10 --> S_6      18 --> C_4v     26 --> D_3h     *
	//       3 --> C_2      11 --> D_3      19 --> D_2d     27 --> D_6h     *
	//       4 --> C_1h     12 --> C_3v     20 --> D_4h     28 --> T        *
	//       5 --> C_2h     13 --> D_3d     21 --> C_6      29 --> T_h      *
	//       6 --> D_2      14 --> C_4      22 --> C_3h     30 --> O        *
	//       7 --> C_2v     15 --> S_4      23 --> C_6h     31 --> T_d      *
	//       8 --> D_2h     16 --> C_4h     24 --> D_6      32 --> O_h      *
	//-------------------------------------------------------------------------

	//there are four trivial cases which could be easily determined
	//because the number of their elements are exclusive
    if (GlobalV::test_symmetry) ModuleBase::TITLE("Symmetry_Basic", "pointgroup");

    std::vector<std::string> pgdict = { "none", "C_1", "S_2", "C_2", "C_1h", "C_2h", "D_2", "C_2v", "D_2h", "C_3", "S_6", "D_3", "C_3v", "D_3d", "C_4", "S_4", "C_4h", "D_4", "C_4v", "D_2d", "D_4h",
    "C_6", "C_3h", "C_6h", "D_6", "C_6v", "D_3h", "D_6h", "T", "T_h", "O", "T_d", "O_h" };

	if(nrot == 1)
	{
		pgnumber = 1;
		pgname="C_1";
		return;
	}
	if(nrot == 3)
	{
		pgnumber = 9;
		pgname="C_3";
		return;
	}
	if(nrot == 16)
	{
		pgnumber = 20;
		pgname="D_4h";
		return;
	}
	if(nrot == 48)
	{
		pgnumber = 32;
		pgname="O_h";
		return;
	}
	
	//-------------------------------------------------------------------------------
	//all other groups need further investigations and detailed analysis
	//first determine the type of elements and count them
	//Possible elements are E, I, C_2, C_3, C_4, C_6 and S_1, S_3, S_4, S_6 (S_1 = m)
	//The type of a symmetry operation can be identified simply by
	//calculating the trace and the determinant of the rotation matrix. The
	//combination of these two quantities is specific for specific elements:
	//-------------------------------------------------------------------------------

	// Element:         E    I  C_2  C_3  C_4  C_6  S_1  S_6  S_4  S_3
	// Trace:          +3   -3   -1    0   +1   +2   +1    0   -1   -2
	// Determinant:    +1   -1   +1   +1   +1   +1   -1   -1   -1   -1

	int trace = 0;
	int det = 0;
	int ninv = 0;

	int nc2 = 0;
	int nc3 = 0;
	int nc4 = 0;
	int nc6 = 0;
	int ns1 = 0;
	int ns3 = 0; //mohan add 2012-01-15
	int ns4 = 0;
	int ns6 = 0; //mohan add 2012-01-15

//	GlobalV::ofs_running << " " << std::setw(5) << "NROT" << std::setw(15) << "TRACE" << std::setw(15) << "DET" << std::endl;
	for(int i = 0; i < nrot; ++i)
	{
		//calculate the trace of a matrix
		trace = int(gmatrix[i].e11+gmatrix[i].e22+gmatrix[i].e33);
		//calculate the determinant of a matrix
		det = int(gmatrix[i].Det());

//		GlobalV::ofs_running << " " << std::setw(5) << i+1 << std::setw(15) << trace << std::setw(15) << det << std::endl;

		if(trace == 3)
		{
			continue;	//found unity operator (trivial)
		}
		//found inversion
		if(trace == -3)
		{
			ninv = 1;
			continue;
		}

		if(trace == -1 && det == 1) ++nc2;
		else if(trace == 0 && det == 1) ++nc3;
		else if(trace == 1 && det == 1) ++nc4;
		else if(trace == 2 && det == 1) ++nc6;
		else if(trace == 1 && det == -1) ++ns1;
		else if(trace == 0 && det == -1) ++ns6; //mohan add 2012-01-15
		else if(trace == -1 && det == -1) ++ns4;
		else if(trace == -2 && det == -1) ++ns3; //mohan add 2012-01-15
	}

        if(test_brav)
	{
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"C2",nc2);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"C3",nc3);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"C4",nc4);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"C6",nc6);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"S1",ns1);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"S3",ns3);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"S4",ns4);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"S6",ns6);
	}
	
	if(nrot == 2)
	{
		if(ninv == 1)
		{
			pgnumber = 2;
			pgname="S_2";
			return;
		}
		if(nc2 == 1)
		{
			pgnumber = 3;
			pgname="C_2";
			return;
		}
		if(ns1 == 1)
		{
			pgnumber = 4;
			pgname="C_1h";
			return;
		}
	}
	if(nrot == 4)
	{
		if(ninv == 1)
		{
			pgnumber = 5;
			pgname="C_2h";
			return;
		}
		if(nc2 == 3)
		{
			pgnumber = 6;
			pgname="D_2";
			return;
		}
		if(ns1 == 2) 
		{
			pgnumber = 7;
			pgname="C_2v";
			return;
		}
		if(nc4 == 2)
		{
			pgnumber = 14;
			pgname="C_4";
			return;
		}
		if(ns4 == 2)
		{
			pgnumber = 15;
			pgname="S_4";
			return;
		}
	}
	if(nrot == 6)
	{
		if(ninv == 1)
		{
			pgnumber = 10;
			pgname="S_6";
			return;
		}
		if(nc2 == 3)
		{
			pgnumber = 11;
			pgname="D_3";
			return;
		}
		if(ns1 == 3)
		{
			pgnumber = 12;
			pgname="C_3v";
			return;
		}
		if(nc2 == 1)
		{
			pgnumber = 21;
			pgname="C_6";
			return;
		}
		if(ns1 == 1)
		{
			pgnumber = 22;
			pgname="C_3h";
			return;
		}
	}
	if(nrot == 8)
	{
		if(ns1 == 3)
		{
			pgnumber = 8;
			pgname="D_2h";
			return;
		}
		if(ns1 == 1)
		{
			pgnumber = 16;
			pgname="C_4h";
			return;
		}
		if(ns1 == 0)
		{
			pgnumber = 17;
			pgname="D_4";
			return;
		}
		if(ns1 == 4)
		{
			pgnumber = 18;
			pgname="C_4v";
			return;
		}
		if(ns1 == 2)
		{
			pgnumber = 19;
			pgname="D_2d";
			return;
		}
	}
	if(nrot == 12)
	{
		if(ns1 == 3)
		{
			pgnumber = 13;
			pgname="D_3d";
			return;
		}
		if(ns1 == 1)
		{
			pgnumber = 23;
			pgname="C_6h";
			return;
		}
		if(nc2 == 7)
		{
			pgnumber = 24;
			pgname="D_6";
			return;
		}
		if(ns1 == 6)
		{
			pgnumber = 25;
			pgname="C_6v";
			return;
		}
		if(ns1 == 4)
		{
			pgnumber = 26;
			pgname="D_3h";
			return;
		}
		if(nc3 == 8)
		{
			pgnumber = 28;
			pgname="T";
			return;
		}
	}
	if(nrot == 24)
	{
		if(nc6 == 2)
		{
			pgnumber = 27;
			pgname="D_6h";
			return;
		}
		if(ninv == 1)
		{
			pgnumber = 29;
			pgname="T_h";
			return;
		}
		if(nc4 == 6)
		{
			pgnumber = 30;
			pgname="O";
			return;
		}
		if(ns4 == 6)
		{
			pgnumber = 31;
			pgname="T_d";
			return;
		}
	}
    GlobalV::ofs_running << "\n WARNING: Symmetry operations cannot completely constitute a point group.\n\
    It'll be better to try another `symmetry_prec`.\n  Now search the subgroups ..." << std::endl;
    pgnumber = this->subgroup(nrot, ninv, nc2, nc3, nc4, nc6, ns1, ns3, ns4, ns6);
    pgname = pgdict[pgnumber];
    this->valid_group = false;
    return;
}


void Symmetry_Basic::rotate( ModuleBase::Matrix3 &gmatrix, ModuleBase::Vector3<double> &gtrans, 
		int i, int j, int k, // FFT grid index.
		const int nr1, const int nr2, const int nr3, // dimension of FFT grid. 
		int &ri, int &rj, int &rk)
{
	static ModuleBase::Matrix3 g;
	g.e11 = gmatrix.e11;
	g.e21 = gmatrix.e21 * (double)nr1 / (double)nr2;
	g.e31 = gmatrix.e31 * (double)nr1 / (double)nr3;
	g.e12 = gmatrix.e12 * (double)nr2 / (double)nr1;
	g.e22 = gmatrix.e22;
	g.e32 = gmatrix.e32 * (double)nr2 / (double)nr3;
	g.e13 = gmatrix.e13 * (double)nr3 / (double)nr1;
	g.e23 = gmatrix.e23 * (double)nr3 / (double)nr2;
	g.e33 = gmatrix.e33;

	ri = int(g.e11 * i + g.e21 * j + g.e31 * k) + (int)(gtrans.x *  nr1);
	if (ri < 0)
	{
		ri += 10 * nr1;
	}
	ri = ri%nr1;
	rj = int(g.e12 * i + g.e22 * j + g.e32 * k) + (int)(gtrans.y  * nr2);
	if (rj < 0)
	{
		rj += 10 * nr2;
	}
	rj = rj%nr2;
	rk = int(g.e13 * i + g.e23 * j + g.e33 * k) + (int)(gtrans.z  * nr3);
	if (rk < 0)
	{
		rk += 10 * nr3;
	}
	rk = rk%nr3;
	return;
}

// atom ordering for each atom type 
// by a "weighted function" f
// (instead of ordering by x, y, z directly)
void Symmetry_Basic::atom_ordering_new(double *posi, const int natom, int *subindex) const
{
	//order the atomic positions inside a supercell by a unique ordering scheme	
	subindex[0] = 0;

	if(natom == 1)
	{
		//if there is only one atom, it is not necessary to order
		return;
	}

	std::vector<double> tmpx(natom);
	std::vector<double> tmpy(natom);
	std::vector<double> tmpz(natom);
	for(int i=0; i<natom; i++)
	{
		tmpx[i] = posi[i*3];
		tmpy[i] = posi[i*3+1];
		tmpz[i] = posi[i*3+2];
	}
	double x_max = *max_element(tmpx.begin(),tmpx.end());
	double x_min = *min_element(tmpx.begin(),tmpx.end());
	double y_max = *max_element(tmpy.begin(),tmpy.end());
	double y_min = *min_element(tmpy.begin(),tmpy.end());
	double z_max = *max_element(tmpz.begin(),tmpz.end());
	double z_min = *min_element(tmpz.begin(),tmpz.end());

	double*  weighted_func = new double[natom];
	
	//the first time: f(x, y, z)
	for(int i=0; i<natom; i++)
	{
		weighted_func[i]=1/epsilon/epsilon*tmpx[i]+1/epsilon*tmpy[i]+tmpz[i];
	}
	ModuleBase::heapsort(natom, weighted_func, subindex);
	this->order_atoms(posi, natom, subindex);
	for(int i=0; i<natom; i++)
	{
		tmpx[i] = posi[i*3];
		tmpy[i] = posi[i*3+1];
		tmpz[i] = posi[i*3+2];
	}

	//the second time: f(y, z) for fixed x	
	for(int i=0; i<natom-1;)
	{
		int ix_right=i+1;	//right bound is no included
		while(ix_right<natom && equal(tmpx[ix_right],tmpx[i])) ++ix_right;
		int nxequal=ix_right-i;
		if(nxequal>1)	//need a new sort
		{
			subindex[0] = 0;
			for(int j=0; j<nxequal; ++j)
			{
				weighted_func[j]=1/epsilon*tmpy[i+j]+tmpz[i+j];
			}
			ModuleBase::heapsort(nxequal, weighted_func, subindex);
			this->order_atoms(&posi[i*3], nxequal, subindex);
		}
		i=ix_right;
	}
	
	delete[] weighted_func;
	return;
}

void Symmetry_Basic::test_atom_ordering(double *posi, const int natom, int *subindex) const
{
	//an interface to test a protected function
	this->atom_ordering_new(posi, natom, subindex);
}
}