#include "../src_io/output.h"
#include "read_pp.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring> // Peize Lin fix bug about strcpy 2016-08-02


using namespace std;

Pseudopot_upf::Pseudopot_upf()
{
	this->els = new string[1];
	this->lchi = new int[1];
	this->oc = new double[1];

	this->r = new double[1];
	this->rab = new double[1];
	this->vloc = new double[1];

	this->kkbeta = new int[1];
	this->lll = new int[1];

	this->rho_at = new double[1];
	this->rho_atc = new double[1];

	this->nn = new int[1];//zhengdy-soc
	this->jchi = new double[1];
	this->jjj = new double[1];

	functional_error = 0;//xiaohui add 2015-03-24
}

Pseudopot_upf::~Pseudopot_upf()
{
	delete [] els; 
	delete [] lchi;
	delete [] oc;

	delete [] r;    //mesh_1
	delete [] rab;  //mesh_2
	delete [] vloc;  //local_1

	delete [] kkbeta; // nl_1
	delete [] lll; // nl_2

	delete [] rho_at;// psrhoatom_1
	delete [] rho_atc;

	delete [] nn;
	delete [] jjj;
	delete [] jchi;
}

int Pseudopot_upf::init_pseudo_reader(const string &fn)
{
    TITLE("Pseudopot_upf","init");
    // First check if this pseudo-potential has spin-orbit information
    ifstream ifs(fn.c_str(), ios::in);

	// can't find the file.
	if (!ifs)
    {
        return 1;
    }

    if(global_pseudo_type=="auto") //zws
	{
		set_pseudo_type(fn);
	}

	// read in the .UPF type of pseudopotentials
	if(global_pseudo_type=="upf")
	{
		int info = read_pseudo_upf(ifs);
		return info;
	}
	// read in the .vwr type of pseudopotentials
	else if(global_pseudo_type=="vwr")
	{
		int info = read_pseudo_vwr(ifs);
		return info;
	}
	else if(global_pseudo_type=="upf201")
	{
		int info = read_pseudo_upf201(ifs);
		return info;
	}

	return 0;
}


//----------------------------------------------------------
// setting the type of the pseudopotential file
//----------------------------------------------------------
int Pseudopot_upf::set_pseudo_type(const string &fn) //zws add
{
    ifstream pptype_ifs(fn.c_str(), ios::in);
    string dummy;
	string strversion;

	if (pptype_ifs.good())
	{
		getline(pptype_ifs,dummy);

		stringstream wdsstream(dummy);
		getline(wdsstream,strversion,'"');
		getline(wdsstream,strversion,'"');

		if ( trim(strversion) == "2.0.1" )
		{
			global_pseudo_type = "upf201";
		}
		else
		{
			global_pseudo_type = "upf";
		}
	}
	return 0;
}

string& Pseudopot_upf::trim(string &in_str)
{
    static const string deltri = " \t" ; // delete tab or space
    string::size_type position = in_str.find_first_of(deltri, position);
    if (position == string::npos)
	{
        return in_str;
	}
    return trim(in_str.erase(position, 1) );
}

string Pseudopot_upf::trimend(string &in_str)
{
    const string &deltri =" \t" ;
    string::size_type position = in_str.find_last_not_of(deltri)+1;
    string tmpstr=in_str.erase(position);
    return tmpstr.erase(0,tmpstr.find_first_not_of(deltri));
} //zws


int Pseudopot_upf::average_p(void)
{
	int error = 0;
	if(!this->has_so && LSPINORB) 
	{
		error++; 
		cout<<"warning_quit! no soc upf used for lspinorb calculation, error!"<<endl; 
		return error;
	}
	//WARNING_QUIT("average_p", "no soc upf used for lspinorb calculation, error!");

	if(!this->has_so || LSPINORB) 
	{
		return error; 
	}

	int new_nbeta = 0; //calculate the new nbeta
	for(int nb=0; nb< this->nbeta; nb++)
	{
		new_nbeta++;
		if(this->lll[nb] != 0 && abs(this->jjj[nb] - this->lll[nb] - 0.5) < 1e-6) //two J = l +- 0.5 average to one
		{
			new_nbeta--;
		}
	}

	this->nbeta = new_nbeta;
	matrix dion_new;
	dion_new.create(this->nbeta, this->nbeta);

	int old_nbeta=-1;
	for(int nb=0; nb<this->nbeta; nb++)
	{
		old_nbeta++;
		int l = this->lll[old_nbeta];
		int ind=0, ind1=0;
		if(l != 0)
		{
			if(abs(this->jjj[old_nbeta] - this->lll[old_nbeta] + 0.5) < 1e-6)
			{
				if(abs(this->jjj[old_nbeta+1]-this->lll[old_nbeta+1]-0.5)>1e-6) 
					WARNING_QUIT("average_p", "error beta function 1 !");
				ind = old_nbeta +1;
				ind1 = old_nbeta;
			}
			else
			{
				if(abs(this->jjj[old_nbeta+1]-this->lll[old_nbeta+1]+0.5)>1e-6)
					WARNING_QUIT("average_p", "error beta function 2 !");
				ind = old_nbeta;
				ind1 = old_nbeta +1;
			}
			double vion1 = ((l+1.0) * this->dion(ind,ind) + l * this->dion(ind1,ind1)) / (2.0*l+1.0);
			//average beta (betar)
			for(int ir = 0; ir<this->mesh;ir++)
			{
				this->beta(nb, ir) = 1.0 / (2.0 * l + 1.0) * 
						( (l + 1.0) * sqrt(this->dion(ind,ind) / vion1) *
						this->beta(ind, ir) + 
						l * sqrt(this->dion(ind1,ind1) / vion1) *
						this->beta(ind1, ir) ) ;
			}
			//average the dion matrix
			this->dion(nb, nb) = vion1;
			old_nbeta++;	
		}
		else
		{
			for(int ir = 0; ir<this->mesh;ir++)
				this->beta(nb, ir) = this->beta(old_nbeta, ir);
			this->dion(nb, nb) = this->dion(old_nbeta, old_nbeta);
		}
		this->lll[nb] = this->lll[old_nbeta]; //reset the lll index, ignore jjj index
	}

	//store the old dion and then recreate dion 
	for(int i=0;i<this->nbeta; i++)
	{
		for(int j=0;j<this->nbeta;j++)
		{
			dion_new(i,j) = this->dion(i,j);
		}
	}

	this->dion = dion_new;
//	this->dion.create(this->nbeta, this->nbeta);
//	for(int i=0;i<this->nbeta; i++)
//		for(int j=0;j<this->nbeta;j++)
//			this->dion(i,j) = dion_new(i,j);
	
	int new_nwfc = 0;
	for(int nb=0; nb<this->nwfc; nb++)
	{
		new_nwfc++;
		if(this->lchi[nb] != 0 && abs(this->jchi[nb] - this->lchi[nb] - 0.5)<1e-6)
		{
			new_nwfc--;
		}
	}

	this->nwfc = new_nwfc;
	int old_nwfc=0;
	for(int nb=0; nb<this->nwfc; nb++)
	{
		old_nwfc++;
		int l = this->lchi[old_nwfc];
		int ind=0, ind1=0;
		if(l!=0)
		{
			if(abs(this->jchi[old_nwfc] - this->lchi[old_nwfc] + 0.5) < 1e-6)
			{
				if(abs(this->jchi[old_nwfc+1]-this->lchi[old_nwfc+1]-0.5)>1e-6) 
				{error++; cout<<"warning_quit! error chi function 1 !"<<endl; return error;}
//					WARNING_QUIT("average_p", "error chi function 1 !");
				ind = old_nwfc +1;
				ind1 = old_nwfc;
			}
			else
			{
				if(abs(this->jchi[old_nwfc+1]-this->lchi[old_nwfc+1]+0.5)>1e-6)
				{error++; cout<<"warning_quit! error chi function 2 !"<<endl; return error;}
//					WARNING_QUIT("average_p", "error chi function 2 !");
				ind = old_nwfc;
				ind1 = old_nwfc +1;
			}
			//average chi
			for(int ir = 0; ir<this->mesh;ir++)
			{
				this->chi(nb, ir) = 1.0 / (2.0 * l + 1.0) *
					( (l+1.0)*this->chi(ind,ir) + (l*this->chi(ind1,ir)) );
				old_nwfc++;
			}
		}
		else{
			for(int ir = 0; ir<this->mesh;ir++)
				this->chi(nb, ir) = this->chi(old_nwfc, ir);
		}
		this->lchi[nb] = this->lchi[old_nwfc]; //reset lchi index
	}
	this->has_so = 0;	
	return error;
}

// Peize Lin add for bsse 2021.04.07
void Pseudopot_upf::set_empty_element(void)
{
	this->zp = 0;
	for(int ir=0; ir<mesh; ++ir)
	{
		this->vloc[ir] = 0;
	}
	for(int i=0; i<nbeta; ++i)
	{
		for(int j=0; j<nbeta; ++j)
		{
			this->dion(i,j) = 0;
		}
	}
	for(int ir=0; ir<mesh; ++ir)
	{
		this->rho_at[ir] = 0;
	}
	return;
}
