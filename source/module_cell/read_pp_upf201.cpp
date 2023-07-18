#include "read_pp.h"
//int Number[2]; // added by zhangwenshuai

//qianrui rewrite it 2021-5-10
int Pseudopot_upf::read_pseudo_upf201(std::ifstream &ifs)
{

    std::string word;
    int ONCVPSP;
	//--------------------------------------
	//-              PP_HEADER             - 
	//--------------------------------------
	if(!ModuleBase::GlobalFunc::SCAN_BEGIN(ifs,"<PP_HEADER"))	ModuleBase::WARNING_QUIT("read_pseudo_upf201","Found no PP_HEADER");
	std::string *name=new std::string[50];
	std::string *val=new std::string[50];
	int nparameter;
	this->getnameval(ifs, nparameter, name, val);
	ONCVPSP = 1;
	
	for(int ip = 0 ; ip < nparameter; ++ip)
	{
		if(name[ip]=="generated")
		{
			//add something//
		}
		else if(name[ip]=="author"){}
		else if(name[ip]=="date"){}
		else if(name[ip]=="comment"){}
		else if(name[ip]=="element"){
			psd = val[ip];
		}
		else if(name[ip]=="pseudo_type"){
			pp_type = val[ip];
			if(pp_type!="NC") 
			{
				ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_header","unknown pseudo type");
			}
		}
		else if(name[ip]=="relativistic"){}
		else if(name[ip]=="is_ultrasoft"){
			if(val[ip]=="T" || val[ip]=="TRUE" || val[ip]=="True" || val[ip]=="true")
			{
				ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_header","ULTRASOFT PSEUDOPOTENTIAL IS NOT SUPPORTED !!!");
			}
		}
		else if(name[ip]=="is_paw"){
			if(val[ip]=="T" || val[ip]=="TRUE" || val[ip]=="True" || val[ip]=="true")
			{
				ModuleBase::WARNING_QUIT("Pseudopot_upf::read_pseudo_header","PAW PSEUDOPOTENTIAL IS NOT SUPPORTED !!!");
			}
		}
		else if(name[ip]=="is_coulomb"){}
		else if(name[ip]=="has_so"){
			if(val[ip]=="T" || val[ip]=="TRUE" || val[ip]=="True" || val[ip]=="true")
				has_so = true;
			else
				has_so = false;
		}
		else if(name[ip]=="has_wfc"){}
		else if(name[ip]=="has_gipaw"){
			ONCVPSP = 0;
		}
		else if(name[ip]=="paw_as_gipaw"){
			ONCVPSP = 0;
		}
		else if(name[ip]=="core_correction"){
			if(val[ip]=="T" || val[ip]=="TRUE" || val[ip]=="True" || val[ip]=="true")
				nlcc = true;
			else
				nlcc = false;
		}
		else if(name[ip]=="functional"){
			xc_func = val[ip];
		}
		else if(name[ip]=="z_valence"){
			zp = atoi(val[ip].c_str());
		}
		else if(name[ip]=="total_psenergy"){
			etotps = atof(val[ip].c_str());
		}
		else if(name[ip]=="wfc_cutoff"){}
		else if(name[ip]=="rho_cutoff"){}
		else if(name[ip]=="l_max"){
			lmax = atoi(val[ip].c_str());
		}
		else if(name[ip]=="l_max_rho"){}
		else if(name[ip]=="l_local"){}
		else if(name[ip]=="mesh_size"){
			mesh = atoi(val[ip].c_str());
		}
		else if(name[ip]=="number_of_wfc"){
			nwfc = atoi(val[ip].c_str());
		}
		else if(name[ip]=="number_of_proj"){
			nbeta = atoi(val[ip].c_str());
		}
		else
		{
			std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
			ModuleBase::WARNING("PP_HEADRER reading", warningstr);
		}
	}
			
	//--------------------------------------
	//-              PP_MESH               - 
	//--------------------------------------
	if(ONCVPSP == 0)
	{
		ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_MESH");
		this->getnameval(ifs, nparameter, name, val);
		for(int ip = 0 ; ip < nparameter; ++ip)
		{
			if(name[ip]=="dx"){}
			else if(name[ip]=="mesh"){}
			else if(name[ip]=="xmin"){}
			else if(name[ip]=="rmax"){}
			else if(name[ip]=="zmesh"){}
			else
			{
				std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
				ModuleBase::WARNING("PP_MESH reading", warningstr);
			}

		}
	}
	else
	{
		ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_MESH>");
	}

	if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_R", true, false))
	{
		ifs.ignore(150, '>');
		// ModuleBase::GlobalFunc::READ_VALUE(ifs, word); // type size columns
		this->read_pseudo_upf201_r(ifs);
	}
	else if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_R>"))
	{
		this->read_pseudo_upf201_r(ifs);
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_R>");

    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RAB", true, false))
	{
		ifs.ignore(150, '>');
		// ModuleBase::GlobalFunc::READ_VALUE(ifs, word); // type size columns
		this->read_pseudo_upf201_rab(ifs);
	}
	else if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RAB>"))
	{
		this->read_pseudo_upf201_rab(ifs);
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RAB>");
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_MESH>");

	//--------------------------------------
	//-              PP_NLCC               - 
	//--------------------------------------
	if (this->nlcc)
	{
		ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NLCC");
		ModuleBase::GlobalFunc::READ_VALUE(ifs, word);    // type size columns
		delete[] rho_atc;
		this->rho_atc = new double[mesh];
		ModuleBase::GlobalFunc::ZEROS(rho_atc, mesh);
		for (int ir = 0;ir < mesh;ir++)
		{
			ifs >> this->rho_atc[ir];
		}
		ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_NLCC>");
	}

	//--------------------------------------
	//-              PP_LOCAL              - 
	//--------------------------------------
	if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_LOCAL"))
	{
		ifs.ignore(150, '>'); // skip type, size, columns and so on.
	}
	else
	{
		ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_LOCAL>");	
	}
	delete[] vloc;
	this->vloc = new double[mesh];
	ModuleBase::GlobalFunc::ZEROS(vloc, mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->vloc[ir];
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_LOCAL>");

	//--------------------------------------
	//-            PP_NONLOCAL             - 
	//--------------------------------------
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_NONLOCAL>");
	delete[] kkbeta;
	delete[] lll;
	this->kkbeta = new int[nbeta];
	this->lll = new int[nbeta];
	this->beta.create(nbeta , mesh);
	this->dion.create(nbeta , nbeta);
	for(int ib=0;ib<nbeta;ib++)
	{
		ifs >> word; //number of beta
		this->getnameval(ifs, nparameter, name, val);
		bool read_cutoff_radius_index = false;
		for(int ip = 0 ; ip < nparameter; ++ip)
		{
			if(name[ip]=="type"){}
			else if(name[ip]=="size"){}
			else if(name[ip]=="columns"){}
			else if(name[ip]=="index"){}
			else if(name[ip]=="label"){}
			else if(name[ip]=="angular_momentum"){
				lll[ib] = atoi(val[ip].c_str());
			}
			else if(name[ip]=="cutoff_radius_index"){
				kkbeta[ib] = atoi (val[ip].c_str());
				read_cutoff_radius_index = true;
			}
			else if(name[ip]=="cutoff_radius"){}
			else if(name[ip]=="ultrasoft_cutoff_radius"){}
			else
			{
				std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
				ModuleBase::WARNING("PP_BETA reading", warningstr);
			}
		}
		if(!read_cutoff_radius_index) // if cutoff_radius_index is not set in pp, we set it as mesh.
		{
			kkbeta[ib] = mesh;
		}
		for (int ir=0;ir<mesh;ir++)
		{
			ifs >> this->beta(ib, ir);
		}
		ifs >> word; //number of beta
	}

	if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_DIJ", true, false))
	{
		ifs.ignore(150, '>');
		// ModuleBase::GlobalFunc::READ_VALUE(ifs, word);  // type size columns
		this->read_pseudo_upf201_dij(ifs);
	}
	else if ( ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_DIJ>"))
	{
		this->read_pseudo_upf201_dij(ifs);
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_DIJ>");
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_NONLOCAL>");

	//--------------------------------------
	//-            PP_PSWFC                - 
	//--------------------------------------
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_PSWFC>");
	delete[] els;
	delete[] lchi;
	delete[] oc;
	this->els = new std::string[nwfc];
	this->lchi = new int[nwfc];
	this->oc = new double[nwfc];
	ModuleBase::GlobalFunc::ZEROS(lchi, nwfc); // angular momentum of each orbital
	ModuleBase::GlobalFunc::ZEROS(oc, nwfc);//occupation of each orbital
	this->chi.create(this->nwfc, this->mesh);
	for(int iw=0;iw<nwfc;iw++)
	{
		ifs >> word; //number of chi
		this->getnameval(ifs, nparameter, name, val);
		for(int ip = 0 ; ip < nparameter; ++ip)
		{
			if(name[ip]=="type"){}
			else if(name[ip]=="size"){}
			else if(name[ip]=="columns"){}
			else if(name[ip]=="index"){}
			else if(name[ip]=="label"){
				els[iw] = val[ip];
			}
			else if(name[ip]=="l"){
				lchi[iw] = atoi(val[ip].c_str());
			}
			else if(name[ip]=="occupation"){
				oc[iw] = atof(val[ip].c_str());
			}
			else if(name[ip]=="n"){}
			else if(name[ip]=="pseudo_energy"){}
			else if(name[ip]=="cutoff_radius"){}
			else if(name[ip]=="ultrasoft_cutoff_radius"){}
			else
			{
				std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
				ModuleBase::WARNING("PP_CHI reading", warningstr);
			}
		}
		for (int ir=0;ir<mesh;ir++)
		{
			ifs >> this->chi(iw, ir);
		}
		ifs >> word; //number of chi
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_PSWFC>");

	//--------------------------------------
	//-          PP_RHOATOM                - 
	//--------------------------------------
	if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RHOATOM", true, false))
	{
    	ifs.ignore(150, '>');
		// ModuleBase::GlobalFunc::READ_VALUE(ifs, word); // type size columns
		this->read_pseudo_upf201_rhoatom(ifs);
	}
	else if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_RHOATOM>"))
	{
		this->read_pseudo_upf201_rhoatom(ifs);
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_RHOATOM>");

	//--------------------------------------
	//-          PP_SPIN_ORB               - 
	//--------------------------------------
	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<PP_SPIN_ORB>", true, false);
	//added by zhengdy-soc
	delete[] this->jchi;
	delete[] this->jjj;
	delete[] this->nn;
	this->jchi = new double [nwfc];
	this->jjj = new double [nbeta];
	this->nn = new int [nwfc];
	ModuleBase::GlobalFunc::ZEROS(jchi,nwfc);
	ModuleBase::GlobalFunc::ZEROS(jjj,nbeta);
	ModuleBase::GlobalFunc::ZEROS(nn,nwfc);

	for(int round=0;round<2;round++)
	{
		ifs>>word;
		if(word=="<PP_RELBETA.1")
		{
			for(int nb = 0;nb<nbeta;nb++)
			{
				if(nb!=0) ifs>>word; //RELBETA
				this->getnameval(ifs, nparameter, name, val);
				for(int ip = 0 ; ip < nparameter; ++ip)
				{
					if(name[ip]=="index"){}
					else if(name[ip]=="lll"){
						lll[nb] = atoi(val[ip].c_str());
					}
					else if(name[ip]=="jjj"){
						jjj[nb] = atof(val[ip].c_str());
					}
					else
					{
						std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
						ModuleBase::WARNING("PP_RELBETA reading", warningstr);
					}
				}
			}
		}
		else if(word=="<PP_RELWFC.1")
		{
			for(int nw = 0;nw<nwfc;nw++)
			{
				if(nw!=0) ifs>>word;     //RELWFC
				this->getnameval(ifs, nparameter, name, val);
				for(int ip = 0 ; ip < nparameter; ++ip)
				{
					if(name[ip]=="index"){}
					else if(name[ip]=="els"){
						els[nw] = val[ip];
					}
					else if(name[ip]=="nn"){
						nn[nw] = atoi(val[ip].c_str());
					}
					else if(name[ip]=="lchi"){
						lchi[nw]=atoi(val[ip].c_str());
					}
					else if(name[ip]=="jchi"){
						jchi[nw]=atof(val[ip].c_str());
					}
					else if(name[ip]=="oc"){
						oc[nw] = atof(val[ip].c_str());
					}
					else
					{
						std::string warningstr = name[ip] + " is not read in. Please add this parameter in read_pp_upf201.cpp if needed.";
						ModuleBase::WARNING("PP_RELWFC reading", warningstr);
					}
				}
			}
		}
		else if(round==0)
		{
			this->has_so = 0;
			//	std::cout<<"ignore SPIN_ORB part!"<<std::endl;
			break;
		}
	}
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</PP_SPIN_ORB>", false);
	if (mesh%2 == 0)
	{
		mesh -= 1;
	}
	
	ModuleBase::GlobalFunc::SCAN_END(ifs, "</UPF>", false);
	delete []name;
	delete []val;
	
	if(GlobalV::DFT_FUNCTIONAL!="default")
	{
		std::string xc_func1 = GlobalV::DFT_FUNCTIONAL;
		transform(xc_func1.begin(), xc_func1.end(), xc_func1.begin(), (::toupper));
		if(xc_func1 != xc_func)
		{
			functional_error = 1;

			std::cout << " dft_functional readin is: " << GlobalV::DFT_FUNCTIONAL << std::endl;
			std::cout << " dft_functional in pseudopot file is: " << xc_func << std::endl;
			GlobalV::ofs_warning << " dft_functional readin is: " << GlobalV::DFT_FUNCTIONAL << std::endl;
			GlobalV::ofs_warning << " dft_functional in pseudopot file is: " << xc_func << std::endl;
		}
	}
	return 0;

}
void Pseudopot_upf:: getnameval(std::ifstream& ifs,int &n, std::string * name, std::string *val)
{
	std::string txt,word;
	//get long txt
	ifs>>txt;
	while(ifs>>word)
	{
		size_t wl= word.length()-1;
		txt=txt+" "+word;
		if(word.substr(wl,1)==">")
		{
			break;
		}
	}

	//count number of parameters according to "="
	size_t pos=0;
	n=0;
	while(1)
	{
		pos = txt.find("=",pos);
		if(pos == std::string::npos) break;
		pos++;
		n++;
	}

	//get name & value
	pos=0;
	size_t pos2, ll;
	for(int i = 0; i < n; ++i)
	{
		pos2 = txt.find("=",pos);
		for(; pos2 > pos ; --pos2)//There may be a space before "=";
		{
			if(txt.substr(pos2-1,1) != " ") break;
		}
		ll=pos2-pos;
		name[i] = txt.substr(pos,ll);
		//std::cout<<i<<" "<<name[i]<<std::endl;
		std::string mark;
		bool findmark=false;
		for(int j = 0; j < 100; ++j)//The mark can be ' or " or .
		{
			mark = txt.substr(pos2,1);
			pos2++;
			if(mark=="\""||mark=="\'"||mark==".")
			{
				findmark = true;
				break;
			}
		}
		if(!findmark) ModuleBase::WARNING_QUIT("read_upf201",
		"The values are not in \' or \". Please improve the program in read_pp_upf201.cpp");
		pos = pos2;
		pos2 = txt.find(mark,pos);
		ll=pos2-pos;
		std::string tmpval = txt.substr(pos,ll);
		tmpval = trim(tmpval);
		val[i]=tmpval;
		pos=pos2+1;
		for(int j = 0; j < 100; ++j)
		{
			if(txt.substr(pos,1)==" " || txt.substr(pos,1)==",")
				pos++;
			else
				break;
		}
		//std::cout<<name[i]<<"=\""<<val[i]<<"\""<<std::endl;
	}
	return;
}

/*void Pseudopot_upf::get_char( std::string ss)
{
    int i, q;
    //char b[1]; //LiuXh 20171109
    char b='\"'; //LiuXh 20171109
    q =0;
    //strcpy(b,"\""); //LiuXh 20171109

    for(i=0;i<200;i++)
    {
        //if(ss[i]== b[0]) //LiuXh 20171109
        if(ss[i]== b) //LiuXh 20171109
        {
           Number[q] = i;
           q++;
        }

    }

    return;
}*/
void Pseudopot_upf::read_pseudo_upf201_r(std::ifstream &ifs)
{
	delete[] r;
	delete[] rab;
	assert(mesh>0);
	this->r = new double[mesh];
	this->rab = new double[mesh];
	ModuleBase::GlobalFunc::ZEROS(r,mesh);
	ModuleBase::GlobalFunc::ZEROS(rab,mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->r[ir];
	}
}
void Pseudopot_upf::read_pseudo_upf201_rab(std::ifstream &ifs)
{
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->rab[ir];
	}
}
void Pseudopot_upf::read_pseudo_upf201_dij(std::ifstream &ifs)
{
	this->nd = nbeta * nbeta;
	for(int i=0;i<nbeta;i++)
	{
		for(int j=0;j<nbeta;j++)
		{
			ifs >> dion(i,j);
			// if ( i != j  && dion(i,j) != 0.0 )
			// {
			// 	ModuleBase::WARNING_QUIT("read_pseudo_upf201","Error: for i != j, Dij of Pseudopotential must be 0.0");
			// }
		}
	}
}
void Pseudopot_upf::read_pseudo_upf201_rhoatom(std::ifstream &ifs)
{
	delete[] rho_at;
	this->rho_at = new double[mesh];
	ModuleBase::GlobalFunc::ZEROS(rho_at, mesh);
	for (int ir = 0;ir < mesh;ir++)
	{
		ifs >> this->rho_at[ir];
	}
}