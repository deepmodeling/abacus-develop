#include "read_pp.h"

int Number[2];

// added by zhangwenshuai
int Pseudopot_upf::read_pseudo_upf201(ifstream &ifs)
{
    string dummy, word;
    int i, j, ir, ONCVPSP;

	while (ifs.good())
	{
		ifs >> dummy;
		// We start from PP_Header
		if(dummy=="<PP_HEADER")
		{
			// Read header
			READ_VALUE(ifs, word);   // generated
			READ_VALUE(ifs, word);   // author
			READ_VALUE(ifs, word);   // date
			READ_VALUE(ifs, word);   // comment

			READ_VALUE(ifs, word);   // element
			stringstream wdsstream(word);
			getline(wdsstream,this->psd,'"'); 
			getline(wdsstream,this->psd,'"'); 

			ifs >> word;   // pseudo_type
			if(word == "pseudo_type=\"")
			{
				ifs >> word;
				get_char(word);
				this->pp_type = word.substr(0,Number[0]);
			}
			else
			{
				get_char(word);
				this->pp_type = word.substr(Number[0]+1,(Number[1]-Number[0]-1));
			}

			if(pp_type!="NC") 
			{
				WARNING_QUIT("Pseudopot_upf::read_pseudo_header","unknown pseudo type");
			}

			READ_VALUE(ifs, word);   // relativistic
			READ_VALUE(ifs, word);   // is_ultrasoft
			if ( word.find("\"T\"") < word.length() ) // zws add 20160108
			{
				cout << "\n WARNING: ULTRASOFT PSEUDOPOTENTIAL IS NOT SUPPORTED !!! \n" << endl;
			}
			READ_VALUE(ifs, word);   // is_paw
			if ( word.find("\"T\"") < word.length() )
			{
				cout << "\n WARNING: PAW PSEUDOPOTENTIAL IS NOT SUPPORTED !!! \n" << endl;
			}

			READ_VALUE(ifs, word);   // is_coulomb
			ifs >> word;   // has_so
			string so;

			if(word == "has_so=\"")
			{
				ifs >> word;
				get_char(word);
				so = word.substr(0,Number[0]);
			}
			else
			{
				get_char(word);
				so = word.substr(Number[0]+1,(Number[1]-Number[0]-1));
			}

			if (so == "T")
			{
				this->has_so = true;
			}
			else
			{
				this->has_so = false;
			}

			READ_VALUE(ifs, word);   // has_wfc
			READ_VALUE(ifs, word);   // has_gipaw

			string nlc;
			//char p[13] = "paw_as_gipaw";
			ifs >> word;             // paw_as_gipaw?
			//cout << "word.substr(0,30) = " << word.substr(0,30) << "."<< endl;
			if( word.substr(0,13) == "paw_as_gipaw" )
			{
				ONCVPSP = 0;
				ifs >> word;     // core_correction
				if(word == "core_correction=\"")
				{
					ifs >> word;
					get_char(word);
					nlc = word.substr(0,Number[0]);
				}
				else
				{
					get_char(word);
					nlc = word.substr(Number[0]+1,(Number[1]-Number[0]-1));
				}

			}
			else
			{
				ONCVPSP = 1; // Generated using ONCVPSP code by D. R. Hamann, SG15 DOJO
				if(word == "core_correction=\"")
				{
					ifs >> word;
					get_char(word);
					nlc = word.substr(0,Number[0]);
				}
				else
				{
					get_char(word);
					nlc = word.substr(Number[0]+1,(Number[1]-Number[0]-1));
				}
			}

			//cout << "nlc = " << nlc << endl;

			if (nlc == "T")
			{
				this->nlcc = true;
			}
			else
			{
				this->nlcc = false;
			}

			READ_VALUE(ifs, word);   // functional
			//cout << "word = " << word << endl;
			//                        this->dft[0]="SLA";
			//                        this->dft[1]="PZ";
			//                        this->dft[2]="NOGX";
			//                        this->dft[3]="NOGC";

			string funcstr;  //{zws 01-06-16
			wdsstream.str("");
			wdsstream.clear();
			wdsstream << word;
			for ( int idft = 0; idft < 2; idft++)
			{
				getline(wdsstream,funcstr,'"');
			}
			wdsstream.str("");
			wdsstream.clear();
			wdsstream << funcstr;

			for( int idft = 0; idft < 4; idft++ )
			{
				getline(wdsstream,dft[idft],'-');
			}

			do 
			{
				getline(ifs, word);
				//cout << "word       = " << word << endl;
				word.erase(0,word.find_first_not_of(" ") );
				word.erase(word.find_last_not_of(" ")+1 );
				//word = trim(word);
				//cout << "trim(word) = " << word << endl;
				get_char(word);
				//cout << " Number = " << Number[0] << ", " << Number[1] << endl;
				//cout << word.substr(0,Number[0])  << "__" << word.substr(Number[0]+1, Number[1]-Number[0]-1) << endl;
				dummy = word.substr(0,Number[0]) ;
				//cout << " dummy = " << dummy << endl;
				if( dummy == "z_valence=" ) 
				{
					this->zp = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
				}
				else if ( dummy == "total_psenergy=" ) 
				{
					this->etotps = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
				}
				else if ( dummy == "rho_cutoff=" )
				{
				}
				else if ( dummy == "wfc_cutoff=" ) 
				{
				}
				else if ( dummy == "l_max=" ) 
				{
					this->lmax = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
				}
				else if ( dummy == "mesh_size=" ) 
				{
					this->mesh = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
				}
				else if ( dummy == "number_of_wfc=" ) 
				{
					this->nwfc = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
				}
				else if ( dummy == "number_of_proj=" ) 
				{
					this->nbeta = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
				}
				//break;

			}while( word.substr(word.length()-1, 1) !=">" ); 
			//cout << "word.substr(word.length()-1, 1)=" <<  word.substr(word.length()-1, 1)  << endl;
			//exit(0);


			//ifs >> word;   // zp
			////cout << "word = " << word << endl;
			//{
			//     if(word == "z_valence=\"")
			//     {
			//        ifs >> word;
			//        get_char(word);
			//        this->zp = atoi(word.substr(0,Number[0]).c_str());
			//     }
			//     else
			//     {
			//        get_char(word);
			//        this->zp = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
			//     }
			//     //cout << "zp = " << this->zp << endl;
			//}

			//ifs >> word;   // total_psenergy
			//{
			//     if(word == "total_psenergy=\"")
			//     {
			//        ifs >> word;
			//        get_char(word);
			//        this->etotps = atof(word.substr(0,Number[0]).c_str());
			//     }
			//     else
			//     {
			//        get_char(word);
			//        this->etotps = atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
			//     }
			//     //cout << "etotps = " << this->etotps << endl;
			//}
			////cout << " word (total_psenergy) = " << word << endl;


			//if(ONCVPSP == 0)    //zws modify 20160108
			//{
			//	READ_VALUE(ifs, word);   // wfc_cutoff
			//	//cout << "word = " << word << endl;
			//}
			//READ_VALUE(ifs, word); // rho_cutoff
			//cout << "word (cutoff) = " << word << endl;


			//ifs >> word;             // lmax
			////cout << "word (lmax) = " << word << endl;
			//{
			//        if(word == "l_max=\"")
			//        {
			//             ifs >> word;
			//             get_char(word);
			//             this->lmax = atoi(word.substr(0,Number[0]).c_str());
			//        }
			//        else
			//        {
			//             get_char(word);
			//             this->lmax = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
			//        }

			//}

			////cout << "lmax = " << this->lmax << endl;

			//if(ONCVPSP == 0)
			//{
			//   READ_VALUE(ifs, word);   // l_max_rho
			//}

			//READ_VALUE(ifs, word);   // l_local

			//ifs >> word;   // mesh_size
			////cout << "word (mesh) = " << word << endl;
			//{
			//     if(word == "mesh_size=\"")
			//     {
			//             ifs >> word;
			//             get_char(word);
			//             this->mesh = atoi(word.substr(0,Number[0]).c_str());
			//     }
			//     else
			//     {
			//             get_char(word);
			//             this->mesh = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
			//     }
			//     //cout << "mesh = " << this->mesh << endl;
			//}



			//ifs >> word;  // number_of_wfc
			////cout << "word = " << word << endl;
			//{
			//     if(word == "number_of_wfc=\"")
			//     {
			//             ifs >> word;
			//             get_char(word);
			//             this->nwfc = atoi(word.substr(0,Number[0]).c_str());

			//     }
			//     else
			//     {
			//             get_char(word);
			//             this->nwfc = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
			//     }
			//     //cout << "nwfc = " << this->nwfc << endl;
			//}
			//     
			//ifs >> word;   // number_of_proj
			////cout << "word = " << word << endl;
			//{
			//     if(word == "number_of_proj=\"")
			//     {
			//             ifs >> word;
			//             get_char(word);
			//             this->nbeta = atoi(word.substr(0,Number[0]).c_str());

			//     }
			//     else
			//     {
			//             get_char(word);
			//             this->nbeta = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
			//     }
			//     //cout << "nbeta = " << this->nbeta << endl;
			//}


			// READ Mesh
			if(ONCVPSP == 0)
			{
				SCAN_BEGIN(ifs, "<PP_MESH");
			}
			else
			{
				SCAN_BEGIN(ifs, "<PP_MESH>");
			}

			assert(mesh>0);
			if(ONCVPSP == 0)
			{
				ifs >> word;             // dx
				ifs >> word;             // mesh
				ifs >> word;             // xmin
				ifs >> word;             // rmax
				ifs >> word;             // zmesh
			}

			SCAN_BEGIN(ifs, "<PP_R"); 
			READ_VALUE(ifs, word);    // type  size  columns

			//                        double  rmesh0 = 1;    //{zws add160108 delete160328
			//                        int 	nmeshdel = 0;
			//                        ifs >> rmesh0;
			//                        if ( abs(rmesh0) < 1.0e-15 )
			//                        {
			//                            mesh       -= 1;
			//                            nmeshdel   += 1;
			//                        }
			//                        cout << " mesh =" << mesh << endl;
			//                    	if (mesh%2 == 0)
			//                    	{
			//                    	    mesh     -= 1;
			//                    	    nmeshdel += 1;
			//                    	}    //}zws add 20160108
			//                    	cout << " nmeshdel =" << nmeshdel << endl;


			delete[] r;
			delete[] rab;
			this->r = new double[mesh];
			this->rab = new double[mesh];
			ZEROS(r,mesh);
			ZEROS(rab,mesh);


			//                        if (nmeshdel == 0)    //{zws add160108 delete160328
			//                        {
			//                            this->r[0] = rmesh0;
			//                            for (ir = 1;ir < mesh;ir++)
			//                            {
			//                                ifs >> this->r[ir];
			//                            }
			//                        }
			//                        else
			//                        {
			//                            for ( int idel=0; idel < nmeshdel-1; idel++)
			//                        	{
			//                            	cout << "skip " << nmeshdel << "grid point(s) in PP mesh" << endl;
			//                        	    double	tmpdel;
			//                        	    ifs >> tmpdel;
			//                        	}
			//                            for (ir = 0;ir < mesh;ir++)
			//                            {
			//                                 ifs >> this->r[ir];
			//                            }
			//                        }    //}zws 20160108
			for (ir = 0;ir < mesh;ir++)
			{
				ifs >> this->r[ir];
			}
			SCAN_END(ifs, "</PP_R>");

			SCAN_BEGIN(ifs, "<PP_RAB");
			READ_VALUE(ifs, word);    // type size columns

			//                        for ( int idel=0; idel < nmeshdel; idel++)    //{zws add 20160108
			//                    	{
			//                    	    double	tmpdel;
			//                    	    ifs >> tmpdel;
			//                    	}    //}zws add 20160108
			for (ir = 0;ir < mesh;ir++)
			{
				ifs >> this->rab[ir];
			}
			SCAN_END(ifs, "</PP_RAB>");
			SCAN_END(ifs, "</PP_MESH>");

			// READ NLCC
			if (this->nlcc)
			{
				SCAN_BEGIN(ifs, "<PP_NLCC");
				READ_VALUE(ifs, word);    // type size columns

				assert(mesh>0);
				delete[] rho_atc;
				this->rho_atc = new double[mesh];
				ZEROS(rho_atc, mesh);

				for (ir = 0;ir < mesh;ir++)
				{
					ifs >> this->rho_atc[ir];
				}
				SCAN_END(ifs, "</PP_NLCC>");

			}

			// READ VLOCAL
			SCAN_BEGIN(ifs, "<PP_LOCAL");
			READ_VALUE(ifs, word);    // type size columns

			assert(mesh>0);
			delete[] vloc;
			this->vloc = new double[mesh];
			ZEROS(vloc, mesh);

			for (ir = 0;ir < mesh;ir++)
			{
				ifs >> this->vloc[ir];
			}

			SCAN_END(ifs, "</PP_LOCAL>");

			// READ NONLOCAL
			SCAN_BEGIN(ifs, "<PP_NONLOCAL>");

			delete[] kkbeta;
			delete[] lll;
			this->kkbeta = new int[nbeta];
			this->lll = new int[nbeta];
			this->beta.create(nbeta , mesh);
			this->dion.create(nbeta , nbeta);

			for(i=0;i<nbeta;i++)
			{
				ifs >> word;  //number
				ifs >> word;  //type
				if(word == "type=\"")
				{
					ifs >> word;
				}
				ifs >> word;  //size
				if(word == "size=\"")
				{
					ifs >> word;
				}

				ifs >> word;  //columns
				if(word == "columns=\"")
				{
					ifs >> word;
				}

				ifs >> word;  //index
				{
					if(word == "index=\"")
					{
						ifs >> word;
						get_char(word);
						//idum = atoi(word.substr(0,Number[0]).c_str());
					}
					else
					{
						get_char(word);
						//idum = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
					}
					//cout << "idum = " << idum << endl;
				}

				if(ONCVPSP == 0)
				{
					ifs >> word;  //label
				}

				ifs >> word;  //angular_momentum
				if(word == "angular_momentum=\"")
				{
					ifs >> word;
					get_char(word);
					this->lll[i] = atoi(word.substr(0,Number[0]).c_str());

				}
				else
				{
					get_char(word);
					this->lll[i] = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
				}

				ifs >> word;  //cutoff_radius_index
				if(word == "cutoff_radius_index=\"")
				{
					ifs >> word;
					get_char(word);
					this->kkbeta[i] = atoi(word.substr(0,Number[0]).c_str());

				}
				else
				{
					get_char(word);
					this->kkbeta[i] = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
				}
				//cout << "kkbeta[i] = " << this->kkbeta[i] << endl;

				if(ONCVPSP ==0) 
				{
					ifs >> word;  //cutoff_radius
					ifs >> word;  //ultrasoft_cutoff_radius
				}
				else
				{
					READ_VALUE(ifs, word); // cutoff_radius
				}

				for (ir=0;ir<mesh;ir++)
				{
					ifs >> this->beta(i, ir);

				}

				ifs >> word;  //number

			}

			// READ DIJ
			SCAN_BEGIN(ifs, "<PP_DIJ");
			READ_VALUE(ifs, word);  // type size columns

			this->nd = nbeta * nbeta;
			for(i=0;i<nbeta;i++)
			{
				for(j=0;j<nbeta;j++)
				{
					ifs >> dion(i,j);
					if ( i != j  && dion(i,j) != 0.0 )
					{
						cout << " error: for i != j, Dij of Pseudopotential must be 0.0 " << endl;
						exit(1);
					}
				}
			}
			SCAN_END(ifs, "</PP_DIJ>");

			SCAN_END(ifs, "</PP_NONLOCAL>");

			// READ PSWFC
			SCAN_BEGIN(ifs, "<PP_PSWFC>");

			delete[] els;
			delete[] lchi;
			delete[] oc;
			this->els = new string[nwfc];
			this->lchi = new int[nwfc];
			this->oc = new double[nwfc];
			ZEROS(lchi, nwfc); // angular momentum of each orbital
			ZEROS(oc, nwfc);//occupation of each orbital

			this->chi.create(this->nwfc, this->mesh);
			for (i=0;i<nwfc;i++)
			{
				ifs >> word;  // number
				ifs >> word;  // type
				ifs >> word;  // size 
				{
					if(word == "size=\"")
					{
						ifs >> word;
						word = "\"" + word ;
					}
				}
				ifs >> word;  // columns
				ifs >> word;  // index
				ifs >> word;  // occupation
				{
					if(word == "occupation=\"")
					{
						ifs >> word;
						word = "\"" + word ;
					}
					get_char(word);
					oc[i] = atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
				}

				ifs >> word;  // pseudo_energy
				if(word == "pseudo_energy=\"")
				{
					ifs >> word;
					word = "\"" + word ;
				}
				get_char(word);

				ifs >> word;  // label
				if(word == "label=\"")
				{
					ifs >> word;
					word = "\"" + word ;
				}
				get_char(word);
				els[i] = word.substr(Number[0]+1,(Number[1]-Number[0]-1));

				ifs >> word;  // l
				if(word == "l=\"")
				{
					ifs >> word;
					word = "\"" + word ;
				}
				get_char(word);
				lchi[i] = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
				//cout << " lchi[i] = " << lchi[i] << endl;

				ifs >> word; // >
				if ( word !=  ">" )
				{
					cout << " error: bad end while reading CHI" << i <<  " of PSWFC" << endl;
					exit(1);
				}

				for (ir = 0;ir < mesh;ir++)
				{
					ifs >> this->chi(i, ir);
				}
				ifs >> word;  // number
			}

			SCAN_END(ifs, "</PP_PSWFC>");

			// READ RHOATOM
			SCAN_BEGIN(ifs, "<PP_RHOATOM");
			READ_VALUE(ifs, word); // type size columns

			delete[] rho_at;
			this->rho_at = new double[mesh];
			ZEROS(rho_at, mesh);

			for (ir = 0;ir < mesh;ir++)
			{
				ifs >> this->rho_at[ir];
			}
			SCAN_END(ifs, "</PP_RHOATOM>");

			SCAN_BEGIN(ifs, "<PP_SPIN_ORB>");
			//added by zhengdy-soc
			delete[] this->jchi;
			delete[] this->jjj;
			delete[] this->nn;
			this->jchi = new double [nwfc];
			this->jjj = new double [nbeta];
			this->nn = new int [nwfc];
			ZEROS(jchi,nwfc);
			ZEROS(jjj,nbeta);
			ZEROS(nn,nwfc);

			for(int round=0;round<2;round++)
			{
				ifs>>word;
				if(word=="<PP_RELBETA.1")
				{
					for(int nb = 0;nb<nbeta;nb++)
					{
						if(nb!=0) ifs>>word; //RELBETA
						ifs>>word;           //index
						ifs>>word;           //lll
						get_char(word);
						this->lll[nb] = atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
						ifs>>word;           //jjj
						get_char(word);
						this->jjj[nb] = atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
					}
				}
				else if(word=="<PP_RELWFC.1")
				{
					if(round==0)
					{
						for(int nw = 0;nw<nwfc;nw++)
						{
							if(nw!=0) ifs>>word;     //RELWFC
							ifs>>word;               //index
							ifs>>word;               //els
							get_char(word);
							this->els[nw]= word.substr(Number[0]+1,(Number[1]-Number[0]-1));
							ifs>>word;               //nn
							get_char(word);
							this->nn[nw]= atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
							ifs>>word;               //lchi
							get_char(word);
							this->lchi[nw]= atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
							ifs>>word;               //jchi
							get_char(word);
							this->jchi[nw]= atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
							ifs>>word;               //oc
							get_char(word);
							this->oc[nw]= atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());

						}
					}
					else
					{
						for(int nw = 0;nw<nwfc;nw++)
						{
							if(nw!=0) ifs>>word;//RELWFC
							ifs>>word;          //index
							ifs>>word;          //lchi
							get_char(word);
							this->lchi[nw]= atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
							ifs>>word;          //jchi
							get_char(word);
							this->jchi[nw]= atof(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());
							ifs>>word;
							get_char(word);
							this->nn[nw]= atoi(word.substr(Number[0]+1,(Number[1]-Number[0]-1)).c_str());

						}
					}
				}
				else if(round==0)
				{
					this->has_so = 0;
					//	cout<<"ignore SPIN_ORB part!"<<endl;
					break;
				}
			}
			SCAN_END(ifs, "</PP_SPIN_ORB>");

			if (mesh%2 == 0)
			{
				mesh -= 1;
			}

			SCAN_END(ifs, "</UPF>");
			break;
		}
	}
	return 0;
}


void Pseudopot_upf::get_char( string ss)
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
}
