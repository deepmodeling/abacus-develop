#include "common.h"
void PRINTCM(const string &s, const ComplexMatrix &m)
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    cout << "\n" << s << " " <<  b1 << " " << b2 ;

    if (b1*b2 == 0) return;

	cout << "\n norm : ";
    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) cout << "\n ";
			double norm = ( conj( m(i, j) ) * m(i, j) ).real();
			if( abs(norm) > 1.0e-9) cout<< setw(12) << norm;
			else cout<<setw(12)<<"0";
        }
    }

	/*
	cout << "\n real part : ";
	for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) cout << "\n ";
			double n = ( m(i, j) ).real();
			if( abs(n) > 1.0e-9) cout<< setw(12) << n;
			else cout<<setw(12)<<"0";
        }
    }


	cout << "\n imag part : ";
	for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) cout << "\n ";
			double n = ( m(i, j) ).imag();
			if( abs(n) > 1.0e-9) cout<< setw(12) << n;
			else cout<<setw(12)<<"0";
        }
    }
	*/

    return;
}


bool SCAN_BEGIN(ifstream &ifs, const string &TargetName, const bool restart, const bool quit)
{
    string SearchName;
    bool find = false;
    if (restart)
    {
        ifs.clear();
        ifs.seekg(0);
    }
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> SearchName;
        if( SearchName == TargetName)
        {
            find = true;
//			cout << " Find the block:" << TargetName << endl;
            break;
        }
    }
    if (!find && quit)
    {
        cout << "\n Can't find : " << TargetName;
        exit(0);
    }
    return find;
}

void SCAN_END(ifstream &ifs, const string &TargetName)
{
    string SearchName;
    ifs >> SearchName;
    if( SearchName != TargetName)
    {
        cout<<"\n"<<"In SCAN_END, can't find: "<<TargetName<<" block."<<endl;
    }
    return;
}

void WARNING_QUIT(const string &file,const string &description)
{
	cout << "\n !!!!!!!!!!!!!!!  WARNING  !!!!!!!!!!!!!!!!!!!!" ;
	cout << "\n "<<file<<"  Warning : "<<description<<endl;
	QUIT();
}

void BLOCK_HERE( const string &description)
{
	//  return;
    cout << "\n********************************************";
    cout << "\n Here is a Block, input anything you like..";
    cout << "\n " << description;
    cout << "\n********************************************" << endl;
    string ok;
    cin >> ok;

}

void QUIT(void)
{
    timer::finish();
#ifdef __MPI
    MPI_Finalize();
#endif
    exit(0);
}

void TITLE(const string &class_name,const string &function_name)
{
//  ofs_running<<"\n\n ==> "<<class_name<<"::"<<function_name<<endl;
	return; //mohan add 2012-06-21
	cout<<"\n -------------------------------------" << endl;
    cout<<" ==> "<<class_name<<"::"<<function_name << endl;
	cout<<" -------------------------------------" << endl;
    return;
}

void TITLE(ofstream &ofs, const string &class_name,const string &function_name)
{
	return;
	cout<<"\n -------------------------------------" << endl;
    cout<<" ==> "<<class_name<<"::"<<function_name << endl;
	cout<<" -------------------------------------" << endl;

	ofs <<" -------------------------------------" << endl;
    ofs <<" ==> "<<class_name<<"::"<<function_name << endl;
	ofs <<" -------------------------------------" << endl;
    return;
}


