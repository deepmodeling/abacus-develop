#include "tool_check.h"
#include "tool_quit.h"

void CHECK_NAME(ifstream &ifs,const string &name_in,bool quit)
{
    string name;
    ifs >> name;
    if ( name != name_in)
    {
		if(quit)
		{
			//ofs_warning << "\n name = " <<name;
			//ofs_warning << "\n should be = " << name_in;
			cout << "\n name = " <<name;
			cout << "\n should be = " << name_in;
        	WARNING_QUIT("CHECK_NAME","Some parameter name is wrong!");
		}
		else
		{
        	cout <<"\n Can not match : "<<name<<"(readin)  "<<name_in<<endl;
		}
    }
    return;
}

void CHECK_INT(ifstream &ifs,const int &v,bool quit)
{
	int v_in;
	ifs >> v_in;
	if( v!= v_in)
	{
		if(quit)
		{
			cout << "\n value = " << v_in;
			cout << "\n should be = " << v;
			WARNING_QUIT("CHECK_INT","Some parameter name is wrong!");
		}
		else
		{
			cout <<"\n Can not match well: "<<v_in<<"(readin)  "<<v<<endl;
		}
	}
	return;
}

void CHECK_DOUBLE(ifstream &ifs,const double &v,bool quit)
{
	const double tiny = 1.0e-5;
	double v_in;
	ifs >> v_in;
	if( (v - v_in) > tiny )
	{
		if(quit)
		{
			cout << " read in value = " << v_in << endl;
			cout << " the value should be = " << v << endl;
			WARNING_QUIT("CHECK_DOUBLE","the name of parameter wrong!");
		}
		else
		{
			cout <<" can not match well (1.0e-5): "<< v_in <<"(readin)  "<<v<<endl;
		}
	}
	return;
}

void CHECK_STRING(ifstream &ifs,const string &v,bool quit)
{
	string v_in;
	ifs >> v_in;
	if( v_in != v )
	{
		if(quit)
		{
			cout << " read in value = " << v_in << endl;
			cout << " the value should be = " << v << endl;
			WARNING_QUIT("CHECK_DOUBLE","the name of parameter wrong!");
		}
		else
		{
			cout <<" can not match well : "<<v_in<<"(readin)  "<<v<<endl;
		}
	}
	return;
}
