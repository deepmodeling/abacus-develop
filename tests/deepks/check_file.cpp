#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

int main(int argc, char **argv)
{
    string ref_file = argv[1];
    string cal_file = argv[2];

    string temp = argv[3];
    istringstream string_temp(temp);
    double threshold;
    string_temp >> threshold;

    ifstream ref(ref_file.c_str(), ios::in);
    ifstream cal(cal_file.c_str(), ios::in);

    if (!ref) 
	{
        cout << "\e[1;31m [  FAILED  ]  \e[0m";
		cout << "Can't find " << ref_file << " !" << endl;
		return 0;
	}

    if (!cal) 
	{
        cout << "\e[1;31m [  FAILED  ]  \e[0m";
		cout << "Can't find " << cal_file << " !" << endl;
		return 0;
	}

    int num_ref=0, num_cal=0;
    string word_ref, word_cal;

    while (ref.good())
    {
        ref >> word_ref;
        num_ref++;
    }

    while (cal.good())
    {
        cal >> word_cal;
        num_cal++;
    }

    if(num_ref != num_cal)
    {
        cout << "\e[1;31m [  FAILED  ]  \e[0m";
        cout << "The num of data is different !" << endl;
        return 0;
    }

    ref.clear();
    ref.seekg(0);
    ref.rdstate();

    cal.clear();
    cal.seekg(0);
    cal.rdstate();

    double data_ref, data_cal;

    for(int index=0; index<num_ref; ++index)
    {
        ref >> word_ref;
        cal >> word_cal;

        if(word_cal != word_ref)
        {
            istringstream string_cal(word_cal);
            istringstream string_ref(word_ref);
            string_cal >> data_cal;
            string_ref >> data_ref;
            if(data_ref==0 || data_cal==0)
            {
                cout << "\e[1;31m [  FAILED  ]  \e[0m";
                cout << "Data No." << index+1 << "  cal=" << word_cal 
                     << "  ref=" << word_ref << endl;
            }
            else if(fabs(data_ref-data_cal) > threshold)
            {
                cout << "\e[1;31m [  FAILED  ]  \e[0m";
                cout << "Data No." << index+1 << "  cal=" << data_cal 
                     << "  ref=" << data_ref << "  deviation=" << fabs(data_ref-data_cal) << endl;
                return 0;
            }
        }
    }
    
    return 1;
}