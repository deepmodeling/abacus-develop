#include "sltk_atom.h"
#include <iostream>

using namespace std;

class AdjacentSet;

//int FAtom::count1 = 0;
//int FAtom::count2 = 0;

/*** Constructors and destructor ***/
FAtom::FAtom()
{
	d_x = 0.0;	
	d_y = 0.0;	
	d_z = 0.0;	
	this->as = NULL;
	type = 0;		
	natom = 0;
	allocate = false;
}

FAtom::~FAtom(){ 
//	if(this->as!=NULL)
//	{
//		delete this->as; 
//		cout << "\n destroy Adjacent!";
//		count2++;
//		cout << "\n count2 = " << count2;
//	}
	if(allocate) 
	{
		delete this->as;
	}
}

void FAtom::delete_vector(void)
{
	as->delete_vector();
}
