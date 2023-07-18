#include "sltk_atom.h"
#include <iostream>

class AdjacentSet;

//int FAtom::count1 = 0;
//int FAtom::count2 = 0;

/*** Constructors and destructor ***/
FAtom::FAtom()
{
	d_x = 0.0;	
	d_y = 0.0;	
	d_z = 0.0;	
	as = nullptr;
	type = 0;		
	natom = 0;
}

FAtom::~FAtom()
{
}

void FAtom::delete_vector(void)
{
	if (as) { as->delete_vector(); }
}
