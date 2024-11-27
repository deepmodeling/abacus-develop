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
	type = 0;		
	natom = 0;
}
