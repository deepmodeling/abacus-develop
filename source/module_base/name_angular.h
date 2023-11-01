#ifndef NAME_ANGULAR_H
#define NAME_ANGULAR_H

namespace ModuleBase
{
	const std::string Name_Angular[5][11] =
	{
    	{"s"},
    	{"pz", "px", "py"},
    	{"dz^2", "dxz", "dyz", "dx^2-y^2", "dxy"},
    	{"fz^3", "fxz^2", "fyz^2", "fzx^2-zy^2", "fxyz", "fx^3-3*xy^2", "f3yx^2-y^3"},
    	{"g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9"}
	};          // name of atomic orbital    jiyy add 2022-05-10
}

#endif
