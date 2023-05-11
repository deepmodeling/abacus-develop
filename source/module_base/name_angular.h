#ifndef NAME_ANGULAR_H
#define NAME_ANGULAR_H

namespace ModuleBase
{
	const std::string Name_Angular[5][11] =
	{
    	{"s"},
    	{"px", "py", "pz"},
    	{"d3z^2-r^2", "dxy", "dxz", "dx^2-y^2", "dyz"},
    	{"f5z^2-3r^2", "f5xz^2-xr^2", "f5yz^2-yr^2", "fzx^2-zy^2", "fxyz", "fx^3-3*xy^2", "f3yx^2-y^3"},
    	{"g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9"}
	};          // name of atomic orbital    jiyy add 2022-05-10
}

#endif
