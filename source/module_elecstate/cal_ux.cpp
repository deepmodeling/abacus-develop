#include "cal_ux.h"
#include "module_parameter/parameter.h"

namespace elecstate {

void cal_ux(UnitCell& ucell) {

    if (PARAM.inp.nspin != 4)
    {
        return;
    }

    const double absolute_mag_thr = 1.0e-6;

    double amag = 0.0;
    double uxmod = 0.0;

    int starting_it = 0;
    int starting_ia = 0;
    bool is_paraller = false;

    // do not sign feature in teh general case
    ucell.magnet.lsign_ = false;
    ModuleBase::GlobalFunc::ZEROS(ucell.magnet.ux_, 3);

	for (int it = 0; it < ucell.ntype; it++) 
	{
		for (int ia = 0; ia < ucell.atoms[it].na; ia++) 
		{
            // m_loc_: local magnetization vector for each atom
            amag = pow(ucell.atoms[it].m_loc_[ia].x, 2)
                   + pow(ucell.atoms[it].m_loc_[ia].y, 2)
                   + pow(ucell.atoms[it].m_loc_[ia].z, 2);

            // find the first atom (it,ia) whose magnetism is not zero
            // compute ux
			if (amag > absolute_mag_thr) 
			{
				ucell.magnet.ux_[0] = ucell.atoms[it].m_loc_[ia].x;
				ucell.magnet.ux_[1] = ucell.atoms[it].m_loc_[ia].y;
				ucell.magnet.ux_[2] = ucell.atoms[it].m_loc_[ia].z;

				starting_it = it;
				starting_ia = ia;

				ucell.magnet.lsign_ = true;
				break;
			}
		}

        // if any atom has magnetism, then break the for iteration
		if (ucell.magnet.lsign_) 
		{
			break;
		}
	}

    // whether the initial magnetizations is parallel
	for (int it = starting_it; it < ucell.ntype; it++) 
	{
		for (int ia = 0; ia < ucell.atoms[it].na; ia++) 
		{
			if (it > starting_it || ia > starting_ia) 
			{
				ucell.magnet.lsign_
					= ucell.magnet.lsign_
					&& judge_parallel(ucell.magnet.ux_, ucell.atoms[it].m_loc_[ia]);
			}
		}
	}

    // if all of the atoms have the same parallel magnetism direction,
    // then set the direction to a unit vector
	if (ucell.magnet.lsign_) 
	{
		uxmod = pow(ucell.magnet.ux_[0], 2) 
			+ pow(ucell.magnet.ux_[1], 2)
			+ pow(ucell.magnet.ux_[2], 2);

		if (uxmod < absolute_mag_thr) 
		{
			ModuleBase::WARNING_QUIT("elecstate::cal_ux", "wrong uxmod");
		}

        // reset the magnetism for each direction
		for (int i = 0; i < 3; i++) 
		{
			ucell.magnet.ux_[i] *= 1 / sqrt(uxmod);
		}
	}
	return;
}

bool judge_parallel(double a[3], ModuleBase::Vector3<double> b) {
    bool jp = false;
    double cross;
    cross = pow((a[1] * b.z - a[2] * b.y), 2)
            + pow((a[2] * b.x - a[0] * b.z), 2)
            + pow((a[0] * b.y - a[1] * b.x), 2);
    jp = (fabs(cross) < 1e-6);
    return jp;
}
}
