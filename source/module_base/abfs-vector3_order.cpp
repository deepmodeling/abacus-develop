#include "abfs-vector3_order.h"

template<>
bool operator< ( const Abfs::Vector3_Order<double> &v1, const Abfs::Vector3_Order<double> &v2 )
{
	// threshold is used when v1 and v2 are approximately equal with double numerical error,
	// but a strange bug occurs at LRI_CV::dVws.
	// Peize Lin change at 2023.04.29
	if      ( v1.x < v2.x ) return true;
	else if ( v1.x > v2.x ) return false;
	if      ( v1.y < v2.y ) return true;
	else if ( v1.y > v2.y ) return false;
	if      ( v1.z < v2.z ) return true;
	else if ( v1.z > v2.z ) return false;
	return false;
	//constexpr double threshold = 1E-8;
	//if      ( v1.x < v2.x - threshold ) return true;
	//else if ( v1.x > v2.x + threshold ) return false;
	//if      ( v1.y < v2.y - threshold ) return true;
	//else if ( v1.y > v2.y + threshold ) return false;
	//if      ( v1.z < v2.z - threshold ) return true;
	//else if ( v1.z > v2.z + threshold ) return false;
	//return false;
}

template<>
bool operator< ( const Abfs::Vector3_Order<int> &v1, const Abfs::Vector3_Order<int> &v2 )
{
	if      ( v1.x < v2.x ) return true;
	else if ( v1.x > v2.x ) return false;
	if      ( v1.y < v2.y ) return true;
	else if ( v1.y > v2.y ) return false;
	if      ( v1.z < v2.z ) return true;
	else if ( v1.z > v2.z ) return false;
	return false;
}