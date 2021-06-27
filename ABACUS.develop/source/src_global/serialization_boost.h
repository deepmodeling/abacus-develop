#ifndef SERIALIZATION_BOOST_H
#define SERIALIZATION_BOOST_H

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/string.hpp>

#include "src_global/vector3.h"
#include "src_lcao/abfs-vector3_order.h"
#include "module_base/matrix.h"
#include "src_global/matrix_wrapper.h"

namespace boost
{
	namespace serialization
	{
		// Vector3<T>								Peize Lin add 2018-07-01
		template<typename Archive, typename T> 
		inline void serialize(Archive & ar, Vector3<T> & v, const unsigned int version)
		{
			ar & v.x;
			ar & v.y;
			ar & v.z;
		}
		
		// Abfs::Vector3_Order<T>					Peize Lin add 2018-07-01
		template<typename Archive, typename T> 
		inline void serialize(Archive & ar, Abfs::Vector3_Order<T> & v, const unsigned int version)
		{
			ar & v.x;
			ar & v.y;
			ar & v.z;
		}
		
		// matrix									Peize Lin add 2018-07-01
		template<typename Archive> 
		inline void save( Archive & ar, const matrix & m, const unsigned int /*file_version*/ )
		{
			const collection_size_type nr(m.nr), nc(m.nc);
			ar << BOOST_SERIALIZATION_NVP(nr) << BOOST_SERIALIZATION_NVP(nc);
			if( nr && nc )
				ar << serialization::make_array( m.c, nr*nc );
		}
		template<class Archive>
		inline void load( Archive & ar, matrix &m,  const unsigned int /*file_version*/ )
		{
			collection_size_type nr, nc;
			ar >> BOOST_SERIALIZATION_NVP(nr) >> BOOST_SERIALIZATION_NVP(nc);
			m.create(nr,nc,false);
			if( nr && nc )
				ar >> serialization::make_array( m.c, nr*nc );
		}
		template<class Archive>
		inline void serialize( Archive & ar, matrix & m, const unsigned int file_version )
		{
			boost::serialization::split_free(ar, m, file_version);
		}
		
		// MatrixWrapper									Peize Lin add 2018-07-31
		template<typename Archive> 
		inline void save( Archive & ar, const Matrix_Wrapper & m, const unsigned int /*file_version*/ )
		{
			const collection_size_type nr(m.nr), nc(m.nc);
			ar << BOOST_SERIALIZATION_NVP(nr) << BOOST_SERIALIZATION_NVP(nc);
			if( nr && nc )
				ar << serialization::make_array( m.c, nr*nc );
		}
		template<class Archive>
		inline void load( Archive & ar, Matrix_Wrapper &m,  const unsigned int /*file_version*/ )
		{
			collection_size_type nr, nc;
			ar >> BOOST_SERIALIZATION_NVP(nr) >> BOOST_SERIALIZATION_NVP(nc);
			m.create(nr,nc,false);
			if( nr && nc )
				ar >> serialization::make_array( m.c, nr*nc );
		}
		template<class Archive>
		inline void serialize( Archive & ar, Matrix_Wrapper & m, const unsigned int file_version )
		{
			boost::serialization::split_free(ar, m, file_version);
		}
	}
}


#include <boost/serialization/collection_traits.hpp>

BOOST_SERIALIZATION_COLLECTION_TRAITS(Vector3)
BOOST_SERIALIZATION_COLLECTION_TRAITS(Abfs::Vector3_Order)

#endif