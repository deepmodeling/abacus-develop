#ifndef SERIALIZATION_BOOST_H
#define SERIALIZATION_BOOST_H

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/std::vector.hpp>
#include <boost/serialization/std::map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/std::string.hpp>

#include "../module_base/vector3.h"
#include "../module_base/abfs-vector3_order.h"
#include "../module_base/matrix.h"
#include "../module_base/matrix_wrapper.h"

namespace boost
{
	namespace serialization
	{
		// ModuleBase::Vector3<T>								Peize Lin add 2018-07-01
		template<typename Archive, typename T>
		inline void serialize(Archive & ar, ModuleBase::Vector3<T> & v, const unsigned int version)
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
		inline void save( Archive & ar, const ModuleBase::matrix & m, const unsigned int /*file_version*/ )
		{
			const collection_size_type nr(m.nr), nc(m.nc);
			ar << BOOST_SERIALIZATION_NVP(nr) << BOOST_SERIALIZATION_NVP(nc);
			if( nr && nc )
				ar << serialization::make_array( m.c, nr*nc );
		}
		template<class Archive>
		inline void load( Archive & ar, ModuleBase::matrix &m,  const unsigned int /*file_version*/ )
		{
			collection_size_type nr, nc;
			ar >> BOOST_SERIALIZATION_NVP(nr) >> BOOST_SERIALIZATION_NVP(nc);
			m.create(nr,nc,false);
			if( nr && nc )
				ar >> serialization::make_array( m.c, nr*nc );
		}
		template<class Archive>
		inline void serialize( Archive & ar, ModuleBase::matrix & m, const unsigned int file_version )
		{
			boost::serialization::split_free(ar, m, file_version);
		}

		// MatrixWrapper									Peize Lin add 2018-07-31
		template<typename Archive>
		inline void save( Archive & ar, const ModuleBase::Matrix_Wrapper & m, const unsigned int /*file_version*/ )
		{
			const collection_size_type nr(m.nr), nc(m.nc);
			ar << BOOST_SERIALIZATION_NVP(nr) << BOOST_SERIALIZATION_NVP(nc);
			if( nr && nc )
				ar << serialization::make_array( m.c, nr*nc );
		}
		template<class Archive>
		inline void load( Archive & ar, ModuleBase::Matrix_Wrapper &m,  const unsigned int /*file_version*/ )
		{
			collection_size_type nr, nc;
			ar >> BOOST_SERIALIZATION_NVP(nr) >> BOOST_SERIALIZATION_NVP(nc);
			m.create(nr,nc,false);
			if( nr && nc )
				ar >> serialization::make_array( m.c, nr*nc );
		}
		template<class Archive>
		inline void serialize( Archive & ar, ModuleBase::Matrix_Wrapper & m, const unsigned int file_version )
		{
			boost::serialization::split_free(ar, m, file_version);
		}
	}
}


#include <boost/serialization/collection_traits.hpp>

BOOST_SERIALIZATION_COLLECTION_TRAITS(ModuleBase::Vector3)
BOOST_SERIALIZATION_COLLECTION_TRAITS(Abfs::Vector3_Order)

#endif
