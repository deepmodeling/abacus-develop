#ifdef __MPI //liyuanbo 2022/2/23
#ifndef EXX_ABFS_INVERSE_MATRIX_DOUBLE_H
#define EXX_ABFS_INVERSE_MATRIX_DOUBLE_H

#include "exx_abfs.h"

#include "../module_base/matrix.h"

#include<memory>

class Exx_Abfs::Inverse_Matrix_Double
{
public:	
	enum class Method{dpotrf,dsyev};

	ModuleBase::matrix A;
	void init( const int &dim_in);
	
	void input( const ModuleBase::matrix &m );
	void input( const std::shared_ptr<ModuleBase::matrix> &pm ) { input(*pm); }
	void input( 
		const ModuleBase::matrix &m_00,
		const ModuleBase::matrix &m_01,
		const ModuleBase::matrix &m_10,
		const ModuleBase::matrix &m_11);
	void input( 
		const std::shared_ptr<ModuleBase::matrix> &pm_00,
		const std::shared_ptr<ModuleBase::matrix> &pm_01,
		const std::shared_ptr<ModuleBase::matrix> &pm_10,
		const std::shared_ptr<ModuleBase::matrix> &pm_11) { input( *pm_00, *pm_01, *pm_10, *pm_11 ); }
		
	void cal_inverse( const Method &method );

	void output( ModuleBase::matrix &m ) const;
	void output( const std::shared_ptr<ModuleBase::matrix> &pm ) const { output(*pm); }
	void output( 
		ModuleBase::matrix &m_00,
		ModuleBase::matrix &m_01,
		ModuleBase::matrix &m_10,
		ModuleBase::matrix &m_11) const;
	void output( 
		const std::shared_ptr<ModuleBase::matrix> &pm_00,
		const std::shared_ptr<ModuleBase::matrix> &pm_01,
		const std::shared_ptr<ModuleBase::matrix> &pm_10,
		const std::shared_ptr<ModuleBase::matrix> &pm_11) const { output( *pm_00, *pm_01, *pm_10, *pm_11 ); }

private:

	void using_dpotrf();
	void using_dsyev( const double &threshold_condition_number = 0 );
	
	void copy_down_triangle();
	
	int dim;
	int info;
};

#endif	// EXX_ABFS_INVERSE_MATRIX_DOUBLE_H
#endif