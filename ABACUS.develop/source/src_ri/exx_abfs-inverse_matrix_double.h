#ifndef EXX_ABFS_INVERSE_MATRIX_DOUBLE_H
#define EXX_ABFS_INVERSE_MATRIX_DOUBLE_H

#include "exx_abfs.h"

#include "module_base/matrix.h"

#include<memory>

class Exx_Abfs::Inverse_Matrix_Double
{
public:	
	enum class Method{dpotrf,dsyev};

	matrix A;
	void init( const int &dim_in);
	
	void input( const matrix &m );
	void input( const shared_ptr<matrix> &pm ) { input(*pm); }
	void input( 
		const matrix &m_00,
		const matrix &m_01,
		const matrix &m_10,
		const matrix &m_11);
	void input( 
		const shared_ptr<matrix> &pm_00,
		const shared_ptr<matrix> &pm_01,
		const shared_ptr<matrix> &pm_10,
		const shared_ptr<matrix> &pm_11) { input( *pm_00, *pm_01, *pm_10, *pm_11 ); }
		
	void cal_inverse( const Method &method );

	void output( matrix &m ) const;
	void output( const shared_ptr<matrix> &pm ) const { output(*pm); }
	void output( 
		matrix &m_00,
		matrix &m_01,
		matrix &m_10,
		matrix &m_11) const;
	void output( 
		const shared_ptr<matrix> &pm_00,
		const shared_ptr<matrix> &pm_01,
		const shared_ptr<matrix> &pm_10,
		const shared_ptr<matrix> &pm_11) const { output( *pm_00, *pm_01, *pm_10, *pm_11 ); }

private:

	void using_dpotrf();
	void using_dsyev( const double &threshold_condition_number = 0 );
	
	void copy_down_triangle();
	
	int dim;
	int info;
};

#endif	// EXX_ABFS_INVERSE_MATRIX_DOUBLE_H