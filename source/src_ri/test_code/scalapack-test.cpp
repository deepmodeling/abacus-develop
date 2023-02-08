#ifdef __MPI
#include<mpi.h>
#endif
#include "module_base/scalapack-connector.h"
#include "module_base/global_function.h"

void test_pblas()
{
	MPI_Init(NULL,NULL);
	int comm_sz;	MPI_Comm_size( MPI_COMM_WORLD, &comm_sz );
	int my_rank;	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

	auto ofs_error_code = []()
	{
		std::cout<<"MPI_SUCCESS          "<<MPI_SUCCESS          <<std::endl;
		std::cout<<"MPI_ERR_BUFFER       "<<MPI_ERR_BUFFER       <<std::endl;
		std::cout<<"MPI_ERR_COUNT        "<<MPI_ERR_COUNT        <<std::endl;
		std::cout<<"MPI_ERR_TYPE         "<<MPI_ERR_TYPE         <<std::endl;
		std::cout<<"MPI_ERR_TAG          "<<MPI_ERR_TAG          <<std::endl;
		std::cout<<"MPI_ERR_COMM         "<<MPI_ERR_COMM         <<std::endl;
		std::cout<<"MPI_ERR_RANK         "<<MPI_ERR_RANK         <<std::endl;
		std::cout<<"MPI_ERR_REQUEST      "<<MPI_ERR_REQUEST      <<std::endl;
		std::cout<<"MPI_ERR_ROOT         "<<MPI_ERR_ROOT         <<std::endl;
		std::cout<<"MPI_ERR_GROUP        "<<MPI_ERR_GROUP        <<std::endl;
		std::cout<<"MPI_ERR_OP           "<<MPI_ERR_OP           <<std::endl;
		std::cout<<"MPI_ERR_TOPOLOGY     "<<MPI_ERR_TOPOLOGY     <<std::endl;
		std::cout<<"MPI_ERR_DIMS         "<<MPI_ERR_DIMS         <<std::endl;
		std::cout<<"MPI_ERR_ARG          "<<MPI_ERR_ARG          <<std::endl;
		std::cout<<"MPI_ERR_UNKNOWN      "<<MPI_ERR_UNKNOWN      <<std::endl;
		std::cout<<"MPI_ERR_TRUNCATE     "<<MPI_ERR_TRUNCATE     <<std::endl;
		std::cout<<"MPI_ERR_OTHER        "<<MPI_ERR_OTHER        <<std::endl;
		std::cout<<"MPI_ERR_INTERN       "<<MPI_ERR_INTERN       <<std::endl;
		std::cout<<"MPI_ERR_IN_STATUS    "<<MPI_ERR_IN_STATUS    <<std::endl;
		std::cout<<"MPI_ERR_PENDING      "<<MPI_ERR_PENDING      <<std::endl;
		std::cout<<"MPI_ERR_KEYVAL       "<<MPI_ERR_KEYVAL       <<std::endl;
		std::cout<<"MPI_ERR_NO_MEM       "<<MPI_ERR_NO_MEM       <<std::endl;
		std::cout<<"MPI_ERR_BASE         "<<MPI_ERR_BASE         <<std::endl;
		std::cout<<"MPI_ERR_INFO_KEY     "<<MPI_ERR_INFO_KEY     <<std::endl;
		std::cout<<"MPI_ERR_INFO_VALUE   "<<MPI_ERR_INFO_VALUE   <<std::endl;
		std::cout<<"MPI_ERR_INFO_NOKEY   "<<MPI_ERR_INFO_NOKEY   <<std::endl;
		std::cout<<"MPI_ERR_SPAWN        "<<MPI_ERR_SPAWN        <<std::endl;
		std::cout<<"MPI_ERR_PORT         "<<MPI_ERR_PORT         <<std::endl;
		std::cout<<"MPI_ERR_SERVICE      "<<MPI_ERR_SERVICE      <<std::endl;
		std::cout<<"MPI_ERR_NAME         "<<MPI_ERR_NAME         <<std::endl;
		std::cout<<"MPI_ERR_WIN          "<<MPI_ERR_WIN          <<std::endl;
		std::cout<<"MPI_ERR_SIZE         "<<MPI_ERR_SIZE         <<std::endl;
		std::cout<<"MPI_ERR_DISP         "<<MPI_ERR_DISP         <<std::endl;
		std::cout<<"MPI_ERR_INFO         "<<MPI_ERR_INFO         <<std::endl;
		std::cout<<"MPI_ERR_LOCKTYPE     "<<MPI_ERR_LOCKTYPE     <<std::endl;
		std::cout<<"MPI_ERR_ASSERT       "<<MPI_ERR_ASSERT       <<std::endl;
		std::cout<<"MPI_ERR_RMA_CONFLICT "<<MPI_ERR_RMA_CONFLICT <<std::endl;
		std::cout<<"MPI_ERR_RMA_SYNC     "<<MPI_ERR_RMA_SYNC     <<std::endl;
	};

//	auto init_ictxt = [&](const char RC, const bool c2f)
//	{
		int dim[2] = {0,0};
		MPI_Dims_create( comm_sz, 2, dim );
//		MPI_Comm ictxt = c2f ? MPI_Comm_c2f(MPI_COMM_WORLD) : MPI_COMM_WORLD;
		MPI_Comm ictxt = MPI_Comm_c2f(MPI_COMM_WORLD);
//		ScalapackConnector::blacs_gridinit( ictxt, RC, dim[0], dim[1] );
		ScalapackConnector::blacs_gridinit( ictxt, 'R', dim[0], dim[1] );
		int nprow, npcol, myprow, mypcol;
		ScalapackConnector::blacs_gridinfo(ictxt,nprow, npcol, myprow, mypcol);
		assert(nprow==3);	assert(npcol==2);

		auto ofs_pro = [&]()
		{
			std::ofstream ofs("ofs_"+ModuleBase::GlobalFunc::TO_STRING(my_rank));
			ofs<<comm_sz<<"\t"<<dim[0]<<"\t"<<dim[1]<<"\t"<<nprow<<"\t"<<npcol<<"\t"<<myprow<<"\t"<<mypcol<<"\t"<<std::endl;
			ofs.close();
		};
		ofs_pro();
		
		auto init_desc = [&]( const int m, const int n, const int mb, const int nb, const int irsrc, const int icsrc )
		{
			const int mlld = ScalapackConnector::numroc( m, mb, myprow, irsrc, nprow );
			const int nlld = ScalapackConnector::numroc( n, nb, mypcol, icsrc, npcol );
			int info;
			std::vector<int> descv(9);
			ScalapackConnector::descinit(ModuleBase::GlobalFunc::VECTOR_TO_PTR(descv), m, n, mb, nb, irsrc, icsrc, ictxt, mlld, info);
			return descv;
		};		
//		return init_desc;
//	};

	auto ofs_matrix = [&](const std::vector<double>& v)
	{
		std::ofstream ofs("matrix_"+ModuleBase::GlobalFunc::TO_STRING(my_rank));
		for(const double i : v )
			ofs<<i<<std::endl;
		ofs.close();
	};
	auto ofs_desc = [&](const std::vector<int> &descv)
	{
		std::ofstream ofs("desc_"+ModuleBase::GlobalFunc::TO_STRING(my_rank));
		for(int i=0; i<9; ++i)
			ofs<<descv[i]<<"\t";
		ofs<<std::endl;
		ofs.close();
	};
	
/*
	std::vector<int> descv_5321_32 = init_ictxt('R',true)(5,3,2,1,0,0);
	int *desc_5321_32 = ModuleBase::GlobalFunc::VECTOR_TO_PTR(descv_5321_32);
	
	std::vector<int> descv_5521_32 = init_ictxt('R',true)(5,5,2,1,0,0);
	int *desc_5521_32 = ModuleBase::GlobalFunc::VECTOR_TO_PTR(descv_5521_32);
	
	std::vector<int> descv_5522_32 = init_ictxt('R',true)(5,5,2,2,0,0);
	int *desc_5522_32 = ModuleBase::GlobalFunc::VECTOR_TO_PTR(descv_5522_32);
	
	std::vector<int> descv_5511_32 = init_ictxt('R',true)(5,5,1,1,0,0);
	int *desc_5511_32 = ModuleBase::GlobalFunc::VECTOR_TO_PTR(descv_5511_32);
	ofs_desc(descv_5511_32);
*/
	std::vector<int> descv_5321_32 = init_desc(5,3,2,1,0,0);
	int *desc_5321_32 = ModuleBase::GlobalFunc::VECTOR_TO_PTR(descv_5321_32);
	
	std::vector<int> descv_5521_32 = init_desc(5,5,2,1,0,0);
	int *desc_5521_32 = ModuleBase::GlobalFunc::VECTOR_TO_PTR(descv_5521_32);
	
	std::vector<int> descv_5522_32 = init_desc(5,5,2,2,0,0);
	int *desc_5522_32 = ModuleBase::GlobalFunc::VECTOR_TO_PTR(descv_5522_32);
	
	std::vector<int> descv_5511_32 = init_desc(5,5,1,1,0,0);
	int *desc_5511_32 = ModuleBase::GlobalFunc::VECTOR_TO_PTR(descv_5511_32);
	
	/*
	auto init_matrix_A5624_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{0,1,2,3,6,7,8,9};
			case 1:	return std::vector<double>{4,5,10,11};
			case 2:	return std::vector<double>{12,13,14,15,18,19,20,21};
			case 3:	return std::vector<double>{16,17,22,23};
			case 4:	return std::vector<double>{24,25,26,27};
			case 5:	return std::vector<double>{28,29};
		}
	};
	auto init_matrix_I21 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,0,0,0};
			case 1:	return std::vector<double>{0,0,0,1,0,0};
			case 2:	return std::vector<double>{0,1,0,0,0,0};
			case 3:	return std::vector<double>{0,0,0,0,1,0};
			case 4:	return std::vector<double>{0,0,1,0,0,0};
			case 5:	return std::vector<double>{0,0,0,0,0,1};
		}
	};
	auto init_matrix_I11 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,0,0,0};
			case 1:	return std::vector<double>{0,0,0,0,1,0};
			case 2:	return std::vector<double>{0,0,0,0,0,1};
			case 3:	return std::vector<double>{1,0,0,0,0,0};
			case 4:	return std::vector<double>{0,1,0,0,0,0};
			case 5:	return std::vector<double>{0,0,0,0,0,1};
		}
	};
	auto init_matrix_I4411 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,1};
			case 1:	return std::vector<double>{0,0,0,0};
			case 2:	return std::vector<double>{0,0,0,0};
			case 3:	return std::vector<double>{1,0,0,1};
		}
	};
	auto init_matrix_I4422 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,1};
			case 1:	return std::vector<double>{0,0,0,0};
			case 2:	return std::vector<double>{0,0,0,0};
			case 3:	return std::vector<double>{1,0,0,1};
		}
	};
	auto init_matrix_I5522 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,0,1,0,0,0,1};
			case 1:	return std::vector<double>{0,0,0,0,0,0};
			case 2:	return std::vector<double>{0,0,0,0,0,0};
			case 3:	return std::vector<double>{1,0,0,1};
		}
	};
	auto init_matrix_I5511_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,0,0};
			case 1:	return std::vector<double>{0,0,0,1};
			case 2:	return std::vector<double>{0,0,0,0,0,1};
			case 3:	return std::vector<double>{1,0,0,0};
			case 4:	return std::vector<double>{0,1,0};
			case 5:	return std::vector<double>{0,0};
		}
	};
	auto init_matrix_I5521_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,0,0,0};
			case 1:	return std::vector<double>{0,0,1,0};
			case 2:	return std::vector<double>{0,1,0,0,0,0};
			case 3:	return std::vector<double>{0,0,0,1};
			case 4:	return std::vector<double>{0,0,1};
			case 5:	return std::vector<double>{0,0};
		}
	};
	auto init_matrix_I5522_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,0,1,0};
			case 1:	return std::vector<double>{0,0,0,0};
			case 2:	return std::vector<double>{0,0,0,0,0,0};
			case 3:	return std::vector<double>{1,0,0,1};
			case 4:	return std::vector<double>{0,0,1};
			case 5:	return std::vector<double>{0,0};
		}
	};
	auto init_matrix_I3321_11 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,0,1,0,0,0,1};
		}
	};
	auto init_matrix_I3321_21 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,0,0,1};
			case 1:	return std::vector<double>{0,1,0};
		}
	};
	auto init_matrix_A5522_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,2,5,6,7,10};
			case 1:	return std::vector<double>{3,4,8,9};
			case 2:	return std::vector<double>{11,12,15,16,17,20};
			case 3:	return std::vector<double>{13,14,18,19};
			case 4:	return std::vector<double>{21,22,25};
			case 5:	return std::vector<double>{23,24};
		}		
	};
	*/
	auto init_matrix_I3311_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0};
			case 1:	return std::vector<double>{0};
			case 2:	return std::vector<double>{0,0};
			case 3:	return std::vector<double>{1};
			case 4:	return std::vector<double>{0,1};
			case 5:	return std::vector<double>{0};
		}
	};
	auto init_matrix_I5511_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,0,0,0};
			case 1:	return std::vector<double>{0,0,0,1};
			case 2:	return std::vector<double>{0,0,0,0,0,1};
			case 3:	return std::vector<double>{1,0,0,0};
			case 4:	return std::vector<double>{0,1,0};
			case 5:	return std::vector<double>{0,0};
		}
	};
	auto init_matrix_I5521_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,0,0,0};
			case 1:	return std::vector<double>{0,1,0,0};
			case 2:	return std::vector<double>{0,0,1,0,0,0};
			case 3:	return std::vector<double>{0,0,0,1};
			case 4:	return std::vector<double>{0,0,1};
			case 5:	return std::vector<double>{0,0};
		}
	};
	auto init_matrix_I5522_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,0,0,1,0,0};
			case 1:	return std::vector<double>{0,0,0,0};
			case 2:	return std::vector<double>{0,0,0,0,0,0};
			case 3:	return std::vector<double>{1,0,0,1};
			case 4:	return std::vector<double>{0,0,1};
			case 5:	return std::vector<double>{0,0};
		}
	};
	auto init_matrix_A5321_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,6,3,8};
			case 1:	return std::vector<double>{2,7};
			case 2:	return std::vector<double>{11,16,13,18};
			case 3:	return std::vector<double>{12,17};
			case 4:	return std::vector<double>{21,23};
			case 5:	return std::vector<double>{22};
		}
	};
	auto init_matrix_A5521_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,6,3,8,5,10};
			case 1:	return std::vector<double>{2,7,4,9};
			case 2:	return std::vector<double>{11,16,13,18,15,20};
			case 3:	return std::vector<double>{12,17,14,19};
			case 4:	return std::vector<double>{21,23,25};
			case 5:	return std::vector<double>{22,24};
		}
	};
	auto init_matrix_A5522_32 = [&]()
	{
		switch(my_rank)
		{
			case 0:	return std::vector<double>{1,6,2,7,5,10};
			case 1:	return std::vector<double>{3,8,4,9};
			case 2:	return std::vector<double>{11,16,12,17,15,20};
			case 3:	return std::vector<double>{13,18,14,19};
			case 4:	return std::vector<double>{21,22,25};
			case 5:	return std::vector<double>{23,24};
		}
	};

	std::vector<double> mA5321_32 = init_matrix_A5321_32();
	std::vector<double> mA5522_32 = init_matrix_A5522_32();
	std::vector<double> mI3311_32 = init_matrix_I3311_32();
	std::vector<double> mI5511_32 = init_matrix_I5511_32();
	std::vector<double> mI5521_32 = init_matrix_I5521_32();
	std::vector<double> mC(100);
	ScalapackConnector::pgemm( 'T','N', 3,5,5, 1, ModuleBase::GlobalFunc::VECTOR_TO_PTR(mA5321_32),1,1,desc_5321_32, ModuleBase::GlobalFunc::VECTOR_TO_PTR(mI5511_32),1,1,desc_5511_32, 0, ModuleBase::GlobalFunc::VECTOR_TO_PTR(mC),1,1,desc_5522_32 );
	
	auto check_matrix = [&](const std::vector<double> &v1, const std::vector<double> &v2)
	{
		const size_t len = min(v1.size(),v2.size());
		for(int i=0; i<len; ++i)
			if(v1[i]!=v2[i])
				std::cout<<my_rank<<std::endl;
	};
	check_matrix(mA5522_32,mC);
	ofs_matrix(mC);

	MPI_Finalize();
}
