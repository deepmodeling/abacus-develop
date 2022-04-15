#include "cmd_neighbor.h"
#include "../src_parallel/parallel_common.h"
#include "../module_base/timer.h"

CMD_neighbor::CMD_neighbor()
{
    dim = 1;
    nlist = new int[dim];
    list = new int*[dim];
    for(int i=0; i<dim; i++)
	{
		list[i] = new int[dim];
	}
}

CMD_neighbor::~CMD_neighbor()
{
    delete[] nlist;
    for(int i=0; i<dim; i++)
	{
		delete[] list[i];
	}
    delete[] list;
}

//------------------------------------------------------
// In big cell, for atom i,
// at most one mirror atom of atom j is its adjent atom
//------------------------------------------------------
ModuleBase::Vector3<double> CMD_neighbor::cell_periodic(const ModuleBase::Vector3<double> a, const ModuleBase::Vector3<double> b)
{
	ModuleBase::Vector3<double> temp = a-b;
	while(temp.x<-0.5) temp.x+=1;
	while(temp.x>0.5) temp.x-=1;
	while(temp.y<-0.5) temp.y+=1;
    while(temp.y>0.5) temp.y-=1;
	while(temp.z<-0.5) temp.z+=1;
    while(temp.z>0.5) temp.z-=1;
	return temp;
}

#include "../module_base/mathzone.h"
//build the neighbor list
void CMD_neighbor::neighbor(const UnitCell_pseudo &ucell_c)
{
    ModuleBase::TITLE("CMD_neighbor", "Neighbor");
    ModuleBase::timer::tick("CMD_neighbor", "Neighbor");

    delete[] nlist;
    for(int i=0; i<dim; i++)
	{
		delete[] list[i];
	}
    delete[] list;

    dim = ucell_c.nat;
    nlist = new int[dim]();
    list = new int*[dim];
    for(int i=0; i<dim; i++)
	{
		list[i] = new int[dim]();
	}

#ifdef __MPI
    int *nlist_local = new int[dim]();
    int **list_local = new int*[dim];
    for(int i=0; i<dim; i++)
	{
		list_local[i] = new int[dim]();
	}
#endif

	for(int i=0; i<dim; i++)
	{
		for(int j=i+1; j<dim; j++)
		{
            int rank = ((2*dim + 1 - i)*i/2 + j - i - 1) % GlobalV::NPROC;
            if(rank == GlobalV::MY_RANK)
            {
                int T1 = ucell_c.iat2it[i];
                int I1 = ucell_c.iat2ia[i];
                ModuleBase::Vector3<double> taud1 = ucell_c.atoms[T1].taud[I1];

                int T2 = ucell_c.iat2it[j];
                int I2 = ucell_c.iat2ia[j];
                ModuleBase::Vector3<double> taud2 = ucell_c.atoms[T2].taud[I2];

			    ModuleBase::Vector3<double> tempd = cell_periodic(taud1,taud2);
                ModuleBase::Vector3<double> temp;
                ModuleBase::Mathzone::Direct_to_Cartesian(
			        tempd.x, tempd.y, tempd.z,
			        ucell_c.latvec.e11, ucell_c.latvec.e12, ucell_c.latvec.e13,
			        ucell_c.latvec.e21, ucell_c.latvec.e22, ucell_c.latvec.e23,
			        ucell_c.latvec.e31, ucell_c.latvec.e32, ucell_c.latvec.e33,
			        temp.x, temp.y, temp.z);

			    if(temp.norm()*ucell_c.lat0 <= GlobalV::SEARCH_RADIUS)
			    {
#ifdef __MPI
                    list_local[i][nlist_local[i]]=j;
                    nlist_local[i]++;
                    list_local[j][nlist_local[j]]=i;
                    nlist_local[j]++;
#else
				    list[i][nlist[i]]=j;
                    nlist[i]++;
                    list[j][nlist[j]]=i;
                    nlist[j]++;
#endif
			    }
            }
        }
	}

#ifdef __MPI
    comm_list(dim, nlist_local, list_local, nlist, list);

    delete[] nlist_local;
    for(int i=0; i<dim; i++)
	{
		delete[] list_local[i];
	}
    delete[] list_local;
#endif

    ModuleBase::timer::tick("CMD_neighbor", "Neighbor");
    return;
}

void CMD_neighbor::comm_list(const int num, int *nlist_in, int **list_in, int *nlist_out, int **list_out)
{
#ifdef __MPI
    MPI_Allreduce(nlist_in, nlist_out, num, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int *send_num = new int[GlobalV::NPROC]();
    int *send_dis = new int[GlobalV::NPROC]();

    for(int i=0; i<num; ++i)
    {
        MPI_Allgather(&nlist_in[i], 1, MPI_INT, send_num, 1, MPI_INT, MPI_COMM_WORLD);
       
        for(int j=1; j<GlobalV::NPROC; ++j)
        {
            send_dis[j] = send_dis[j-1] + send_num[j-1];
        }

        MPI_Allgatherv(list_in[i], nlist_in[i], MPI_INT, 
                       list_out[i], send_num, send_dis, MPI_INT, MPI_COMM_WORLD);
    }

    delete[] send_num;
    delete[] send_dis;
#endif
    return;
}