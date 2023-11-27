#include "./hcontainer_funcs.h"
#include "./transfer.h"

#ifdef __MPI
#include <mpi.h>

//#include <chrono>

namespace hamilt
{
// transfer the HContainer from serial object to parallel object
template <typename TR>
void transferSerial2Parallels(const hamilt::HContainer<TR>& hR_s, hamilt::HContainer<TR>* hR_p, const int serial_rank)
{
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hamilt::HTransSerial<TR>* trans_s = nullptr;
    if (my_rank == serial_rank)
    {
        trans_s = new hamilt::HTransSerial<TR>(size, const_cast<hamilt::HContainer<TR>*>(&hR_s));
    }
    hamilt::HTransPara<TR> trans_p(size, hR_p);
    // plan indexes
    //std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    if (my_rank == serial_rank)
    {
        // send indexes to other ranks
        for (int i = 0; i < size; ++i)
        {
            if (i == serial_rank)
                continue;
            trans_s->send_ap_indexes(i);
        }

        { // transfer serial_rank to serial_rank
            std::vector<int> tmp_indexes;
            trans_s->cal_ap_indexes(serial_rank, &tmp_indexes);
            trans_p.receive_ap_indexes(serial_rank, tmp_indexes.data(), tmp_indexes.size());
            trans_p.cal_orb_indexes(serial_rank, &tmp_indexes);
            trans_s->receive_orb_indexes(serial_rank, tmp_indexes.data(), tmp_indexes.size());
        }
        
        // receive indexes from other ranks
        for (int i = 0; i < size; ++i)
        {
            if (i == serial_rank)
                continue;
            trans_s->receive_orb_indexes(i);
        }
    }
    else // my_rank != serial_rank
    {
        // receive ap_indexes from serial_rank and then send orb_indexes to serial_rank
        trans_p.receive_ap_indexes(serial_rank);
        trans_p.send_orb_indexes(serial_rank);
    }
    std::vector<TR> all_values;
    long max_size = 0;
    if (my_rank == serial_rank)
    {
        // calculate max size of values
        max_size = trans_s->get_max_size();
        all_values.resize(max_size * size);
    }
    MPI_Bcast(&max_size, 1, MPI_LONG, serial_rank, MPI_COMM_WORLD);
    std::vector<TR> receive_values(max_size);
    /*std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // send data
    if (my_rank == serial_rank)
    {
        for (int i = 0; i < size; ++i)
        {
            trans_s->pack_data(i, (all_values.data() + i * max_size));
        }
    }

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> pre_scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // MPI_scatter to send values
    MPI_Scatter(all_values.data(),
                max_size,
                MPITraits<TR>::datatype(),
                receive_values.data(),
                max_size,
                MPITraits<TR>::datatype(),
                serial_rank,
                MPI_COMM_WORLD);
    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // receive data
    trans_p.receive_data(serial_rank, receive_values.data());

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> post_scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << " S2P: my_rank = " << my_rank << " indexes_time = " << elapsed_time0.count()
              << " data_trans_time = " << pre_scatter_time.count()<<" "<<scatter_time.count()
              <<" "<<post_scatter_time.count() << std::endl;*/

    if (my_rank == serial_rank)
    {
        delete trans_s;
    }
}

// transfer the HContainer from parallel object to serial object
template <typename TR>
void transferParallels2Serial(const hamilt::HContainer<TR>& hR_p, hamilt::HContainer<TR>* hR_s, const int serial_rank)
{
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hamilt::HTransSerial<TR>* trans_s = nullptr;
    if (my_rank == serial_rank)
    {
        trans_s = new hamilt::HTransSerial<TR>(size, hR_s);
    }
    hamilt::HTransPara<TR> trans_p(size, const_cast<hamilt::HContainer<TR>*>(&hR_p));

    // plan indexes
    //std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    if (my_rank == serial_rank)
    {
        { // transfer serial_rank to serial_rank
            std::vector<int> tmp_indexes;
            trans_s->cal_ap_indexes(serial_rank, &tmp_indexes);
            trans_p.receive_ap_indexes(serial_rank, tmp_indexes.data(), tmp_indexes.size());
            trans_p.cal_orb_indexes(serial_rank, &tmp_indexes);
            trans_s->receive_orb_indexes(serial_rank, tmp_indexes.data(), tmp_indexes.size());
        }

        // send indexes to other ranks
        for (int i = 0; i < size; ++i)
        {
            if (i == serial_rank)
                continue;
            trans_s->send_ap_indexes(i);
        }
        // receive indexes from other ranks
        for (int i = 0; i < size; ++i)
        {
            if (i == serial_rank)
                continue;
            trans_s->receive_orb_indexes(i);
        }
    }
    else // my_rank != serial_rank
    {
        trans_p.receive_ap_indexes(serial_rank);
        trans_p.send_orb_indexes(serial_rank);
    }
    /*std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // send data
    std::vector<TR> receive_values;
    long max_size;
    if (my_rank == serial_rank)
    {
        max_size = trans_s->get_max_size();
        receive_values.resize(max_size * size);
    }
    MPI_Bcast(&max_size, 1, MPI_LONG, serial_rank, MPI_COMM_WORLD);

    // MPI_gather to receive values
    std::vector<TR> send_values;
    send_values.resize(max_size);
    trans_p.pack_data(serial_rank, send_values.data());
    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> pre_gather_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/
    MPI_Gather(send_values.data(),
               max_size,
               MPITraits<TR>::datatype(),
               receive_values.data(),
               max_size,
               MPITraits<TR>::datatype(),
               serial_rank,
               MPI_COMM_WORLD);
    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> gather_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/
    if (my_rank == serial_rank)
    {
        for (int i = 0; i < size; ++i)
        {
            trans_s->receive_data(i, &receive_values[i * max_size]);
        }
    }

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> post_gather_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << " P2S: my_rank = " << my_rank << " elapsed_time0 = " << elapsed_time0.count()
              << " data_trans_time = " << pre_gather_time.count()<<" "<<gather_time.count()
              <<" "<<post_gather_time.count() << std::endl;*/

    if (my_rank == serial_rank)
    {
        delete trans_s;
    }
}

// transferSerials2Parallels
template <typename TR>
void transferSerials2Parallels(const hamilt::HContainer<TR>& hR_s, hamilt::HContainer<TR>* hR_p)
{
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hamilt::HTransSerial<TR> trans_s(size, const_cast<hamilt::HContainer<TR>*>(&hR_s));
    hamilt::HTransPara<TR> trans_p(size, hR_p);
    // plan indexes
    //std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    // transfer indexes with other ranks


    { // begin of indexes_transfer
    // -----------------------------------
    // int tools for MPI_alltoallv
    std::vector<int> sendbuf, receivebuf;
    std::vector<int> sendcounts(size), recvcounts(size), sdispls(size), rdispls(size);
    // -----------------------------------
    // prepare sendbuf and sendcounts and sdispls and size of receivebuf   
    for (int i = 0; i < size; ++i)
    { // transfer in same process
        std::vector<int> tmp_indexes;
        trans_s.cal_ap_indexes(i, &tmp_indexes);
        sendcounts[i] = tmp_indexes.size();
        sdispls[i] = sendbuf.size();
        sendbuf.insert(sendbuf.end(), tmp_indexes.begin(), tmp_indexes.end());
    }

    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // resize the receivebuf
    long recvbuf_size = 0;
    for (int i = 0; i < size; ++i)
    {
        recvbuf_size += recvcounts[i];
    }
    receivebuf.resize(recvbuf_size);
    rdispls[0] = 0;
    for (int i = 1; i < size; ++i)
    {
        rdispls[i] = rdispls[i - 1] + recvcounts[i - 1];
    }

    // MPI_Alltoallv to send indexes
    MPI_Alltoallv(sendbuf.data(),
                  sendcounts.data(),
                  sdispls.data(),
                  MPI_INT,
                  receivebuf.data(),
                  recvcounts.data(),
                  rdispls.data(),
                  MPI_INT,
                  MPI_COMM_WORLD);

    // receive indexes from other ranks
    sendbuf.clear();
    for (int i = 0; i < size; ++i)
    {
        trans_p.receive_ap_indexes(i, &receivebuf[rdispls[i]], recvcounts[i]);
        std::vector<int> tmp_indexes;
        trans_p.cal_orb_indexes(i, &tmp_indexes);
        sendcounts[i] = tmp_indexes.size();
        sdispls[i] = sendbuf.size();
        sendbuf.insert(sendbuf.end(), tmp_indexes.begin(), tmp_indexes.end());
    }

    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // resize the receivebuf
    recvbuf_size = 0;
    for (int i = 0; i < size; ++i)
    {
        recvbuf_size += recvcounts[i];
    }
    receivebuf.resize(recvbuf_size);
    rdispls[0] = 0;
    for (int i = 1; i < size; ++i)
    {
        rdispls[i] = rdispls[i - 1] + recvcounts[i - 1];
    }

    // MPI_Alltoallv to send indexes
    MPI_Alltoallv(sendbuf.data(),
                  sendcounts.data(),
                  sdispls.data(),
                  MPI_INT,
                  receivebuf.data(),
                  recvcounts.data(),
                  rdispls.data(),
                  MPI_INT,
                  MPI_COMM_WORLD);

    // receive indexes from other ranks
    for (int i = 0; i < size; ++i)
    {
        trans_s.receive_orb_indexes(i, &receivebuf[rdispls[i]], recvcounts[i]);
    }

    }//end of indexes_transfer

    { // begin of data_transfer
    // -----------------------------------
    // TR tools for MPI_alltoallv
    std::vector<TR> sendbuf, receivebuf;
    std::vector<int> sendcounts(size), recvcounts(size), sdispls(size), rdispls(size);
    // -----------------------------------
    // prepare sendbuf and sendcounts and sdispls and size of receivebuf

    trans_s.get_value_size(sendcounts.data());
    sdispls[0] = 0;
    long sendbuf_size = sendcounts[0];
    for (int i = 1; i < size; ++i)
    {
        sdispls[i] = sdispls[i - 1] + sendcounts[i - 1];
        sendbuf_size += sendcounts[i];
    }
    sendbuf.resize(sendbuf_size);
    trans_p.get_value_size(recvcounts.data());

    long recvbuf_size = 0;
    for (int i = 0; i < size; ++i)
    {
        rdispls[i] = recvbuf_size;
        recvbuf_size += recvcounts[i];
    }
    receivebuf.resize(recvbuf_size);

    /*std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // send data
    for (int i = 0; i < size; ++i)
    {
        if(sendcounts[i] > 0)
        {
            trans_s.pack_data(i, (sendbuf.data() + sdispls[i]));
        }
    }

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> pre_scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // MPI_Alltoallv to send values
    MPI_Alltoallv(sendbuf.data(),
                  sendcounts.data(),
                  sdispls.data(),
                  MPITraits<TR>::datatype(),
                  receivebuf.data(),
                  recvcounts.data(),
                  rdispls.data(),
                  MPITraits<TR>::datatype(),
                  MPI_COMM_WORLD);

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // receive data
    for (int i = 0; i < size; ++i)
    {
        if(recvcounts[i] > 0)
        {
            trans_p.receive_data(i, (receivebuf.data() + rdispls[i]));
        }
    }
    } // end of data_transfer

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> post_scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << " S2P: my_rank = " << my_rank << " indexes_time = " << elapsed_time0.count()
              << " data_trans_time = " << pre_scatter_time.count()<<" "<<scatter_time.count()
              <<" "<<post_scatter_time.count() << std::endl;*/

}

// transferParallels2Serials
template <typename TR>
void transferParallels2Serials(const hamilt::HContainer<TR>& hR_p, hamilt::HContainer<TR>* hR_s)
{
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    hamilt::HTransPara<TR> trans_p(size, const_cast<hamilt::HContainer<TR>*>(&hR_p));
    hamilt::HTransSerial<TR> trans_s(size, hR_s);
    // plan indexes
    //std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    // transfer indexes with other ranks


    { // begin of indexes_transfer
    // -----------------------------------
    // int tools for MPI_alltoallv
    std::vector<int> sendbuf, receivebuf;
    std::vector<int> sendcounts(size), recvcounts(size), sdispls(size), rdispls(size);
    // -----------------------------------
    
    for (int i = 0; i < size; ++i)
    { // transfer in same process
        std::vector<int> tmp_indexes;
        trans_s.cal_ap_indexes(i, &tmp_indexes);
        sendcounts[i] = tmp_indexes.size();
        sdispls[i] = sendbuf.size();
        sendbuf.insert(sendbuf.end(), tmp_indexes.begin(), tmp_indexes.end());
    }

    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // resize the receivebuf
    long recvbuf_size = 0;
    for (int i = 0; i < size; ++i)
    {
        recvbuf_size += recvcounts[i];
    }
    receivebuf.resize(recvbuf_size);
    rdispls[0] = 0;
    for (int i = 1; i < size; ++i)
    {
        rdispls[i] = rdispls[i - 1] + recvcounts[i - 1];
    }

    // MPI_Alltoallv to send indexes
    MPI_Alltoallv(sendbuf.data(),
                  sendcounts.data(),
                  sdispls.data(),
                  MPI_INT,
                  receivebuf.data(),
                  recvcounts.data(),
                  rdispls.data(),
                  MPI_INT,
                  MPI_COMM_WORLD);

    // receive indexes from other ranks
    sendbuf.clear();
    for (int i = 0; i < size; ++i)
    {
        trans_p.receive_ap_indexes(i, &receivebuf[rdispls[i]], recvcounts[i]);
        std::vector<int> tmp_indexes;
        trans_p.cal_orb_indexes(i, &tmp_indexes);
        sendcounts[i] = tmp_indexes.size();
        sdispls[i] = sendbuf.size();
        sendbuf.insert(sendbuf.end(), tmp_indexes.begin(), tmp_indexes.end());
    }

    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // resize the receivebuf
    recvbuf_size = 0;
    for (int i = 0; i < size; ++i)
    {
        recvbuf_size += recvcounts[i];
    }
    receivebuf.resize(recvbuf_size);
    rdispls[0] = 0;
    for (int i = 1; i < size; ++i)
    {
        rdispls[i] = rdispls[i - 1] + recvcounts[i - 1];
    }

    // MPI_Alltoallv to send indexes
    MPI_Alltoallv(sendbuf.data(),
                  sendcounts.data(),
                  sdispls.data(),
                  MPI_INT,
                  receivebuf.data(),
                  recvcounts.data(),
                  rdispls.data(),
                  MPI_INT,
                  MPI_COMM_WORLD);

    // receive indexes from other ranks
    for (int i = 0; i < size; ++i)
    {
        trans_s.receive_orb_indexes(i, &receivebuf[rdispls[i]], recvcounts[i]);
    }

    }//end of indexes_transfer

    { // begin of data_transfer
    // -----------------------------------
    // TR tools for MPI_alltoallv
    std::vector<TR> sendbuf, receivebuf;
    std::vector<int> sendcounts(size), recvcounts(size), sdispls(size), rdispls(size);
    // -----------------------------------
    // prepare sendbuf and sendcounts and sdispls and size of receivebuf

    trans_p.get_value_size(sendcounts.data());
    sdispls[0] = 0;
    long sendbuf_size = sendcounts[0];
    for (int i = 1; i < size; ++i)
    {
        sdispls[i] = sdispls[i - 1] + sendcounts[i - 1];
        sendbuf_size += sendcounts[i];
    }
    sendbuf.resize(sendbuf_size);
    trans_s.get_value_size(recvcounts.data());

    long recvbuf_size = 0;
    for (int i = 0; i < size; ++i)
    {
        rdispls[i] = recvbuf_size;
        recvbuf_size += recvcounts[i];
    }
    receivebuf.resize(recvbuf_size);

    /*std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // send data
    for (int i = 0; i < size; ++i)
    {
        if(sendcounts[i] > 0)
        {
            trans_p.pack_data(i, (sendbuf.data() + sdispls[i]));
        }
    }

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> pre_scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // MPI_Alltoallv to send values
    MPI_Alltoallv(sendbuf.data(),
                  sendcounts.data(),
                  sdispls.data(),
                  MPITraits<TR>::datatype(),
                  receivebuf.data(),
                  recvcounts.data(),
                  rdispls.data(),
                  MPITraits<TR>::datatype(),
                  MPI_COMM_WORLD);

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();*/

    // receive data
    for (int i = 0; i < size; ++i)
    {
        if(recvcounts[i] > 0)
        {
            trans_s.receive_data(i, (receivebuf.data() + rdispls[i]));
        }
    }
    } // end of data_transfer

    /*end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> post_scatter_time
        = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << " S2P: my_rank = " << my_rank << " indexes_time = " << elapsed_time0.count()
              << " data_trans_time = " << pre_scatter_time.count()<<" "<<scatter_time.count()
              <<" "<<post_scatter_time.count() << std::endl;*/

}

template<typename TR>
void gatherParallels(const hamilt::HContainer<TR>& hR_p,
                     hamilt::HContainer<TR>* hR_s,
                     const int serial_rank)
{
    // gather <IJR>s from all ranks to serial_rank
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> para_ijrs = hR_p.get_ijr_info();
    if (my_rank == serial_rank)
    {
        hR_s->insert_ijrs(&para_ijrs);
        for (int i = 0; i < size; ++i)
        {
            if (i == serial_rank)
                continue;
            std::vector<int> tmp_ijrs;
            MPI_Status status;
            int tmp_size = 0;
            MPI_Recv(&tmp_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            tmp_ijrs.resize(tmp_size);
            MPI_Recv(tmp_ijrs.data(),
                     tmp_ijrs.size(),
                     MPI_INT,
                     i,
                     0,
                     MPI_COMM_WORLD,
                     &status);
            hR_s->insert_ijrs(&tmp_ijrs);
        }
        hR_s->allocate();
    }
    else
    {
        int tmp_size = para_ijrs.size();
        MPI_Send(&tmp_size, 1, MPI_INT, serial_rank, 0, MPI_COMM_WORLD);
        MPI_Send(para_ijrs.data(), para_ijrs.size(), MPI_INT, serial_rank, 0, MPI_COMM_WORLD);
    }
    // gather values from Parallels to target serial_rank
    transferParallels2Serial(hR_p, hR_s, serial_rank);
}

// specialize for double and std::complex<double>
template void transferSerial2Parallels(const hamilt::HContainer<double>& hR_s,
                                      hamilt::HContainer<double>* hR_p,
                                      const int serial_rank);
template void transferSerial2Parallels(const hamilt::HContainer<std::complex<double>>& hR_s,
                                      hamilt::HContainer<std::complex<double>>* hR_p,
                                      const int serial_rank);
template void transferParallels2Serial(const hamilt::HContainer<double>& hR_p,
                                      hamilt::HContainer<double>* hR_s,
                                      const int serial_rank);
template void transferParallels2Serial(const hamilt::HContainer<std::complex<double>>& hR_p,
                                      hamilt::HContainer<std::complex<double>>* hR_s,
                                      const int serial_rank);
template void transferSerials2Parallels(const hamilt::HContainer<double>& hR_s,
                                        hamilt::HContainer<double>* hR_p);
template void transferSerials2Parallels(const hamilt::HContainer<std::complex<double>>& hR_s,
                                        hamilt::HContainer<std::complex<double>>* hR_p);

template void transferParallels2Serials(const hamilt::HContainer<double>& hR_p,
                                        hamilt::HContainer<double>* hR_s);
template void transferParallels2Serials(const hamilt::HContainer<std::complex<double>>& hR_p,
                                        hamilt::HContainer<std::complex<double>>* hR_s);
template void gatherParallels(const hamilt::HContainer<double>& hR_p,
                              hamilt::HContainer<double>* hR_s,
                              const int serial_rank);
template void gatherParallels(const hamilt::HContainer<std::complex<double>>& hR_p,
                                hamilt::HContainer<std::complex<double>>* hR_s,
                                const int serial_rank);

} // namespace hamilt

#endif // __MPI