#include "custom_collectives.h"

#include <mpi.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <cstring>
#include <vector>

void Custom_Scatter(int* sendbuf, int sendcount, MPI_Datatype sendtype,
                    int* recvbuf, int recvcount, MPI_Datatype recvtype,
                    int root, MPI_Comm comm) {
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    MPI_Status status;
    int d = ceil(log2(size));
    bool ab = false;
    if(log2(size) != d) {
        ab = true;
    }
    int flip = 1 << (d - 1);
    int mask = flip - 1;

    int ab_limit = (1 << d) / 2;
    int ab_step = ab_limit;
    int send_size = (recvcount) * (1 << d) / 2;
    int* tempbuf = (int*)malloc(size * recvcount * sizeof(int));
    if(rank == root) {
        memcpy(tempbuf, sendbuf, size * sendcount * sizeof(int));
    }

    for(int i = d - 1; i > -1; i--) {
        if(((rank ^ root) & mask) == 0) {
            int rec_send_size = send_size;
            if((rank ^ flip) < size) {
                if(((rank ^ flip) >= ab_limit || (rank ^ root) >= ab_limit) && ab) {
                    rec_send_size = (size - ab_limit) * (recvcount);
                }
                if(((rank ^ root) & flip) == 0) {
                    MPI_Send(tempbuf + ((rank ^ flip) - (rank ^ root)) * recvcount, rec_send_size, recvtype, rank ^ flip, 1, comm); 
                } else {
                    MPI_Recv(tempbuf, rec_send_size, recvtype, rank ^ flip, 1, comm, &status);
                }
            }
        }
        mask = mask >> 1;
        flip = flip >> 1;
        send_size = send_size >> 1;   
        ab_limit = ab_limit + (ab_step >> 2);
        ab_step = ab_step >> 2;
    }
    
    memcpy(recvbuf, tempbuf, recvcount * sizeof(int));
    free(tempbuf);
    return;
}


void ktoAll(int* buff, int count, int offset, MPI_Datatype datatype, std::vector<int> sends, std::vector<int> recvs, MPI_Comm comm) {

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    while(!recvs.empty()) {
        int index;
        int sender = false, receiver = false;
        auto position = std::find(sends.begin(), sends.end(), rank);
        if(position != sends.end()) {
            index = position - sends.begin();
            sender = true;
        } else {
            position = std::find(recvs.begin(), recvs.end(), rank);
            if(position == recvs.end())
                return;
            index = position - recvs.begin();
            receiver = true;
        }

        if(sender && index < (int)recvs.size()) {
            
            MPI_Send(buff + offset, count, datatype, recvs[index], 0, comm);
        } else if(receiver && index < (int)sends.size()) {
            
            MPI_Recv(buff + offset, count, datatype, sends[index], 0, comm, MPI_STATUS_IGNORE);
        }

        int min_size = std::min(sends.size(), recvs.size());
        sends.insert(sends.end(), recvs.begin(), recvs.begin() + min_size);
        recvs.erase(recvs.begin(), recvs.begin() + min_size);
    }
    return;
}

void Custom_Allgather(int* sendbuf, int sendcount, MPI_Datatype sendtype,
                      int* recvbuf, int recvcount, MPI_Datatype recvtype,
                      MPI_Comm comm) {
     /* get rank and size */
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    /* get number of steps */
    int steps = (int)std::ceil(std::log2(size));

    /* copy sendbuf's data at correct location in recvbuf */
    memcpy(recvbuf + rank * recvcount, sendbuf, sendcount * sizeof(int));

    int flip = 1, step_size = 1;

    /* perform log(p) steps of all_gather */
    for (int i = 0; i < steps; i++) {
        int neighbor = rank ^ flip;
        int r_count = step_size * recvcount;
        int s_count = step_size * sendcount;
        int send_offset = (rank/step_size) * step_size * sendcount;
        int recv_offset = (neighbor/step_size) * step_size * recvcount;

        if(send_offset + s_count > size * sendcount)
            s_count = size * sendcount - send_offset;
        if(recv_offset + r_count > size * recvcount)
            r_count = size * recvcount - recv_offset;
       if(neighbor < size) {
            MPI_Sendrecv(recvbuf + send_offset, s_count, sendtype, neighbor, 0,
                     recvbuf + recv_offset, r_count, recvtype, neighbor, 0,
                     comm, MPI_STATUS_IGNORE);
       }
        /* update mask for next step */
        flip <<= 1;
        step_size *= 2;
        /*MPI_Barrier(comm);*/
    }

    std::vector<int> sends;
    std::vector<int> recvs;
    int l = ceil(log2(size));
    int cap = 1 << l;
    int half = cap / 2;
    for(int i = 0; i < size; i++) {
        if(i < size - half || i >= half) 
            sends.push_back(i);
        else 
            recvs.push_back(i);
    }
    if(!sends.empty() && !recvs.empty())
        ktoAll(recvbuf, recvcount * (size - half), half * recvcount, recvtype, sends, recvs, comm);

    return;
}


void Custom_Allreduce(int* sendbuf, int* recvbuf, int count,
                      MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
    //Implement the code below
    //op will be MPI_SUM
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int depth = log2(size);
    int tempSize = size;

    if ((size & (size - 1)) != 0) {
        tempSize = 1 << (depth + 1);  
        depth = log2(tempSize);    
    }

    for (int i = 0; i < count; i++) {
        recvbuf[i] = sendbuf[i];
    }
    memcpy(recvbuf, sendbuf, sizeof(int) * count);

    int* tempbuf = (int*) malloc(count * sizeof(int));

    for (int j = depth - 1; j >= 0; j--) {
        int flip = 1 << j;
        int neighbor = rank ^ flip;

        if (neighbor < size && rank < size) {

            if (rank < neighbor) {
                // printf("Rank %d sends to Rank %d\n", rank, neighbor);
                MPI_Send(recvbuf, count, datatype, neighbor, 1, comm);

                // printf("Rank %d waits to receive from Rank %d\n", rank, neighbor);
                MPI_Recv(tempbuf, count, datatype, neighbor, 1, comm, MPI_STATUS_IGNORE);
            } else {
                // printf("Rank %d waits to receive from Rank %d\n", rank, neighbor);
                MPI_Recv(tempbuf, count, datatype, neighbor, 1, comm, MPI_STATUS_IGNORE);

                // printf("Rank %d sends to Rank %d\n", rank, neighbor);
                MPI_Send(recvbuf, count, datatype, neighbor, 1, comm);
            }

            if (op == MPI_SUM) {
                for (int k = 0; k < count; k++) {
                    recvbuf[k] += tempbuf[k];
                }
            }
        }
    }
    free(tempbuf);

    std::vector<int> sends;
    std::vector<int> recvs;
    int l = ceil(log2(size));
    int cap = 1 << l;
    int half = cap / 2;
    for(int i = 0; i < size; i++) {
        if(i < half) 
            sends.push_back(i);
        else 
            recvs.push_back(i);
    }
    if(!sends.empty() && !recvs.empty())
        ktoAll(recvbuf, count, 0, datatype, sends, recvs, comm);
    return;
}   

void Custom_Alltoall_Hypercube(int* sendbuf, int sendcount, MPI_Datatype sendtype,
                               int* recvbuf, int recvcount, MPI_Datatype recvtype,
                               MPI_Comm comm) {
   //Implement the code below

    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    MPI_Status status;

    int depth = log2(size);

    memcpy(recvbuf, sendbuf, sendcount * size * sizeof(int)); //using memcpy for speed since it's faster than the for loops I had earlier

    int partition = size >> 1;
    int* tempbuf = (int*) malloc(recvcount * size * sizeof(int));

    for (int j = depth - 1; j >= 0; j--) {
        int flip = 1 << j;
        int neighbor = rank ^ flip;
        // printf("Rank %d communicating with neighbor %d\n", rank, neighbor);

        if (neighbor < size) {
            if (rank < neighbor) {
                // printf("FROM LEFT Rank %d sends to Rank %d\n", rank, neighbor);
                MPI_Sendrecv(recvbuf + (recvcount * partition), recvcount * partition, recvtype, neighbor, 0,
                             tempbuf + (recvcount * partition), recvcount * partition, recvtype, neighbor, 0, comm, &status);

                memcpy(tempbuf, recvbuf, recvcount * partition * sizeof(int));

            } else {
                // printf("FROM RIGHT Rank %d sends to Rank %d\n", rank, neighbor);
                MPI_Sendrecv(recvbuf, recvcount * partition, recvtype, neighbor, 0,
                             tempbuf, recvcount * partition, recvtype, neighbor, 0, comm, &status);
                memcpy(tempbuf + (recvcount * partition), recvbuf + (recvcount * partition), recvcount * partition * sizeof(int));
            }
        }
        for (int i = 0; i < partition; i++) {
            memcpy(recvbuf + i * 2 * recvcount, tempbuf + i * recvcount, recvcount * sizeof(int));
            memcpy(recvbuf + (i * 2 + 1) * recvcount, tempbuf + (partition + i) * recvcount, recvcount * sizeof(int));
        }
    }

    free(tempbuf);
    return; 
}


void Custom_Alltoall_Arbitrary(int* sendbuf, int sendcount, MPI_Datatype sendtype,
                               int* recvbuf, int recvcount, MPI_Datatype recvtype,
                               MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    memcpy(recvbuf + (rank * recvcount), sendbuf + (rank * sendcount), recvcount * sizeof(int));
    for (int j = 1; j < size; j++) {
        int send_to = (rank + j) % size;
        int recv_from = ((rank - j) + size) % size;
        MPI_Sendrecv(sendbuf + (send_to * sendcount), sendcount, sendtype, send_to, 0, 
            recvbuf + (recv_from * recvcount), recvcount, recvtype, recv_from, 0, comm, MPI_STATUS_IGNORE);
    }
    return;
}