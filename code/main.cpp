#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <cstring>
#include "custom_collectives.h"

int find_arg_idx(int argc, char **argv, const char *option)
{
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], option) == 0)
        {
            return i;
        }
    }
    return -1;
}

void test_scatter(int n, int size, int root, MPI_Comm comm);
void test_allgather(int n, int size, MPI_Comm comm);
void test_allreduce(int n, int size, MPI_Comm comm);
void test_alltoall_hypercubic(int n, int size, MPI_Comm comm);
void test_alltoall_arbitrary(int n, int size, MPI_Comm comm);

void test_scatter(int n, int size, int root, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<int> sendbuf(n);
    std::vector<int> recvbuf(n / size);
    std::vector<int> expected(n / size);

    if (rank == 0)
    {
        for (int i = 0; i < n; ++i)
        {
            sendbuf[i] = i;
        }

        std::cout << "Problem Size: " << n << std::endl;
        std::cout << "Number of Processors: " << size << std::endl;
    }

    double start, end;
    if (rank == 0)
    {
        start = MPI_Wtime();
    }
    MPI_Scatter(sendbuf.data(), n / size, MPI_INT, expected.data(), n / size, MPI_INT, 0, comm);
    if (rank == 0)
    {
        end = MPI_Wtime();
        std::cout << "MPI_Scatter time: " << end - start << std::endl;
    }

    if (rank == 0)
    {
        start = MPI_Wtime();
    }
    Custom_Scatter(sendbuf.data(), n / size, MPI_INT, recvbuf.data(), n / size, MPI_INT, 0, comm);
    if (rank == 0)
    {
        end = MPI_Wtime();
        std::cout << "Custom_Scatter time: " << end - start << std::endl;
    }

    int correct = std::equal(recvbuf.begin(), recvbuf.end(), expected.begin());
    int result;
    MPI_Reduce(&correct, &result, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank == 0)
    {
        if (result == size)
        {
            std::cout << "Implementation is correct!" << std::endl;
        }
        else
        {
            std::cout << "Implementation failed" << std::endl;
        }
    }
}

void test_allgather(int n, int size, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<int> sendbuf(n / size);
    std::vector<int> recvbuf(n);
    std::vector<int> expected(n);

    for (int i = 0; i < n / size; ++i)
    {
        sendbuf[i] = rank * n / size + i;
    }

    if (rank == 0)
    {
        std::cout << "Problem Size: " << n << std::endl;
        std::cout << "Number of Processors: " << size << std::endl;
    }

    double start, end;
    if (rank == 0)
    {
        start = MPI_Wtime();
    }
    MPI_Allgather(sendbuf.data(), n / size, MPI_INT, expected.data(), n / size, MPI_INT, comm);
    if (rank == 0)
    {
        end = MPI_Wtime();
        std::cout << "MPI_Allgather time: " << end - start << std::endl;
//        for (int i = 0; i < n; ++i)
//        {
//            std::cout << expected[i] << " ";
//        }
//        std::cout << std::endl;
    }

    if (rank == 0)
    {
        start = MPI_Wtime();
    }
    Custom_Allgather(sendbuf.data(), n / size, MPI_INT, recvbuf.data(), n / size, MPI_INT, comm);
    if (rank == 0)
    {
        end = MPI_Wtime();
        std::cout << "Custom_Allgather time: " << end - start << std::endl;
//        for (int i = 0; i < n; ++i)
//        {
//            std::cout << recvbuf[i] << " ";
//        }
//        std::cout << std::endl;
    }

    int correct = std::equal(recvbuf.begin(), recvbuf.end(), expected.begin());
    int result;
    MPI_Reduce(&correct, &result, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank == 0)
    {
        if (result == size)
        {
            std::cout << "Implementation is correct!" << std::endl;
        }
        else
        {
            std::cout << "Implementation failed" << std::endl;
            printf("Expected: %d and Received: %d\n", size, result);
        }
    }
}

void test_allreduce(int n, int size, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<int> sendbuf(n);
    std::vector<int> recvbuf(n);
    std::vector<int> expected(n);

    for (int i = 0; i < n; ++i)
    {
        sendbuf[i] = rank * n + i;
    }

    if (rank == 0)
    {
        std::cout << "Problem Size: " << n << std::endl;
        std::cout << "Number of Processors: " << size << std::endl;
    }

    double start, end;
    if (rank == 0)
    {
        start = MPI_Wtime();
    }
    MPI_Allreduce(sendbuf.data(), expected.data(), n, MPI_INT, MPI_SUM, comm);
    if (rank == 0)
    {
        end = MPI_Wtime();
        std::cout << "MPI_Allreduce time: " << end - start << std::endl;
    }

    if (rank == 0)
    {
        start = MPI_Wtime();
    }
    Custom_Allreduce(sendbuf.data(), recvbuf.data(), n, MPI_INT, MPI_SUM, comm);
    if (rank == 0)
    {
        end = MPI_Wtime();
        std::cout << "Custom_Allreduce time: " << end - start << std::endl;
    }

    int correct = std::equal(recvbuf.begin(), recvbuf.end(), expected.begin());
    int result;
    MPI_Reduce(&correct, &result, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank == 0)
    {
        if (result == size)
        {
            std::cout << "Implementation is correct!" << std::endl;
        }
        else
        {
            std::cout << "Implementation failed" << std::endl;
        }
    }
}

void test_alltoall_hypercubic(int n, int size, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<int> sendbuf(n);
    std::vector<int> recvbuf_hypercube(n);
    std::vector<int> expected(n);

    for (int i = 0; i < n; ++i)
    {
        sendbuf[i] = rank * n + i;
    }

    if (rank == 0)
    {
        std::cout << "Problem Size: " << n << std::endl;
        std::cout << "Number of Processors: " << size << std::endl;
    }

    double start, end;
    if (rank == 0)
    {
        start = MPI_Wtime();
    }
    MPI_Alltoall(sendbuf.data(), n / size, MPI_INT, expected.data(), n / size, MPI_INT, comm);

    if (rank == 0)
    {
        end = MPI_Wtime();
        std::cout << "MPI_Alltoall time: " << end - start << std::endl;
    }

    if (rank == 0)
    {
        start = MPI_Wtime();
    }
    Custom_Alltoall_Hypercube(sendbuf.data(), n / size, MPI_INT, recvbuf_hypercube.data(), n / size, MPI_INT, comm);
    if (rank == 0)
    {
        end = MPI_Wtime();
        std::cout << "Custom_Alltoall_Hypercubic time: " << end - start << std::endl;
    }

    int correct_hypercube = std::equal(recvbuf_hypercube.begin(), recvbuf_hypercube.end(), expected.begin());
    int result_hypercube;
    MPI_Reduce(&correct_hypercube, &result_hypercube, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank == 0)
    {
        if (result_hypercube == size)
        {
            std::cout << "Implementation is correct!" << std::endl;
        }
        else
        {
            std::cout << "Implementation failed" << std::endl;
        }
    }
}

void test_alltoall_arbitrary(int n, int size, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<int> sendbuf(n);
    std::vector<int> recvbuf_arbitrary(n);
    std::vector<int> expected(n);

    for (int i = 0; i < n; ++i)
    {
        sendbuf[i] = rank * n + i;
    }

    if (rank == 0)
    {
        std::cout << "Problem Size: " << n << std::endl;
        std::cout << "Number of Processors: " << size << std::endl;
    }

    double start, end;
    if (rank == 0)
    {
        start = MPI_Wtime();
    }
    MPI_Alltoall(sendbuf.data(), n / size, MPI_INT, expected.data(), n / size, MPI_INT, comm);
    if (rank == 0)
    {
        end = MPI_Wtime();
        std::cout << "MPI_Alltoall time: " << end - start << std::endl;
    }

    if (rank == 0)
    {
        start = MPI_Wtime();
    }
    Custom_Alltoall_Arbitrary(sendbuf.data(), n / size, MPI_INT, recvbuf_arbitrary.data(), n / size, MPI_INT, comm);
    if (rank == 0)
    {
        end = MPI_Wtime();
        std::cout << "Custom_Alltoall_Arbitrary time: " << end - start << std::endl;
    }

    int correct_arbitrary = std::equal(recvbuf_arbitrary.begin(), recvbuf_arbitrary.end(), expected.begin());
    int result_arbitrary;
    MPI_Reduce(&correct_arbitrary, &result_arbitrary, 1, MPI_INT, MPI_SUM, 0, comm);

    if (rank == 0)
    {
        if (result_arbitrary == size)
        {
            std::cout << "Implementation is correct!" << std::endl;
        }
        else
        {
            std::cout << "Implementation failed" << std::endl;
        }
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int root = std::stoi(argv[3]);
    int prob_size = std::stoi(argv[2]);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(NULL) + rank);

    if (find_arg_idx(argc, argv, "-s") >= 0)
    {

        test_scatter(prob_size, size, root, MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    }

    if (find_arg_idx(argc, argv, "-g") >= 0)
    {

        test_allgather(prob_size, size, MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    }

    if (find_arg_idx(argc, argv, "-r") >= 0)
    {

        test_allreduce(prob_size, size, MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    }

    if (find_arg_idx(argc, argv, "-a") >= 0)
    {

        test_alltoall_arbitrary(prob_size, size, MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    }

    if (find_arg_idx(argc, argv, "-h") >= 0)
    {

        test_alltoall_hypercubic(prob_size, size, MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    }

    MPI_Finalize();
    return 0;
}