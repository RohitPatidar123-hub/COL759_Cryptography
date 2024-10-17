#include <gmpxx.h>
#include <mpi.h>
#include <iostream>
#include <cstring>

int main(int argc, char *argv[]) {
    // Initialize MPI environment
    MPI_Init(&argc, &argv);

    int world_size; // Number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank; // Rank of the current process
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Define the target number 'n' and compute 'm'
    mpz_class n, m;
    std::string n_str, m_str;
    int n_len, m_len;
    int root_process = 0;

    if (world_rank == root_process) {
        n = 1828844741;                // Assign number to n
        m = sqrt(n);                   // Compute square root of n

        n_str = n.get_str();           // Convert mpz_class to string
        m_str = m.get_str();           // Convert mpz_class to string

        n_len = n_str.size() + 1;      // Get string lengths (+1 for null terminator)
        m_len = m_str.size() + 1;
    }

    // Broadcast the lengths of the strings first
    MPI_Bcast(&n_len, 1, MPI_INT, root_process, MPI_COMM_WORLD);
    MPI_Bcast(&m_len, 1, MPI_INT, root_process, MPI_COMM_WORLD);

    // Allocate memory for strings on non-root processes
    if (world_rank != root_process) {
        n_str.resize(n_len);
        m_str.resize(m_len);
    }

    // Broadcast the actual string data
    MPI_Bcast(&n_str[0], n_len, MPI_CHAR, root_process, MPI_COMM_WORLD);
    MPI_Bcast(&m_str[0], m_len, MPI_CHAR, root_process, MPI_COMM_WORLD);

    // Convert strings back to mpz_class on non-root processes
    if (world_rank != root_process) {
        n.set_str(n_str, 10);
        m.set_str(m_str, 10);
    }

    // Print the values to verify the broadcast
    std::cout << "Process " << world_rank << ": n = " << n << ", m = " << m << std::endl;

    // Finalize MPI environment
    MPI_Finalize();
    return 0;
}
