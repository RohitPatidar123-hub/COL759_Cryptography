#include <gmpxx.h>
#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    std::string n_str = ""; // Example input
    mpz_class n(n_str);
    
    // Check if we are the root process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    if (world_rank == 0) {
        mpz_class square;
        mpz_sqrt(square.get_mpz_t(), n.get_mpz_t()); // Get the integer square root of n
        
        unsigned long int exponent = 2;
        mpz_class result;
        mpz_mul(square.get_mpz_t(), n.get_mpz_t(), n.get_mpz_t());
        if (result == n) {
            std::cout << "Non-trivial factor found: " << square << std::endl;
            std::cout << "Complementary factor: " << square << std::endl;
            MPI_Finalize();
            return 0;
        }
    } else {
        mpz_class square_root;
        mpz_sqrt(square_root.get_mpz_t(), n.get_mpz_t()); // Get the integer square root of n
        
        unsigned long int exponent = 2;
        mpz_class result;
        mpz_pow_ui(result.get_mpz_t(), square_root.get_mpz_t(), exponent);
        
        if (result == n) {
            // cout << "Non-trivial factor found: " << result << endl;
            // cout << "Complementary factor: " << result << endl;
            MPI_Finalize();
            return 0;
        }
    }
    
    MPI_Finalize();
    return 0;
}

