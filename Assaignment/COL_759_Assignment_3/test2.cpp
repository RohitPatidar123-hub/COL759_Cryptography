#include <algorithm>
#include <cmath>
#include <gmpxx.h>
#include <iostream>
#include <mpi.h>
#include <vector>

// Function to generate all prime numbers up to 'limit' using Sieve of
// Eratosthenes
std::vector<unsigned long> generatePrimes(unsigned long limit) {
  std::vector<bool> is_prime(limit + 1, true);
  is_prime[0] = is_prime[1] = false;

  unsigned long sqrt_limit = static_cast<int>(std::sqrt(limit));
  for (unsigned long p = 2; p <= sqrt_limit; ++p) {
    if (is_prime[p]) {
      for (unsigned long multiple = p * p; multiple <= limit; multiple += p) {
        is_prime[multiple] = false;
      }
    }
  }

  std::vector<unsigned long> primes;
  for (unsigned long p = 2; p <= limit; ++p) {
    if (is_prime[p]) {
      primes.push_back(p);
    }
  }
  return primes;
}

// Function to compute the Legendre symbol (a/p) using GMP
// Returns 1 if 'a' is a quadratic residue modulo 'p'
// Returns -1 if 'a' is a non-residue modulo 'p'
// Returns 0 if 'p' divides 'a'
int legendreSymbol(mpz_class n, int p) {
   mpz_t a; 
  mpz_init_set(a, n.get_mpz_t());
  if (mpz_divisible_ui_p(a, p)) {
    return 0;
  }
  // Compute a^((p-1)/2) mod p
  mpz_t exponent, result, mod;
  mpz_inits(exponent, result, mod, NULL);

  mpz_set_ui(mod, p);
  mpz_set_ui(exponent, (p - 1) / 2);
  mpz_powm(result, a, exponent, mod);

  // Convert result to integer
  unsigned long res = mpz_get_ui(result);
  mpz_clears(exponent, result, mod, NULL);

  if (res == 1) {
    return 1;
  } else {
    return -1;
  }
}

// Function to compute Q(x) = (x + m)^2 - n
void compute_Qx(mpz_class &result,int a, mpz_class m, mpz_class n) { 
            mpz_class x=a;
            result = (x + m) * (x + m) - n; 
           // std::cout<<"Qx : "<<result<<"  ,";
            
     }

 mpz_class multiplyByUInt(mpz_class a, unsigned int multiplier) {
    mpz_class result;
    mpz_mul_ui(result.get_mpz_t(), a.get_mpz_t(), multiplier);
    return result;
}    

// Function to factorize Q(x) over the Factor Base
// Returns a vector of exponents (including -1) if Q(x) is smooth; otherwise,
// returns an empty vector
std::vector<int> factorize_Qx(int Qx, const std::vector<int> &factor_base) {
  std::vector<int> exponents(factor_base.size(), 0);
  int original_Qx = Qx;

  if (Qx < 0) {
    exponents[0] = 1; // Exponent of -1 is 1
    Qx = -Qx;
  }

  for (size_t i = 1; i < factor_base.size(); ++i) {
    int p = factor_base[i];
    while (Qx % p == 0) {
      exponents[i]++;
      Qx /= p;
    }
  }

  if (Qx == 1) { // Successfully factorized over Factor Base
    return exponents;
  } else {
    return {}; // Not smooth
  }
}

// Function to print the exponent matrix
void printMatrix(const std::vector<std::vector<int>> &matrix,
                 const std::vector<int> &factor_base) {
  std::cout << "Exponent Matrix (mod 2):\n";
  // Header
  std::cout << "Row\t";
  for (const auto &p : factor_base) {
    std::cout << p << "\t";
  }
  std::cout << "\n";

  // Rows
  for (size_t i = 0; i < matrix.size(); ++i) {
    std::cout << i + 1 << "\t";
    for (const auto &val : matrix[i]) {
      std::cout << val % 2 << "\t";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

// Function to perform Gaussian Elimination over GF(2) and find dependencies
std::vector<std::vector<int>>
findDependencies(std::vector<std::vector<int>> matrix_mod2) {
  int num_rows = matrix_mod2.size();
  if (num_rows == 0)
    return {};

  int num_cols = matrix_mod2[0].size();
  std::vector<int> pivot_col(num_cols, -1);
  // Initialize an identity matrix to track dependencies
  std::vector<std::vector<int>> identity(num_rows,
                                         std::vector<int>(num_rows, 0));
  for (int i = 0; i < num_rows; ++i) {
    identity[i][i] = 1;
  }

  // Perform Gaussian elimination
  for (int col = 0; col < num_cols; ++col) {
    // Find a pivot row
    int pivot_row = -1;
    for (int row = col; row < num_rows; ++row) {
      if (matrix_mod2[row][col] == 1) {
        pivot_row = row;
        break;
      }
    }

    if (pivot_row == -1) {
      continue; // No pivot in this column
    }

    // Swap current row with pivot_row if necessary
    if (pivot_row != col) {
      std::swap(matrix_mod2[col], matrix_mod2[pivot_row]);
      std::swap(identity[col], identity[pivot_row]);
    }

    pivot_col[col] = col;

    // Eliminate all other 1's in this column
    for (int row = 0; row < num_rows; ++row) {
      if (row != col && matrix_mod2[row][col] == 1) {
        for (int c = 0; c < num_cols; ++c) {
          matrix_mod2[row][c] ^= matrix_mod2[col][c];
        }
        for (int c = 0; c < num_rows; ++c) {
          identity[row][c] ^= identity[col][c];
        }
      }
    }
  }

  // Identify dependencies (nullspace vectors)
  std::vector<std::vector<int>> dependencies;

  // Rows without a pivot correspond to dependencies
  for (int row = 0; row < num_rows; ++row) {
    bool is_zero = true;
    for (int col = 0; col < num_cols; ++col) {
      if (matrix_mod2[row][col] != 0) {
        is_zero = false;
        break;
      }
    }
    if (is_zero) {
      // The corresponding row in the identity matrix represents the dependency
      dependencies.push_back(identity[row]);
    }
  }

  return dependencies;
}

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
       n.set_str("3000000000000000003", 10); 
      // n.set_str("21", 10);                // Assign number to n
        m = sqrt(n) +1;                   // Compute square root of n

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

  std::vector<int> factor_base_primes; // Primes where (n/p) = 1
  std::vector<int> factor_base;        // Including -1

  if (world_rank == 0) {
    std::cout << "Quadratic Sieve (QS) Implementation\n";
    std::cout << "===================================\n";
    std::cout << "Target number (n): " << n << std::endl;
    std::cout << "Computed m (ceil(sqrt(n))): " << m << "\n" << std::endl;

    // Step 1: Generate primes up to 'm'
    std::vector<unsigned long> primes;
    //unsigned long m_1 = m.get_ui();
    primes = generatePrimes(100000);  

    // Step 2: Exclude primes that divide 'n'
    std::vector<int> primes_filtered;
    for (const auto &p : primes) {
      if (n % p != 0) { // Exclude if p divides n
        primes_filtered.push_back(p);
      }
    }

    // Initialize GMP variable for 'n'
     mpz_class mpz_n = mpz_class(n);
    

    std::cout << "Primes up to " << m << ":\n";
    for (const auto &p : primes) {
      std::cout << p << " ";
    }
    std::cout << "\n\nPrimes after excluding those that divide " << n << ":\n";
    for (const auto &p : primes_filtered) {
      std::cout << p << " ";
    }
    std::cout << "\n\n";

    for (const auto &p : primes_filtered) {
      int ls = legendreSymbol(mpz_n, p);
      std::cout << "Legendre symbol (" << n << "/" << p << ") = " << ls
                << std::endl;
      if (ls == 1) { // Include in Factor Base if n is a quadratic residue modulo p
        factor_base_primes.push_back(p);
        std::cout << p << " is a quadratic residue modulo " << p
                  << ". Included in Factor Base.\n"
                  << std::endl;
      } else {
        std::cout << p << " is NOT a quadratic residue modulo " << p
                  << ". Excluded from Factor Base.\n"
                  << std::endl;
      }
    }
    // Step 4: Construct the Factor Base by adding -1
    factor_base.push_back(-1); // Always include -1
    for (const auto &p : factor_base_primes) {
      factor_base.push_back(p);
    }

    // Display the Factor Base
    std::cout << "Final Factor Base (including -1):\n";
    for (const auto &p : factor_base) {
      std::cout << p << " ";
    }
    std::cout << "\n\n";

    // Send the size of the Factor Base to other processes
    int fb_size = factor_base.size();
    MPI_Bcast(&fb_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Send the Factor Base to other processes
    MPI_Bcast(factor_base.data(), fb_size, MPI_INT, 0, MPI_COMM_WORLD);

  } else {
    // Receive the size of the Factor Base
    int fb_size;
    MPI_Bcast(&fb_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Receive the Factor Base
    factor_base.resize(fb_size);
    MPI_Bcast(factor_base.data(), fb_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Extract factor_base_primes (exclude -1)
    factor_base_primes.assign(factor_base.begin() + 1, factor_base.end());
  }

  // Ensure all processes have received the Factor Base
  MPI_Barrier(MPI_COMM_WORLD);

  // Step 5: Sieving Process - Compute Q(x) for x = 0 to x = 99
  const int x_min = 0;
  const int x_max = 99;
  const int total_x = x_max - x_min + 1; // 100

  // Determine the number of x's per process
  int x_per_process = total_x / world_size;
  int remainder = total_x % world_size;

  // Determine the start and end x for each process
  int local_x_start, local_x_end;

  if (world_rank < remainder) {
    // Processes with rank < remainder get (x_per_process + 1) x's
    local_x_start = world_rank * (x_per_process + 1);
    local_x_end = local_x_start + x_per_process;
  } else {
    // Processes with rank >= remainder get x_per_process x's
    local_x_start = world_rank * x_per_process + remainder;
    local_x_end = local_x_start + x_per_process - 1;
  }

  // Handle edge cases where x_end might exceed x_max
  if (local_x_end > x_max) {
    local_x_end = x_max;
  }

  // Each process computes Q(x) for its assigned x's
  std::vector<int> local_Qx;
  std::vector<int> local_x;
  std :: cout<<"world rank :"<<world_rank<<"local_x_start :"<<local_x_start<<"local_x_end"<<local_x_end;
  for (int x = local_x_start; x <= local_x_end; ++x) {
    mpz_class Qx;
    compute_Qx(Qx,x, m, n);
    int Qx_int=Qx.get_si();
    std::cout<<"Qx : "<<Qx_int<<"  ,";
    local_Qx.push_back(Qx_int);
    local_x.push_back(x);
  }
  MPI_Barrier(MPI_COMM_WORLD);  //Ensure all process calculate Qx for their respective class
  // Now, gather all Q(x) values to the root process
  // First, gather the counts from each process
  int  local_count = local_Qx.size();
  std::vector<int> recv_counts(world_size, 0);

  MPI_Gather(&local_count, 1,MPI_INT, recv_counts.data(), 1,MPI_INT, 0,
             MPI_COMM_WORLD);        // send sixe of Qx
  std::cout<<"\n\n";
  // Now, prepare for Gatherv
  std::vector<int> displs;      // It tells MPI where each process's data starts in the final collected array.
  std::vector<int> all_Qx;      // will hold all corresponding all the gathered Q(x)
  std::vector<int> all_x;       // will hold all corresponding x value
  if (world_rank == 0) {
    displs.resize(world_size, 0);
    int total_recv = recv_counts[0];       //received from the root process itself Q(x)
    for (int i = 1; i < world_size; ++i) {
      displs[i] = displs[i - 1] + recv_counts[i - 1];
      std::cout<<" displs["<<i<<"] = "<<displs[i]<<recv_counts[i] <<" ";
      total_recv += recv_counts[i];
    }
    std::cout<<"\ntotal_recv : "<<total_recv<<"\n";
    all_Qx.resize(total_recv, 0); // Resize to hold all received Q(x)
    all_x.resize(total_recv, 0);  // Resize to hold all received x's
  }

  // Gather all Q(x) values
  MPI_Gatherv(local_Qx.data(), local_count, MPI_INT, all_Qx.data(),
              recv_counts.data(), displs.data(),MPI_INT, 0, MPI_COMM_WORLD);

  // Gather all x values
  //send_mpz_class_vector(local_Qx, 0, MPI_COMM_WORLD);
  MPI_Gatherv(local_x.data(), local_count, MPI_INT, all_x.data(),
              recv_counts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  // Root process assembles the final Q(x) array
  if (world_rank == 0) {
    // Initialize an array of size 100 with default values
    std::vector<int> final_Qx(total_x, 0);

    // Iterate through all gathered x's and Q(x)'s to populate the array
    for (size_t i = 0; i < all_x.size(); ++i) {
      int x = all_x[i];
      int Qx = all_Qx[i];
      if (x >= x_min && x <= x_max) {
        final_Qx[x] = Qx;
      }
    }

    // Display the final Q(x) array
    std::cout << "Final Q(x) Array:\n";
    for (int x = x_min; x <= x_max; ++x) {
      std::cout << "Q(" << x << ") = " << final_Qx[x]<<" yes or no " << std::endl;
    }
    std::cout << "\n";

    // Step 6: Factorize Q(x) and identify smooth relations
    std::vector<std::vector<int>> smooth_relations; // Exponent vectors
    std::vector<int> smooth_x;                      // Corresponding x values
    std::vector<int> smooth_Qx;                     // Corresponding Q(x) values

    std::cout << "Factoring Q(x) over the Factor Base and identifying smooth "
                 "relations...\n"
              << std::endl;

    for (int x = x_min; x <= x_max; ++x) {
      int Qx = final_Qx[x];
      std::vector<int> exponents = factorize_Qx(Qx, factor_base);
      if (!exponents.empty()) { // Q(x) is smooth
        smooth_relations.emplace_back(exponents);
        smooth_x.push_back(x);
        smooth_Qx.push_back(Qx);
        std::cout << "Q(" << x << ") = " << Qx
                  << " is smooth over the Factor Base." << std::endl;
        std::cout << "Exponents: ";
        for (const auto &exp : exponents) {
          std::cout << exp << " ";
        }
        std::cout << "\nExponents (mod 2): ";
        for (const auto &exp : exponents) {
          std::cout << exp % 2 << " ";
        }
        std::cout << "\n" << std::endl;
      } else {
        std::cout << "Q(" << x << ") = " << Qx
                  << " is NOT smooth over the Factor Base." << std::endl;
        std::cout << "Discarding this relation.\n" << std::endl;
      }
    }

    std::cout << "Total smooth relations found: " << smooth_relations.size()
              << "\n"
              << std::endl;

    if (smooth_relations.empty()) {
      std::cout << "No smooth relations found. Increase the sieving range or "
                   "adjust the Factor Base."
                << std::endl;
      MPI_Finalize();
      return 0;
    }

    // Step 7: Construct the Exponent Matrix
    // Each row corresponds to a smooth relation
    // Each column corresponds to a prime in the Factor Base (including -1)
    // The entries are the exponents modulo 2

    std::cout << "Constructing the Exponent Matrix...\n" << std::endl;
    printMatrix(smooth_relations, factor_base);

    // Step 8: Perform Gaussian Elimination to Find Dependencies
    std::cout << "Performing Gaussian Elimination over GF(2) to find "
                 "dependencies...\n"
              << std::endl;
    // Create a copy of the matrix with exponents modulo 2
    std::vector<std::vector<int>> matrix_mod2 = smooth_relations;
    for (auto &row : matrix_mod2) {
      for (auto &val : row) {
        val = val % 2;
      }
    }
    std::vector<std::vector<int>> dependencies = findDependencies(matrix_mod2);

    // Display Dependencies
    if (dependencies.empty()) {
      std::cout << "No dependencies found.\n" << std::endl;
    } else {
      std::cout << "Dependencies Found:\n";
      for (size_t i = 0; i < dependencies.size(); ++i) {
        std::cout << "Dependency " << i + 1 << ": ";
        for (size_t j = 0; j < dependencies[i].size(); ++j) {
          if (dependencies[i][j] == 1) {
            std::cout << "Relation " << j + 1 << " ";
          }
        }
        std::cout << "\n";
      }
      std::cout << "\n";

      // Step 9: Use dependencies to compute 'a' and 'b', and find factors
      std::cout << "Attempting to find factors using dependencies...\n"
                << std::endl;

    mpz_t mpz_n;             // Declare the mpz_t variable
    mpz_init(mpz_n);         // Initialize it
    mpz_set(mpz_n, n.get_mpz_t()); // Set its value from mpz_class

      for (size_t i = 0; i < dependencies.size(); ++i) {
        // Initialize 'a' and 'b'
        mpz_t a, b;
        mpz_inits(a, b, NULL);
        mpz_set_ui(a, 1);
        mpz_set_ui(b, 1);

        // Exponent vector for 'b'
        std::vector<int> total_exponents(factor_base.size(), 0);

        // Multiply corresponding x + m for 'a' and collect exponents for 'b'
        for (size_t j = 0; j < dependencies[i].size(); ++j) {
          if (dependencies[i][j] == 1) {
            mpz_class x = smooth_x[j];
            mpz_class x_plus_m = x + m;
            //mpz_mul_ui(a, a, x_plus_m);
            mpz_class a = 1; // Example initialization of a
            a = multiplyByUInt(a, x_plus_m.get_ui());

            // Sum exponents
            for (size_t k = 0; k < factor_base.size(); ++k) {
              total_exponents[k] += smooth_relations[j][k];
            }
          }
        }

        // Divide exponents by 2 for 'b' (since exponents are even)
        for (size_t k = 0; k < total_exponents.size(); ++k) {
          total_exponents[k] /= 2;
        }

        // Compute 'b' as the product of primes raised to the total_exponents
        for (size_t k = 0; k < factor_base.size(); ++k) {
          if (total_exponents[k] > 0) {
            mpz_t temp;
            mpz_init(temp);
            mpz_set_si(temp, factor_base[k]);
            mpz_pow_ui(temp, temp, total_exponents[k]);
            mpz_mul(b, b, temp);
            mpz_clear(temp);
          }
        }

        // Compute 'a mod n' and 'b mod n'
        mpz_mod(a, a, mpz_n);
        mpz_mod(b, b, mpz_n);

        // Compute gcd(a - b, n)
        mpz_t gcd;
        mpz_init(gcd);
        mpz_t diff;
        mpz_init(diff);
        mpz_sub(diff, a, b);
        mpz_abs(diff, diff); // Ensure positive
        mpz_gcd(gcd, diff, mpz_n);

        // Check if gcd is a non-trivial factor
        if (mpz_cmp_ui(gcd, 1) > 0 && mpz_cmp(gcd, mpz_n) < 0) {
          // Compute the complementary factor
          mpz_t complementary_factor;
          mpz_init(complementary_factor);
          mpz_divexact(complementary_factor, mpz_n, gcd);

          std::cout << "Non-trivial factor found: ";
          mpz_out_str(stdout, 10, gcd);
          std::cout << "\nComplementary factor: ";
          mpz_out_str(stdout, 10, complementary_factor);
          std::cout << "\n" << std::endl;

          mpz_clears(a, b, gcd, diff, complementary_factor, NULL);
          mpz_clear(mpz_n);
          MPI_Finalize();
          return 0;
        }

        // Compute gcd(a + b, n)
        mpz_add(diff, a, b);
        mpz_mod(diff, diff, mpz_n); // Ensure within n
        mpz_gcd(gcd, diff, mpz_n);

        if (mpz_cmp_ui(gcd, 1) > 0 && mpz_cmp(gcd, mpz_n) < 0) {
          // Compute the complementary factor
          mpz_t complementary_factor;
          mpz_init(complementary_factor);
          mpz_divexact(complementary_factor, mpz_n, gcd);

          std::cout << "Non-trivial factor found: ";
          mpz_out_str(stdout, 10, gcd);
          std::cout << "\nComplementary factor: ";
          mpz_out_str(stdout, 10, complementary_factor);
          std::cout << "\n" << std::endl;

          mpz_clears(a, b, gcd, diff, complementary_factor, NULL);
          mpz_clear(mpz_n);
          MPI_Finalize();
          return 0;
        }

        mpz_clears(a, b, gcd, diff, NULL);
      }

      mpz_clear(mpz_n);

      std::cout
          << "No non-trivial factors found with the current dependencies.\n"
          << std::endl;
    }
  }

  // Finalize MPI environment
  MPI_Finalize();

  return 0;
}
