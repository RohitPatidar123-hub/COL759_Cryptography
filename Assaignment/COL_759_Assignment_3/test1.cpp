#include <iostream>
#include <gmp.h>

int main() {
    // Initialize a GMP arbitrary precision integer and floating point
    mpz_t num;
    mpf_t sqrt_result;

    // Set the number (40-digit number example)
    mpz_init_set_str(num, "1234567890123456789012345678901234567890", 10); // Base 10
    mpf_init(sqrt_result);

    // Set precision (can adjust this for better precision)
    mpf_set_default_prec(1000);

    // Convert integer to floating point for sqrt
    mpf_set_z(sqrt_result, num);

    // Calculate square root
    mpf_sqrt(sqrt_result, sqrt_result);

    // Print result
    std::cout << "Square root: ";
    mpf_out_str(stdout, 10, 0, sqrt_result); // Print in base 10
    std::cout << std::endl;

    // Clear memory
    mpz_clear(num);
    mpf_clear(sqrt_result);

    return 0;
}
