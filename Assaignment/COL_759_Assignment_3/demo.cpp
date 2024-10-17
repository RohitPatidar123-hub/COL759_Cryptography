#include <gmp.h>
#include <stdio.h>
#include <gmpxx.h>
void foo(mpz_t result, const mpz_t param, unsigned long n) {
    unsigned long i;

    // Multiply param by n and store the result in "result"
    mpz_mul_ui(result, param, n);

    // Loop to perform additional operations on "result"
    for (i = 1; i < n; i++) {
        mpz_add_ui(result, result, i * 7);  // Add i*7 to the result
    }
}

int main(void) {
    mpz_t r, n;  // Declare GMP variables

    // Initialize r (result) and n (input parameter)
    mpz_init(r);
    mpz_init_set_str(n, "1000000000000000000500023", 10);  // Initialize n with a string value

    // Call foo function with r as the result, n as the input, and 20 as the multiplier
    foo(r, n, 20L);
    std::vector<long long> prime;
    generate(n,prime);

    // Print the result
    gmp_printf("Result: %Zd\n", r);

    // Clear the variables (free the memory)
    mpz_clear(r);
    mpz_clear(n);

    return 0;
}
