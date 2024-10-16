/*
GMP's conversion functions facilitate the interaction between multi-precision integers and standard C data types,
 making it easier to integrate GMP with other parts of your program. 
 These functions provide flexible ways to handle conversions while considering the limitations of the target types.
*/

#include<iostream>
#include<gmp.h>

int main()
{
    mpz_t big_number;
    mpz_init_set_str(big_number, "12345678901234567890", 10);
    signed long int small = mpz_get_si(big_number);
    double result = mpz_get_d(big_number);
    signed long int exponent;
    double mantissa = mpz_get_d_2exp(&exponent, big_number);
    char *string_rep = mpz_get_str(NULL, 10, big_number);

    gmp_printf("big number = %Zd\n \n",big_number);

    std :: cout<<"small : "<<small<<", result : "<<result<<", Exponent : "<<exponent<<", Mantisa :"<<mantissa<<", string_rep="<<string_rep;

    mpz_clear(big_number);

}
