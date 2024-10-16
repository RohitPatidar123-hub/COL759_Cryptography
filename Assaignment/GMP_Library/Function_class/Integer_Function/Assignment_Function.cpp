 /*
mpz_set(): Copy one mpz_t value to another.
mpz_set_ui(): Set mpz_t from an unsigned long integer.
mpz_set_si(): Set mpz_t from a signed long integer.
mpz_set_d(): Set mpz_t from a double (integer part only).
mpz_set_q(): Set mpz_t from a rational number (mpq_t).
mpz_set_f(): Set mpz_t from a floating-point number (mpf_t).
mpz_set_str(): Set mpz_t from a string representation in a specified base.
mpz_swap(): Efficiently swap values between two mpz_t variables.

*/

#include<iostream>
#include<gmp.h>
 
 int main()
 {
       mpz_t a,b,c,d,e,f,g,h,large;    //declare varible of int type
       mpq_t  rational;
       mpf_t fl ;
       mpz_init(b);
       mpz_init(c);
       mpz_init(d);
       mpz_init(e);
       mpz_init(f);
       mpz_init(g);
       mpz_init(f);
       mpz_init(large);
       mpq_init(rational);
       mpf_init(fl);
       mpq_set_ui(rational,5,2);
       mpf_set_d(fl, 123.456);
       mpz_init_set_ui(a,15100);          //set value of mpz_t into unsigned long int 

       mpz_set(b,a);                      //copy one mpz_t variable into another mpz_t variable
       mpz_set_si(c,-1024);               // use this function to store signed long int    
       mpz_set_d(d,33333.333);            //Use this function when you want to convert a floating-point number to an mpz_t integer.
       mpz_set_q(e,rational);             //Use this function when you want to convert a rational number (mpq_t) to an integer (mpz_t).
       mpz_set_f(f,fl);                   //Use this function when you want to convert a floating-point number (mpf_t) to an integer (mpz_t).
       mpz_set_str(large, "3000000000000000000", 10);   //This function is useful for converting string representations of numbers to mpz_t.
       

       
       
       gmp_printf("a = %Zd , b = %Zd, c = %Zd , d = %Zd,e = %Zd , f = %Zd, large =%Zd \n \n ",a,b,c,d,e,f,large);
       
       mpz_swap(f,large);
       gmp_printf("After Swap : f = %Zd , large = %Zd \n \n",f,large);
        mpz_clear(a);
        mpz_clear(b);
        mpz_clear(c);
        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpf_clear(fl);
        mpq_clear(rational);
    

       

 }



