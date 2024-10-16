#include<iostream>
#include<gmp.h>
 
 int main()
 {
        mpz_t n,m;
        mpz_init(m);
        if(mpz_init_set_str(n,"100000000000000000000000000000000000000000000000000",10)==0)
          {
                 std ::cout<<"Successfully set value of n to 1000000000000000001000048000000000000000007000287\n \n \n";
          }
        else {
                  printf("error occur \n \n");
             }
        mpz_sqrt(m,n);
        gmp_printf("Square root of %Zd is %Zd ",n,m);     

        mpz_clears(n,m,NULL);     

 }
