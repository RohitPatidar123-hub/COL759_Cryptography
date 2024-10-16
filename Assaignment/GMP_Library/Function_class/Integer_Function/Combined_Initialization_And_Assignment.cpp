#include<iostream>
#include<gmp.h>
 
 int main()
 {
       mpz_t a,b,c,d,e;    //declare varible of int type
       mpz_init_set_ui(a,1234);
       mpz_init_set(b,a);

       mpz_init_set_si(c,-1234);

       mpz_init_set_d(d,1234.3333);
       if(mpz_init_set_str(e,"1000000000000000001000048000000000000000007000287",10)==0)
          {
                 std ::cout<<"Successfully set value of e to 1000000000000000001000048000000000000000007000287\n \n \n";
          }
        else {
                  printf("error occur \n \n");
             }
        gmp_printf("a =%Zd , b= %Zd ,c=%Zd ,d=%Zd,e = %Zd\n \n",a,b,c,d,e);      


       mpz_clears(a,b,c,d,e,NULL);


 }
