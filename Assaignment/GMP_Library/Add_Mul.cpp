#include<iostream>
#include<gmp.h>
 
 int main()
 {
       mpz_t a,b,c,x,y,z;    //declare varible of int type
       
       mpz_init(a);    // initialize a,b,c by 0
       mpz_init(b);
       mpz_init(c);
       mpz_init(x);
       mpz_init(y);
       mpz_init(z);

     
       mpz_set_ui(a,3000000000000000);
       mpz_set_ui(b,5000000000000000);    // assaignment function 
       mpz_set_str(x, "30000000000000000000000000000000000000000000000000", 10);  // Base 10
       mpz_set_str(y, "50000000000000000000000000000000000000000000000000", 10);  // Base 10
      


       gmp_printf("Large no. is : %Zd \nLarge no. is : %Zd \n \n ",x,y);
       mpz_add(z,x,y);                //arithmatic addition operation 
       gmp_printf("%Zd + %Zd = %Zd\n \n",x,y,z);

       mpz_mul(z,x,y);              //arithmetic mul operation
       gmp_printf("%Zd * %Zd = %Zd\n \n",x,y,z);
       mpz_init(c);                   //to clear previous value   

       mpz_add(c,a,b);                //arithmatic addition operation 
       gmp_printf("%Zd + %Zd = %Zd\n \n",a,b,c);
       mpz_init(c);                   //to clear previous value 

       mpz_sub(c,a,b);               // arithmetic subtraction operation 
       gmp_printf("%Zd - %Zd = %Zd\n \n",a,b,c);
       mpz_init(c);

       mpz_mul(c,a,b);              //arithmetic mul operation
       gmp_printf("%Zd * %Zd = %Zd\n \n",a,b,c);
       

        mpz_clear(a);  //clear memory when done
        mpz_clear(b);
        mpz_clear(c);
        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(z);


 }

