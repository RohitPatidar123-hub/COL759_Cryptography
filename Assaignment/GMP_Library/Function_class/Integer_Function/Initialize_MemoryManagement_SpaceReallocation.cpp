#include<iostream>
#include<gmp.h>
 
 int main()
 {
       mpz_t a,b,c,d,e,f,x,y,z,large;    //declare varible of int type
      
       mpz_init(a);                       //initialize single variable by 0 
       mpz_inits(b,c,d,e,f,x,y,z,NULL);   //initialize multiple varible by zero  , necessary to use NULL otherwise come error
       
      // GMP will automatically grow the variable as needed if it exceeds n bits. 
      //This function helps avoid frequent reallocation when the maximum size is known
       mpz_init2(large,1024);             //large no is initialize with the space of 1024 bits.save time
       gmp_printf("a =%Zd , b = %Zd,c = %Zd ,d=%Zd ,e=%Zd ,f=%Zd ,x=%Zd ,y=%Zd ,z=%Zd\n \n",a,b,c,d,e,f,x,y,z);
       gmp_printf("Large = %Zd \n \n",large);

       // mpz_realloc2 This function is typically used to avoid repeated automatic reallocations when you know the maximum size of the number in advance
       mpz_realloc2(2024);   
       mpz_clear(z);    //clear one variable at the time

       mpz_clears(a,b,c,d,e,f,x,y,NULL);    // clear multiple variable at the time 

       mpz_clear(large);     

 }

