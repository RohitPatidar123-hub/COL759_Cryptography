#include<iostream>
#include<gmp.h>
 
 int main()
 {
        mpz_t n,d,r,q,c;
        mpz_inits(r,q,NULL);
        mpz_init_set_ui(n,10);
        mpz_init_set_ui(d,5);
         mpz_init_set_ui(c,5);
//Ceiling Division Functions
       std :: cout<<"Ceiling Division Functions\n ";
        mpz_cdiv_q(q,n,d);                                            // Computes the quotient  q=⌈n/d⌉
        gmp_printf(" %Zd =⌈%Zd / %Zd⌉ \n",q,n,d);

        mpz_cdiv_r(r,n,d);                                            //compute the remainder r = n mod d with oppoosite sign of d 
        gmp_printf(" %Zd = %Zd mod %Zd \n",r,n,d); 

        mpz_cdiv_qr(q,r,n,d);                                        // Computes both quotient and remainder.
        gmp_printf("q = %Zd ,r = %Zd ,n = %Zd ,d = %Zd \n\n",q,r,n,d) ;


//Floor Division Functions
        std ::cout<<"Floor Division Functions \n" ;
        mpz_inits(r,q,NULL);
        mpz_fdiv_q(q,n,d);                                            
        gmp_printf(" %Zd = ⌊%Zd / %Zd ⌋\n",q,n,d);                   // Computes the quotient q=⌊n/d⌋.

        mpz_fdiv_r(r,n,d);                                            
        gmp_printf(" %Zd = %Zd mod %Zd \n",r,n,d);                 // Computes the remainder  r = nmod d with the same sign as d

        mpz_fdiv_qr(q,r,n,d);                                        
        gmp_printf("q = %Zd ,r = %Zd ,n = %Zd ,d = %Zd \n\n",q,r,n,d) ;  // Computes both quotient and remainder.
               
//Truncate Division Functions
       std :: cout<<"Truncate Division Functions\n" ;
        mpz_inits(r,q,NULL);
        mpz_tdiv_q(q,n,d);                                            
        gmp_printf(" %Zd = ⌊%Zd / %Zd ⌋\n",q,n,d);                   // Computes the quotient q=⌊n/d⌋.

        mpz_tdiv_r(r,n,d);                                            
        gmp_printf(" %Zd = %Zd mod %Zd \n",r,n,d);                 // Computes the remainder  r = nmod d with the same sign as d

        mpz_tdiv_qr(q,r,n,d);                                        
        gmp_printf("q = %Zd ,r = %Zd ,n = %Zd ,d = %Zd \n\n",q,r,n,d) ;  // Computes both quotient and remainder.
               
// Unsigned Long Division Functions
    std :: cout<<"Unsigned Long Division Functions\n";
          unsigned long int Q,R;
          mpz_inits(r,q,NULL);
          Q= mpz_cdiv_q_ui(q,n,3);
          gmp_printf(" %Zd =⌈%Zd / 3⌉ \n",q,n);
          std ::cout<<"Q = "<<Q<<"\n";
          R= mpz_cdiv_r_ui(r,n,3);
          gmp_printf(" %Zd = %Zd mod 3\n",r,n);
           std ::cout<<"R = "<<R<<"\n";
           R=mpz_cdiv_qr_ui(q,r,n,3);
           gmp_printf("q = %Zd ,r = %Zd ,n = %Zd ,d = 3 \n",q,r,n) ;  // Computes both quotient and remainder.
            std ::cout<<"R = "<<R<<"\n";


//Exact Division Functions 
       std:: cout<<"Exact Division Functions\n";                  
       mpz_divexact(q,n,d);                                       //Computes q=n/d only if d divides n exactly.
       gmp_printf(" %Zd = %Zd / %Zd \n",q,n,d);

        mpz_divexact_ui(q,n,2);
        gmp_printf(" %Zd = %Zd / 2 \n",q,n);                      //Same as above, using an unsigned long divisor.

 //Modulo Functions
       std:: cout<<"Modulo Functions\n"; 
       mpz_mod(r,n,d);
       gmp_printf(" %Zd = %Zd mod %Zd\n",r,n,d);                      //Computes r = n mod d, ensuring the result is non-negative.

       mpz_mod_ui(r,n,30000000000000000);                                             //same as above, using an unsigned long divisor
       gmp_printf(" %Zd = %Zd mod 3\n",r,n);

//Divisibility Check Functions  
       std:: cout<<"Divisibility Check Functions \n";   
       
       if(mpz_divisible_p(n,d))                                   //Returns non-zero if n is divisible by d
          {
            gmp_printf("%Zd is divisible by %Zd \n",n,d);
          } 
        else  gmp_printf("%Zd is not divisible by %Zd \n",n,d);  

        if(mpz_divisible_ui_p(n,3))                                   //Checks divisibility with an unsigned long divisor
          {
            gmp_printf("%Zd is divisible by 3 \n",n);
          } 
        else  gmp_printf("%Zd is not divisible by 3 \n",n);


// Congruence Check Functions
       std:: cout<<"Congruence  Check Functions \n";                              //Checks if n is congruent to c mod d.

       if(mpz_congruent_p(n,c,d))
         {
            gmp_printf("%Zd is congruent to %Zd modulo %Zd \n",n,c,d);
         }
       else {
                gmp_printf("%Zd is not  congruent to %Zd modulo %Zd \n",n,c,d);
           }  


        if(mpz_congruent_ui_p(n,2,2))                                          //Congruence check with unsigned long integers
         {
            gmp_printf("%Zd is congruent to 2 modulo 2 \n",n);
         }
       else {
                gmp_printf("%Zd is not  congruent to 2 modulo 2 \n",n);
           }  
   
       



             





       

           


        mpz_clears(n,d,r,c,q,NULL);


 }
