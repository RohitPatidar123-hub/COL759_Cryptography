#include<iostream>
#include<gmp.h>
 
 int main()
 {
       mpz_t a,b,c,add,sub,mul,mul1,mul2,addmul1,addmul2,submul1,submul2,result,neg,abs;    //declare varible of int type
       mpz_init_set_ui(a,2000000000000);
       mpz_init_set_ui(b,3000000000000);
       mpz_init_set_ui(c,8000000000000);
       mpz_init(result);
       mpz_init(add);
       mpz_init(neg);
       mpz_init(abs);
       mpz_init(sub);
       mpz_init(mul);
       mpz_init(mul1);
       mpz_init(mul2);
       mpz_init_set_ui(addmul1,500000);
       mpz_init_set_ui(addmul2,500000);
       mpz_init_set_ui(submul1,50000000000000);
       mpz_init_set_ui(submul2,50000000000000);




       
       
       mpz_add(add,a,b);
       mpz_sub(sub,a,b);
       mpz_mul(mul,a,b);
       mpz_mul_si(mul1,a,-10);                     // Computes rop = op1 * op2 where op2 is a signed long.
       mpz_mul_ui(mul2,a,10);                      //Computes rop = op1 * op2 where op2 is an unsigned long.
       mpz_mul_2exp(result,c,10);                   //result = c*2^10;






       gmp_printf("%Zd  + %Zd = %Zd\n",a,b,add);
       gmp_printf("%Zd  - %Zd = %Zd\n",a,b,sub);
       gmp_printf("%Zd  * %Zd = %Zd\n",a,b,mul);
       gmp_printf("%Zd  * -10 = %Zd\n",a,mul1);
       gmp_printf("%Zd  * 10= %Zd\n",a,mul2);
       gmp_printf("%Zd ",addmul1);
       mpz_addmul(addmul1,a,b) ;                    //compute addmul1=addmul1 + a*b
       gmp_printf(" + %Zd  * %Zd = %Zd\n",a,b,addmul1);

       gmp_printf("%Zd ",addmul2);
       mpz_addmul_ui(addmul2,a,10) ;                    //compute addmul2=addmul2 + a*10
       gmp_printf(" + %Zd  * 10 = %Zd\n",a,addmul2);

       gmp_printf("%Zd ",submul1);
       mpz_submul(submul1,a,b) ;                    //compute submul1=submul1 - a*b
       gmp_printf(" - %Zd  * %Zd = %Zd\n",a,b,submul1);

       gmp_printf("%Zd ",submul2);
       mpz_submul_ui(submul2,a,10) ;                    //compute submul2=submul2 + a*10
       gmp_printf(" - %Zd  * 10 = %Zd\n",a,submul2);
  
       gmp_printf("%Zd * 2^10 = %Zd\n",c,result );

       mpz_neg(neg,submul1);
       mpz_abs(abs,submul1);

       gmp_printf("negation of %Zd = %Zd \n",submul1,neg);
       gmp_printf("absolutevalue  of %Zd = %Zd \n",submul1,abs);









       mpz_clears(a,b,add,sub,mul,mul1,mul2,addmul1,addmul2,submul1,submul2,result,neg,abs,NULL);
 }
