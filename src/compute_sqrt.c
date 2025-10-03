#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>

#include "structures.h"

void build_sqrt(mpz_t n,unsigned long dim, bool vector[dim], dyn_array_classic primes, mpz_t x, mpz_t y, dyn_array relations, dyn_array smooth)
{
    mpz_t tmp,tmp2;
    mpz_init(tmp);
    mpz_init(tmp2);
    mpz_set_ui(y,1);
    mpz_set_ui(x,1);
    unsigned long vector_y[primes.len];
    for (unsigned long i = 0 ; i < primes.len ; i++) vector_y[i] = 0;
    for (unsigned long e = 0 ; e < dim ; e++)
    {
        if (vector[e])
        {
            mpz_mul(x,x,*(smooth.start+e));
            mpz_mod(x,x,n);
            for (unsigned long i = 0 ; i < primes.len ; i++)
            {
                mpz_set_ui(tmp,*(primes.start+i));
                while(mpz_divisible_p(*(relations.start+e),tmp))
                {
                    vector_y[i]++;
                    mpz_mul_ui(tmp,tmp,*(primes.start+i));
                }
            }
        }
    }
    for (unsigned long k = 0 ; k < primes.len ; k++)
    {
        mpz_set_ui(tmp2,vector_y[k]);
        mpz_div_ui(tmp2,tmp2,2);
        mpz_set_ui(tmp,*(primes.start+k));
        mpz_powm(tmp,tmp,tmp2,n);
        mpz_mul(y,y,tmp);
        mpz_mod(y,y,n);
    }
    mpz_clears(tmp,tmp2,NULL);
}