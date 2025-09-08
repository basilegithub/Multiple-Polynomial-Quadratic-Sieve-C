#include <gmp.h>

#include "structures.h"
#include "generate_primes.h"

void smoothB(mpz_t B, dyn_array_classic* primes)
{
    mpz_t m;
    mpz_init(m);
    unsigned long i = 0;

    mpz_sqrt(m,B);
    mpz_add_ui(m,m,1);
    dyn_array_classic array;
    init_len_classic(&array,mpz_get_ui(B));
    while (mpz_cmp_ui(B,i) >= 0)
    {
        *(array.start+i) = 1;
        i++;
    }
    i = 1;
    while (mpz_cmp_ui(B,i) >= 0)
    {
        *(array.start+i) = 0;
        i += 2;
    }
    unsigned long index = 2;
    append_classic(primes,index);
    index = 3;
    unsigned long tmp = 2,square_index;
    while (index < mpz_get_ui(m))
    {
        if (*(array.start+tmp) == 1)
        {
            square_index = index*index;
            append_classic(primes,index);
            i = square_index;
            tmp = i-1;
            while(mpz_cmp_ui(B,i) >= 0)
            {
                *(array.start+tmp) = 0;
                i += index;
                tmp = i-1;
            }
        }
        index += 2;
        tmp = index-1;
    }
    i = mpz_get_ui(m);
    while (mpz_cmp_ui(B,i) >= 0)
    {
        tmp = i-1;
        if (*(array.start+tmp) == 1)
        {
            append_classic(primes,i);
        }
        i++;
    }
    mpz_clear(m);
}