#include <gmp.h>
#include <stdlib.h>

#include "structures.h"

void sieve(dyn_array_small* sieve,
    const unsigned long length,
    const unsigned long L,
    const dyn_array_classic* primes,
    const unsigned long * logs,
    const dyn_array* roots,
    const mpz_t n,
    const mpz_t poly_a,
    const mpz_t poly_b,
    const mpz_t poly_c,
    const dyn_array* inverse_a,
    const dyn_array_classic* way_to_root,
    const dyn_array_classic* locations,
    const unsigned long skipped,
    const unsigned long prime_start,
    dyn_array* smooth, 
    const unsigned long smooth_bound)
{
    mpz_t z,tmp,tmpb;
    mpz_init(z);
    mpz_init(tmp);
    mpz_init(tmpb);
    int valid;

    unsigned long root, increment, p;
    unsigned long r;
    unsigned long b;

    mpz_mul_2exp(tmpb, poly_b, 1);

    mpz_neg(tmp, poly_c);
    mpz_invert(tmpb,tmpb,n);

    mpz_mul(tmp,tmp,tmpb);

    const unsigned short * restrict sieve_array_start = sieve->start;
    const unsigned long * restrict way_to_root_start = way_to_root->start;
    const unsigned long * restrict locations_start = locations->start;
    const unsigned long * restrict log_start = logs;
    const unsigned long * restrict primes_start = primes->start;
    const mpz_t * restrict roots_start = roots->start;
    const mpz_t * restrict inversa_a_start = inverse_a->start;

    for (unsigned long i = prime_start ; i < primes->len ; i++)
    {
        p = *(primes_start+i);
        increment = log_start[i];
        valid = 1;
        for (unsigned long j = 0 ; j < locations->len ; j++)
        {
            if (*(locations_start+j) == i)
            {
                valid = 0;
                break;
            }
        }
        switch (valid)
        {
            case 1:
                mpz_sub(z,*(roots_start+i),poly_b);
                root = mpz_mod_ui(z,z,p);
                root = (root * mpz_get_ui(*(inversa_a_start+i)))%p;
                root = (root+L)%p;

                for (unsigned short* ptr = sieve_array_start+root ; ptr < sieve_array_start+length ; ptr += p)
                {
                    *ptr += increment;
                }
                root = (root + p - *(way_to_root_start+i))%p;
                
                for (unsigned short* ptr = sieve_array_start+root ; ptr < sieve_array_start+length ; ptr += p)
                {
                    *ptr += increment;
                }
                break;
            case 0:
                root = mpz_mod_ui(z,tmp,p);
                root = (root+L)%p;
                
                for (unsigned short* ptr = sieve_array_start+root; ptr < sieve_array_start+length ; ptr += p)
                {
                    *ptr += increment;
                }
                break;
        }
    }

    mpz_t k;
    mpz_init(k);

    signed long tmp_array[length];
    size_t tmp_count = 0;

    for (unsigned long x = L-1 ; --x ;)
    {
        if (*(sieve_array_start+x)+skipped >= smooth_bound)
        {
            tmp_array[tmp_count++] = (long)x - (long)L;
        }
        if (*(sieve_array_start+length-1-x)+skipped >= smooth_bound)
        {
            tmp_array[tmp_count++] = (long)L - (long)x;
        }
    }

    if (*(sieve_array_start+L)+skipped >= smooth_bound)
    {
        tmp_array[tmp_count++] = 0;
    }

    for (size_t i = 0 ; i < tmp_count ; i++)
    {
        append_only_si(smooth, tmp_array[i]);
    }

    mpz_clears(tmp,tmpb,z,k,NULL);
}