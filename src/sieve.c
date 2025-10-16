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
    signed long *tmp_array,
    const unsigned long skipped,
    const unsigned long prime_start,
    dyn_array* smooth, 
    const unsigned long smooth_bound)
{

    mpz_t z, tmp, inv_2b_mod_n;
    mpz_inits(z, tmp, inv_2b_mod_n, NULL);
    bool excluded;

    unsigned long root, log_weight, p;
    unsigned long r;
    unsigned long b;

    mpz_mul_2exp(inv_2b_mod_n, poly_b, 1);
    mpz_invert(inv_2b_mod_n, inv_2b_mod_n, n);

    mpz_neg(tmp, poly_c);

    mpz_mul(tmp, tmp, inv_2b_mod_n);

    unsigned short * restrict sieve_array_start = sieve->start;
    const unsigned long * restrict root_offsets = way_to_root->start;
    const unsigned long * restrict excluded_primes = locations->start;
    const unsigned long * restrict log_start = logs;
    const unsigned long * restrict primes_start = primes->start;
    const mpz_t * restrict roots_start = roots->start;
    const mpz_t * restrict inverse_a_start = inverse_a->start;

    for (size_t i = prime_start ; i < primes->len ; i++)
    {
        p = primes_start[i];
        log_weight = log_start[i];
        excluded = false;
        for (size_t j = 0 ; j < locations->len ; j++)
        {
            if (excluded_primes[j] == i)
            {
                excluded = true;
                break;
            }
        } 

        if (!excluded)
        {
            mpz_sub(z, roots_start[i], poly_b);
            root = mpz_fdiv_ui(z, p);
            root = (root * mpz_get_ui(inverse_a_start[i]))%p;

            root = (root+L)%p;
            
            #pragma omp simd
            for (unsigned short* ptr = sieve_array_start+root ; ptr < sieve_array_start+length ; ptr += p)
            {
                *ptr += log_weight;
            }

            root = (root + p - root_offsets[i])%p;
            
            #pragma omp simd
            for (unsigned short* ptr = sieve_array_start+root ; ptr < sieve_array_start+length ; ptr += p)
            {
                *ptr += log_weight;
            }

        }
        else
        {
            root = mpz_fdiv_ui(z, p);
            root = (root+L)%p;
            
            #pragma omp simd
            for (unsigned short* ptr = sieve_array_start+root; ptr < sieve_array_start+length ; ptr += p)
            {
                *ptr += log_weight;
            }

        }
    }


    mpz_t k;
    mpz_init(k);

    size_t tmp_count = 0;

    for (size_t x = L-1 ; --x ;)
    {
        if (sieve_array_start[x]+skipped >= smooth_bound)
        {
            tmp_array[tmp_count++] = (long)x - (long)L;
        }
        if (sieve_array_start[length-1-x]+skipped >= smooth_bound)
        {
            tmp_array[tmp_count++] = (long)L - (long)x;
        }
    }

    if (sieve_array_start[L]+skipped >= smooth_bound)
    {
        tmp_array[tmp_count++] = 0;
    }

    for (size_t i = 0 ; i < tmp_count ; i++)
    {
        append_only_si(smooth, tmp_array[i]);
    }

    mpz_clears(tmp, inv_2b_mod_n, z, k, NULL);
}