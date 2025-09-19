#include <gmp.h>
#include <stdlib.h>

#include "structures.h"

void batch_smooth(dyn_array* reported, dyn_array* large_primes, dyn_array_classic* smooth, dyn_array* tmp_array, mpz_t prod_primes, mpz_t limit, unsigned long prime)
{
    for (mpz_t* ptr = reported->start ; ptr < reported->start+reported->len ; ptr++)
    {
        while (mpz_divisible_ui_p(*ptr,prime)) mpz_divexact_ui(*ptr,*ptr,prime);
    }
    mpz_t tmp;
    mpz_init_set_ui(tmp,2);
    unsigned long e = 0;
    for (mpz_t* ptr = reported->start ; ptr < reported->start+reported->len ; ptr++)
    {
        while (mpz_cmp(tmp,*ptr) < 0)
        {
            mpz_mul(tmp,tmp,tmp);
            e++;
        }
    }
    // computed log(log(max[reported]))

    // computing prod tree
    unsigned long tmp_long = reported->len;
    unsigned long index1 = tmp_long, index2 = tmp_long;
    mpz_t boring;
    mpz_init(boring);
    mpz_add_ui(boring,prod_primes,1);
    tmp_long >>= 1;
    for (unsigned long i = 0 ; i < tmp_long ; i++)
    {
        mpz_mul(tmp,*(reported->start+2*i),*(reported->start+2*i+1));
        if (mpz_cmp(tmp,prod_primes) < 1) mpz_set(*(tmp_array->start+index1+i),tmp);
        else mpz_set(*(tmp_array->start+index1+i),boring);
    }
    index2 += tmp_long;
    while (tmp_long > 1)
    {
        for (unsigned long i = 0 ; i < tmp_long ; i++)
        {
            if (!mpz_cmp(*(tmp_array->start+index1+2*i),boring) || !mpz_cmp(*(tmp_array->start+index1+2*i+1),boring)) mpz_set(*(tmp_array->start+index2+i),boring);
            else
            {
                mpz_mul(tmp,*(tmp_array->start+index1+2*i),*(tmp_array->start+index1+2*i+1));
                if (mpz_cmp(tmp,prod_primes) < 1) mpz_set(*(tmp_array->start+index2+i),tmp);
                else mpz_set(*(tmp_array->start+index2+i),boring);
            }
        }
        index1 = index2;
        tmp_long >>= 1;
        index2 += tmp_long;
    }
    // prod tree created

    //computing remainder tree
    index1 = 0;
    index2 = 1;
    unsigned long stored_len = tmp_array->len;
    mpz_mod(*(tmp_array->start+stored_len-1),prod_primes,*(tmp_array->start+stored_len-1));
    while (tmp_long < reported->len>>1)
    {
        for (unsigned long i = 0 ; i < tmp_long ; i++)
        {
            if (mpz_cmp(*(tmp_array->start+stored_len-1-index2-2*i),boring) != 0) mpz_mod(*(tmp_array->start+stored_len-1-index2-2*i),*(tmp_array->start+stored_len-1-index1-i),*(tmp_array->start+stored_len-1-index2-2*i));
            else mpz_set(*(tmp_array->start+stored_len-1-index2-2*i),prod_primes);
            if (mpz_cmp(*(tmp_array->start+stored_len-1-index2-2*i),boring) != 0) mpz_mod(*(tmp_array->start+stored_len-2-index2-2*i),*(tmp_array->start+stored_len-1-index1-i),*(tmp_array->start+stored_len-2-index2-2*i));
            else mpz_set(*(tmp_array->start+stored_len-2-index2-2*i),prod_primes);
        }
        index1 = index2;
        tmp_long <<= 1;
        index2 += tmp_long;
    }
    for (unsigned long i = 0 ; i < tmp_long ; i++)
    {
        mpz_mod(*(tmp_array->start+2*i),*(tmp_array->start+reported->len+i),*(reported->start+2*i));
        mpz_mod(*(tmp_array->start+2*i+1),*(tmp_array->start+reported->len+i),*(reported->start+2*i+1));
    }

    for (unsigned long i = 0 ; i < reported->len ; i++)
    {
        tmp_long = 0;
        mpz_set(tmp,*(tmp_array->start+i));
        while (mpz_cmp_ui(tmp,0) && tmp_long < e)
        {
            mpz_mul(tmp,tmp,tmp);
            mpz_mod(tmp,tmp,*(reported->start+i));
            tmp_long++;
        }
        if (!mpz_cmp_ui(tmp,0)) *(smooth->start+i) = 1;
        else
        {
            mpz_gcd(boring,*(reported->start+i),tmp);
            mpz_divexact(boring,*(reported->start+i),boring);
            if (mpz_cmp(boring,limit) < 1)
            {
                *(smooth->start+i) = 1;
                mpz_set(*(large_primes->start+i),boring);
            }
        }
    }
    mpz_clears(tmp,boring,NULL);
}

void naive_smooth(dyn_array* reported, dyn_array* large_primes, dyn_array_classic* smooth, dyn_array_classic primes, mpz_t limit)
{
    mpz_t tmp, tmp2;
    mpz_init(tmp);
    mpz_init(tmp2);
    unsigned long i = 0;

    for (mpz_t* ptr = reported->start ; ptr < reported->start+reported->len ; ptr++)
    {
        mpz_set(tmp, *ptr);
        for (unsigned long i = 0 ; i < primes.len ; i++)
        {
            mpz_mod_ui(tmp2, tmp, *(primes.start+i));

            while (!mpz_cmp_ui(tmp2, 0))
            {
                mpz_div_ui(tmp, tmp, *(primes.start+i));
                mpz_mod_ui(tmp2, tmp, *(primes.start+i));
            }
        }

        if (!mpz_cmp_ui(tmp, 1))
        {
            *(smooth->start+i) = 1;
        }
        else if (mpz_cmp(tmp, limit) < 1)
        {
            *(smooth->start+i) = 1;
            mpz_set(*(large_primes->start+i), tmp);
        }

        i++;
    }
    mpz_clears(tmp, tmp2, NULL);
}