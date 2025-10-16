#include <gmp.h>
#include <stdlib.h>
#include <time.h>

#include "structures.h"
#include "utils.h"

void pollard_rho(const mpz_t m, mpz_t p1, mpz_t p2, gmp_randstate_t state)
{
    mpz_t a, b, d, e, x, y, X, Y, tmp;
    mpz_inits(a, b, d, e, x, y, X, Y, tmp, NULL);

    while (1)
    {
        mpz_sub_ui(tmp, m, 4);
        mpz_urandomm(a, state, tmp);
        mpz_add_ui(a, a, 1);

        mpz_sub_ui(tmp, m, 1);
        mpz_urandomm(b, state, tmp);

        mpz_set(x, b);
        mpz_set(y, b);
        mpz_set_ui(d, 1);

        while (mpz_cmp_ui(d, 1) == 0)
        {
            mpz_set_ui(e, 1);
            mpz_set(X, x);
            mpz_set(Y, y);

            for (unsigned int k = 0 ; k < 100 ; k++)
            {
                mpz_mul(tmp, x, x);
                mpz_add(tmp, tmp, a);
                mpz_mod(x, tmp, m); // x = (x*x+a)%m

                mpz_mul(tmp, y, y);
                mpz_add(tmp, tmp, a);
                mpz_mod(y, tmp, m); // y = (y*y+a)%m

                mpz_mul(tmp, y, y);
                mpz_add(tmp, tmp, a);
                mpz_mod(y, tmp, m); // y = (y*y+a)%m

                mpz_sub(tmp, x, y);
                mpz_mod(tmp, tmp, m);
                mpz_mul(b, e, tmp);
                mpz_mod(e, b, m); // e = e*(x-y)%m
            }

            mpz_gcd(d, e, m);
        }

        if (mpz_cmp(d, m) == 0)
        {
            mpz_set(x, X);
            mpz_set(y, Y);
            mpz_set_ui(d, 1);

            while (mpz_cmp_ui(d, 1) == 0)
            {
                mpz_mul(tmp, x, x);
                mpz_add(tmp, tmp, a);
                mpz_mod(x, tmp, m); // x = (x*x+a)%m

                mpz_mul(tmp, y, y);
                mpz_add(tmp, tmp, a);
                mpz_mod(y, tmp, m); // y = (y*y+a)%m

                mpz_mul(tmp, y, y);
                mpz_add(tmp, tmp, a);
                mpz_mod(y, tmp, m); // y = (y*y+a)%m

                mpz_sub(tmp, x, y);
                mpz_mod(tmp, tmp, m);
                mpz_gcd(d, tmp, m);
            }
        }

        if (mpz_cmp(d, m)) // If d != m
        {
            mpz_divexact(tmp, m, d);
            if (mpz_cmp(tmp, d) < 0)
            {
                mpz_set(p1, tmp);
                mpz_set(p2, d);
            }
            else
            {
                mpz_set(p1, d);
                mpz_set(p2, tmp);
            }
            mpz_clears(a, b, d, e, x, y, X, Y, tmp, NULL);
            return;
        }
    }
}

unsigned long build_product_tree(dyn_array* reported, dyn_array* tmp_array, mpz_t prod_primes, mpz_t prod_primes_p1, mpz_t tmp)
{
    unsigned long tmp_long = reported->len;
    unsigned long index1 = tmp_long, index2 = tmp_long;

    mpz_add_ui(prod_primes_p1, prod_primes, 1);

    tmp_long >>= 1;
    for (unsigned long i = 0 ; i < tmp_long ; i++)
    {
        mpz_mul(tmp, reported->start[i<<1], reported->start[(i<<1)+1]);
        if (mpz_cmp(tmp, prod_primes) < 1) mpz_set(tmp_array->start[index1+i], tmp);
        else mpz_set(tmp_array->start[index1+i], prod_primes_p1);
    }

    index2 += tmp_long;
    
    while (tmp_long > 1)
    {
        for (unsigned long i = 0 ; i < tmp_long ; i++)
        {
            if (!mpz_cmp(tmp_array->start[index1+(i<<1)], prod_primes_p1) || !mpz_cmp(tmp_array->start[index1+(i<<1)+1], prod_primes_p1))
            {
                mpz_set(tmp_array->start[index2+i], prod_primes_p1);
            }
            else
            {
                mpz_mul(tmp, tmp_array->start[index1+(i<<1)], tmp_array->start[index1+(i<<1)+1]);

                if (mpz_cmp(tmp, prod_primes) < 1) mpz_set(tmp_array->start[index2+i], tmp);
                else mpz_set(tmp_array->start[index2+i], prod_primes_p1);
            }
        }
        index1 = index2;
        tmp_long >>= 1;
        index2 += tmp_long;
    }
    return tmp_long;
}

void build_remainder_tree(dyn_array* reported, dyn_array* tmp_array, mpz_t prod_primes, mpz_t prod_primes_p1, unsigned long tmp_long)
{
    unsigned long index1 = 0;
    unsigned long index2 = 1;
    unsigned long stored_len = tmp_array->len;
    mpz_mod(tmp_array->start[stored_len-1], prod_primes, tmp_array->start[stored_len-1]);

    bool flag_overflow;

    while (tmp_long < reported->len>>1)
    {
        for (unsigned long i = 0 ; i < tmp_long ; i++)
        {
            flag_overflow = mpz_cmp(tmp_array->start[stored_len-1-index2-(i<<1)], prod_primes_p1)
            if (mpz_cmp(tmp_array->start[stored_len-1-index2-(i<<1)], prod_primes_p1))
            {
                mpz_mod(tmp_array->start[stored_len-1-index2-(i<<1)], tmp_array->start[stored_len-1-index1-i], tmp_array->start[stored_len-1-index2-(i<<1)]);
            }
            else mpz_set(tmp_array->start[stored_len-1-index2-(i<<1)], prod_primes);

            if (flag_overflow)
            {
                mpz_mod(tmp_array->start[stored_len-2-index2-(i<<1)], tmp_array->start[stored_len-1-index1-i], tmp_array->start[stored_len-2-index2-(i<<1)]);
            }
            else mpz_set(tmp_array->start[stored_len-2-index2-(i<<1)], prod_primes);
        }
        index1 = index2;
        tmp_long <<= 1;
        index2 += tmp_long;
    }


    for (unsigned long i = 0 ; i < tmp_long ; i++)
    {
        mpz_mod(tmp_array->start[i<<1], tmp_array->start[reported->len+i], reported->start[i<<1]);
        mpz_mod(tmp_array->start[(i<<1)+1], tmp_array->start[reported->len+i], reported->start[(i<<1)+1]);
    }
}

void batch_smooth(dyn_array* reported, dyn_array* tmp_array, PartialRelation *tmp_array2, mpz_t prod_primes, mpz_t prod_primes_p1, mpz_t limit, mpz_t limit_2, unsigned long prime, gmp_randstate_t state)
{
    mpz_t tmp, tmp2;
    mpz_inits(tmp, tmp2, NULL);
    mpz_set_ui(tmp, prime);

    for (size_t i = 0 ; i < reported->len ; i++)
    {
        mpz_remove(reported->start[i], reported->start[i], tmp);
    } // Remove occurences of 2

    mpz_set_ui(tmp, 2);

    unsigned long squaring_steps = 0;
    
    for (mpz_t* ptr = reported->start ; ptr < reported->start+reported->len ; ptr++)
    {
        while (mpz_cmp(tmp, *ptr) < 0)
        {
            mpz_mul(tmp2, tmp, tmp);
            mpz_set(tmp, tmp2);
            squaring_steps++;
        }
    }
    // computed log(log(max[reported]))

    // computing prod tree
    unsigned long tmp_long = build_product_tree(reported, tmp_array, prod_primes, prod_primes_p1, tmp);

    //computing remainder tree
    build_remainder_tree(reported, tmp_array, prod_primes, prod_primes_p1, tmp_long);

    for (unsigned long i = 0 ; i < reported->len ; i++)
    {
        tmp_long = 0;
        mpz_set(tmp, tmp_array->start[i]);
        while (mpz_cmp_ui(tmp, 0) && tmp_long < squaring_steps)
        {
            mpz_mul(tmp, tmp, tmp);
            mpz_mod(tmp, tmp, reported->start[i]);
            tmp_long++; 
        }

        if (!mpz_cmp_ui(tmp, 0))
        {
            mpz_set_ui((tmp_array2+i)->small_p, 1);
            mpz_set_ui((tmp_array2+i)->big_p, 1);
        }
        else
        {
            mpz_gcd(tmp2, reported->start[i], tmp);
            mpz_divexact(tmp2,reported->start[i], tmp2);

            if (mpz_cmp(tmp2, limit) < 1)
            {
                mpz_set_ui((tmp_array2+i)->small_p, 1);
                mpz_set((tmp_array2+i)->big_p, tmp2);
            }
            else if (mpz_cmp(tmp2, limit_2) < 1 && !fermat_primality(tmp2))
            {
                pollard_rho(tmp2, tmp, tmp2, state);
                mpz_set((tmp_array2+i)->small_p, tmp);
                mpz_set((tmp_array2+i)->big_p, tmp2);
            }
            else mpz_set_ui((tmp_array2+i)->small_p, 0);
        }
    }
    mpz_clears(tmp, tmp2, NULL);
}

void naive_smooth(dyn_array* reported, PartialRelation *tmp_array, dyn_array_classic primes, mpz_t limit, mpz_t limit_2, gmp_randstate_t state)
{
    mpz_t tmp, tmp2, tmp3;
    mpz_init(tmp);
    mpz_inits(tmp2, tmp3, NULL);
    unsigned long i = 0;

    for (mpz_t* ptr = reported->start ; ptr < reported->start+reported->len ; ptr++)
    {
        mpz_set(tmp, *ptr);
        for (unsigned long i = 0 ; i < primes.len ; i++)
        {
            mpz_mod_ui(tmp2, tmp, primes.start[i]);

            while (!mpz_cmp_ui(tmp2, 0))
            {
                mpz_div_ui(tmp, tmp, primes.start[i]);
                mpz_mod_ui(tmp2, tmp, primes.start[i]);
            }
        }

        if (!mpz_cmp_ui(tmp, 1))
        {
            mpz_set_ui((tmp_array+i)->small_p, 1);
            mpz_set_ui((tmp_array+i)->big_p, 1);
        }
        else if (mpz_cmp(tmp, limit) < 1)
        {
            mpz_set_ui((tmp_array+i)->small_p, 1);
            mpz_set((tmp_array+i)->big_p, tmp);
        }
        else if (mpz_cmp(tmp, limit_2) < 1 && !fermat_primality(tmp))
        {
            pollard_rho(tmp, tmp2, tmp3, state);
            mpz_set((tmp_array+i)->small_p, tmp2);
            mpz_set((tmp_array+i)->big_p, tmp3);
        }
        else mpz_set_ui((tmp_array+i)->small_p, 0);

        i++;
    }
    mpz_clears(tmp, tmp2, NULL);
}