#include <gmp.h>
#include <stdlib.h>

#include "structures.h"

void sieve(dyn_array_small* sieve,unsigned long length, unsigned long L, dyn_array_classic* primes, unsigned long* logs, dyn_array* roots, mpz_t n, mpz_t poly_a, mpz_t poly_b, mpz_t poly_c, dyn_array* inverse_a, dyn_array_classic* way_to_root, dyn_array_classic* locations, unsigned long skipped, unsigned long prime_start, dyn_array* smooth, unsigned long smooth_bound)
{
    mpz_t z,tmp,tmpb;
    mpz_init(z);
    mpz_init(tmp);
    mpz_init(tmpb);
    int valid;
    unsigned long root;
    mpz_mul_ui(tmpb,poly_b,2);
    unsigned long tmplolol, p;
    for (unsigned long i = prime_start ; i < primes->len ; i++)
    {
        p = *(primes->start+i);
        tmplolol = logs[i];
        valid = 1;
        for (unsigned long j = 0 ; j < locations->len && valid; j++)
        {
            if (*(locations->start+j) == i) valid = 0;
        }
        switch (valid)
        {
            case 1:
                mpz_sub(z,*(roots->start+i),poly_b);
                mpz_mul(z,z,*(inverse_a->start+i));
                mpz_mod_ui(z,z,p);
                root = (mpz_get_ui(z)+L)%p;
                for (unsigned short* ptr = sieve->start+root ; ptr < sieve->start+length ; ptr += p)
                {
                    *ptr += tmplolol;
                }
                root = (root-*(way_to_root->start+i))%p;
                for (unsigned short* ptr = sieve->start+root ; ptr < sieve->start+length ; ptr += p)
                {
                    *ptr += tmplolol;
                }
                break;
            case 0:
                mpz_neg(z,poly_c);
                mpz_mul_ui(tmp,poly_b,2);
                mpz_invert(tmp,tmp,n);
                mpz_mul(z,z,tmp);
                mpz_add_ui(z,z,L);
                mpz_mod_ui(z,z,p);
                for (unsigned short* ptr = sieve->start+mpz_get_ui(z) ; ptr < sieve->start+length ; ptr += p)
                {
                    *ptr += tmplolol;
                }
                break;
        }
    }
    mpz_t k;
    mpz_init(k);
    for (unsigned long x = L-1 ; --x ;)
    {
        mpz_set_ui(k,x);
        mpz_sub_ui(k,k,L);

        if (*(sieve->start+x)+skipped >= smooth_bound) append_only(smooth,k);
        if (*(sieve->start+length-1-x)+skipped >= smooth_bound)
        {
            mpz_neg(k,k);
            append_only(smooth,k);
        }
    }
    if (*(sieve->start+L)+skipped >= smooth_bound)
    {
        mpz_set_ui(k,0);
        append_only(smooth,k);
    }
    mpz_clears(tmp,tmpb,z,k,NULL);
}