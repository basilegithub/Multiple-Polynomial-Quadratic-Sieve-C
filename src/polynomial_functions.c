#include <gmp.h>
#include <stdlib.h>

#include "structures.h"
#include "utils.h"

void create_polynomial(mpz_t a, dyn_array* sol_needed, dyn_array* second_part, dyn_array_classic* locations, mpz_t n, dyn_array_classic* primes, dyn_array* roots, unsigned long bounds[4], mpf_t target, mpf_t ln2, mpf_t ln10, dyn_array* inverse_a, dyn_array_classic* way_to_root, mpf_t best_bound, dyn_array_classic* tmp_where, mpf_t e, unsigned long mult)
{
    mpf_t bound;
    mpf_init(bound);
    mpz_set_ui(a,1);
    sol_needed->len = 0;
    second_part->len = 0;
    inverse_a->len = 0;
    way_to_root->len = 0;
    mpf_set_ui(best_bound,10);
    dyn_array selected;
    init(&selected);
    mpf_t tmpf,tmpf2;
    mpf_init(tmpf);
    mpf_init(tmpf2);
    unsigned long location;
    mpz_t tmp;
    mpz_init(tmp);
    unsigned long tmp_long;
    int valid;
    mpz_t tmp2;
    mpz_init(tmp2);
    for (unsigned long i = 0 ; i < 6 ; i++)
    {
        tmp_where->len = 0;
        mpf_set(bound,target);
        while (mpf_cmp_ui(bound,1) > 0)
        {
            if (mpf_cmp_ui(bound,7) > -1) location = rand()%(bounds[3]-bounds[2])+bounds[2];
            else if (mpf_cmp_ui(bound,4) > -1) location = rand()%(bounds[2]-bounds[1])+bounds[1];
            else location = rand()%(bounds[1]-bounds[0])+bounds[0];

            valid = 1;
            if (*(primes->start+location) == mult) valid = 0;
            for (unsigned long j = 0 ; j < tmp_where->len && valid; j++)
            {
                if (*(tmp_where->start+j) == location) valid = 0;
            }
            if (valid == 1)
            {
                append_classic(tmp_where,location);
                mpz_set_ui(tmp,*(primes->start+location));
                mpf_set_z(tmpf,tmp);
                my_log10(tmpf,tmpf,ln2,ln10,e);
                mpf_sub(bound,bound,tmpf);
            }
        }
        mpf_abs(tmpf,bound);
        if (mpf_cmp_ui(best_bound,10) == 0 || mpf_cmp(tmpf,tmpf2) < 0)
        {
            locations->len = 0;
            for (unsigned long j = 0 ; j < tmp_where->len ; j++) append_classic(locations,*(tmp_where->start+j));
            mpf_set(best_bound,bound);
            mpf_abs(tmpf2,best_bound);
        }
    }
    for (unsigned long i = 0 ; i < locations->len ; i++)
    {
        tmp_long = *(locations->start+i);
        mpz_set_ui(tmp2,*(primes->start+tmp_long));
        mpz_mul(a,a,tmp2);
        append(&selected,tmp2);
        append_only(sol_needed,*(roots->start+tmp_long));
    }
    for (unsigned long i = 0 ; i < selected.len ; i++)
    {
        mpz_divexact(tmp,a,*(selected.start+i));
        mpz_invert(tmp2,tmp,*(selected.start+i));
        mpz_mul(tmp2,tmp2,tmp);
        mpz_mod(tmp2,tmp2,a);
        append_only(second_part,tmp2);
    }
    for (unsigned long i = 0 ; i < primes->len ; i++)
    {
        valid = 1;
        for (unsigned long j = 0 ; j < locations->len ; j++)
        {
            if (*(locations->start+j) == i) valid = 0;
        }
        if (valid == 1)
        {
            mpz_set_ui(tmp2,*(primes->start+i));
            mpz_invert(tmp,a,tmp2);
            append_only(inverse_a,tmp);
            mpz_mul(tmp,tmp,*(roots->start+i));
            mpz_mul_2exp(tmp, tmp, 1);
            mpz_mod_ui(tmp,tmp,*(primes->start+i));
            append_classic(way_to_root,mpz_get_ui(tmp));
        }
        else
        {
            mpz_set_ui(tmp,0);
            append_only(inverse_a,tmp);
            append_classic(way_to_root,0);
        }
    }
    for (unsigned long i = 0 ; i < selected.len ; i++) mpz_clear(*(selected.start+i));
    free(selected.start);
    selected.start = NULL;
    mpz_clears(tmp,tmp2,NULL);
    mpf_clears(bound,tmpf,tmpf2,NULL);
}

void CRT(mpz_t res, dyn_array* moduli, mpz_t a, dyn_array* second_part)
{
    mpz_t total;
    mpz_init_set_ui(total, 0);

    for (unsigned long i = 0 ; i < moduli->len ; i++)
    {
        mpz_addmul(total, moduli->start[i], second_part->start[i]);
    }
    mpz_mod(res, total, a);
    mpz_clear(total);
}