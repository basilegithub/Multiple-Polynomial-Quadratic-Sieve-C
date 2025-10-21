#include <gmp.h>
#include <stdlib.h>

#include "structures.h"
#include "utils.h"

void create_polynomial(
    dyn_array* sol_needed,
    dyn_array* second_part,
    dyn_array* roots,
    dyn_array* inverse_a,
    dyn_array_classic* locations,
    dyn_array_classic* primes,
    dyn_array_classic* way_to_root,
    dyn_array_classic* tmp_where,
    mpz_t a,
    mpz_t n,
    mpf_t target,
    mpf_t ln2,
    mpf_t ln10,
    mpf_t best_bound,
    mpf_t e,
    unsigned long bounds[4],
    unsigned long mult
)
{
    mpz_set_ui(a, 1);
    sol_needed->len = 0;
    second_part->len = 0;
    inverse_a->len = 0;
    way_to_root->len = 0;

    mpf_set_ui(best_bound, 10);

    dyn_array selected;
    init(&selected);

    mpf_t bound, tmp_mpf, tmp_mpf2;
    mpf_inits(bound, tmp_mpf, tmp_mpf2, NULL);

    unsigned long location, tmp_long;

    mpz_t tmp_mpz, tmp_mpz2;
    mpz_inits(tmp_mpz, tmp_mpz2, NULL);

    bool valid;

    for (unsigned int attempt = 0 ; attempt < 6 ; attempt++)
    {
        tmp_where->len = 0;
        mpf_set(bound, target);
        while (mpf_cmp_ui(bound, 1) > 0)
        {
            if (mpf_cmp_ui(bound, 7) > -1) location = rand()%(bounds[3]-bounds[2])+bounds[2];
            else if (mpf_cmp_ui(bound, 4) > -1) location = rand()%(bounds[2]-bounds[1])+bounds[1];
            else location = rand()%(bounds[1]-bounds[0])+bounds[0];

            valid = true;
            if (primes->start[location] == mult) valid = false;
            for (size_t j = 0 ; j < tmp_where->len && valid; j++)
            {
                if (tmp_where->start[j] == location) valid = false;
            }
            if (valid)
            {
                append_classic(tmp_where, location);
                mpz_set_ui(tmp_mpz, primes->start[location]);
                mpf_set_z(tmp_mpf, tmp_mpz);
                my_log10(tmp_mpf, tmp_mpf, ln2, ln10, e);
                mpf_sub(bound, bound, tmp_mpf);
            }
        }
        mpf_abs(tmp_mpf, bound);
        if (!mpf_cmp_ui(best_bound, 10) || mpf_cmp(tmp_mpf, tmp_mpf2) < 0)
        {
            locations->len = 0;
            for (size_t j = 0 ; j < tmp_where->len ; j++) append_classic(locations, tmp_where->start[j]);
            mpf_set(best_bound, bound);
            mpf_abs(tmp_mpf2, best_bound);
        }
    }

    for (size_t i = 0 ; i < locations->len ; i++)
    {
        tmp_long = locations->start[i];
        mpz_set_ui(tmp_mpz2, primes->start[tmp_long]);
        mpz_mul(a, a, tmp_mpz2);
        append(&selected, tmp_mpz2);
        append_only(sol_needed, roots->start[tmp_long]);
    }

    for (size_t i = 0 ; i < selected.len ; i++)
    {
        mpz_divexact(tmp_mpz, a, selected.start[i]);
        mpz_invert(tmp_mpz2, tmp_mpz, selected.start[i]);
        mpz_mul(tmp_mpz2, tmp_mpz2, tmp_mpz);
        mpz_mod(tmp_mpz2, tmp_mpz2, a);
        append_only(second_part, tmp_mpz2);
    }

    for (size_t i = 0 ; i < primes->len ; i++)
    {
        valid = true;
        for (unsigned long j = 0 ; j < locations->len ; j++)
        {
            if (locations->start[j] == i) valid = false;
        }
        if (valid)
        {
            mpz_set_ui(tmp_mpz2, primes->start[i]);
            mpz_invert(tmp_mpz, a, tmp_mpz2);
            append_only(inverse_a, tmp_mpz);
            mpz_mul(tmp_mpz, tmp_mpz, roots->start[i]);
            mpz_mul_2exp(tmp_mpz, tmp_mpz, 1);
            mpz_mod_ui(tmp_mpz, tmp_mpz, primes->start[i]);
            append_classic(way_to_root, mpz_get_ui(tmp_mpz));
        }
        else
        {
            mpz_set_ui(tmp_mpz, 0);
            append_only(inverse_a, tmp_mpz);
            append_classic(way_to_root, 0);
        }
    }
    
    for (size_t i = 0 ; i < selected.len ; i++) mpz_clear(selected.start[i]);
    free(selected.start);
    selected.start = NULL;
    mpz_clears(tmp_mpz, tmp_mpz2, NULL);
    mpf_clears(bound, tmp_mpf, tmp_mpf2, NULL);
}

void CRT(dyn_array* moduli, dyn_array* second_part, mpz_t a, mpz_t res)
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