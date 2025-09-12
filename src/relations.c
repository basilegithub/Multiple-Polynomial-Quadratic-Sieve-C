#include <gmp.h>

#include "structures.h"

void handle_relations(
    dyn_array* relations,
    dyn_array* smooth_numbers,
    dyn_array* store_partial,
    dyn_array* psmooth,
    dyn_array* partial,
    dyn_array large_primes,
    dyn_array block,
    dyn_array coefficient,
    dyn_array_classic is_smooth,
    mpz_t n,
    mpz_t value,
    mpz_t tmp_bin,
    mpz_t* tmp_vec2,
    unsigned long k,
    unsigned long tmp_a,
    unsigned long tmp_b,
    unsigned long tmplol,
    unsigned long* full_found,
    unsigned long* partial_found,
    unsigned long* indexp,
    int* need_append
)
{
    if (*(is_smooth.start+k))
    {
        *(is_smooth.start+k) = 0;
        mpz_set(value,*(block.start+k));
        mpz_mul(value,value,*(coefficient.start+2*k));
        mpz_add(value,value,*(coefficient.start+2*k+1));
        if (!mpz_cmp_ui(*(large_primes.start+k),1))
        {
            append(smooth_numbers,value);
            mpz_mul(value,value,value);
            mpz_sub(value,value,n);
            append(relations,value);
            (*full_found)++;
        }
        else
        {
            if (*need_append == 1)
            {
                *need_append = 0;
                append(psmooth,value);
                mpz_mul(value,value,value);
                mpz_sub(value,value,n);
                append(partial,value);
                mpz_set(tmp_vec2[1],*(large_primes.start+k));
                append_block(store_partial,2,tmp_vec2);
                (*indexp)++;
                mpz_add_ui(tmp_vec2[0],tmp_vec2[0],1);
            } else {
                mpz_set(tmp_vec2[1],*(large_primes.start+k));
                tmp_a = 0;
                tmp_b = *indexp - 1;
                tmplol = (tmp_a+tmp_b)/2;
                mpz_set(tmp_bin,*(store_partial->start+tmplol*2+1));
                while (mpz_cmp(tmp_bin,tmp_vec2[1]) != 0 && tmp_a <= tmp_b)
                {
                    if (mpz_cmp(tmp_bin,tmp_vec2[1]) < 0)
                    {
                        tmp_a = tmplol+1;
                    } else {
                        if (tmplol == 0) break;
                        tmp_b = tmplol-1;
                    }
                    tmplol = (tmp_a+tmp_b)/2;
                    mpz_set(tmp_bin,*(store_partial->start+tmplol*2+1));
                }
                if (mpz_cmp(tmp_bin,tmp_vec2[1]) != 0)
                {
                    append(psmooth,value);
                    mpz_mul(value,value,value);
                    mpz_sub(value,value,n);
                    append(partial,value);
                    insert_block(store_partial,tmp_a*2,2,tmp_vec2);
                    (*indexp)++;
                    mpz_add_ui(tmp_vec2[0],tmp_vec2[0],1);
                } else {
                    tmplol = mpz_get_ui(*(store_partial->start+2*tmplol));
                    mpz_mul(tmp_bin,value,value);
                    mpz_sub(tmp_bin,tmp_bin,n);
                    mpz_mul(tmp_bin,tmp_bin,*(partial->start+tmplol));
                    mpz_divexact(tmp_bin,tmp_bin,tmp_vec2[1]);
                    mpz_divexact(tmp_bin,tmp_bin,tmp_vec2[1]);
                    append(relations,tmp_bin);

                    mpz_mul(value,value,*(psmooth->start+tmplol));
                    mpz_invert(tmp_bin,tmp_vec2[1],n);
                    mpz_mul(value,value,tmp_bin);
                    mpz_mod(value,value,n);
                    append(smooth_numbers,value);
                    (*partial_found)++;
                }
            }
            mpz_set_ui(*(large_primes.start+k),1);
        }
    }
}