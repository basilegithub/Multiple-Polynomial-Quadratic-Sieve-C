#include <gmp.h>
#include <stdio.h>

#include "structures.h"

void handle_relations(
    dyn_array* relations,
    dyn_array* smooth_numbers,
    PartialRelation *tmp_array,
    Hashmap_PartialRelation *partial_relations,
    dyn_array block,
    dyn_array coefficient,
    mpz_t n,
    mpz_t value,
    mpz_t tmp_bin,
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
    mpz_t tmp;
    mpz_init(tmp);
    if (mpz_cmp_ui((tmp_array+k)->small_p, 0)) // If small_p is not 0, ie relation has at most two large primes
    {
        mpz_set(value, *(block.start+k));
        mpz_mul(tmp, value, *(coefficient.start+2*k));
        mpz_add(value, tmp, *(coefficient.start+2*k+1));
        if (!mpz_cmp_ui((tmp_array+k)->big_p, 1)) // If big_p is 1, ie relation is full
        {
            append(smooth_numbers, value);
            mpz_mul(tmp, value, value);
            mpz_sub(value, tmp, n);
            append(relations, value);
            (*full_found)++;
        }
        else
        {
            if (*need_append == 1)
            {
                *need_append = 0;
                mpz_set((tmp_array+k)->x, value);
                mpz_mul(tmp,value,value);
                mpz_sub(value,tmp,n);
                mpz_set((tmp_array+k)->y, value);
                hashmap_2d_put_node(partial_relations, *(tmp_array+k));
                (*indexp)++;
            } else {
                bool flag = hashmap_2d_is_present_mpz(partial_relations, (tmp_array+k)->small_p, (tmp_array+k)->big_p);
                if (!flag)
                {
                    mpz_set((tmp_array+k)->x, value);
                    mpz_mul(tmp,value,value);
                    mpz_sub(value,tmp,n);
                    mpz_set((tmp_array+k)->y, value);
                    hashmap_2d_put_node(partial_relations, *(tmp_array+k));
                    (*indexp)++;
                } else {
                    PartialRelation to_combine_node;
                    mpz_inits(to_combine_node.x, to_combine_node.y, to_combine_node.small_p, to_combine_node.big_p, NULL);
                    hashmap_2d_get_from_mpz(partial_relations, (tmp_array+k)->small_p, (tmp_array+k)->big_p, &to_combine_node);

                    mpz_mul(tmp_bin, value, value);
                    mpz_sub(tmp_bin, tmp_bin, n);
                    mpz_mul(tmp, tmp_bin, to_combine_node.y);
                    mpz_divexact(tmp_bin, tmp, to_combine_node.big_p);
                    mpz_divexact(tmp, tmp_bin, to_combine_node.big_p);
                    append(relations,tmp);

                    mpz_mul(value, value, to_combine_node.x);
                    mpz_invert(tmp_bin, to_combine_node.big_p, n);
                    mpz_mul(tmp, value, tmp_bin);
                    mpz_mod(value, tmp, n);
                    append(smooth_numbers,value);
                    (*partial_found)++;

                    mpz_clears(to_combine_node.x, to_combine_node.y, to_combine_node.small_p, to_combine_node.big_p, NULL);
                }
            }
        }
    }
    mpz_clear(tmp);
}