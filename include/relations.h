#ifndef RELATIONS_H
#define RELATIONS_H

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
);

#endif // RELATIONS_H