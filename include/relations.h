#ifndef RELATIONS_H
#define RELATIONS_H

#include <gmp.h>

#include "structures.h"

void handle_relations(
    dyn_array* relations,
    dyn_array* smooth_numbers,
    PartialRelation *large_primes,
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
);

#endif // RELATIONS_H