#ifndef POLYNOMIAL_FUNCTIONS_H
#define POLYNOMIAL_FUNCTIONS_H

#include <gmp.h>

#include "structures.h"

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
    mpf_t target,
    mpf_t ln2,
    mpf_t ln10,
    mpf_t best_bound,
    mpf_t e,
    unsigned long bounds[4],
    unsigned long mult
);

void CRT(dyn_array* moduli, dyn_array* second_part, mpz_t a, mpz_t res);

#endif // POLYNOMIAL_FUNCTIONS_H