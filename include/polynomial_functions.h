#ifndef POLYNOMIAL_FUNCTIONS_H
#define POLYNOMIAL_FUNCTIONS_H

#include <gmp.h>

#include "structures.h"

void create_polynomial(mpz_t a, dyn_array* sol_needed, dyn_array* second_part, dyn_array_classic* locations, mpz_t n, dyn_array_classic* primes, dyn_array* roots, unsigned long bounds[4], mpf_t target, mpf_t ln2, mpf_t ln10, dyn_array* inverse_a, dyn_array_classic* way_to_root, mpf_t best_bound, dyn_array_classic* tmp_where, mpf_t e, unsigned long mult);
void CRT(mpz_t res, dyn_array* moduli, mpz_t a, dyn_array* second_part);

#endif // POLYNOMIAL_FUNCTIONS_H