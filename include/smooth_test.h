#ifndef SMOOTH_TEST_H
#define SMOOTH_TEST_H

#include <gmp.h>

#include "structures.h"

void batch_smooth(dyn_array* reported, dyn_array *tmp_array, PartialRelation* tmp_array2, mpz_t limit, mpz_t limit_2, mpz_t prod_primes, unsigned long prime);
void naive_smooth(dyn_array* reported, PartialRelation *tmp_array, dyn_array_classic primes, mpz_t limit);
void pollard_rho(const mpz_t m, mpz_t p1, mpz_t p2);

#endif // SMOOTH_TEST_H