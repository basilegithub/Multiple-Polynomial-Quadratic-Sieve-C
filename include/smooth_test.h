#ifndef SMOOTH_TEST_H
#define SMOOTH_TEST_H

#include <gmp.h>

#include "structures.h"

void pollard_rho(const mpz_t m, mpz_t p1, mpz_t p2, gmp_randstate_t state);
unsigned long build_product_tree(dyn_array* reported, dyn_array* tmp_array, mpz_t prod_primes, mpz_t prod_primes_p1, mpz_t tmp);
void build_remainder_tree(dyn_array* reported, dyn_array* tmp_array, mpz_t prod_primes, mpz_t prod_primes_p1, unsigned long tmp_long);
void batch_smooth(PartialRelation* tmp_array2, dyn_array* reported, dyn_array *tmp_array, mpz_t limit, mpz_t limit_2, mpz_t prod_primes, mpz_t prod_primes_p1, unsigned long prime, gmp_randstate_t state);
void naive_smooth(PartialRelation *tmp_array, dyn_array* reported, dyn_array_classic primes, mpz_t limit, mpz_t limit_2, gmp_randstate_t state);

#endif // SMOOTH_TEST_H