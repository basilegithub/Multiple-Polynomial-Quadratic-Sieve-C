#ifndef SIEVE_H
#define SIEVE_H

#include <gmp.h>

#include "structures.h"

void sieve(
    dyn_array* smooth,
    dyn_array* roots,
    dyn_array* inverse_a,
    dyn_array_classic *primes,
    dyn_array_classic *way_to_root,
    dyn_array_classic *locations,
    dyn_array_small *sieve,
    mpz_t n,
    mpz_t poly_a,
    mpz_t poly_b,
    mpz_t poly_c,
    unsigned long *logs,
    unsigned long length,
    unsigned long L,
    unsigned long skipped,
    unsigned long prime_start,
    unsigned long smooth_bound,
    signed long *tmp_array
);

#endif // SIEVE_H