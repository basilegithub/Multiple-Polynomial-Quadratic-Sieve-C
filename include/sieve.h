#ifndef SIEVE_H
#define SIEVE_H

#include <gmp.h>

#include "structures.h"

void sieve(dyn_array_small* sieve,
    unsigned long length,
    unsigned long L,
    dyn_array_classic* primes,
    unsigned long * logs,
    dyn_array* roots,
    mpz_t n,
    mpz_t poly_a,
    mpz_t poly_b,
    mpz_t poly_c,
    dyn_array* inverse_a,
    dyn_array_classic* way_to_root,
    dyn_array_classic* locations,
    unsigned long skipped,
    unsigned long prime_start,
    dyn_array* smooth,
    unsigned long smooth_bound);

#endif // SIEVE_H