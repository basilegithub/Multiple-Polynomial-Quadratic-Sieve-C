#ifndef MONO_CPU_SIEVING_H
#define MONO_CPU_SIEVING_H

#include <gmp.h>
#include <time.h>

#include "structures.h"

void mono_cpu_sieve(
    dyn_array* relations,
    dyn_array* smooth_numbers,
    dyn_array a,
    dyn_array_classic primes,
    mpz_t n,
    mpz_t prod_primes,
    mpz_t cst,
    mpz_t cst2,
    mpf_t target,
    mpf_t ln2,
    mpf_t ln10,
    mpf_t e,
    unsigned long* full_found,
    unsigned long* partial_found,
    unsigned long* indexp,
    unsigned long* bounds,
    unsigned long* logs,
    unsigned long best_mult,
    unsigned long sieve_len,
    unsigned long batch_size,
    unsigned long half,
    unsigned long dim,
    unsigned long addup,
    unsigned long skipped,
    unsigned long prime_start,
    unsigned long smooth_bound,
    unsigned long prime,
    unsigned long time_diff,
    unsigned long objective,
    unsigned long seconds,
    int* need_append,
    int flag_batch_smooth,
    double nb_large,
    time_t second1,
    time_t second2
);

#endif // MONO_CPU_SIEVING_H