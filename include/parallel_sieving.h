#ifndef PARALLEL_SIEVING_H
#define PARALLEL_SIEVING_H

#include <gmp.h>
#include <time.h>

#include "structures.h"

void parallel_sieve(
    dyn_array* relations,
    dyn_array* smooth_numbers,
    dyn_array* store_partial,
    dyn_array* psmooth,
    dyn_array* partial,
    dyn_array_classic primes,
    dyn_array a,
    mpz_t n,
    mpz_t prod_primes,
    mpz_t cst,
    mpz_t tmp_bin,
    mpz_t* tmp_vec2,
    mpf_t nb_large,
    mpf_t target,
    mpf_t ln2,
    mpf_t ln10,
    mpf_t e,
    mpf_t var1,
    mpf_t var2,
    mpf_t var3,
    mpf_t tmpf2,
    unsigned long best_mult,
    unsigned long time_seed,
    unsigned long sieve_len,
    unsigned long batch_size,
    unsigned long half,
    unsigned long dim,
    unsigned long addup,
    unsigned long skipped,
    unsigned long prime_start,
    unsigned long smooth_bound,
    unsigned long prime,
    unsigned long tmp_a,
    unsigned long tmp_b,
    unsigned long tmplol,
    unsigned long time_diff,
    unsigned long objective,
    unsigned long seconds,
    unsigned long* full_found,
    unsigned long* partial_found,
    unsigned long* indexp,
    unsigned long* bounds,
    unsigned long* logs,
    int* need_append,
    int nb_cpu_sieve,
    int flag_batch_smooth,
    time_t second1,
    time_t second2
);

#endif // PARALLEL_SIEVING_H