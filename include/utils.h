#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>
#include <gmp.h>

#include "structures.h"

// Functions declaration

// Precomputation functions

void compute_e(mpf_t e);

// Exponential functions

void recursive_exp(mpf_t res, mpz_t pow, mpf_t e);
void myexp(mpf_t res, mpf_t x, mpf_t e);

// Logarithm functions

void my_int_log2(mpz_t n);
unsigned long log2_ui(unsigned long n);
void natural_log(mpf_t res, mpf_t x, mpf_t ln2, mpf_t e);
void my_log10(mpf_t res, mpf_t x, mpf_t ln2, mpf_t ln10, mpf_t e);

// Number theory functions

int my_legendre(mpz_t n, unsigned long p);
void sqrt_mod(mpz_t n, unsigned long p, gmp_randstate_t state);

// Linear algebra functions

void multiply(unsigned long n, unsigned long index, dyn_array_classic A, bool b[n], bool res[n]);
bool dot_prod(unsigned long n, bool lbd[n], bool x[n]);

// Wiedemann polynomial functions

void poly_prod(mpz_t res, mpz_t poly_a, mpz_t poly_b);
void div_poly(mpz_t quotient, mpz_t remainder, mpz_t poly_a, mpz_t poly_b);
void poly_eval(unsigned long n, mpz_t poly, bool x[n], bool res[n], dyn_array_classic A, unsigned long limit);

#endif // UTILS_H