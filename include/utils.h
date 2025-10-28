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
bool fermat_primality(mpz_t n);

// Linear algebra functions

void multiply(dyn_array_classic A, unsigned long n, unsigned long index, bool b[n], bool res[n]);
void multiply_sparse(const dyn_array_classic A, const unsigned long dim_out, const unsigned long index, const size_t *b, size_t *res);
bool dot_prod(unsigned long n, bool lbd[n], size_t x[n]);
void add_vectors(size_t *output, const size_t *vec_a, const size_t *vec_b, const size_t N);
void identity(size_t *output, const size_t N);
void concatenate(size_t *output, const size_t *matrix_A, const size_t *matrix_B, const size_t N, const size_t dim_out);
void dense_multiply(size_t *output, const size_t *matrix_A, const size_t *matrix_B, const size_t len_A, const size_t len_B);
void sparse_multiply_transpose(const dyn_array_classic sparse_matrix, const size_t *vector, size_t *output, const unsigned long limit, const unsigned long dim);
void dense_multiply_transpose(size_t *output, size_t *matrix, size_t *vector, size_t dim1, size_t dim2);
void transpose_dense(mpz_t *output, size_t *matrix, size_t dim1, size_t dim2);

// Wiedemann polynomial functions

void poly_prod(mpz_t res, mpz_t poly_a, mpz_t poly_b);
void div_poly(mpz_t quotient, mpz_t remainder, mpz_t poly_a, mpz_t poly_b);
void poly_eval(dyn_array_classic A, mpz_t poly, unsigned long n, bool x[n], bool res[n], unsigned long limit);

#endif // UTILS_H