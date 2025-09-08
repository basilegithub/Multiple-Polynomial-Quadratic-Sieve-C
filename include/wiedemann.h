#ifndef WIEDEMANN_H
#define WIEDEMANN_H

#include <gmp.h>

#include "structures.h"

void poly_anul(mpz_t D, mpz_t B, unsigned long m);
void wiedemann(dyn_array_classic A, unsigned long n, unsigned char vec[n], unsigned long limit, mpz_t poly_res, unsigned long degree);

#endif // WIEDEMANN_H