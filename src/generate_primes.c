#include <gmp.h>
#include <stdlib.h>

#include "structures.h"
#include "generate_primes.h"

void erasthotenes_sieve(dyn_array_classic* primes, mpz_t bound)
{
    mpz_t sqrtB;
    mpz_init(sqrtB);

    mpz_sqrt(sqrtB, bound);
    mpz_add_ui(sqrtB, sqrtB, 1);

    unsigned long bound_ui = mpz_get_ui(bound);
    unsigned long sqrtB_ui = mpz_get_ui(sqrtB);

    bool *array = calloc(bound_ui, sizeof(bool));

    for (size_t i = 0 ; i < bound_ui ; i += 2) array[i] = true;

    for (size_t i = 1 ; i < bound_ui ; i += 2) array[i] = false;

    append_classic(primes, 2);

    unsigned long index = 3;
    unsigned long tmp = 2, square_index, i;

    while (index < sqrtB_ui)
    {
        if (array[tmp])
        {
            square_index = index*index;
            append_classic(primes, index);
            i = square_index;
            while(i <= bound_ui)
            {
                array[i-1] = false;
                i += index;
            }
        }
        index += 2;
        tmp = index-1;
    }

    i = sqrtB_ui;
    while (i <= bound_ui)
    {
        if (array[i-1] == 1)
        {
            append_classic(primes, i);
        }
        i++;
    }
    mpz_clear(sqrtB);
    free(array);
}