#include <gmp.h>
#include <stdlib.h>

#include "structures.h"
#include "utils.h"

void poly_anul(mpz_t D, mpz_t B, unsigned long m)
{
    mpz_t tmp;
    mpz_init(tmp);
    mpz_t A, C, E, Q, R;
    mpz_init_set_ui(A,1);
    mpz_mul_2exp(A,A,2*m);
    mpz_init_set_ui(C,0);
    mpz_set_ui(D,1);
    mpz_init(E);
    mpz_init(Q);
    mpz_init(R);
    mpz_set(tmp,B);
    my_int_log2(tmp);
    while (mpz_cmp_ui(tmp,m) > -1)
    {
        div_poly(Q,R,A,B);
        poly_prod(tmp,Q,D);
        mpz_xor(E,C,tmp);
        mpz_set(C,D);
        mpz_set(D,E);
        mpz_set(A,B);
        mpz_set(B,R);
        mpz_set(tmp,B);
        my_int_log2(tmp);
    }
    mpz_clears(A, C, E, Q, R, tmp, NULL);
}

void wiedemann(dyn_array_classic A, unsigned long n, unsigned char vec[n], unsigned long limit, mpz_t poly_res, unsigned long degree)
{
    unsigned char tmp[n];
    unsigned char block[n];
    for (unsigned long i = 0 ; i < n ; i++)
    {
        tmp[i] = 0;
        block[i] = 0;
    }
    multiply(n,limit,A,vec,block);
    mpz_t P;
    unsigned long d = degree;
    mpz_init_set_ui(P,1);
    mpz_t sequence;
    mpz_init(sequence);
    unsigned char lbd[n];
    int michel = 1;
    mpz_t pi;
    mpz_init(pi);
    while (michel)
    {
        for (unsigned long i = 0 ; i < n ; i++)
        {
            lbd[i] = rand()%2;
            tmp[i] = block[i];
        }
        mpz_set_ui(sequence,0);
        for (unsigned long i = 0 ; i < 2*(n-d)-1 ; i++)
        {
            mpz_add_ui(sequence,sequence,dot_prod(n,lbd,tmp));
            mpz_mul_ui(sequence,sequence,2);
            multiply(n,limit,A,tmp,tmp);
        }
        mpz_add_ui(sequence,sequence,dot_prod(n,lbd,tmp));
        mpz_set_ui(pi,1);

        poly_anul(pi,sequence,n-d);
        poly_prod(P,P,pi);
        poly_eval(n,pi,block,block,A,limit);
        my_int_log2(pi);
        d += mpz_get_ui(pi);

        michel = 0;
        for (unsigned long i = 0 ; i < n && !michel; i++)
        {
            if (block[i]) michel = 1;
        }
    }
    mpz_set(poly_res,P);
    mpz_clears(P,sequence,pi,NULL);
}