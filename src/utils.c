#include <gmp.h>
#include <stdlib.h>

#include "structures.h"

void compute_e(mpf_t e)
{
    mpf_set_ui(e, 1);
    mpf_set_prec(e, 96);
    mpf_t fac, tmp;
    mpf_init_set_ui(fac, 1);
    mpf_init(tmp);
    for (unsigned long i = 1 ; i < 50 ; i++)
    {
        mpf_mul_ui(fac, fac, i);
        mpf_ui_div(tmp, 1, fac);
        mpf_add(e, e, tmp);
    }
    mpf_clears(fac, tmp, NULL);
}

void recursive_exp(mpf_t res, mpz_t pow, mpf_t e)
{
    if (!mpz_cmp_ui(pow, 0)) mpf_set_ui(res, 1);
    else if (!mpz_cmp_ui(pow, 1)) mpf_set(res, e);
    else
    {
        mpf_t tmpf;
        mpf_init(tmpf);
        mpz_t tmp;
        mpz_init(tmp);
        if (mpz_divisible_ui_p(pow, 2))
        {
            mpz_div_ui(tmp, pow, 2);
            recursive_exp(tmpf, tmp, e);
            mpf_mul(res, tmpf, tmpf);
        }
        else
        {
            mpz_div_ui(tmp, pow, 2);
            recursive_exp(tmpf, tmp, e);
            mpf_mul(tmpf, tmpf, tmpf);
            mpf_mul(res, e, tmpf);
        }
        mpf_clear(tmpf);
        mpz_clear(tmp);
    }
}

void myexp(mpf_t res, mpf_t x, mpf_t e)
{
    int flag = 0;
    if (mpf_cmp_ui(x, 0) < 0)
    {
        flag = 1;
        mpf_neg(x, x);
    }
    mpf_t tmp_res, rest;
    mpz_t tmp;

    mpf_init(rest);
    mpf_init(tmp_res);
    mpz_init(tmp);

    mpz_set_f(tmp, x);
    mpf_set_z(rest, tmp);
    mpf_sub(rest, x, rest);

    recursive_exp(tmp_res, tmp, e);

    mpf_t tmp_rest, fact, tmpf, tmpf2;

    mpf_init_set(tmp_rest, rest);
    mpf_init_set_ui(fact, 1);
    mpf_init_set_ui(tmpf, 1);
    mpf_init(tmpf2);

    for (unsigned long i = 2 ; i <= 10 ; i++)
    {
        mpf_div(tmpf2, tmp_rest, fact);
        mpf_add(tmpf, tmpf, tmpf2);
        mpf_mul_ui(fact, fact, i);
        mpf_mul(tmp_rest, tmp_rest, rest);
    }
    mpf_mul(res, tmp_res, tmpf);
    if (flag)
    {
        mpf_set_ui(tmpf, 1);
        mpf_div(res, tmpf, res);
    }
    mpz_clear(tmp);
    mpf_clears(tmpf, tmpf2, rest, tmp_res, tmp_rest, fact, NULL);
}

void my_int_log2(mpz_t n)
{
    mpz_set_ui(n, mpz_sizeinbase(n, 2) - 1);
}

void natural_log(mpf_t res, mpf_t x, mpf_t ln2, mpf_t e)
{
    mpz_t tmp;

    mpz_init(tmp);
    mpz_set_f(tmp, x);
    mpz_set_ui(tmp, mpz_sizeinbase(tmp, 2));

    mpf_t a, tmpf;

    mpf_init(tmpf);
    mpf_init(a);
    mpf_set_prec(a, 32);
    mpf_set_z(a, tmp);
    mpf_mul(a, a, ln2);

    for (unsigned long i = 0 ; i < 5 ; i++)
    {
        myexp(tmpf, a, e);
        mpf_div(tmpf, x, tmpf);
        mpf_sub_ui(a, a, 1);
        mpf_add(a, a, tmpf);
    }
    mpf_set(res, a);
    mpz_clear(tmp);
    mpf_clears(tmpf, a, NULL);
}

void my_log10(mpf_t res, mpf_t x, mpf_t ln2, mpf_t ln10, mpf_t e)
{
    mpz_t tmp;

    mpz_init(tmp);
    mpz_set_f(tmp, x);
    mpz_set_ui(tmp,mpz_sizeinbase(tmp, 2));

    mpf_t a, tmpf;

    mpf_init(tmpf);
    mpf_init(a);
    mpf_set_prec(a, 32);
    mpf_set_z(a, tmp);
    mpf_mul(a, a, ln2);

    for (unsigned long i = 0 ; i < 6 ; i++)
    {
        myexp(tmpf, a, e);
        mpf_div(tmpf, x, tmpf);
        mpf_sub_ui(a, a, 1);
        mpf_add(a, a, tmpf);
    }
    mpf_set(res, a);
    mpf_div(res, res, ln10);
    mpz_clear(tmp);
    mpf_clears(tmpf, a, NULL);
}

int my_legendre(mpz_t n, unsigned long p)
{
    mpz_t tmp;
    mpz_init(tmp);
    unsigned long tmpl;
    mpz_mod_ui(tmp,n,p);
    tmpl = mpz_get_ui(tmp);
    int t = 1;
    unsigned long tmps;
    while (tmpl)
    {
        while (tmpl%2 == 0)
        {
            tmpl >>= 1;
            if (p%8 == 3 || p%8 == 5) t = -t;
        }
        tmps = tmpl;
        tmpl = p;
        p = tmps;
        if (tmpl%4 == p%4 && p%4 == 3) t = -t;
        tmpl %= p;
    }
    mpz_clear(tmp);
    if (p == 1) return t;
    return 0;
}

void sqrt_mod(mpz_t n, unsigned long p, gmp_randstate_t state)
{
    mpz_t z,tmp,tmp2,P_value;
    unsigned long P;
    unsigned long tmp3;
    unsigned long r = 0;
    mpz_init(z);
    mpz_init(tmp);

    mpz_mod_ui(n,n,p);
    P = p-1;
    mpz_init_set_ui(P_value,P);
    mpz_set_ui(tmp,p);
    mpz_urandomm(z,state,tmp);
    if (mpz_cmp_ui(z,1) != 1) mpz_set_ui(z,2);
    while (my_legendre(z,p) != -1)
    {
        mpz_urandomm(z,state,tmp);
        if (mpz_cmp_ui(z,1) != 1) mpz_set_ui(z,2);
    }
    while (P%2 == 0)
    {
        P >>= 1;
        r++;
    }
    mpz_t generator, lambda, omega, res, m, two_mpz;
    mpz_init(generator);
    mpz_init(lambda);
    mpz_init(omega);
    mpz_powm_ui(generator,z,P,tmp);
    mpz_powm_ui(lambda,n,P,tmp);
    tmp3 = (P+1)>>1;
    mpz_powm_ui(omega,n,tmp3,tmp);
    mpz_init_set_ui(res,0);
    mpz_init(m);
    mpz_init_set_ui(two_mpz,2);
    mpz_t tmp_l;
    mpz_init(tmp_l);
    mpz_init(tmp2);
    while (1)
    {
        if (mpz_cmp_ui(lambda,0) == 0)
        {
            mpz_set_ui(n,0);
            break;
        }
        if (mpz_cmp_ui(lambda,1) == 0)
        {
            mpz_set(n,omega);
            break;
        }
        mpz_set_ui(m,1);
        while (mpz_cmp_ui(m,r) < 0)
        {
            mpz_powm(tmp_l,two_mpz,m,P_value);
            mpz_powm(tmp_l,lambda,tmp_l,tmp);
            if (mpz_cmp_ui(tmp_l,1) == 0) break;
            mpz_add_ui(m,m,1);
        }
        mpz_ui_sub(tmp2,r,m);
        mpz_sub_ui(tmp2,tmp2,1);
        mpz_powm(tmp_l,two_mpz,tmp2,P_value);
        mpz_mul_ui(tmp2,tmp_l,2);

        mpz_powm(tmp2,generator,tmp2,tmp);
        mpz_mul(lambda,lambda,tmp2);
        mpz_mod(lambda,lambda,tmp);

        mpz_powm(tmp2,generator,tmp_l,tmp);
        mpz_mul(omega,omega,tmp2);
        mpz_mod(omega,omega,tmp);
    }
    mpz_clears(z, tmp, tmp2, P_value, generator, lambda, omega, res, m, two_mpz, NULL);
}

void multiply(unsigned long n, unsigned long index, dyn_array_classic A, unsigned char b[n], unsigned char res[n])
{
    unsigned char tmp[n];
    unsigned long i = 0;
    unsigned char tmp2 = 0;
    for (unsigned long k = 0 ; k < A.len ; k++)
    {
        if (*(A.start+k) == index)
        {
            tmp[i] = tmp2;
            i++;
            tmp2 = 0;
        }
        else tmp2 ^= b[*(A.start+k)];
    }
    for (unsigned long j = 0 ; j < i ; j++) res[j] = tmp[j];
    for (unsigned long j = i ; j < n ; j++) res[j] = 0;
}

unsigned char dot_prod(unsigned long n, unsigned char lbd[n], unsigned char x[n])
{
    unsigned char tmp = 0;
    for (unsigned long i = 0 ; i < n ; i++)
    {
        if (lbd[i]) tmp ^= x[i];
    }
    return tmp;
}

void poly_prod(mpz_t res, mpz_t poly_a, mpz_t poly_b)
{
    mpz_t tmp_poly;
    mpz_init_set_ui(tmp_poly,0);
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set(tmp,poly_a);
    my_int_log2(tmp);
    for (unsigned long i = mpz_get_ui(tmp) ; i > 0 ; i--)
    {
        mpz_div_2exp(tmp,poly_a,i);
        if (!mpz_divisible_ui_p(tmp,2)) mpz_xor(tmp_poly,tmp_poly,poly_b);
        mpz_mul_2exp(tmp_poly,tmp_poly,1);
    }
    if (!mpz_divisible_ui_p(poly_a,2)) mpz_xor(tmp_poly,tmp_poly,poly_b);
    mpz_set(res,tmp_poly);
    mpz_clears(tmp,tmp_poly,NULL);
}

void div_poly(mpz_t quotient, mpz_t remainder, mpz_t poly_a, mpz_t poly_b)
{
    mpz_t tmp, tmp2;
    mpz_init_set(tmp,poly_a); // problem is there
    mpz_init_set(tmp2,poly_b); // or at this line : GNU MP cannot reallocate memory (new_size = 40000; old_size = 4)
    my_int_log2(tmp);
    my_int_log2(tmp2);
    if (mpz_cmp(tmp,tmp2) < 0)
    {
        mpz_set(remainder,poly_a);
        mpz_set_ui(quotient,0);
    }
    else
    {
        mpz_set(remainder,poly_a);
        mpz_set_ui(quotient,0);
        mpz_set(tmp,poly_a);
        mpz_set(tmp2,poly_b);
        my_int_log2(tmp);
        my_int_log2(tmp2);
        signed long dif = mpz_get_ui(tmp)-mpz_get_ui(tmp2);
        for (signed long j = dif ; j > 0 ; j--)
        {
            if (mpz_cmp_ui(remainder,0))
            {
                mpz_add_ui(quotient,quotient,1);
                mpz_mul_2exp(tmp,poly_b,j);
                mpz_xor(remainder,remainder,tmp);
            }
            mpz_mul_ui(quotient,quotient,2);
        }
        mpz_set(tmp,remainder);
        my_int_log2(tmp);
        if (mpz_cmp_ui(remainder,0) && mpz_cmp(tmp,tmp2) > -1)
        {
            mpz_add_ui(quotient,quotient,1);
            mpz_xor(remainder,remainder,poly_b);
        }
    }
    mpz_clears(tmp,tmp2,NULL);
}

void poly_eval(unsigned long n, mpz_t poly, unsigned char x[n], unsigned char res[n], dyn_array_classic A, unsigned long limit)
{
    unsigned char tmp2[n];
    for (unsigned long i = 0 ; i < n ; i++) tmp2[i] = 0;
    mpz_t tmp;
    mpz_init_set(tmp,poly);
    my_int_log2(tmp);
    unsigned long degree = mpz_get_ui(tmp);
    for (unsigned long i = 0 ; i < degree ; i++)
    {
        mpz_div_2exp(tmp,poly,degree-i);
        if (!mpz_divisible_ui_p(tmp,2))
        {
             for (unsigned long j = 0 ; j < n ; j++) tmp2[j] ^= x[j];
        }
        multiply(n,limit,A,tmp2,tmp2);
    }
    if (!mpz_divisible_ui_p(poly,2))
    {
        for (unsigned long j = 0 ; j < n ; j++) tmp2[j] ^= x[j];
    }
    for (unsigned long i = 0 ; i < n ; i++) res[i] = tmp2[i];
    mpz_clear(tmp);
}