#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gmp.h>

#include "structures.h"
#include "utils.h"
#include "generate_primes.h"
#include "mono_cpu_sieving.h"
#include "parallel_sieving.h"
#include "wiedemann.h"
#include "build_matrix.h"
#include "compute_sqrt.h"
#include "reduce_matrix.h"

void count(dyn_array_classic matrix, unsigned long limit,unsigned long dim)
{
    unsigned long nb_lines = 0;
    unsigned long nonzero = 0;
    for (unsigned long i = 0 ; i < matrix.len ; i++)
    {
        if (*(matrix.start+i) == limit) nb_lines++;
        else nonzero++;
    }
    double density = (double) nonzero/nb_lines;
    struct tm tm = *localtime(&(time_t){time(NULL)});
    printf("%d:%d:%d  matrix reduced to %lux%lu ; %lu nonzero values, density = %.2f ",tm.tm_hour,tm.tm_min,tm.tm_sec,dim,nb_lines,nonzero,density);
    mpz_t tmp,tmp2;
    mpz_init(tmp2);
    mpz_init_set_d(tmp,density);
    mpz_mul_ui(tmp,tmp,dim);
    mpz_mul_ui(tmp,tmp,dim);
    mpz_mul_ui(tmp,tmp,13);
    mpz_div_ui(tmp,tmp,220433617);
    mpz_div_ui(tmp2,tmp,3600);
    gmp_printf("%Zdh ",tmp2);
    mpz_mod_ui(tmp2,tmp,3600);
    mpz_div_ui(tmp2,tmp2,60);
    gmp_printf("%Zdm ",tmp2);
    mpz_mod_ui(tmp2,tmp,60);
    gmp_printf("%Zds needed\n",tmp2);
    mpz_clears(tmp,tmp2,NULL);
}

void initialize_1(gmp_randstate_t state, mpf_t ln10, mpf_t ln2, mpf_t e)
{
    mpz_t n, b;
    mpz_init_set_ui(n, 3);
    mpz_init_set_ui(b, 82939);
    mpz_powm_ui(n, n, rand()*23%103, b);
    gmp_randseed(state, n);

    mpf_init_set_d(ln10,2.302585092994046);
    mpf_init_set_d(ln2,0.6931471805599453);

    compute_e(e);

    mpz_clear(n);
    mpz_clear(b);
}

unsigned long compute_best_mult(mpf_t best_time, mpf_t work, mpf_t time1, mpf_t time2, mpf_t e, mpf_t ln2, mpf_t tmpf, mpf_t tmpf2, mpz_t N, mpz_t n, mpz_t tmp, mpz_t tmp2, dyn_array_classic tmp_primes)
{
    int res;
    unsigned long best_mult;

    for (unsigned long k = 1 ; k < 40 ; k++)
    {
        mpz_mul_ui(n,N,k);
        mpf_set_z(tmpf,n);
        natural_log(tmpf,tmpf,ln2,e);
        natural_log(tmpf2,tmpf,ln2,e);
        mpf_mul(tmpf,tmpf,tmpf2);
        mpf_sqrt(tmpf,tmpf);
        mpf_div_ui(work,tmpf,2);

        mpz_mod_ui(tmp2,n,8);
        if (!mpz_cmp_ui(tmp2,1))
        {
            mpf_mul_ui(tmpf2,ln2,2);
            mpf_sub(work,work,tmpf2);
        }
        else if (!mpz_cmp_ui(tmp2,5)) mpf_sub(work,work,ln2);
        else
        {
            mpz_mod_ui(tmp2,n,4);
            if (!mpz_cmp_ui(tmp2,3))
            {
                mpf_div_ui(tmpf2,ln2,2);
                mpf_sub(work,work,tmpf2);
            }
        }
        for (unsigned long i = 1 ; i < tmp_primes.len ; i++)
        {
            mpz_set_ui(tmp,*(tmp_primes.start+i));
            res = mpz_legendre(n,tmp);
            if (res == 1)
            {
                mpf_set_z(tmpf,tmp);
                mpf_set(tmpf2,tmpf);
                natural_log(tmpf,tmpf,ln2,e);
                mpf_sub_ui(tmpf2,tmpf2,1);
                mpf_mul_ui(tmpf,tmpf,2);
                mpf_div(tmpf2,tmpf,tmpf2);
                mpf_sub(work,work,tmpf2);
            }
            else if (!res)
            {
                mpf_set_z(tmpf,tmp);
                mpf_set(tmpf2,tmpf);
                natural_log(tmpf,tmpf,ln2,e);
                mpf_sub_ui(tmpf2,tmpf2,1);
                mpf_div(tmpf2,tmpf,tmpf2);
                mpf_sub(work,work,tmpf2);
            }
        }
        mpf_mul_ui(work,work,2);
        myexp(work,work,e);
        mpf_set_z(tmpf,n);
        mpf_mul(work,work,tmpf);
        natural_log(time2,work,ln2,e);
        natural_log(work,time2,ln2,e);
        mpf_mul(time2,time2,work);
        mpf_sqrt(time2,time2);
        if (mpf_cmp_si(best_time,-1) == 0 || mpf_cmp(time2,best_time) < 0)
        {
            if (mpf_cmp_si(best_time,-1) == 0)
            {
                mpf_set(time1,time2);
            }
            mpf_set(best_time,time2);
            best_mult = k;
        }
    }

    return best_mult;
}

void compute_b(mpf_t bound, mpf_t e, mpf_t ln2, mpf_t tmpf, mpf_t tmpf2, mpz_t b)
{
    natural_log(tmpf,bound,ln2,e);
    natural_log(tmpf2,tmpf,ln2,e);
    mpf_mul(tmpf,tmpf,tmpf2);
    mpf_sqrt(tmpf,tmpf);
    mpf_div_ui(tmpf,tmpf,2);
    myexp(tmpf,tmpf,e);
    mpf_mul_ui(tmpf,tmpf,10);
    mpf_div_ui(tmpf,tmpf,32);
    mpz_set_f(b,tmpf);
}

void print_bits_and_digits(mpz_t n, mpz_t tmp)
{
    mpz_set(tmp,n);
    my_int_log2(tmp);
    mpz_add_ui(tmp, tmp, 1);
    gmp_printf("%Zd bits ; ", tmp);

    unsigned long digits = 0;
    mpz_set(tmp, n);
    while (mpz_cmp_ui(tmp, 0))
    {
        mpz_div_ui(tmp, tmp, 10);
        digits++;
    }
    printf("%lu digits\n", digits);
}

void compute_factor_base(dyn_array_classic* primes, dyn_array* a, dyn_array_classic B, mpz_t* prod_primes, mpz_t n, mpz_t tmp2, gmp_randstate_t state)
{
    mpz_init_set_ui(*prod_primes,2);
    for (unsigned long i = 1 ; i < B.len ; i++)
    {
        if (my_legendre(n,*(B.start+i)) == 1)
        {
            append_classic(primes,*(B.start+i));
            mpz_set(tmp2,n);
            sqrt_mod(tmp2,*(B.start+i),state);
            append(a,tmp2);
            mpz_mul_ui(*prod_primes,*prod_primes,*(B.start+i));
        }
        else if (!my_legendre(n,*(B.start+i)))
        {
            append_classic(primes,*(B.start+i));
            mpz_set_ui(tmp2,0);
            append(a,tmp2);
            mpz_mul_ui(*prod_primes,*prod_primes,*(B.start+i));
        }
    }
}

void compute_bounds(unsigned long* bounds, dyn_array_classic primes)
{
    for (unsigned long i = 1 ; i < primes.len ; i++)
    {
        if (*(primes.start+i) >= 10 && bounds[0] == 0) bounds[0] = i;
        else if (*(primes.start+i) >= 100 && bounds[1] == 0) bounds[1] = i;
        else if (*(primes.start+i) >= 1000 && bounds[2] == 0) bounds[2] = i;
        else if (*(primes.start+i) >= 3000 && bounds[3] == 0)
        {
            bounds[3] = i-1;
            break;
        }
    }
    if (bounds[2] != 0 && bounds[3] == 0 && primes.len-1-bounds[2] > 20)
    {
        bounds[3] = primes.len-1;
    }
}

void compute_logs(dyn_array_classic primes, mpz_t tmp, unsigned long* logs)
{
    mpz_t last,last_log;
    mpz_init_set_ui(last,2);
    mpz_init_set_ui(last_log,1);
    for (unsigned long i = 1 ; i < primes.len ; i++)
    {
        mpz_set(tmp,last);
        mpz_mul_ui(tmp,tmp,1414213562);
        mpz_div_ui(tmp,tmp,1000000000);
        if (mpz_cmp_ui(tmp,*(primes.start+i)) >= 0)
        {
            logs[i] = mpz_get_ui(last_log);
            continue;
        }
        mpz_mul_ui(tmp,last,2);
        if (mpz_cmp_ui(tmp,*(primes.start+i)) >= 0)
        {
            mpz_add_ui(tmp,last_log,1);
            logs[i] = mpz_get_ui(tmp);
            continue;
        }
        while (mpz_cmp_ui(tmp,*(primes.start+i)) == -1)
        {
            mpz_add_ui(last_log,last_log,1);
            mpz_mul_ui(last,last,2);
            mpz_mul_ui(tmp,last,2);
        }
        logs[i] = mpz_get_ui(last_log);
    }
}

int compute_factors(dyn_array relations, dyn_array smooth_numbers, dyn_array_classic bin_matrix, dyn_array_classic primes, mpz_t N, mpz_t tmp, mpz_t tmp2, unsigned long len)
{
    struct tm tm;

    mpz_t x,y;
    mpz_init(x);
    mpz_init(y);
    unsigned char null_space[relations.len];
    mpz_t poly;
    mpz_init_set_ui(poly,1);
    unsigned long tmp_long = 0;
    int flag;
    unsigned long k;
    mpz_t poly_res;
    mpz_init(poly_res);
    mpz_t poly_shift;
    mpz_init_set_ui(poly_shift,1);

    int block_len = 8; // less or equal to 8 for unsigned char to work
    unsigned char michel;
    unsigned char tmp_array[relations.len];
    unsigned char tmp_vec[relations.len];
    unsigned long degree = 0;
    while (1)
    {
        tm = *localtime(&(time_t){time(NULL)});
        printf("\n%d:%d:%d  new block\n",tm.tm_hour,tm.tm_min,tm.tm_sec);
        for (unsigned long i = 0 ; i < relations.len ; i++)
        {
            null_space[i] = 0;
            for (unsigned long j = 1 ; j < block_len ; j++)
            {
                null_space[i] += rand()%2;
                null_space[i] <<= 1;
            }
            null_space[i] += rand()%2;
        }
        if (mpz_cmp_ui(poly,1) > 0) poly_eval(relations.len,poly,null_space,null_space,bin_matrix,len);
        for (unsigned long i = 0 ; i < block_len ; i++)
        {
            michel = 1;
            for (unsigned long j = 0 ; j < relations.len ; j++)
            {
                tmp_array[j] = (null_space[j]>>i)&1;
                tmp_vec[j] = tmp_array[j];
            }
            if (tmp_long > 0)
            {
                k = 0;
                flag = 0;
                for (unsigned long j = 0 ; j < relations.len && !flag ; j++)
                {
                    if (tmp_array[j]) flag = 1;
                }
                while (k <= tmp_long && flag)
                {
                    for (unsigned long j = 0 ; j < relations.len ; j++) tmp_vec[j] = tmp_array[j];
                    multiply(relations.len,len,bin_matrix,tmp_array,tmp_array);
                    k++;
                    flag = 0;
                    for (unsigned long j = 0 ; j < relations.len && !flag ; j++)
                    {
                        if (tmp_array[j]) flag = 1;
                    }
                }
                flag = 1;
                for (unsigned long j = 0 ; j < relations.len && flag ; j++)
                {
                    if (tmp_array[j]) flag = 0;
                }
                if (flag) michel = 0;
                for (unsigned long j = 0 ; j < relations.len ; j++) tmp_array[j] = tmp_vec[j];
            }
            if (michel)
            {
                wiedemann(bin_matrix,relations.len,tmp_array,len,poly_res,degree);
                mpz_set(tmp,poly_res);
                my_int_log2(tmp);
                degree += mpz_get_ui(tmp);
                while (mpz_divisible_ui_p(poly_res,2))
                {
                    mpz_div_ui(poly_res,poly_res,2);
                    tmp_long++;
                    mpz_mul_ui(poly_shift,poly_shift,2);
                }
                poly_prod(poly,poly,poly_res);
                poly_eval(relations.len,poly_res,null_space,null_space,bin_matrix,len);
                for (unsigned long j = 0 ; j < relations.len ; j++) tmp_array[j] = (null_space[j]>>i)&1;
                if (mpz_cmp_ui(poly_shift,1) > 0) poly_eval(relations.len,poly_shift,tmp_array,tmp_array,bin_matrix,len);
            }
            tm = *localtime(&(time_t){time(NULL)});
            printf("%d:%d:%d  ",tm.tm_hour,tm.tm_min,tm.tm_sec);
            if (michel) printf("minimal polynomial updated, ");
            printf("kernel vector found\n");
            mpz_set_ui(x,1);
            mpz_set_ui(y,1);
            build_sqrt(N,relations.len,tmp_array,primes,x,y,relations,smooth_numbers);
            mpz_sub(tmp,x,y);
            mpz_add(tmp2,x,y);
            mpz_gcd(tmp,tmp,N);
            mpz_gcd(tmp2,tmp2,N);
            char array1;
            char array2;
            //gmp_printf("%Zd %Zd\n",tmp,tmp2);
            if (mpz_cmp_ui(tmp,1) != 0 && mpz_cmp(tmp,N) != 0)
            {
                if (mpz_probab_prime_p(tmp,100) > 0)
                {
                    array1 = 'p';
                } else {array1 = 'C';}
                if (mpz_probab_prime_p(tmp2,100) > 0)
                {
                    array2 = 'p';
                } else {array2 = 'C';}
                tm = *localtime(&(time_t){time(NULL)});
                gmp_printf("\n%d:%d:%d  %Zd = %Zd (%c) x %Zd (%c)\n",tm.tm_hour,tm.tm_min,tm.tm_sec,N,tmp,array1,tmp2,array2);
                return 0;
            }
        }
    }
    return 0;
}

int main()
{
    int flag_parallel_sieve = 1;
    int flag_batch_smooth = 1;

    printf("program started\n\n");
    srand(time(NULL));
    mpz_t N,n,b,tmp,tmp2;
    mpz_init_set_ui(n,3);
    mpz_init_set_ui(b,82939);
    mpz_powm_ui(n,n,rand()*23%103,b);
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed(state,n);
    mpf_t ln10,ln2,e;
    mpf_init(e);
    
    initialize_1(state, ln10, ln2, e);

    mpz_init(tmp);
    mpz_init(tmp2);
    mpf_t bound,tmpf,tmpf2;
    char number[100];
    printf("enter the number you wish to factor : ");
    fgets(number,100,stdin);
    mpf_init(bound);
    mpf_set_prec(bound,256);
    mpf_set_str(bound,number,10);
    mpf_init(tmpf);
    mpf_init(tmpf2);
    mpz_init_set_ui(b,10000);
    mpz_init(n);

    mpz_init_set_str(N,number,10);
    dyn_array_classic tmp_primes;
    init_classic(&tmp_primes);
    smoothB(b,&tmp_primes);
    mpf_t best_time, work, time2,time1;
    mpf_init_set_si(best_time,-1);
    mpf_init(work);
    mpf_init(time2);
    mpf_init(time1);
    mpf_set_prec(tmpf,10);
    mpf_set_prec(tmpf2,10);
    mpf_set_prec(work,10);

    unsigned long best_mult = compute_best_mult(best_time, work, time1, time2, e, ln2, tmpf, tmpf2, N, n, tmp, tmp2, tmp_primes);

    mpz_mul_ui(n,N,best_mult);

    mpf_set_z(bound,n);

    compute_b(bound, e, ln2, tmpf, tmpf2, b);

    print_bits_and_digits(n, tmp);

    mpf_sub(tmpf,best_time,time1);
    myexp(tmpf,tmpf,e);
    mpf_neg(tmpf,tmpf);
    mpf_add_ui(tmpf,tmpf,1);
    mpf_mul_ui(tmpf,tmpf,100);
    gmp_printf("multiplier used : %lu ; expected time gain = %.4Ff %%\n",best_mult,tmpf);
    printf("\n");
    mpz_t m,M;
    mpz_init(m);
    mpz_init(M);
    mpz_sqrt(m,n);
    mpz_add_ui(m,m,1);
    mpz_div_ui(M,b,1);
    unsigned long sieve_len = 2*mpz_get_ui(M)+1;

    dyn_array_classic B;
    init_classic(&B);
    smoothB(b,&B);

    dyn_array_classic primes;
    init_classic(&primes);
    append_classic(&primes,2);
    mpz_mod_ui(tmp2,n,2);
    dyn_array a;
    init(&a);
    append(&a,tmp2);
    mpz_t prod_primes;

    compute_factor_base(&primes, &a, B, &prod_primes, n, tmp2, state);

    struct tm tm = *localtime(&(time_t){time(NULL)});
    printf("%d:%d:%d\n\tfactor base created : %lu primes\n",tm.tm_hour,tm.tm_min,tm.tm_sec,primes.len);
    printf("\tp_max = %lu\n",*(primes.start+primes.len-1));

    unsigned long bounds[4] = {0,0,0,0};
    compute_bounds(bounds, primes);
    if (bounds[3] == 0) return -1;

    unsigned long logs[primes.len];
    logs[0] = 1;

    compute_logs(primes, tmp, logs);

    unsigned long dim = primes.len+1;
    dyn_array relations,smooth_numbers;
    init(&relations);
	init(&smooth_numbers);
	dyn_array partial;
    dyn_array psmooth;
    dyn_array store_partial;
    init(&partial);
    init(&psmooth);
    init(&store_partial);
	{
        unsigned long full_found = 0, partial_found = 0;
        unsigned long indexp = 0;
        int need_append = 1;
        unsigned long batch_size = 4096;

        mpf_t nb_primes1, nb_primes2, nb_large;
        mpf_init(nb_primes1);
        mpf_init(nb_primes2);
        mpf_init(nb_large);
        mpf_t var1, var2, var3;
        mpf_init(var1);
        mpf_init(var2);
        mpf_init(var3);
        unsigned long objective, seconds = 0;
        mpz_t multiplier, cst;
        
        mpz_init_set_ui(cst,*(primes.start+primes.len-1));
        mpf_set_z(tmpf,cst);
        natural_log(tmpf2,tmpf,ln2,e);
        mpf_sub_ui(tmpf2,tmpf2,1);
        mpf_div(nb_primes1,tmpf,tmpf2);
        mpz_init_set_ui(multiplier,*(primes.start+primes.len-1));
        mpz_sqrt(multiplier,multiplier);
        mpz_mul(cst,cst,multiplier);
        mpf_set_z(tmpf,cst);
        natural_log(tmpf2,tmpf,ln2,e);
        mpf_sub_ui(tmpf2,tmpf2,1);
        mpf_div(nb_primes2,tmpf,tmpf2);
        mpf_sub(nb_large,nb_primes2,nb_primes1);
        mpf_div_ui(nb_large,nb_large,2);
        gmp_printf("\tlarge prime bound = %Zd = %Zd*p_max\n",cst,multiplier);
        mpz_t tmp_bin;
        mpz_init(tmp_bin);
        unsigned long tmp_a, tmp_b, tmplol;
        mpf_t target;
        mpf_init(target);
        mpf_set_prec(target,256);
        mpf_mul_ui(target,bound,2);
        natural_log(target,target,ln2,e);
        mpf_div_ui(target,target,2);
        mpf_set_z(tmpf,M);
        natural_log(tmpf,tmpf,ln2,e);
        mpf_sub(target,target,tmpf);
        mpf_div(target,target,ln10);
        mpz_t tmp_vec2[2];
        for (unsigned long i = 0 ; i < 2 ; i++) mpz_init_set_ui(tmp_vec2[i],0);
        printf("\n");
        unsigned long prime_start = 20;
        mpf_t tmp_long1, tmp_long2;
        mpf_init(tmp_long1);
        mpf_init(tmp_long2);
        unsigned long half = mpz_get_ui(M);
        mpz_set(tmp,n);
        my_int_log2(tmp);
        mpz_add_ui(tmp,tmp,1);
        mpz_div_ui(tmp,tmp,2);
        mpz_mul_ui(tmp2,M,5);
        mpz_div_ui(tmp2,tmp2,8);
        mpz_neg(tmp2,tmp2);

        long double tmp_skipped = 10;
        for (unsigned long i = 0 ; i < prime_start ; i++)
        {
            mpf_set_ui(tmp_long1,logs[i]);
            mpf_set_ui(tmp_long2,*(primes.start+i));
            mpf_div(tmp_long1,tmp_long1,tmp_long2);
            if (my_legendre(n,*(primes.start+i)) == 1) tmp_skipped += 2*mpf_get_d(tmp_long1);
            else tmp_skipped += mpf_get_d(tmp_long1);
        }
        mpz_set(tmp,cst);
        my_int_log2(tmp);
        tmp_skipped += mpz_get_d(tmp);

        unsigned long skipped = (unsigned long) tmp_skipped;
        unsigned long addup = dim/100;
        unsigned long prime = *(primes.start);
        unsigned long smooth_bound;
        mpz_div_ui(tmp,n,2);
        mpz_sqrt(tmp,tmp);
        mpz_mul_ui(tmp,tmp,half);
        my_int_log2(tmp);
        smooth_bound = mpz_get_ui(tmp);

        time_t second1 = time(NULL);
        time_t second2;
        unsigned long time_diff;

        unsigned long time_seed = time(NULL);
        
        if (flag_parallel_sieve)
        {
            //omp_set_num_threads(4);
            printf("sieving using %d threads...\n",omp_get_max_threads());

            parallel_sieve(
                        &relations,
                        &smooth_numbers,
                        &store_partial,
                        &psmooth,
                        &partial,
                        primes,
                        a,
                        n,
                        prod_primes,
                        cst,
                        tmp_bin,
                        tmp_vec2,
                        nb_large,
                        target,
                        ln2,
                        ln10,
                        e,
                        var1,
                        var2,
                        var3,
                        tmpf2,
                        best_mult,
                        time_seed,
                        sieve_len,
                        batch_size,
                        half,
                        dim,
                        addup,
                        skipped,
                        prime_start,
                        smooth_bound,
                        prime,
                        tmp_a,
                        tmp_b,
                        tmplol,
                        time_diff,
                        objective,
                        seconds,
                        &full_found,
                        &partial_found,
                        &indexp,
                        bounds,
                        logs,
                        &need_append,
                        flag_batch_smooth,
                        second1,
                        second2
                    );
        }
        else
        {
            printf("sieving using 1 thread...\n");

            mono_cpu_sieve(
                        &relations,
                        &smooth_numbers,
                        &store_partial,
                        &psmooth,
                        &partial,
                        primes,
                        a,
                        n,
                        prod_primes,
                        cst,
                        tmp_bin,
                        tmp_vec2,
                        nb_large,
                        target,
                        ln2,
                        ln10,
                        e,
                        var1,
                        var2,
                        var3,
                        tmpf2,
                        best_mult,
                        time_seed,
                        sieve_len,
                        batch_size,
                        half,
                        dim,
                        addup,
                        skipped,
                        prime_start,
                        smooth_bound,
                        prime,
                        tmp_a,
                        tmp_b,
                        tmplol,
                        time_diff,
                        objective,
                        seconds,
                        &full_found,
                        &partial_found,
                        &indexp,
                        bounds,
                        logs,
                        &need_append,
                        flag_batch_smooth,
                        second1,
                        second2
                    );
        }
        

        printf("\r%lu/(%lu+20+%lu) relations found : full = %lu ; partial = %lu (%lu)",relations.len,dim,addup,full_found,partial_found,indexp);
        fflush(stdout);
	}
	mpz_init(tmp);
	for (unsigned long i = 0 ; i < psmooth.len ; i++) mpz_clear(*(psmooth.start+i));
    free(psmooth.start);
    psmooth.start = NULL;
	for (unsigned long i = 0 ; i < partial.len ; i++) mpz_clear(*(partial.start+i));
    free(partial.start);
    partial.start = NULL;
	for (unsigned long i = 0 ; i < store_partial.len ; i++) mpz_clear(*(store_partial.start+i));
    free(store_partial.start);
    store_partial.start = NULL;
	mpz_clear(prod_primes);
	tm = *localtime(&(time_t){time(NULL)});
    printf("\n\n%d:%d:%d  sieving done, reducing set of relations from %lu ",tm.tm_hour,tm.tm_min,tm.tm_sec,relations.len);
    reduce_relations(&relations,&smooth_numbers,&primes,n);
    printf("to %lu, building matrix...\n",relations.len);
    dyn_array_classic bin_matrix,rel_weight;
    init_classic(&bin_matrix);
    init_classic(&rel_weight);
    unsigned long nonzero;
    double density;
    unsigned long nb_lines;
    build_matrix(&bin_matrix,relations,primes,&nonzero,&density,&nb_lines,&rel_weight);
    unsigned long len = relations.len;
    tm = *localtime(&(time_t){time(NULL)});
    printf("%d:%d:%d  matrix built %lux%lu ; %lu nonzero values, density = %.2f ",tm.tm_hour,tm.tm_min,tm.tm_sec,len,nb_lines,nonzero,density);
    mpz_set_d(tmp,density);
    mpz_mul_ui(tmp,tmp,len);
    mpz_mul_ui(tmp,tmp,len);
    mpz_mul_ui(tmp,tmp,13);
    mpz_div_ui(tmp,tmp,220433617);
    mpz_div_ui(tmp2,tmp,3600);
    gmp_printf("%Zdh ",tmp2);
    mpz_mod_ui(tmp2,tmp,3600);
    mpz_div_ui(tmp2,tmp2,60);
    gmp_printf("%Zdm ",tmp2);
    mpz_mod_ui(tmp2,tmp,60);
    gmp_printf("%Zds needed\n",tmp2);
    unsigned long merge_bound = 5;
    reduce_matrix(&bin_matrix,&relations,&smooth_numbers,len,N,&rel_weight,merge_bound);
    count(bin_matrix,len,relations.len);
    bin_matrix.size = bin_matrix.len;
    bin_matrix.start = realloc(bin_matrix.start,bin_matrix.len*sizeof(unsigned long));

    return compute_factors(relations, smooth_numbers, bin_matrix, primes, N, tmp, tmp2, len);
}