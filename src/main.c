#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gmp.h>

#include "logs.h"
#include "config.h"
#include "structures.h"
#include "utils.h"
#include "generate_primes.h"
#include "mono_cpu_sieving.h"
#include "parallel_sieving.h"
#include "build_matrix.h"
#include "reduce_matrix.h"
#include "gaussian_elimination.h"
#include "wiedemann.h"
#include "block_lanczos.h"
#include "compute_sqrt.h"

void count(FILE *logfile, dyn_array_classic matrix, unsigned long limit, unsigned long dim)
{
    unsigned long nb_lines = 0;
    unsigned long nonzero = 0;

    for (unsigned long i = 0 ; i < matrix.len ; i++)
    {
        if (matrix.start[i] == limit) nb_lines++;
        else nonzero++;
    }

    double density = (double) nonzero/nb_lines;
    struct tm tm = *localtime(&(time_t){time(NULL)});
    log_msg(logfile, "matrix reduced to %lux%lu ; %lu nonzero values, density = %.2f", dim, nb_lines, nonzero, density);
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

    mpz_clears(n, b, NULL);
}

unsigned long compute_best_mult(mpf_t best_time, mpf_t work, mpf_t time1, mpf_t time2, mpf_t e, mpf_t ln2, mpf_t tmpf, mpf_t tmpf2, mpz_t N, mpz_t n, mpz_t tmp, mpz_t tmp2, dyn_array_classic tmp_primes)
{
    int res;
    unsigned long best_mult;

    for (unsigned long k = 1 ; k < 40 ; k++)
    {
        mpz_mul_ui(n, N, k);
        mpf_set_z(tmpf, n);
        natural_log(tmpf, tmpf, ln2, e);
        natural_log(tmpf2, tmpf, ln2, e);
        mpf_mul(tmpf, tmpf, tmpf2);
        mpf_sqrt(tmpf, tmpf);
        mpf_div_2exp(work, tmpf, 1);

        mpz_mod_ui(tmp2, n, 8); // tackling n%(2^k)
        if (!mpz_cmp_ui(tmp2, 1))
        {
            mpf_mul_2exp(tmpf2, ln2, 1);
            mpf_sub(work, work, tmpf2);
        }
        else if (!mpz_cmp_ui(tmp2, 5)) mpf_sub(work, work, ln2);
        else
        {
            mpz_mod_ui(tmp2, n, 4);
            if (!mpz_cmp_ui(tmp2, 3))
            {
                mpf_div_2exp(tmpf2, ln2, 1);
                mpf_sub(work, work, tmpf2);
            }
        }

        for (size_t i = 1 ; i < tmp_primes.len ; i++)
        {
            mpz_set_ui(tmp, tmp_primes.start[i]);
            res = mpz_legendre(n, tmp);
            if (res == 1)
            {
                mpf_set_z(tmpf, tmp);
                mpf_set(tmpf2, tmpf);
                natural_log(tmpf, tmpf, ln2, e);
                mpf_sub_ui(tmpf2, tmpf2, 1);
                mpf_mul_2exp(tmpf, tmpf, 1);
                mpf_div(tmpf2, tmpf, tmpf2);
                mpf_sub(work, work, tmpf2);
            }
            else if (!res)
            {
                mpf_set_z(tmpf, tmp);
                mpf_set(tmpf2, tmpf);
                natural_log(tmpf, tmpf, ln2, e);
                mpf_sub_ui(tmpf2, tmpf2, 1);
                mpf_div(tmpf2, tmpf, tmpf2);
                mpf_sub(work, work, tmpf2);
            }
        }

        // Compute estimated time
        mpf_mul_2exp(work, work, 1);
        myexp(work, work, e);
        mpf_set_z(tmpf, n);
        mpf_mul(work, work, tmpf);
        natural_log(time2, work, ln2, e);
        natural_log(work, time2, ln2, e);
        mpf_mul(time2, time2, work);
        mpf_sqrt(time2, time2);

        if (mpf_cmp_si(best_time, -1) == 0 || mpf_cmp(time2, best_time) < 0)
        {
            if (mpf_cmp_si(best_time, -1) == 0)
            {
                mpf_set(time1, time2);
            }
            mpf_set(best_time, time2);
            best_mult = k;
        }
    }

    return best_mult;
}

void compute_b(mpf_t bound, mpf_t e, mpf_t ln2, mpf_t tmpf, mpf_t tmpf2, mpz_t b)
{
    natural_log(tmpf, bound, ln2, e);
    natural_log(tmpf2, tmpf, ln2, e);
    mpf_mul(tmpf, tmpf, tmpf2);
    mpf_sqrt(tmpf, tmpf);
    mpf_div_2exp(tmpf, tmpf, 1);
    myexp(tmpf, tmpf, e);
    mpf_mul_ui(tmpf, tmpf, 10);
    mpf_div_ui(tmpf, tmpf, 80);
    mpz_set_f(b, tmpf);
}

void print_bits_and_digits(FILE *logfile, mpz_t n)
{
    log_msg(logfile, "%lu bits ; %lu digits", mpz_sizeinbase(n, 2), mpz_sizeinbase(n, 10));
}

void compute_factor_base(dyn_array_classic* primes, dyn_array* roots, dyn_array_classic B, mpz_t* prod_primes, mpz_t n, mpz_t tmp2, gmp_randstate_t state)
{
    mpz_init_set_ui(*prod_primes, 2);

    for (size_t i = 1 ; i < B.len ; i++)
    {
        if (my_legendre(n, B.start[i]) == 1)
        {
            append_classic(primes, B.start[i]);
            mpz_set(tmp2, n);

            sqrt_mod(tmp2, B.start[i], state);
            append(roots, tmp2);
            mpz_mul_ui(*prod_primes, *prod_primes, B.start[i]);
        }
        else if (!my_legendre(n, B.start[i]))
        {
            append_classic(primes, B.start[i]);
            mpz_set_ui(tmp2, 0);
            append(roots, tmp2);
            mpz_mul_ui(*prod_primes, *prod_primes, B.start[i]);
        }
    }
}

void compute_bounds(unsigned long* bounds, dyn_array_classic primes)
{
    for (size_t i = 1 ; i < primes.len ; i++)
    {
        if (primes.start[i] >= 10 && bounds[0] == 0) bounds[0] = i;
        else if (primes.start[i] >= 100 && bounds[1] == 0) bounds[1] = i;
        else if (primes.start[i] >= 1000 && bounds[2] == 0) bounds[2] = i;
        else if (primes.start[i] >= 3000 && bounds[3] == 0)
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
    mpz_t last, last_log;
    mpz_init_set_ui(last, 2);
    mpz_init_set_ui(last_log, 1);

    for (size_t i = 1 ; i < primes.len ; i++)
    {
        mpz_set(tmp, last);
        mpz_mul_ui(tmp, tmp, 1414213562);
        mpz_div_ui(tmp, tmp, 1000000000); // tmp ~ last*sqrt2
        if (mpz_cmp_ui(tmp, primes.start[i]) >= 0)
        {
            logs[i] = mpz_get_ui(last_log);
            continue;
        }
        mpz_mul_2exp(tmp, last, 1);
        if (mpz_cmp_ui(tmp, primes.start[i]) >= 0)
        {
            mpz_add_ui(tmp, last_log, 1);
            logs[i] = mpz_get_ui(tmp);
            continue;
        }
        while (mpz_cmp_ui(tmp, primes.start[i]) == -1)
        {
            mpz_add_ui(last_log, last_log, 1);
            mpz_mul_2exp(last, last, 1);
            mpz_mul_2exp(tmp, last, 1);
        }
        logs[i] = mpz_get_ui(last_log);
    }

    mpz_clears(last, last_log, NULL);
}

void convert_to_vec(mpz_t embedding, unsigned long relations_len, bool tmp_vec[relations_len])
{
    mpz_t tmp, tmp2;
    mpz_inits(tmp, tmp2, NULL);
    mpz_set_ui(tmp2, 1);

    for (size_t i = 0 ; i < relations_len ; i++)
    {
        mpz_and(tmp, embedding, tmp2);
        tmp_vec[relations_len-i-1] = (bool) mpz_get_ui(tmp);

        mpz_div_2exp(embedding, embedding, 1);
    }

    mpz_clears(tmp, tmp2, NULL);
}

void compute_factors(FILE *logfile, dyn_array relations, dyn_array smooth_numbers, dyn_array_classic bin_matrix, dyn_array_classic primes, mpz_t N, mpz_t tmp, mpz_t tmp2, unsigned long len, size_t block_size)
{
    struct tm tm;

    mpz_t x, y, minimal_polynomial_estimation;
    mpz_inits(x, y, NULL);
    mpz_init_set_ui(minimal_polynomial_estimation, 1);

    unsigned long k;

    size_t block_len;

    if (block_size <= 16) block_len = block_size;
    else block_len = 16;

    bool *kernel_vec = calloc(relations.len, sizeof(bool));

    while (1)
    {
        dyn_array kernel_vectors;
        init(&kernel_vectors);

        log_msg(logfile, "new block");

        wiedemann(&kernel_vectors, bin_matrix, minimal_polynomial_estimation, block_len, relations.len, len);

        log_msg(logfile, "kernel_size = %lu", kernel_vectors.len);

        for (size_t i = 0 ; i < kernel_vectors.len ; i++)
        {
            convert_to_vec(kernel_vectors.start[i], relations.len, kernel_vec);

            mpz_set_ui(x, 1);
            mpz_set_ui(y, 1);

            build_sqrt(relations, smooth_numbers, primes, N, x, y, relations.len, kernel_vec);

            mpz_sub(tmp, x, y);
            mpz_add(tmp2, x, y);
            mpz_gcd(tmp, tmp, N);
            mpz_gcd(tmp2, tmp2, N);
            char array1;
            char array2;

            if (mpz_cmp_ui(tmp, 1) != 0 && mpz_cmp(tmp, N) != 0)
            {
                if (mpz_probab_prime_p(tmp, 100) > 0)
                {
                    array1 = 'p';
                } else {array1 = 'C';}

                if (mpz_probab_prime_p(tmp2, 100) > 0)
                {
                    array2 = 'p';
                } else {array2 = 'C';}

                log_blank_line(logfile);
                log_gmp_msg(logfile, "%Zd = %Zd (%c) x %Zd (%c)", N, tmp, array1, tmp2, array2);

                return;
            }
        }
    }
}

int main()
{
    printf("program started\n\n");

    FILE *logfile = NULL;

    init_log(&logfile);

    char* config_path = "./config/config.txt";
    int nb_cpu_sieve;
    int flag_batch_smooth;
    int flag_gaussian_elimination;
    int flag_block_lanczos;
    size_t block_size;

    parse_config(config_path, &nb_cpu_sieve, &flag_batch_smooth, &flag_gaussian_elimination, &flag_block_lanczos, &block_size);

    srand(time(NULL));
    mpz_t N, n, b, tmp, tmp2;
    mpz_init_set_ui(n, 3);
    mpz_init_set_ui(b, 82939);
    mpz_powm_ui(n, n, rand()*23%103,b);

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed(state, n);

    mpf_t ln10, ln2, e, bound, tmpf, tmpf2;
    mpf_inits(e, bound, tmpf, tmpf2, NULL);

    mpf_set_prec(bound, 256);
    mpf_set_prec(tmpf, 10);
    mpf_set_prec(tmpf2, 10);
    
    initialize_1(state, ln10, ln2, e);

    char number[100];
    printf("enter the number you wish to factor : ");
    fgets(number, 100, stdin);

    mpz_inits(tmp, tmp2, n, NULL);
    mpz_init_set_ui(b,10000);

    mpz_init_set_str(N, number, 10);
    mpf_set_str(bound, number, 10);

    dyn_array_classic tmp_primes;
    init_classic(&tmp_primes);

    erasthotenes_sieve(&tmp_primes, b);

    mpf_t best_time, work, time2, time1;
    mpf_init_set_si(best_time, -1);
    mpf_inits(work, time1, time2, NULL);
    mpf_set_prec(work, 10);

    unsigned long best_mult = compute_best_mult(best_time, work, time1, time2, e, ln2, tmpf, tmpf2, N, n, tmp, tmp2, tmp_primes);

    mpz_mul_ui(n, N, best_mult);

    mpf_set_z(bound, n);

    compute_b(bound, e, ln2, tmpf, tmpf2, b);

    print_bits_and_digits(logfile, n);

    mpf_sub(tmpf, best_time, time1);
    myexp(tmpf, tmpf, e);
    mpf_neg(tmpf, tmpf);
    mpf_add_ui(tmpf, tmpf, 1);
    mpf_mul_ui(tmpf, tmpf, 100);
    log_gmp_msg(logfile, "multiplier used : %lu ; expected time gain = %.4Ff %%", best_mult, tmpf);
    log_blank_line(logfile);

    mpz_t m, M;
    mpz_inits(m, M, NULL);

    mpz_sqrt(m, n);
    mpz_add_ui(m, m, 1);
    mpz_div_ui(M, b, 1);
    unsigned long sieve_len = 2*mpz_get_ui(M)+1;

    dyn_array_classic B;
    init_classic(&B);
    erasthotenes_sieve(&B, b);

    dyn_array_classic primes;
    init_classic(&primes);
    append_classic(&primes, 2);

    mpz_mod_ui(tmp2, n, 2);
    dyn_array roots;
    init(&roots);
    append(&roots, tmp2);

    mpz_t prod_primes;

    compute_factor_base(&primes, &roots, B, &prod_primes, n, tmp2, state);

    struct tm tm = *localtime(&(time_t){time(NULL)});
    log_msg(logfile, "Factor base created : %lu primes", primes.len);
    log_msg(logfile, "p_max = %lu", primes.start[primes.len-1]);

    unsigned long bounds[4] = {0, 0, 0, 0};
    compute_bounds(bounds, primes);
    if (bounds[3] == 0) return -1;

    unsigned long logs[primes.len];
    logs[0] = 1;

    compute_logs(primes, tmp, logs);

    unsigned long dim = primes.len+1;
    dyn_array relations, smooth_numbers;
    init(&relations);
	init(&smooth_numbers);

	{
        unsigned long full_found = 0, partial_found = 0;
        unsigned long indexp = 0;
        int need_append = 1;
        unsigned long batch_size = 4096;

        mpf_t nb_primes1, nb_primes2, nb_large, var1, var2, var3;
        mpf_inits(nb_primes1, nb_primes2, nb_large, var1, var2, var3, NULL);

        unsigned long objective, seconds = 0;
        mpz_t multiplier, cst, cst2;
        
        mpz_init_set_ui(cst, primes.start[primes.len-1]);
        mpf_set_z(tmpf, cst);
        natural_log(tmpf2, tmpf, ln2, e);
        mpf_sub_ui(tmpf2, tmpf2, 1);
        mpf_div(nb_primes1, tmpf, tmpf2);

        mpz_init_set_ui(multiplier, primes.start[primes.len-1]);
        mpz_sqrt(multiplier, multiplier);
      
        mpz_mul(cst, cst, multiplier);
        mpz_init_set(cst2, cst);
        mpz_mul_ui(tmp2, cst2, primes.start[primes.len-1]);
        mpz_set(cst2, tmp2);
        
        mpf_set_z(tmpf, cst);
        natural_log(tmpf2, tmpf, ln2, e);
        mpf_sub_ui(tmpf2, tmpf2, 1);
        mpf_div(nb_primes2, tmpf, tmpf2);

        mpf_sub(nb_large, nb_primes2, nb_primes1);
        mpf_div_2exp(nb_large, nb_large, 1);
        mpf_mul(nb_large, nb_large, nb_large);

        log_gmp_msg(logfile, "Large prime bound 1 = %Zd = %Zd*p_max", cst, multiplier);
        log_gmp_msg(logfile, "Large prime bound 2 = %Zd = %Zd*p_max^2", cst2, multiplier);
        log_blank_line(logfile);

        mpz_t tmp_bin;
        mpz_init(tmp_bin);

        unsigned long tmp_a, tmp_b, tmplol;

        mpf_t target;
        mpf_init(target);
        mpf_set_prec(target, 256);

        mpf_mul_2exp(target, bound, 1);
        natural_log(target, target, ln2, e);
        mpf_div_2exp(target, target, 1);
        mpf_set_z(tmpf, M);
        natural_log(tmpf, tmpf, ln2, e);
        mpf_sub(target, target, tmpf);
        mpf_div(target, target, ln10);

        unsigned long prime_start = 20;

        mpf_t tmp_long1, tmp_long2;
        mpf_inits(tmp_long1, tmp_long2, NULL);

        unsigned long half = mpz_get_ui(M);
        mpz_set(tmp, n);
        my_int_log2(tmp);
        mpz_add_ui(tmp, tmp, 1);
        mpz_div_2exp(tmp, tmp, 1);
        mpz_mul_ui(tmp2, M, 5);
        mpz_div_2exp(tmp2, tmp2, 3);
        mpz_neg(tmp2, tmp2);

        long double tmp_skipped = 0;
        for (unsigned long i = 0 ; i < prime_start ; i++) // We add the contribution of the skipped small primes
        {
            mpf_set_ui(tmp_long1, logs[i]);
            mpf_set_ui(tmp_long2, primes.start[i]);
            mpf_div(tmp_long1, tmp_long1, tmp_long2);
            if (my_legendre(n, primes.start[i]) == 1) tmp_skipped += 2*mpf_get_d(tmp_long1);
            else tmp_skipped += mpf_get_d(tmp_long1);
        }
        mpz_set(tmp, cst2);
        my_int_log2(tmp);
        tmp_skipped += mpz_get_d(tmp);

        unsigned long skipped = (unsigned long) tmp_skipped;
        unsigned long addup = dim/100;
        unsigned long prime = *primes.start;
        unsigned long smooth_bound;

        smooth_bound = (mpz_sizeinbase(n, 2) - 2) / 2 + log2_ui(half);

        time_t second1 = time(NULL);
        time_t second2;
        unsigned long time_diff;

        unsigned long time_seed = time(NULL);
        
        if (nb_cpu_sieve > 1)
        {
            parallel_sieve(
                        logfile,
                        &relations,
                        &smooth_numbers,
                        roots,
                        primes,
                        n,
                        prod_primes,
                        cst,
                        cst2,
                        tmp_bin,
                        target,
                        ln2,
                        ln10,
                        e,
                        &full_found,
                        &partial_found,
                        &indexp,
                        bounds,
                        logs,
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
                        &need_append,
                        nb_cpu_sieve,
                        flag_batch_smooth,
                        mpf_get_d(nb_large),
                        second1,
                        second2
                    );
        }
        else
        {
            log_msg(logfile, "Sieving using 1 cpu...");

            mono_cpu_sieve(
                        &relations,
                        &smooth_numbers,
                        roots,
                        primes,
                        n,
                        prod_primes,
                        cst,
                        cst2,
                        tmp_bin,
                        target,
                        ln2,
                        ln10,
                        e,
                        &full_found,
                        &partial_found,
                        &indexp,
                        bounds,
                        logs,
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
                        &need_append,
                        flag_batch_smooth,
                        mpf_get_d(nb_large),
                        second1,
                        second2
                    );
        }
        
        printf("\r%lu/(%lu+20+%lu) relations found : full = %lu ; partial = %lu (%lu)",
            relations.len,
            dim,
            addup,
            full_found,
            partial_found,
            indexp
        );
        fflush(stdout);

        printf("\r\033[K");
        fflush(stdout);

        printf("\r");
        fflush(stdout);

        log_msg(logfile, "%lu/(%lu+20+%lu) relations found : full = %lu ; partial = %lu (%lu)",
            relations.len,
            dim,
            addup,
            full_found,
            partial_found,
            indexp
        );
	}

	mpz_init(tmp);

	mpz_clear(prod_primes);
	tm = *localtime(&(time_t){time(NULL)});

    unsigned long relations_len = relations.len;
    reduce_relations(&relations, &smooth_numbers, &primes,n);
    log_blank_line(logfile);
    log_msg(logfile, "Sieving done, reduced set of relations from %lu to %lu, building matrix...", relations_len, relations.len);

    if (flag_gaussian_elimination)
    {

        mpz_t x, y;
        mpz_inits(x, y, NULL);

        unsigned long base_size = primes.len + 1;
        unsigned long relations_len = relations.len;

        mpz_t *dense_matrix = calloc(relations_len, sizeof(mpz_t));
        for (size_t i = 0 ; i < relations_len ; i++) mpz_init(dense_matrix[i]);

        mpz_t *res = calloc(relations_len, sizeof(mpz_t));
        for (size_t i = 0 ; i < relations_len ; i++) mpz_init(res[i]);

        bool tmp_vec[relations_len];

        build_dense_matrix(relations, primes, dense_matrix, relations.len, primes.len+1);

        log_msg(logfile, "matrix built %lux%lu ; performing gaussian elimination...", relations.len, primes.len+1);

        gaussian_elimination(dense_matrix, res, relations.len, primes.len+1);

        for (unsigned long i = 0 ; i < relations.len ; i++)
        {
            if (row_is_zero(dense_matrix, i, relations.len, primes.len+1))
            {
                mpz_set_ui(x, 1);
                mpz_set_ui(y, 1);

                convert_to_vec(res[i], relations_len, tmp_vec);

                build_sqrt(relations, smooth_numbers, primes, N, x, y, relations.len, tmp_vec);

                mpz_sub(tmp, x, y);
                mpz_add(tmp2, x, y);
                mpz_gcd(tmp, tmp, N);
                mpz_gcd(tmp2, tmp2, N);
                char array1;
                char array2;

                if (mpz_cmp_ui(tmp, 1) != 0 && mpz_cmp(tmp, N) != 0)
                {
                    if (mpz_probab_prime_p(tmp, 100) > 0)
                    {
                        array1 = 'p';
                    } else {array1 = 'C';}

                    if (mpz_probab_prime_p(tmp2, 100) > 0)
                    {
                        array2 = 'p';
                    } else {array2 = 'C';}

                    tm = *localtime(&(time_t){time(NULL)});
                    log_blank_line(logfile);
                    log_gmp_msg(logfile, "%Zd = %Zd (%c) x %Zd (%c)", N, tmp, array1, tmp2, array2);
                    if (logfile) fclose(logfile);
                    return 1;
                }
            }
        }
    }

    dyn_array_classic bin_matrix, rel_weight;
    init_classic(&bin_matrix);
    init_classic(&rel_weight);

    unsigned long nonzero;
    double density;
    unsigned long nb_lines;

    build_sparse_matrix(relations, &bin_matrix, &rel_weight, primes, &nonzero, &nb_lines, &density);

    unsigned long len = relations.len;

    tm = *localtime(&(time_t){time(NULL)});
    log_msg(logfile, "matrix built %lux%lu ; %lu nonzero values, density = %.2f", len, nb_lines, nonzero, density);

    unsigned long merge_bound = 5;

    reduce_matrix(&relations, &smooth_numbers, &bin_matrix, &rel_weight, N, len, merge_bound);

    count(logfile, bin_matrix, len, relations.len);

    log_blank_line(logfile);

    bin_matrix.size = bin_matrix.len;
    bin_matrix.start = realloc(bin_matrix.start, bin_matrix.len*sizeof(unsigned long));

    if (!flag_block_lanczos)
    {
        compute_factors(logfile, relations, smooth_numbers, bin_matrix, primes, N, tmp, tmp2, len, block_size);
    }

    else
    {
        mpz_t x, y;
        mpz_inits(x, y, NULL);

        while (true)
        {
            dyn_array kernel_vectors;
            init(&kernel_vectors);

            block_lanczos(&kernel_vectors, bin_matrix, relations.len, block_size, len, logfile);
            if (block_size < 16) block_size <<= 1;

            bool *kernel_vec = calloc(relations.len, sizeof(bool));

            for (size_t i = 0 ; i < kernel_vectors.len ; i++)
            {
                convert_to_vec(kernel_vectors.start[i], relations.len, kernel_vec);

                build_sqrt(relations, smooth_numbers, primes, N, x, y, relations.len, kernel_vec);

                mpz_sub(tmp, x, y);
                mpz_add(tmp2, x, y);
                mpz_gcd(tmp, tmp, N);
                mpz_gcd(tmp2, tmp2, N);
                char array1;
                char array2;

                if (mpz_cmp_ui(tmp, 1) != 0 && mpz_cmp(tmp, N) != 0)
                {
                    if (mpz_probab_prime_p(tmp, 100) > 0)
                    {
                        array1 = 'p';
                    } else {array1 = 'C';}

                    if (mpz_probab_prime_p(tmp2, 100) > 0)
                    {
                        array2 = 'p';
                    } else {array2 = 'C';}

                    tm = *localtime(&(time_t){time(NULL)});
                    log_blank_line(logfile);
                    log_gmp_msg(logfile, "%Zd = %Zd (%c) x %Zd (%c)", N, tmp, array1, tmp2, array2);
                    return 1;
                }
            }

            free_dyn_array(&kernel_vectors);
        }
    }

    if (logfile) fclose(logfile);

    // dyn_array_classic A;
    // init_classic(&A);

    // append_classic(&A, 0);
    // append_classic(&A, 1);
    // append_classic(&A, 2);
    // append_classic(&A, 5);
    // append_classic(&A, 0);
    // append_classic(&A, 1);
    // append_classic(&A, 2);
    // append_classic(&A, 5);

    // dyn_array kernel_vectors;
    // init(&kernel_vectors);

    // mpz_t poly;
    // mpz_init_set_ui(poly, 1);

    // size_t block_size = 4;

    // unsigned long n = 3;

    // unsigned long limit = 5;

    // wiedemann(&kernel_vectors, A, poly, block_size, n, limit);

    // for (unsigned long i = 0 ; i < kernel_vectors.len ; i++)
    // {
    //     gmp_printf("kernel vector %lu = %Zd\n", i, kernel_vectors.start[i]);
    // }
}