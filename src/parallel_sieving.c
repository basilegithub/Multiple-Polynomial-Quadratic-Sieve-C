#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

#include "structures.h"
#include "utils.h"
#include "polynomial_functions.h"
#include "sieve.h"
#include "smooth_test.h"
#include "relations.h"

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
)
{
    if (nb_cpu_sieve > omp_get_max_threads())
    {
        printf("You don't have %d cpu available, sieving with %d cpu instead...\n", nb_cpu_sieve,  omp_get_max_threads());
    }
    else
    {
        omp_set_num_threads(nb_cpu_sieve);
        printf("Sieving with %d cpu ...\n", nb_cpu_sieve);
    }
    #pragma omp parallel
        {
            srand(time_seed^ omp_get_thread_num());
            mpz_t poly_a, poly_b, poly_c;
            mpz_init(poly_a);
            mpz_init(poly_b);
            mpz_init(poly_c);
            mpz_t tmp;
            mpz_init(tmp);

            dyn_array_small sieve_array;
            init_len_small(&sieve_array,sieve_len);
            for (unsigned long i = 0 ; i < sieve_array.len ; i++) *(sieve_array.start+i) = 0;

            dyn_array solutions_needed, second_part, inverse_a;
            dyn_array_classic locations, tmp_where, way_to_root;
            init_classic(&locations);
            init_classic(&tmp_where);
            init_classic(&way_to_root);
            init2_len(&solutions_needed,mpf_get_ui(target));
            init2_len(&second_part,mpf_get_ui(target));
            init2_len(&inverse_a,primes.len);

            mpf_t best_bound;
            mpf_init(best_bound);
            mpz_t value;
            mpz_init(value);

            mpz_t poly_index, tmppolyindex;
            mpz_init_set_ui(poly_index,0);
            mpz_init_set_ui(tmppolyindex,1);

            dyn_array coefficient;
            init2_len(&coefficient,2*batch_size);
            reset(&coefficient);

            unsigned long threshold = 0;
            unsigned long ind;

            dyn_array to_batch, large_primes;
            dyn_array_classic is_smooth;
            init2_len(&to_batch,batch_size);
            reset(&to_batch);
            init2_len(&large_primes,batch_size);
            for (unsigned long i = 0 ; i < batch_size ; i++) mpz_set_ui(*(large_primes.start+i),1);
            init_len_classic(&is_smooth,batch_size);
            for (unsigned long i = 0 ; i < batch_size ; i++) *(is_smooth.start+i) = 0;
            dyn_array batch_array;
            init2_len(&batch_array,2*batch_size-1);
            dyn_array block;
            init2_len(&block,batch_size);
            reset(&block);
            double rate1, rate2;
            mpf_t tmpf;
            mpf_init(tmpf);
            dyn_array tmp_block;
            init2_len(&tmp_block,half);
            reset(&tmp_block);
            mpz_t tmp_poly,tmp_poly2;
            mpz_init(tmp_poly);
            mpz_init(tmp_poly2);
            while(1)
            {
                if (relations->len >= dim+20+addup) break;

                if (mpz_sizeinbase(tmppolyindex,2)-1 == threshold)
                {
                    create_polynomial(poly_a,&solutions_needed,&second_part,&locations,n,&primes,&a,bounds,target,ln2,ln10,&inverse_a,&way_to_root,best_bound,&tmp_where,e,best_mult);
                    mpz_set_ui(poly_index,0);
                    mpz_set_ui(tmppolyindex,1);
                    threshold = locations.len-1;
                } else
                {
                    mpz_xor(tmp_poly,poly_index,tmppolyindex);
                    ind = 0;
                    while (mpz_cmp_ui(tmp_poly,0) > 0)
                    {
                        mpz_neg(*(solutions_needed.start+solutions_needed.len-1-ind),*(solutions_needed.start+solutions_needed.len-1-ind));
                        mpz_div_ui(tmp_poly,tmp_poly,2);
                        ind++;
                    }
                    mpz_add_ui(poly_index,poly_index,1);
                    mpz_add_ui(tmppolyindex,tmppolyindex,1);
                }
                CRT(tmp_poly,&solutions_needed,poly_a,&second_part);
                mpz_div_ui(tmp_poly2,poly_a,2);
                if (mpz_cmp(tmp_poly,tmp_poly2) > 0) mpz_sub(tmp_poly,poly_a,tmp_poly);
                mpz_set(poly_b,tmp_poly);
                mpz_mul(tmp_poly,tmp_poly,tmp_poly);
                mpz_sub(tmp_poly,tmp_poly,n);
                mpz_divexact(tmp_poly,tmp_poly,poly_a);
                mpz_set(poly_c,tmp_poly);
                sieve(&sieve_array,sieve_len,half,&primes,logs,&a,n,poly_a,poly_b,poly_c,&inverse_a,&way_to_root,&locations,skipped,prime_start,&tmp_block,smooth_bound);
                for (unsigned long i = 0 ; i < tmp_block.len ; i++)
                {
                    mpz_set(value,*(tmp_block.start+i));
                    append_only(&block,value);
                    mpz_mul(value,value,poly_a);
                    mpz_add(value,value,poly_b);
                    mpz_add(tmp,value,poly_b);
                    mpz_mul(tmp,tmp,*(tmp_block.start+i));
                    mpz_add(tmp,tmp,poly_c);

                    mpz_abs(tmp,tmp);
                    append_only(&to_batch,tmp);
                    append_only(&coefficient,poly_a);
                    append_only(&coefficient,poly_b);
                    if (to_batch.len == batch_size)
                    {
                        if (flag_batch_smooth)
                        {
                            batch_smooth(&to_batch,&large_primes,&is_smooth,&batch_array,prod_primes,cst,prime);
                        }
                        else
                        {
                            naive_smooth(&to_batch, &large_primes, &is_smooth, primes, cst);
                        }
                        
                        #pragma omp critical
                        {
                            for (unsigned long k = 0 ; k < batch_size ; k++)
                            {
                                handle_relations(
                                    relations,
                                    smooth_numbers,
                                    store_partial,
                                    psmooth,
                                    partial,
                                    large_primes,
                                    block,
                                    coefficient,
                                    is_smooth,
                                    n,
                                    value,
                                    tmp_bin,
                                    tmp_vec2,
                                    k,
                                    tmp_a,
                                    tmp_b,
                                    tmplol,
                                    full_found,
                                    partial_found,
                                    indexp,
                                    need_append
                                );
                            }

                            reset(&block);
                            reset(&to_batch);
                            reset(&coefficient);
                            second2 = time(NULL);
                            time_diff = second2-second1;
                            if (time_diff == 0)
                            {
                                rate1 = (double)(*full_found);
                                rate2 = (double)(*partial_found + *indexp);
                            }
                            else
                            {
                                rate1 = (double)(*full_found)/time_diff;
                                rate2 = (double)(*partial_found + *indexp)/time_diff;
                            }
                            if (dim+20 >= relations->len)
                            {
                                objective = dim+20-relations->len;
                                mpf_set_d(tmpf,rate2);
                                mpf_set_ui(var1,objective);
                                mpf_mul(var1,var1,nb_large);
                                mpf_mul(var1,var1,tmpf);
                                mpf_mul_ui(var1,var1,2);

                                mpf_set_d(var2,rate1);
                                mpf_mul(var2,var2,nb_large);

                                mpf_set_ui(tmpf,partial->len);
                                mpf_div(tmpf,tmpf,nb_large);
                                mpf_set_d(tmpf2,rate2);
                                mpf_mul(tmpf,tmpf,tmpf2);
                                mpf_add(var2,var2,tmpf);
                                mpf_mul(tmpf,var2,var2);
                                mpf_add(tmpf,tmpf,var1);
                                mpf_sqrt(tmpf,tmpf);
                                mpf_sub(tmpf,tmpf,var2);
                                mpf_set_d(var3,rate2);
                                mpf_div(tmpf,tmpf,var3);
                                seconds = mpf_get_ui(tmpf);
                            }
                            #pragma omp master
                            {
                                printf("\r%lu/(%lu+20+%lu) relations found : full = %lu ; partial = %lu (%lu) (r1 = %.2f ; r2 = %.2f ; %luh %lum %lus left)",relations->len,dim,addup,*full_found,*partial_found,*indexp,rate1,rate2,seconds/3600,(seconds%3600)/60,(seconds%60));
                                fflush(stdout);
                            }
                        }
                    }
                }
                for (unsigned long i = 0 ; i < sieve_array.len ; i++) *(sieve_array.start+i) = 0;
                reset(&tmp_block);
            }
            for (unsigned long i = 0 ; i < batch_array.len ; i++) mpz_clear(*(batch_array.start+i));
            free(batch_array.start);
            batch_array.start = NULL;
            for (unsigned long i = 0 ; i < block.len ; i++) mpz_clear(*(block.start+i));
            free(block.start);
            block.start = NULL;
            free(sieve_array.start);
            sieve_array.start = NULL;
        }
}