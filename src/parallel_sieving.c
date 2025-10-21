#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include <math.h>

#include "logs.h"
#include "structures.h"
#include "utils.h"
#include "polynomial_functions.h"
#include "sieve.h"
#include "smooth_test.h"
#include "relations.h"

void parallel_sieve(
    FILE *logfile,
    dyn_array* relations,
    dyn_array* smooth_numbers,
    dyn_array roots,
    dyn_array_classic primes,
    mpz_t n,
    mpz_t prod_primes,
    mpz_t cst,
    mpz_t cst2,
    mpz_t tmp_bin,
    mpf_t target,
    mpf_t ln2,
    mpf_t ln10,
    mpf_t e,
    unsigned long* full_found,
    unsigned long* partial_found,
    unsigned long* indexp,
    unsigned long* bounds,
    unsigned long* logs,
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
    int* need_append,
    int nb_cpu_sieve,
    int flag_batch_smooth,
    double nb_large,
    time_t second1,
    time_t second2
)
{
    if (nb_cpu_sieve > omp_get_max_threads())
    {
        log_msg(logfile, "You don't have %d cpu available, sieving with %d cpu instead...\n", nb_cpu_sieve,  omp_get_max_threads());
        omp_set_num_threads(omp_get_max_threads());
    }
    else
    {
        log_msg(logfile, "Sieving with %d cpu...\n", nb_cpu_sieve);
        omp_set_num_threads(nb_cpu_sieve);
    }
    #pragma omp parallel
        {
            srand(time_seed^ omp_get_thread_num());
            mpz_t poly_a, poly_b, poly_c, tmp, value, poly_index, tmppolyindex;
            mpz_inits(poly_a, poly_b, poly_c, tmp, value, NULL);
            mpz_init_set_ui(poly_index,0);
            mpz_init_set_ui(tmppolyindex,1);

            dyn_array_small sieve_array;
            init_len_small(&sieve_array,sieve_len);
            memset(sieve_array.start, 0, sieve_len*sizeof(unsigned short));

            dyn_array solutions_needed, second_part, inverse_a;
            init2_len(&solutions_needed,mpf_get_ui(target));
            init2_len(&second_part,mpf_get_ui(target));
            init2_len(&inverse_a,primes.len);

            dyn_array_classic locations, tmp_where, way_to_root;
            init_classic(&locations);
            init_classic(&tmp_where);
            init_classic(&way_to_root);

            mpf_t best_bound;
            mpf_init(best_bound);

            dyn_array coefficient;
            init2_len(&coefficient, 2*batch_size);
            reset(&coefficient);

            unsigned long threshold = 0;
            unsigned long ind;

            Hashmap_PartialRelation partial_relations;
            hashmap_2d_create(&partial_relations, 16*4096);

            PartialRelation *tmp_array = calloc(batch_size, sizeof(PartialRelation));
            for (size_t i = 0 ; i < batch_size ; i++)
            {
                mpz_inits(
                    tmp_array[i].x,
                    tmp_array[i].y,
                    tmp_array[i].small_p,
                    tmp_array[i].big_p
                );
            }

            PartialRelation to_combine_node;
            mpz_inits(to_combine_node.x, to_combine_node.y, to_combine_node.small_p, to_combine_node.big_p, NULL);

            dyn_array_stack stack;
            stack_init(&stack);

            Hashmap_graph graph;
            hashmap_graph_create(&graph, 16*4096);

            Hashmap_1D parent;
            hashmap_1d_create(&parent, 16*4096);

            dyn_array to_batch;
            init2_len(&to_batch,batch_size);
            reset(&to_batch);

            dyn_array_classic is_smooth;
            init_len_classic(&is_smooth,batch_size);

            dyn_array batch_array;
            init2_len(&batch_array,2*batch_size-1);

            dyn_array block;
            init2_len(&block,batch_size);
            reset(&block);

            double rate1, rate2;

            mpf_t tmpf;
            mpf_init(tmpf);
            
            dyn_array sieved_candidates;
            init2_len(&sieved_candidates, half);
            reset(&sieved_candidates);

            mpz_t tmp_poly,tmp_poly2;
            mpz_inits(tmp_poly, tmp_poly2, NULL);

            gmp_randstate_t state;
            gmp_randinit_default(state);
            gmp_randseed_ui(state, (unsigned long)time(NULL));

            mpz_t prod_primes_p1;
            mpz_init(prod_primes_p1);
            mpz_add_ui(prod_primes_p1, prod_primes_p1, 1);

            signed long *tmp_array_sieve = calloc(sieve_len, sizeof (signed long));

            while(1)
            {
                if (relations->len >= dim+20+addup) break;

                if (mpz_sizeinbase(tmppolyindex, 2)-1 == threshold)
                {
                    create_polynomial(
                        &solutions_needed,
                        &second_part,
                        &roots,
                        &inverse_a,
                        &locations,
                        &primes,
                        &way_to_root,
                        &tmp_where,
                        poly_a,
                        n,
                        target,
                        ln2,
                        ln10,
                        best_bound,
                        e,
                        bounds,
                        best_mult
                    );

                    mpz_set_ui(poly_index, 0);
                    mpz_set_ui(tmppolyindex, 1);
                    threshold = locations.len-1;
                }
                else
                {
                    mpz_xor(tmp_poly, poly_index, tmppolyindex);
                    ind = 0;
                    while (mpz_cmp_ui(tmp_poly, 0) > 0)
                    {
                        mpz_neg(solutions_needed.start[solutions_needed.len-1-ind], solutions_needed.start[solutions_needed.len-1-ind]);
                        mpz_div_2exp(tmp_poly, tmp_poly, 1);
                        ind++;
                    }
                    mpz_add_ui(poly_index, poly_index, 1);
                    mpz_add_ui(tmppolyindex,tmppolyindex, 1);
                }
                
                CRT(&solutions_needed, &second_part, poly_a, tmp_poly);

                mpz_div_2exp(tmp_poly2, poly_a, 1);
                if (mpz_cmp(tmp_poly,tmp_poly2) > 0) mpz_sub(tmp_poly,poly_a,tmp_poly);
                mpz_set(poly_b,tmp_poly);
                mpz_mul(tmp_poly,tmp_poly,tmp_poly);
                mpz_sub(tmp_poly,tmp_poly,n);
                mpz_divexact(tmp_poly,tmp_poly,poly_a);
                mpz_set(poly_c,tmp_poly);

                sieve(
                    &sieved_candidates,
                    &roots,
                    &inverse_a,
                    &primes,
                    &way_to_root,
                    &locations,
                    &sieve_array,
                    n,
                    poly_a,
                    poly_b,
                    poly_c,
                    logs,
                    sieve_len,
                    half,
                    skipped,
                    prime_start,
                    smooth_bound,
                    tmp_array_sieve
                );

                for (size_t i = 0 ; i < sieved_candidates.len ; i++)
                {
                    mpz_set(value, sieved_candidates.start[i]);
                    append_only(&block, value);
                    mpz_mul(value, value, poly_a);
                    mpz_add(value, value, poly_b);
                    mpz_add(tmp, value, poly_b);
                    mpz_mul(tmp, tmp, sieved_candidates.start[i]);
                    mpz_add(tmp, tmp, poly_c);

                    mpz_abs(tmp, tmp);
                    append_only(&to_batch, tmp);
                    append_only(&coefficient, poly_a);
                    append_only(&coefficient, poly_b);

                    if (to_batch.len == batch_size)
                    {
                        if (flag_batch_smooth)
                        {
                            batch_smooth(
                                &to_batch,
                                &batch_array,
                                tmp_array,
                                prod_primes,
                                prod_primes_p1,
                                cst,
                                cst2,
                                prime,
                                state
                            );
                        }
                        else
                        {
                            naive_smooth(
                                &to_batch,
                                tmp_array,
                                primes,
                                cst,
                                cst2,
                                state
                            );
                        }
                        
                        #pragma omp critical
                        {
                            for (size_t k = 0 ; k < batch_size ; k++)
                            {
                                handle_relations(
                                    &partial_relations,
                                    &graph,
                                    &parent,
                                    tmp_array,
                                    &to_combine_node,
                                    &stack,
                                    relations,
                                    smooth_numbers,
                                    block,
                                    coefficient,
                                    n,
                                    value,
                                    tmp_bin,
                                    full_found,
                                    partial_found,
                                    indexp,
                                    k,
                                    tmp_a,
                                    tmp_b,
                                    tmplol,
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

                                seconds = (double)objective / (rate1 + rate2*rate2/nb_large);
                            }
                            #pragma omp master
                            {
                                printf("\r%lu/(%lu+20+%lu) relations found : full = %lu ; partial = %lu (%lu) (r1 = %.2f ; r2 = %.2f ; %luh %lum %lus left)",relations->len,dim,addup,*full_found,*partial_found,*indexp,rate1,rate2,seconds/3600,(seconds%3600)/60,(seconds%60));
                                fflush(stdout);
                            }
                        }
                    }
                }
                memset(sieve_array.start, 0, sieve_len*sizeof(unsigned short));
                reset(&sieved_candidates);
            }
            gmp_randclear(state);

            for (size_t i = 0 ; i < batch_array.len ; i++) mpz_clear(batch_array.start[i]);
            free(batch_array.start);
            batch_array.start = NULL;

            for (size_t i = 0 ; i < block.len ; i++) mpz_clear(block.start[i]);
            free(block.start);
            block.start = NULL;

            free(sieve_array.start);
            sieve_array.start = NULL;

            for (size_t i = 0 ; i < batch_size ; i++)
            {
                mpz_clears(
                    tmp_array[i].x,
                    tmp_array[i].y,
                    tmp_array[i].small_p,
                    tmp_array[i].big_p,
                    NULL
                );
            }
            free(tmp_array);

            mpz_clears(to_combine_node.x, to_combine_node.y, to_combine_node.small_p, to_combine_node.big_p, NULL);
        }
}