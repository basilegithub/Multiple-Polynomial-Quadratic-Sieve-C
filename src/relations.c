#include <gmp.h>
#include <stdio.h>
#include <stdint.h>

#include "structures.h"

bool DFS(Hashmap_graph graph, mpz_t small_p, mpz_t big_p, dyn_array *path, dyn_array_stack *stack)
{
    append_eco(path, big_p);

    dyn_array *tmp_array = hashmap_graph_get_ptr(&graph, big_p);

    for (size_t i = 0 ; i < tmp_array->len ; i++)
    {
        if (!mpz_cmp(small_p, *(tmp_array->start+i)))
        {
            append_eco(path, small_p);
            return true; // path is found
        }
    }

    stack_reset(stack);

    mpz_t next_node, parent;
    mpz_inits(next_node, parent, NULL);

    mpz_set_si(parent, -1);

    stack_push(stack, tmp_array, parent);

    while (stack->size)
    {
        while (stack->size && stack->data[stack->size-1].next_index == stack->data[stack->size-1].neighbors->len)
        {
            stack_pop(stack);
            delete_dyn_eco(path, path->len-1);
        }
        if (!stack->size) break;

        mpz_set(next_node, *(stack->data[stack->size-1].neighbors->start + stack->data[stack->size-1].next_index)); // Last element in stack
        stack->data[stack->size-1].next_index++;

        if (stack->size == 1) mpz_set(parent, big_p);
        else
        {
            mpz_set(parent, *(path->start + path->len - 2));
        }

        if (!mpz_cmp(next_node, parent)) continue;

        append_eco(path, next_node);

        tmp_array = hashmap_graph_get_ptr(&graph, next_node);
        for (size_t i = 0 ; i < tmp_array->len ; i++)
        {
            if (!mpz_cmp(small_p, *(tmp_array->start+i)))
            {
                append_eco(path, small_p);
                mpz_clears(next_node, parent, NULL);
                return true; // path is found
            }
        }

        stack_push(stack, tmp_array, parent);
    }
    mpz_clears(next_node, parent, NULL);
    return false;
}

bool find_path(Hashmap_graph graph, mpz_t small_p, mpz_t big_p, dyn_array *path, dyn_array_stack *stack)
{
    if (!mpz_cmp_ui(small_p, 1))
    {
        return DFS(graph, small_p, big_p, path, stack);
    }
    else
    {
        mpz_t tmp;
        mpz_init_set_ui(tmp, 1);

        dyn_array path_small_p_to_one, path_big_p_to_one;
        init(&path_small_p_to_one);
        init(&path_big_p_to_one);
        // We excpect the node one to be a hub, thus it is more efficient to try to connect each prime to one
        // Otherwise we would need to "get out" of the node one, which has many neighbors

        if (DFS(graph, tmp, small_p, &path_small_p_to_one, stack))
        {
            if (DFS(graph, tmp, big_p, &path_big_p_to_one, stack))
            {
                size_t index1 = path_small_p_to_one.len-1, index2 = path_big_p_to_one.len-1;

                while (!mpz_cmp(*(path_small_p_to_one.start+index1), *(path_big_p_to_one.start+index2)))
                {
                    index1--;
                    index2--;
                }

                for (size_t i = 0 ; i < index2+2 ; i++)
                {
                    append(path, *(path_big_p_to_one.start+i));
                }
                for (size_t i = index1+1 ; i-- > 0 ; )
                {
                    append(path, *(path_small_p_to_one.start+i));
                }
                mpz_clear(tmp);
                free_dyn_array(&path_small_p_to_one);
                free_dyn_array(&path_big_p_to_one);
                return true;
            }
        }
    }

    return DFS(graph, small_p, big_p, path, stack);
}

void combine_path(Hashmap_PartialRelation *partial_relations, dyn_array path, mpz_t x, mpz_t n, mpz_t res_x, mpz_t res_y, PartialRelation *to_combine_node)
{
    mpz_t tmp, tmp2, tmp3, tmp4;
    mpz_inits(tmp, tmp2, tmp3, tmp4, NULL);

    mpz_set(res_x, x);

    mpz_mul(tmp, x, x);
    mpz_sub(res_y, tmp, n);

    for (size_t i = 0 ; i < path.len-1 ; i++)
    {
        mpz_set(tmp, *(path.start+i));
        mpz_set(tmp2, *(path.start+i+1));
        if (mpz_cmp(tmp, tmp2) > 0)
        {
            mpz_swap(tmp, tmp2);        
        }

        hashmap_2d_get_from_mpz(partial_relations, tmp, tmp2, to_combine_node);

        mpz_mul(tmp3, res_y, to_combine_node->y);
        mpz_set(res_y, tmp3);

        mpz_mul(tmp3, res_x, to_combine_node->x);
        mpz_mod(res_x, tmp3, n);

        mpz_mul(tmp3,to_combine_node->x, to_combine_node->x);
        mpz_sub(tmp3, tmp3, n);
    }


    mpz_mul(tmp3, res_x, res_x);
    mpz_mod(tmp3, tmp3, n);
    mpz_mod(tmp4, res_y, n);

    for (size_t i = 0 ; i < path.len ; i++)
    {
        mpz_set(tmp3, *(path.start+i));
        if (mpz_cmp_ui(tmp3, 1))
        {
            mpz_divexact(res_y, res_y, tmp3);
            mpz_divexact(res_y, res_y, tmp3);

            mpz_invert(tmp3, tmp3, n);
            mpz_mul(tmp4, tmp3, res_x);
            mpz_mod(res_x, tmp4, n);
        }
    }

    mpz_mul(tmp3, res_x, res_x);
    mpz_mod(tmp3, tmp3, n);
    mpz_mod(tmp4, res_y, n);
    
    mpz_clears(tmp, tmp2, tmp3, tmp4, NULL);
}

void handle_relations(
    dyn_array* relations,
    dyn_array* smooth_numbers,
    PartialRelation *tmp_array,
    Hashmap_PartialRelation *partial_relations,
    PartialRelation *to_combine_node,
    Hashmap_graph *graph,
    Hashmap_1D *parent,
    dyn_array_stack *stack,
    dyn_array block,
    dyn_array coefficient,
    mpz_t n,
    mpz_t value,
    mpz_t tmp_bin,
    unsigned long k,
    unsigned long tmp_a,
    unsigned long tmp_b,
    unsigned long tmplol,
    unsigned long* full_found,
    unsigned long* partial_found,
    unsigned long* indexp,
    int* need_append
)
{
    mpz_t tmp;
    mpz_init(tmp);
    if (mpz_cmp_ui((tmp_array+k)->small_p, 0)) // If small_p is not 0, ie relation has at most two large primes
    {
        mpz_set(value, *(block.start+k));
        mpz_mul(tmp, value, *(coefficient.start+2*k));
        mpz_add(value, tmp, *(coefficient.start+2*k+1));
        if (!mpz_cmp_ui((tmp_array+k)->big_p, 1)) // If big_p is 1, ie relation is full
        {
            append(smooth_numbers, value);
            mpz_mul(tmp, value, value);
            mpz_sub(value, tmp, n);
            append(relations, value);
            (*full_found)++;
        }
        else if (mpz_cmp((tmp_array+k)->small_p, (tmp_array+k)->big_p))
        {
            if (*need_append == 1)
            {
                *need_append = 0;
                mpz_set((tmp_array+k)->x, value);
                mpz_mul(tmp,value,value);
                mpz_sub(value,tmp,n);
                mpz_set((tmp_array+k)->y, value);
                hashmap_2d_put_node(partial_relations, *(tmp_array+k));
                (*indexp)++;

                hashmap_graph_put(graph, (tmp_array+k)->small_p, (tmp_array+k)->big_p); // Add edge to graph
                hashmap_graph_put(graph, (tmp_array+k)->big_p, (tmp_array+k)->small_p); // Add edge to graph

                hashmap_1d_put(parent, (tmp_array+k)->small_p, (tmp_array+k)->small_p); // Set the parent node of small_p to itself
                hashmap_1d_put(parent, (tmp_array+k)->big_p, (tmp_array+k)->small_p); // Set the parent node of big_p to small_p
                // We expect a lot of small_p equal to one, thus it seems logical to set small_p (ie 1) as the parent node

                hashmap_1d_get(parent, (tmp_array+k)->small_p, tmp);

            } else {
                bool flag_small_prime = hashmap_graph_is_present(graph, (tmp_array+k)->small_p);
                if (!flag_small_prime) // We have never seen small_p
                {
                    mpz_set((tmp_array+k)->x, value);
                    mpz_mul(tmp,value,value);
                    mpz_sub(value,tmp,n);
                    mpz_set((tmp_array+k)->y, value);
                    hashmap_2d_put_node(partial_relations, *(tmp_array+k));
                    (*indexp)++;

                    bool flag_big_prime = hashmap_graph_is_present(graph, (tmp_array+k)->big_p);

                    hashmap_graph_put(graph, (tmp_array+k)->small_p, (tmp_array+k)->big_p); // Add edge to graph
                    hashmap_graph_put(graph, (tmp_array+k)->big_p, (tmp_array+k)->small_p); // Add edge do graph, if we have seen big_p it will append, if we have not seen it it will create a array

                    if (flag_big_prime)
                    {
                        hashmap_1d_put(parent, (tmp_array+k)->small_p, (tmp_array+k)->big_p); // Set the parent node of small_p to big_p
                    }
                    else
                    {
                        hashmap_1d_put(parent, (tmp_array+k)->small_p, (tmp_array+k)->small_p); // Set the parent node of small_p to itself
                        hashmap_1d_put(parent, (tmp_array+k)->big_p, (tmp_array+k)->small_p); // Set the parent node of big_p to small_p
                    }
                } else { // We have seen small_p
                    bool flag_big_prime = hashmap_graph_is_present(graph, (tmp_array+k)->big_p);

                    if (!flag_big_prime)
                    {
                        hashmap_graph_put(graph, (tmp_array+k)->small_p, (tmp_array+k)->big_p); // Add edge to graph
                        hashmap_graph_put(graph, (tmp_array+k)->big_p, (tmp_array+k)->small_p); // Add edge do graph

                        hashmap_1d_put(parent, (tmp_array+k)->big_p, (tmp_array+k)->small_p); // Set the parent node of big_p to small_p

                        mpz_set((tmp_array+k)->x, value);
                        mpz_mul(tmp,value,value);
                        mpz_sub(value,tmp,n);
                        mpz_set((tmp_array+k)->y, value);
                        hashmap_2d_put_node(partial_relations, *(tmp_array+k));
                        (*indexp)++;
                    }
                    else
                    {
                        mpz_t tmp2, parent_small_p, parent_big_p;
                        mpz_inits(tmp2, parent_small_p, parent_big_p, NULL);
                        mpz_set(parent_small_p, (tmp_array+k)->small_p);
                        mpz_set(parent_big_p, (tmp_array+k)->big_p);

                        hashmap_1d_get(parent, parent_small_p, tmp);

                        while (mpz_cmp(tmp, parent_small_p))
                        {
                            hashmap_1d_get(parent, tmp, tmp2); // tmp2 = parent[tmp] = parent[parent[parent_small_p]]
                            hashmap_1d_put(parent, parent_small_p, tmp2); // parent[parent_small_p] = tmp2 = parent[parent[parent_small_p]]
                            mpz_set(parent_small_p, tmp); // parent_small_p = parent[parent_small_p]
                            hashmap_1d_get(parent, parent_small_p, tmp);
                        }

                        hashmap_1d_get(parent, parent_big_p, tmp);

                        while (mpz_cmp(tmp, parent_big_p))
                        {
                            hashmap_1d_get(parent, tmp, tmp2);
                            hashmap_1d_put(parent, parent_big_p, tmp2);
                            mpz_set(parent_big_p, tmp);
                            hashmap_1d_get(parent, parent_big_p, tmp);
                        }

                        mpz_clear(tmp2);

                        if (mpz_cmp(parent_small_p, parent_big_p)) // If both primes are not part of the same connected component, ie the new edge adds no cycle
                        {
                            hashmap_graph_put(graph, (tmp_array+k)->small_p, (tmp_array+k)->big_p); // Add edge to graph
                            hashmap_graph_put(graph, (tmp_array+k)->big_p, (tmp_array+k)->small_p); // Add edge do graph
                            
                            if (mpz_cmp(parent_small_p, parent_big_p) < 0) hashmap_1d_put(parent, parent_big_p, parent_small_p);
                            else hashmap_1d_put(parent, parent_small_p, parent_big_p);

                            mpz_set((tmp_array+k)->x, value);
                            mpz_mul(tmp,value,value);
                            mpz_sub(value,tmp,n);
                            mpz_set((tmp_array+k)->y, value);
                            hashmap_2d_put_node(partial_relations, *(tmp_array+k));
                            (*indexp)++;
                        }
                        else // The two primes are part of the same connected component, one new cycle is detected, we now have to find the exact path
                        {
                            dyn_array path;
                            init(&path);

                            mpz_t res_x, res_y;
                            mpz_inits(res_x, res_y, NULL);
                            
                            find_path(*graph, (tmp_array+k)->small_p, (tmp_array+k)->big_p, &path, stack);

                            combine_path(partial_relations, path, value, n, res_x, res_y, to_combine_node);

                            append(relations, res_y);
                            append(smooth_numbers, res_x);

                            (*partial_found)++;

                            mpz_clears(res_x, res_y, NULL);
                            free_dyn_array(&path);
                        }

                        mpz_clears(parent_small_p, parent_big_p, NULL);
                    }
                }
            }
        }
    }
    mpz_clear(tmp);
}