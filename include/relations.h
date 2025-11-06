#ifndef RELATIONS_H
#define RELATIONS_H

#include <gmp.h>

#include "structures.h"

bool DFS(
    Hashmap_graph graph,
    dyn_array_stack *stack,
    dyn_array *path,
    mpz_t small_p,
    mpz_t big_p
);

void find_path(
    Hashmap_graph graph,
    dyn_array_stack *stack,
    dyn_array *path,
    mpz_t small_p, mpz_t big_p
);

void combine_path(
    Hashmap_PartialRelation *partial_relations,
    PartialRelation *to_combine_node,
    dyn_array path,
    mpz_t x,
    mpz_t n,
    mpz_t res_x,
    mpz_t res_y
);

void handle_relations(
    Hashmap_PartialRelation *partial_relations,
    Hashmap_graph *graph,
    Hashmap_1D *parent,
    PartialRelation *tmp_array,
    PartialRelation *to_combine_node,
    dyn_array_stack *stack,
    dyn_array* relations,
    dyn_array* smooth_numbers,
    dyn_array block,
    dyn_array coefficient,
    mpz_t n,
    mpz_t value,
    unsigned long* full_found,
    unsigned long* partial_found,
    unsigned long* indexp,
    unsigned long k,
    int* need_append
);

#endif // RELATIONS_H