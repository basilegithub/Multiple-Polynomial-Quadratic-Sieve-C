#ifndef RELATIONS_H
#define RELATIONS_H

#include <gmp.h>

#include "structures.h"

bool DFS(Hashmap_graph graph, mpz_t small_p, mpz_t big_p, dyn_array *path, dyn_array_stack *stack);
void find_path(Hashmap_graph graph, mpz_t small_p, mpz_t big_p, dyn_array *path, dyn_array_stack *stack);
void combine_path(Hashmap_PartialRelation *partial_relations, dyn_array path, mpz_t x, mpz_t n, mpz_t res_x, mpz_t res_y, PartialRelation *to_combine_node);

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
);

#endif // RELATIONS_H