#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <gmp.h>
#include <stdbool.h>

typedef struct
{
    mpz_t* start;
    unsigned long len;
    unsigned long size;
    unsigned long initialized;
} dyn_array; // dynamic array of large integers

typedef struct
{
    unsigned long* start;
    unsigned long len;
    unsigned long size;
} dyn_array_classic; // dynamic array of unsigned long integers

typedef struct
{
    unsigned short int* start;
    unsigned long len;
    unsigned long size;
} dyn_array_small; // dynamic array of unsigned short integers

typedef struct HashNode1D {
    mpz_t key;
    mpz_t value;
    struct HashNode1D *next;
} HashNode1D;

typedef struct Hashmap_1D {
    size_t buckets;
    HashNode1D **table;
} Hashmap_1D; // For the parent hashmap

typedef struct HashNodeGraph {
    mpz_t key;
    dyn_array value;
    struct HashNodeGraph *next;
} HashNodeGraph;

typedef struct Hashmap_graph {
    size_t buckets;
    HashNodeGraph **table;
} Hashmap_graph; // For the graph hashmap

typedef struct PartialRelation {
    mpz_t x; // x (mod n)
    mpz_t y; // x*x - n
    mpz_t small_p;
    mpz_t big_p; // y = smooth_part*small_p*big_p, small_p < big_p
} PartialRelation;

typedef struct HashNodePartialRelation {
    mpz_t key;
    PartialRelation *value;
    struct HashNodePartialRelation *next;
} HashNodePartialRelation;

typedef struct Hashmap_PartialRelation {
    size_t buckets;
    HashNodePartialRelation **table;
} Hashmap_PartialRelation;

typedef struct {
    dyn_array *neighbors;   // pointer to adjacency list in the graph
    mpz_t parent;
    size_t next_index;      // index of the next neighbor to process
} stack_node;

typedef struct {
    stack_node *data;
    size_t size;
    size_t capacity;
} dyn_array_stack; // For the stack in DFS for cycle finding

// Functions declaration

// Initialize functions

void init(dyn_array* array);
void init_classic(dyn_array_classic* array);
void init_small(dyn_array_small* array);
void init_len(dyn_array* array, unsigned long length);
void init_len_classic(dyn_array_classic* array, unsigned long length);
void init_len_small(dyn_array_small* array, unsigned long length);
void init2_len(dyn_array* array, unsigned long length);

// Append functions

void append(dyn_array* array, mpz_t element);
void append_eco(dyn_array* array, mpz_t element);
void append_only(dyn_array* array, mpz_t element);
void append_only_si(dyn_array* array, signed long element);
void append_block(dyn_array* array, unsigned long block_len, mpz_t tmp_vec[block_len]);
void append_classic(dyn_array_classic* array, unsigned long element);
void append_block_classic(dyn_array_classic* array, unsigned long block_len, unsigned long tmp_vec[block_len]);

// Delete functions

void delete_classic(dyn_array_classic* array, unsigned long index);
void delete_classic_first(dyn_array_classic* array);
void delete_dyn(dyn_array* array, unsigned long index);
void delete_dyn_eco(dyn_array* array, unsigned long index);
void delete_dyn_unsorted(dyn_array* array, unsigned long index);

// Insert functions

void insert(dyn_array* array, mpz_t element, unsigned long index);
void insert_classic(dyn_array_classic* array, unsigned long element, unsigned long index);
void insert_block(dyn_array* array, unsigned long index, unsigned long block_len, mpz_t element[block_len]);

// Liberating arrays

void reset(dyn_array* array);
void free_dyn_array(dyn_array* array);

// Query functions

int is_present(dyn_array* array, mpz_t element);
int is_present_ui(dyn_array* array, unsigned long param);

// 1D Hashmap functions

void hashmap_1d_create(Hashmap_1D *graph, const size_t buckets);
size_t hash_1d_mpz_strong(const Hashmap_1D *graph, const mpz_t key);
void hashmap_1d_put(Hashmap_1D *graph, const mpz_t key, const mpz_t value);
bool hashmap_1d_get(Hashmap_1D *graph, const mpz_t key, mpz_t output);
bool hashmap_1d_is_present(Hashmap_1D *hashmap, const mpz_t key);
void hashmap_1d_free(Hashmap_1D *graph);

// Graph hashmap functions

void hashmap_graph_create(Hashmap_graph *hashmap, const size_t buckets);
size_t hash_graph_mpz_strong(const Hashmap_graph *hashmap, const mpz_t key);
void hashmap_graph_put(Hashmap_graph *hashmap, const mpz_t key, mpz_t value);
bool hashmap_graph_is_present(Hashmap_graph *hashmap, const mpz_t key);
void hashmap_graph_get(Hashmap_graph *hashmap, const mpz_t key, dyn_array *output);
dyn_array* hashmap_graph_get_ptr(Hashmap_graph *hashmap, const mpz_t key);
void hashmap_graph_free(Hashmap_graph *hashmap);

// 2D partial_relations Hashmap functions

void hashmap_2d_create(Hashmap_PartialRelation *partial_relations, const size_t buckets);
void hashmap_2d_compute_key(const unsigned long small_p, const unsigned long big_p, mpz_t key);
void hashmap_2d_compute_key_from_mpz(const mpz_t small_p, const mpz_t big_p, mpz_t key);
size_t hash_2d_mpz_strong(const Hashmap_PartialRelation *partial_relations, const mpz_t key);
void hashmap_2d_put(Hashmap_PartialRelation *partial_relations, const unsigned long small_p, const unsigned long big_p, const PartialRelation value);
void hashmap_2d_put_node(Hashmap_PartialRelation *partial_relations, PartialRelation new_node);
bool hashmap_2d_get_from_mpz(Hashmap_PartialRelation *partial_relations, const mpz_t small_p, const mpz_t big_p, PartialRelation *output);
bool hashmap_2d_is_present_mpz(Hashmap_PartialRelation *partial_relations, const mpz_t small_p, const mpz_t big_p);
void hashmap_2d_free(Hashmap_PartialRelation *partial_relations);

// Stack functions

void stack_init(dyn_array_stack *stack);
void stack_push(dyn_array_stack *stack, const dyn_array *neighbors, const mpz_t parent_index);
void stack_pop(dyn_array_stack *stack);
void stack_reset(dyn_array_stack *stack);
void stack_free(dyn_array_stack *stack);

#endif // STRUCTURE_H