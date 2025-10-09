#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "structures.h"
#include "utils.h"

void init(dyn_array* array)
{
	array->start = calloc(1, sizeof(mpz_t));
	array->len = 0;
	array->size = 1;
    array->initialized = 0;
}

void init_classic(dyn_array_classic* array)
{
	array->start = calloc(1, sizeof(unsigned long));
	array->len = 0;
	array->size = 1;
}

void init_small(dyn_array_small* array)
{
	array->start = calloc(1, sizeof(unsigned short int));
	array->len = 0;
	array->size = 1;
}

void init_len(dyn_array* array, unsigned long length)
{
    array->size = 1;
	while (array->size <= length) array->size <<= 1;
	array->start = calloc(array->size, sizeof(mpz_t));
	array->len = length;
    array->initialized = 0;
}

void init_len_classic(dyn_array_classic* array, unsigned long length)
{
    array->size = 1;
	while (array->size <= length) array->size <<= 1;
	array->start = calloc(array->size, sizeof(unsigned long));
	array->len = length;
}

void init_len_small(dyn_array_small* array, unsigned long length)
{
    array->size = 1;
	while (array->size <= length) array->size <<= 1;
	array->start = calloc(array->size, sizeof(unsigned short int));
	array->len = length;
}

void init2_len(dyn_array* array, unsigned long length)
{
    array->size = 1;
	while (array->size <= length) array->size <<= 1;
	array->start = calloc(array->size, sizeof(mpz_t));
	array->len = length;
	for (unsigned long i = 0 ; i < length ; i++) mpz_init(*(array->start+i));
    array->initialized = length;
}

void append(dyn_array* array, mpz_t element)
{
    if (array->size <= array->len)
    {
        array->size <<= 1;
        array->start = realloc(array->start, sizeof(mpz_t)*(array->size));
    }
    mpz_init_set(*(array->start+array->len), element);
    array->len++;
    array->initialized++;
}

void append_eco(dyn_array* array, mpz_t element)
{
    if (array->size <= array->len)
    {
        array->size <<= 1;
        array->start = realloc(array->start, sizeof(mpz_t)*(array->size));
    }

    if (array->len == array->initialized)
    {
        mpz_init_set(*(array->start+array->len), element);
        array->initialized++;
    }
    else
    {
        mpz_set(*(array->start+array->len), element);
    }
    array->len++;
}

void append_only(dyn_array* array, mpz_t element)
{
    mpz_set(*(array->start+array->len), element);
    array->len++;
}

void append_only_si(dyn_array* array, signed long element)
{
    mpz_set_si(*(array->start+array->len), element);
    array->len++;
}

void append_block(dyn_array* array, unsigned long block_len, mpz_t tmp_vec[block_len])
{
    if (array->size <= array->len + block_len)
    {
        while (array->size <= array->len + block_len)
        {
            array->size <<= 1;
        }
        array->start = realloc(array->start, sizeof(mpz_t)*array->size);
    }
    for (unsigned long i = 0 ; i < block_len ; i++)
    {
        mpz_init_set(*(array->start+array->len+i), tmp_vec[i]);
    }
    array->len += block_len;
}

void append_classic(dyn_array_classic* array, unsigned long element)
{
    if (array->size <= array->len)
    {
        array->size <<= 1;
        array->start = realloc(array->start, sizeof(unsigned long)*(array->size));
    }
    *(array->start+array->len) = element;
    array->len++;
}

void append_block_classic(dyn_array_classic* array, unsigned long block_len, unsigned long tmp_vec[block_len])
{
    if (array->size <= array->len+block_len)
    {
        while (array->size <= array->len+block_len)
        {
            array->size <<= 1;
        }
        array->start = realloc(array->start, sizeof(unsigned long)*array->size);
    }
    for (unsigned long i = 0 ; i < block_len ; i++)
    {
        *(array->start+array->len+i) = tmp_vec[i];
    }
    array->len += block_len;
}

void delete_classic(dyn_array_classic* array, unsigned long index)
{
    for (unsigned long i = index ; i < array->len-1 ; i++) *(array->start+i) = *(array->start+i+1);
    array->len--;
}

void delete_dyn(dyn_array* array, unsigned long index)
{
    for (unsigned long i = index ; i < array->len-1 ; i++) mpz_set(*(array->start+i), *(array->start+i+1));
    array->len--;
    mpz_clear(*(array->start + array->len));
    array->initialized--;
}

void delete_dyn_eco(dyn_array* array, unsigned long index)
{
    for (unsigned long i = index ; i < array->len-1 ; i++) mpz_set(*(array->start+i), *(array->start+i+1));
    array->len--;
}

void delete_dyn_unsorted(dyn_array* array, unsigned long index)
{
    array->len--;
    mpz_set(*(array->start+index), *(array->start+array->len));
    mpz_clear(*(array->start + array->len));
    array->initialized--;
}

void insert(dyn_array* array, mpz_t element, unsigned long index)
{
    if (array->size <= array->len)
    {
        array->size <<= 1;
        array->start = realloc(array->start, sizeof(mpz_t)*(array->size));
    }
    mpz_init(*(array->start + array->len));
    for (unsigned long i = array->len ; i > index ; i--)
    {
        mpz_set(*(array->start+i), *(array->start+i-1));
    }
    mpz_set(*(array->start+index), element);
    array->len++;
    array->initialized++;
}

void insert_classic(dyn_array_classic* array, unsigned long element, unsigned long index)
{
    if (array->size <= array->len)
    {
        array->size <<= 1;
        array->start = realloc(array->start, sizeof(unsigned long)*(array->size));
    }
    for (unsigned long i = array->len ; i > index ; i--)
    {
        *(array->start + i) = *(array->start+i-1);
    }
    *(array->start + index) = element;
    array->len++;
}

void insert_block(dyn_array* array, unsigned long index, unsigned long block_len, mpz_t element[block_len])
{
    if (array->size <= array->len + block_len)
    {
        while (array->size <= array->len + block_len)
        {
            array->size <<= 1;
        }
        array->start = realloc(array->start, sizeof(mpz_t)*array->size);
    }
    for (unsigned long i = array->len+block_len-1 ; i >= array->len ; i--)
    {
        mpz_init_set(*(array->start + i), *(array->start+i-block_len));
    }
    for (unsigned long i = array->len-1 ; i >= index + block_len ; i--)
    {
        mpz_set(*(array->start+i), *(array->start+i-block_len));
    }
    for (unsigned long i = 0 ; i < block_len ; i++)
    {
        mpz_set(*(array->start+index+i), element[i]);
    }
    array->len += block_len;
    array->initialized++;
}

void reset(dyn_array* array)
{
    array->len = 0;
}

void free_dyn_array(dyn_array* array) {
    if (!array || !array->start) return;  // safety check

    // Clear each GMP integer
    for (unsigned long i = 0; i < array->initialized; i++) {
        mpz_clear(array->start[i]);
    }

    // Free the memory block holding the mpz_t structs
    free(array->start);

    // Reset the struct fields
    array->start = NULL;
    array->len = 0;
    array->size = 0;
    array->initialized = 0;
}

int is_present(dyn_array* array, mpz_t element)
{
    if (array->len == 0) return 0;
    unsigned long a = 0;
    unsigned long b = array->len-1;
    unsigned long tmp = (a+b)/2;
    while (mpz_cmp(*(array->start + tmp), element) != 0 && a <= b)
    {
        if (mpz_cmp(*(array->start + tmp), element) < 0) a = tmp + 1;
        else
        {
            if (tmp == 0) return 0;
            b = tmp - 1;
        }
        tmp = (a+b)/2;
    }
    if (mpz_cmp(*(array->start + tmp), element) == 0) return 1;
    return 0;
}

int is_present_ui(dyn_array* array, unsigned long param)
{
    mpz_t element;
    mpz_init_set_ui(element, param);
    if (array->len == 0) return 0;
    unsigned long a = 0;
    unsigned long b = array->len-1;
    unsigned long tmp = (a+b)/2;
    while (mpz_cmp(*(array->start + tmp), element) != 0 && a <= b)
    {
        if (mpz_cmp(*(array->start + tmp), element) < 0) a = tmp + 1;
        else
        {
            if (tmp == 0) return 0;
            b = tmp - 1;
        }
        tmp = (a+b)/2;
    }
    if (mpz_cmp(*(array->start + tmp), element) == 0) return 1;
    return 0;
}

// 1D Hashmap functions

void hashmap_1d_create(Hashmap_1D *hashmap, const size_t buckets)
{
    hashmap->buckets = buckets;
    hashmap->table = calloc(hashmap->buckets, sizeof(HashNode1D*));
    if (!hashmap->table) {
        exit(EXIT_FAILURE);
    }
}

size_t hash_1d_mpz_strong(const Hashmap_1D *hashmap, const mpz_t key) {
    return (size_t) mpz_fdiv_ui(key, hashmap->buckets);
}

void hashmap_1d_put(Hashmap_1D *hashmap, const mpz_t key, const mpz_t value) {
    size_t index = hash_1d_mpz_strong(hashmap, key);
    HashNode1D *node = hashmap->table[index];

    while (node) {
        if (!mpz_cmp(node->key, key)) {
            mpz_set(node->value, value); // update existing
            return;
        }
        node = node->next;
    }

    // not found → insert
    node = malloc(sizeof(HashNode1D));
    mpz_init_set(node->key, key);
    mpz_init_set(node->value, value);
    node->next = hashmap->table[index];
    hashmap->table[index] = node;
}

bool hashmap_1d_get(Hashmap_1D *hashmap, const mpz_t key, mpz_t output) {
    size_t index = hash_1d_mpz_strong(hashmap, key);
    HashNode1D *node = hashmap->table[index];

    while (node) {
        if (!mpz_cmp(node->key, key)) {
            mpz_set(output, node->value);
            return true;  // found
        }
        node = node->next;
    }

    return false; // not found
}

bool hashmap_1d_is_present(Hashmap_1D *hashmap, const mpz_t key) {
    size_t index = hash_1d_mpz_strong(hashmap, key);
    HashNode1D *node = hashmap->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            return true;  // found
        }
        node = node->next;
    }

    return false; // not found
}

void hashmap_1d_free(Hashmap_1D *hashmap) {
    for (size_t i = 0; i < hashmap->buckets; i++) {
        HashNode1D *node = hashmap->table[i];
        while (node) {
            HashNode1D *next = node->next;
            mpz_clear(node->key);
            mpz_clear(node->value);
            free(node);
            node = next;
        }
    }
    free(hashmap->table);
    hashmap->table = NULL;
    hashmap->buckets = 0;
}

// Graph hashmap functions

void hashmap_graph_create(Hashmap_graph *hashmap, const size_t buckets)
{
    hashmap->buckets = buckets;
    hashmap->table = calloc(hashmap->buckets, sizeof(HashNodeGraph*));
    if (!hashmap->table) {
        exit(EXIT_FAILURE);
    }
}

size_t hash_graph_mpz_strong(const Hashmap_graph *hashmap, const mpz_t key)
{
    return (size_t) mpz_fdiv_ui(key, hashmap->buckets);
}

void hashmap_graph_put(Hashmap_graph *hashmap, const mpz_t key, mpz_t value)
{
    size_t index = hash_graph_mpz_strong(hashmap, key);
    HashNodeGraph *node = hashmap->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            append(&node->value, value); // Append to existing dynamic array
            return;
        }
        node = node->next;
    }

    // not found → create a new dynamic array and put the value in it
    node = malloc(sizeof(HashNodeGraph));
    init(&node->value);
    mpz_init_set(node->key, key);
    append(&node->value, value);
    node->next = hashmap->table[index];
    hashmap->table[index] = node;
}

bool hashmap_graph_is_present(Hashmap_graph *hashmap, const mpz_t key)
{
    size_t index = hash_graph_mpz_strong(hashmap, key);
    HashNodeGraph *node = hashmap->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            return true;  // found
        }
        node = node->next;
    }

    return false; // not found
}

void hashmap_graph_get(Hashmap_graph *hashmap, const mpz_t key, dyn_array *output)
{
    size_t index = hash_graph_mpz_strong(hashmap, key);
    HashNodeGraph *node = hashmap->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            output = &node->value;
            return;
        }
        node = node->next;
    }
    return;  // not found
}

dyn_array* hashmap_graph_get_ptr(Hashmap_graph *hashmap, const mpz_t key)
{
    size_t index = hash_graph_mpz_strong(hashmap, key);
    HashNodeGraph *node = hashmap->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            return &node->value;  // return pointer to internal array
        }
        node = node->next;
    }
    return NULL;  // not found
}

void hashmap_graph_free(Hashmap_graph *hashmap) {
    for (size_t i = 0; i < hashmap->buckets; i++) {
        HashNodeGraph *node = hashmap->table[i];
        while (node) {
            HashNodeGraph *next = node->next;
            mpz_clear(node->key);
            free_dyn_array(&node->value);
            free(node);
            node = next;
        }
    }
    free(hashmap->table);
    hashmap->table = NULL;
    hashmap->buckets = 0;
}

// 2D Hashmap functions

void hashmap_2d_create(Hashmap_PartialRelation *partial_relations, const size_t buckets)
{
    partial_relations->buckets = buckets;
    partial_relations->table = calloc(partial_relations->buckets, sizeof(HashNodePartialRelation*));
    if (!partial_relations->table) {
        exit(EXIT_FAILURE);
    }
}

void hashmap_2d_compute_key(const unsigned long small_p, const unsigned long big_p, mpz_t key)
{
    mpz_t tmp;
    mpz_init_set_ui(tmp, big_p);
    size_t bits = mpz_sizeinbase(tmp, 2);    
    mpz_set_ui(tmp, small_p);  // bits required
    mpz_mul_2exp(key, tmp, bits);            // key = small_p << bits
    mpz_add_ui(key, key, big_p);                    // key += big_p
}

void hashmap_2d_compute_key_from_mpz(const mpz_t small_p, const mpz_t big_p, mpz_t key)
{
    size_t bits = mpz_sizeinbase(big_p, 2);      // bits required
    mpz_mul_2exp(key, small_p, bits);            // key = small_p << bits
    mpz_add(key, key, big_p);                    // key += big_p
}

size_t hash_2d_mpz_strong(const Hashmap_PartialRelation *partial_relations, const mpz_t key) {
    return (size_t) mpz_fdiv_ui(key, partial_relations->buckets);
}

void hashmap_2d_put(Hashmap_PartialRelation *partial_relations, const unsigned long small_p, const unsigned long big_p, const PartialRelation value)
{
    mpz_t key;
    mpz_init(key);
    hashmap_2d_compute_key(small_p, big_p, key);
    size_t index = hash_2d_mpz_strong(partial_relations, key);

    HashNodePartialRelation *node = partial_relations->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            mpz_set(node->value->x, value.x); // update existing
            mpz_set(node->value->y, value.y); // update existing
            mpz_set(node->value->small_p, value.small_p); // update existing
            mpz_set(node->value->big_p, value.big_p); // update existing
            mpz_clear(key);
            return;
        }
        node = node->next;
    }

    // not found → insert
    node = malloc(sizeof(HashNodePartialRelation));
    node->value = malloc(sizeof(PartialRelation));
    mpz_init_set(node->key, key);

    mpz_init_set(node->value->x, value.x);
    mpz_init_set(node->value->y, value.y);
    mpz_init_set(node->value->small_p, value.small_p);
    mpz_init_set(node->value->big_p, value.big_p);

    node->next = partial_relations->table[index];
    partial_relations->table[index] = node;

    mpz_clear(key);
}

void hashmap_2d_put_node(Hashmap_PartialRelation *partial_relations, PartialRelation new_node)
{
    mpz_t key;
    mpz_init(key);
    hashmap_2d_compute_key_from_mpz(new_node.small_p, new_node.big_p, key);
    size_t index = hash_2d_mpz_strong(partial_relations, key);

    HashNodePartialRelation *node = partial_relations->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            mpz_set(node->value->x, new_node.x); // update existing
            mpz_set(node->value->y, new_node.y); // update existing
            mpz_set(node->value->small_p, new_node.small_p); // update existing
            mpz_set(node->value->big_p, new_node.big_p); // update existing
            mpz_clear(key);
            return;
        }
        node = node->next;
    }

    // not found → insert
    node = malloc(sizeof(HashNodePartialRelation));
    node->value = malloc(sizeof(PartialRelation));
    mpz_init_set(node->key, key);

    mpz_init_set(node->value->x, new_node.x);
    mpz_init_set(node->value->y, new_node.y);
    mpz_init_set(node->value->small_p, new_node.small_p);
    mpz_init_set(node->value->big_p, new_node.big_p);

    node->next = partial_relations->table[index];
    partial_relations->table[index] = node;

    mpz_clear(key);
}

bool hashmap_2d_get_from_mpz(Hashmap_PartialRelation *partial_relations, const mpz_t small_p, const mpz_t big_p, PartialRelation *output)
{
    mpz_t key;
    mpz_init(key);
    hashmap_2d_compute_key_from_mpz(small_p, big_p, key);
    size_t index = hash_2d_mpz_strong(partial_relations, key);
    
    HashNodePartialRelation *node = partial_relations->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            mpz_set(output->x, node->value->x);
            mpz_set(output->y, node->value->y);
            mpz_set(output->small_p, node->value->small_p);
            mpz_set(output->big_p, node->value->big_p);
            mpz_clear(key);
            return true;  // found
        }
        node = node->next;
    }

    mpz_clear(key);
    return false; // not found
}

bool hashmap_2d_is_present_mpz(Hashmap_PartialRelation *partial_relations, const mpz_t small_p, const mpz_t big_p)
{
     mpz_t key;
    mpz_init(key);
    hashmap_2d_compute_key_from_mpz(small_p, big_p, key);
    size_t index = hash_2d_mpz_strong(partial_relations, key);
    
    HashNodePartialRelation *node = partial_relations->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            mpz_clear(key);
            return true;  // found
        }
        node = node->next;
    }

    mpz_clear(key);

    return false; // not found
}

void hashmap_2d_free(Hashmap_PartialRelation *partial_relations) {
    for (size_t i = 0; i < partial_relations->buckets; i++) {
        HashNodePartialRelation *node = partial_relations->table[i];
        while (node) {
            HashNodePartialRelation *next = node->next;
            mpz_clear(node->key);
            mpz_clear(node->value->x);
            mpz_clear(node->value->y);
            mpz_clear(node->value->small_p);
            mpz_clear(node->value->big_p);
            free(node);
            node = next;
        }
    }
    free(partial_relations->table);
    partial_relations->table = NULL;
    partial_relations->buckets = 0;
}

// Stack functions

void stack_init(dyn_array_stack *stack)
{
    stack->data = malloc(sizeof(stack_node) * 4);  // initial capacity
    stack->size = 0;
    stack->capacity = 4;
}

void stack_push(dyn_array_stack *stack, const dyn_array *neighbors, const mpz_t parent)
{
    if (stack->size == stack->capacity) {
        stack->capacity *= 2;
        stack->data = realloc(stack->data, sizeof(stack_node) * stack->capacity);
    }

    stack->data[stack->size].neighbors = neighbors;
    stack->data[stack->size].next_index = 0;  // start at first neighbor
    mpz_init_set(stack->data[stack->size].parent, parent); // store parent index to skip
    stack->size++;
}

void stack_pop(dyn_array_stack *stack)
{
    mpz_clear(stack->data[--stack->size].parent);
}

void stack_reset(dyn_array_stack *stack)
{
    for (size_t i = 0; i < stack->size; i++)
        mpz_clear(stack->data[i].parent);

    stack->size = 0;
}

void stack_free(dyn_array_stack *stack)
{
    for (size_t i = 0; i < stack->size; i++)
        mpz_clear(stack->data[i].parent);

    free(stack->data);
    stack->data = NULL;
    stack->size = 0;
    stack->capacity = 0;
}