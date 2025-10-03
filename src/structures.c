#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>

#include "structures.h"
#include "utils.h"

void init(dyn_array* array)
{
	array->start = calloc(1, sizeof(mpz_t));
	array->len = 0;
	array->size = 1;
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
}

void append(dyn_array* array, mpz_t element)
{
    if (array->size <= array->len)
    {
        array->size <<= 1;
        array->start = realloc(array->start, sizeof(mpz_t)*(array->size));
    }
    mpz_init_set(*(array->start+array->len), element);
    array->len++;int mpz_equal(const mpz_t a, const mpz_t b);
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

void delete_classic_first(dyn_array_classic* array)
{
    array->start++;
    array->len--;
}

void delete_dyn(dyn_array* array, unsigned long index)
{
    for (unsigned long i = index ; i < array->len-1 ; i++) mpz_set(*(array->start+i), *(array->start+i+1));
    array->len--;
    mpz_clear(*(array->start + array->len));
}

void delete_dyn_unsorted(dyn_array* array, unsigned long index)
{
    array->len--;
    mpz_set(*(array->start+index), *(array->start+array->len));
    mpz_clear(*(array->start + array->len));
}

void delete_dyn_first(dyn_array* array)
{
    array->start++;
    array->len--;
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
}

void reset(dyn_array* array)
{
    array->len = 0;
}

void free_dyn_array(dyn_array* array) {
    if (!array || !array->start) return;  // safety check

    // Clear each GMP integer
    for (unsigned long i = 0; i < array->len; i++) {
        mpz_clear(array->start[i]);
    }

    // Free the memory block holding the mpz_t structs
    free(array->start);

    // Reset the struct fields
    array->start = NULL;
    array->len = 0;
    array->size = 0;
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

void hashmap_1d_create(Hashmap_1D *graph, const size_t buckets) {
    graph->buckets = buckets;
    graph->table = calloc(graph->buckets, sizeof(HashNode1D*));
    if (!graph->table) {
        exit(EXIT_FAILURE);
    }
}

size_t hash_1d_mpz_strong(const Hashmap_1D *graph, const mpz_t key) {
    return (size_t) mpz_fdiv_ui(key, graph->buckets);
}

void hashmap_1d_put(Hashmap_1D *graph, const mpz_t key, const mpz_t value) {
    size_t index = hash_mpz_strong(graph, key);
    HashNode1D *node = graph->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            mpz_set(node->value, value); // update existing
            return;
        }
        node = node->next;
    }

    // not found → insert
    node = malloc(sizeof(HashNode1D));
    mpz_init_set(node->key, key);
    mpz_init_set(node->value, value);
    node->next = graph->table[index];
    graph->table[index] = node;
}

bool hashmap_1d_get(Hashmap_1D *graph, const mpz_t key, mpz_t output) {
    size_t index = hash_mpz_strong(graph, key);
    HashNode1D *node = graph->table[index];

    while (node) {
        if (mpz_cmp(node->key, key) == 0) {
            mpz_set(output, node->value);
            return true;  // found
        }
        node = node->next;
    }

    return false; // not found
}

void hashmap_1d_free(Hashmap_1D *graph) {
    for (size_t i = 0; i < graph->buckets; i++) {
        HashNode1D *node = graph->table[i];
        while (node) {
            HashNode1D *next = node->next;
            mpz_clear(node->key);
            mpz_clear(node->value);
            free(node);
            node = next;
        }
    }
    free(graph->table);
    graph->table = NULL;
    graph->buckets = 0;
}