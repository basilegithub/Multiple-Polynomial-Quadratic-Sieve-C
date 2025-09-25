#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <gmp.h>

typedef struct
{
    mpz_t* start;
    unsigned long len;
    unsigned long size;
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
void append_only(dyn_array* array, mpz_t element);
void append_only_si(dyn_array* array, signed long element);
void append_block(dyn_array* array, unsigned long block_len, mpz_t tmp_vec[block_len]);
void append_classic(dyn_array_classic* array, unsigned long element);
void append_block_classic(dyn_array_classic* array, unsigned long block_len, unsigned long tmp_vec[block_len]);

// Delete functions

void delete_classic(dyn_array_classic* array, unsigned long index);
void delete_classic_first(dyn_array_classic* array);
void delete_dyn(dyn_array* array, unsigned long index);
void delete_dyn_unsorted(dyn_array* array, unsigned long index);
void delete_dyn_first(dyn_array* array);

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

#endif // STRUCTURE_H