#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>

#include "structures.h"

void reduce_relations(dyn_array* relations, dyn_array* smooth, dyn_array_classic* primes, mpz_t n)
{
    unsigned long count;
    unsigned long bit_count;
    size_t index = 0, index2 = 0;
    bool modified = true;

    mpz_t tmp_mpz, tmp_mpz2;
    mpz_inits(tmp_mpz, tmp_mpz2, NULL);

    while (modified)
    {
        modified = false;
        for (size_t i = 0 ; i < primes->len ; i++)
        {

            count = 0;
            for (size_t j = 0 ; j < relations->len && count < 3 ; j++)
            {
                mpz_set_ui(tmp_mpz, primes->start[i]);
                bit_count = mpz_remove(tmp_mpz2, relations->start[j], tmp_mpz);
                if (bit_count&1)
                {
                    count++;
                    if (count == 1) index = j;
                    else index2 = j;
                }
            }

            if (count == 1)
            {
                delete_dyn_unsorted(relations, index);
                delete_dyn_unsorted(smooth, index);
                modified = true;
            }
            else if (count == 2)
            {
                mpz_mul(relations->start[index], relations->start[index], relations->start[index2]);
                mpz_mul(smooth->start[index], smooth->start[index], smooth->start[index2]);
                mpz_mod(smooth->start[index], smooth->start[index], n);

                delete_dyn_unsorted(relations, index2);
                delete_dyn_unsorted(smooth, index2);
                modified = true;
            }
        }
    }

    mpz_clears(tmp_mpz, tmp_mpz2, NULL);
}

void bubble_sort_down(dyn_array_classic *sorted, dyn_array_classic weights, unsigned long tmp3)
{
    unsigned long tmp4;
    while (tmp3 > 0 && weights.start[sorted->start[tmp3]] < weights.start[sorted->start[tmp3-1]])
    {
        tmp4 = sorted->start[tmp3-1];
        sorted->start[tmp3-1] = sorted->start[tmp3];
        sorted->start[tmp3] = tmp4;
        tmp3--;
    }
}

void bubble_sort_up(dyn_array_classic *sorted, dyn_array_classic weights, unsigned long tmp3)
{
    unsigned long tmp4;
    while (tmp3 < sorted->len-1 && weights.start[sorted->start[tmp3]] > weights.start[sorted->start[tmp3+1]])
    {
        tmp4 = sorted->start[tmp3+1];
        sorted->start[tmp3+1] = sorted->start[tmp3];
        sorted->start[tmp3] = tmp4;
        tmp3++;
    }
}

size_t get_index(dyn_array_classic sorted, unsigned long line_index)
{
    for (size_t z = 0 ; z < sorted.len ; z++)
    {
        if (sorted.start[z] == line_index)
        {
            return z;
        }
    }
}

bool delete_empty_row(dyn_array_classic *matrix, dyn_array_classic *weights, dyn_array_classic *sorted, size_t i, unsigned long row_index)
{
    delete_classic(matrix, i);
    delete_classic(weights, row_index);
    size_t k = 0;
    while (k < sorted->len)
    {
        if (sorted->start[k] == row_index) delete_classic(sorted, k);
        else
        {
            if (sorted->start[k] > row_index) sorted->start[k]--;
            k++;
        }
    }
    return true;
}

void delete_singleton(
    dyn_array *relations,
    dyn_array *smooth,
    dyn_array_classic *matrix,
    dyn_array_classic *weights, 
    dyn_array_classic *rel_weight,
    dyn_array_classic sorted,
    size_t i,
    unsigned long row_index,
    unsigned long row_delimiter
)
{
    unsigned long tmp = matrix->start[i-1];
    delete_classic(matrix, i);
    delete_classic(matrix, i-1);
    delete_dyn(relations, tmp);
    delete_dyn(smooth, tmp);
    delete_classic(rel_weight, tmp);
    delete_classic(weights, row_index);
    unsigned long tmp2 = 0;
    size_t j = 0;

    while (j < matrix->len)
    {
        if (matrix->start[j] == row_delimiter)
        {
            tmp2++;
            j++;
        }
        else if (matrix->start[j] == tmp)
        {
            delete_classic(matrix, j);
            weights->start[tmp2]--;

            unsigned long tmp3 = get_index(sorted, tmp2);

            bubble_sort_down(&sorted,*weights, tmp3); // update the sorted list of line weights
        }
        else
        {
            if (matrix->start[j] > tmp) matrix->start[j]--;
            j++;
        }
    }

    size_t k = 0;
    while (k < sorted.len)
    {
        if (sorted.start[k] == row_index) delete_classic(&sorted, k);
        else
        {
            if (sorted.start[k] > row_index) sorted.start[k]--;
            k++;
        }
    }
}

void delete_two_elements_row(
    dyn_array *relations,
    dyn_array *smooth,
    dyn_array_classic *matrix,
    dyn_array_classic *rel_weight,
    dyn_array_classic *weights,
    dyn_array_classic sorted,
    mpz_t n,
    size_t i,
    unsigned long row_index,
    unsigned long row_delimiter
)
{
    unsigned long pivot, remain, tmp2, tmp3;
    size_t j;

    pivot = matrix->start[i-2];
    remain = matrix->start[i-1];

    mpz_mul(relations->start[remain], relations->start[remain], relations->start[pivot]);
    mpz_mul(smooth->start[remain], smooth->start[remain], smooth->start[pivot]);
    mpz_mod(smooth->start[remain], smooth->start[remain], n);

    delete_classic(matrix, i-2);
    delete_classic(matrix, i-2);
    delete_classic(matrix, i-2);
    rel_weight->start[remain]--;

    delete_dyn(relations, pivot);
    delete_dyn(smooth, pivot);
    delete_classic(rel_weight, pivot);
    delete_classic(weights, row_index);

    tmp2 = 0;
    j = 0;
    bool has_pivot = false;
    while (j < matrix->len)
    {
        if (matrix->start[j] == row_delimiter)
        {
            if (has_pivot)
            {
                has_pivot = false;
                insert_classic(matrix, remain-1, j);
                rel_weight->start[remain-1]++;
                weights->start[tmp2]++;

                tmp3 = get_index(sorted, tmp2);

                bubble_sort_up(&sorted, *weights, tmp3);
                j++;
            }
            tmp2++;
            j++;
        }

        else if (matrix->start[j] == pivot)
        {
            has_pivot = true;
            delete_classic(matrix, j);
            weights->start[tmp2]--;

            tmp3 = get_index(sorted, tmp2);

            bubble_sort_down(&sorted,*weights, tmp3); // update the sorted list of line weights
        }

        else if (matrix->start[j] == remain && has_pivot)
        {
            has_pivot = false;
            delete_classic(matrix, j);
            rel_weight->start[remain-1]--;
            weights->start[tmp2]--;

            tmp3 = get_index(sorted, tmp2);

            bubble_sort_down(&sorted, *weights, tmp3);
        }
        else if (matrix->start[j] > remain && has_pivot)
        {
            has_pivot = false;
            insert_classic(matrix, remain-1, j);
            rel_weight->start[remain-1]++;
            weights->start[tmp2]++;

            tmp3 = get_index(sorted, tmp2);

            bubble_sort_up(&sorted, *weights, tmp3); // update the sorted list of line weights
            j++;
        }
        else
        {
            if (matrix->start[j] > pivot) matrix->start[j]--;
            j++;
        }
    }

    size_t k = 0;
    while (k < sorted.len)
    {
        if (sorted.start[k] == row_index) delete_classic(&sorted, k);
        else
        {
            if (sorted.start[k] > row_index) sorted.start[k]--;
            k++;
        }
    }
}

/*
 * reduce_matrix:
 *   Simplifies a sparse matrix representing relations.
 *   - Removes empty and singleton rows.
 *   - Fuses 2-element rows into other rows.
 *   - Maintains a sorted list of row weights.
 *   - Merges lightweight rows below 'merge_bound'.
 */
void reduce_matrix(
    dyn_array* relations,
    dyn_array* smooth,
    dyn_array_classic* matrix,
    dyn_array_classic* rel_weight,
    mpz_t n,
    unsigned long row_delimiter,
    unsigned long merge_bound
)
{
    dyn_array_classic weights;
    init_classic(&weights);

    dyn_array_classic sorted;
    init_classic(&sorted);

    unsigned long line_len = 0;
    unsigned long row_index = 0;
    unsigned long tmp, tmp2, tmp3 = 0;
    signed long tmp_index;

    for (size_t i = 0 ; i < matrix->len ; i++)
    {
        if (matrix->start[i] == row_delimiter)
        {
            append_classic(&weights, line_len);
            if (weights.len == 1) append_classic(&sorted, row_index);
            else
            {
                signed long a = 0;
                signed long b = row_index-1;
                tmp_index = (a+b)/2;
                tmp2 = sorted.start[tmp_index];
                while (a <= b && weights.start[tmp2] != line_len)
                {
                    if (weights.start[tmp2] < line_len) a = tmp_index+1;
                    else b = tmp_index-1;
                    tmp_index = (a+b)/2;
                    tmp2 = sorted.start[tmp_index];
                }
                if (weights.start[tmp2] == line_len) insert_classic(&sorted, row_index, tmp_index);
                else insert_classic(&sorted, row_index, a);
            }
            line_len = 0;
            row_index++;
        }
        else line_len++;
    }

    // weights are sorted
    size_t i, j;
    unsigned long k;
    unsigned long nb_lines;
    bool changed = true;
    unsigned long pivot;

    while(changed)
    {
        row_index = 0;
        nb_lines = 0;
        changed = false;
        i = 0;
        line_len = 0;
        while (i < matrix->len)
        {
            if (matrix->start[i] == row_delimiter)
            {
                if (line_len == 0) // delete empty row
                {
                    changed = delete_empty_row(matrix, &weights, &sorted, i, row_index);
                }

                else if (line_len == 1) // delete singleton
                {
                    delete_singleton(relations, smooth, matrix, &weights, rel_weight, sorted, i, row_index, row_delimiter);

                    line_len = 0;
                    i = 0;
                    nb_lines = 0;
                    changed = true;
                    row_index = 0;
                }

                else if (line_len == 2) // fuse row with only two elements
                {
                    delete_two_elements_row(relations, smooth, matrix, rel_weight, &weights, sorted, n, i, row_index, row_delimiter);

                    line_len = 0;
                    i = 0;
                    nb_lines = 0;
                    changed = true;
                    row_index = 0;
                }
                else
                {
                    i++;
                    line_len = 0;
                    nb_lines++;
                    row_index++;
                }
            }
            else
            {
                i++;
                line_len++;
            }
        }

        while (relations->len > nb_lines+15)
        {
            changed = true;
            row_index = 0;
            i = 0;

            while (row_index != sorted.start[0])
            {
                if (matrix->start[i] == row_delimiter) row_index++;
                i++;
            }

            if (matrix->start[i] == row_delimiter)
            {
                delete_classic(matrix, i);
                nb_lines--;
                delete_classic(&weights, row_index);

                for (size_t z = 1 ; z < sorted.len ; z++)
                {
                    if (sorted.start[z] > sorted.start[0]) sorted.start[z]--;
                }
                delete_classic(&sorted, 0);
            }

            else
            {
                weights.start[row_index]--;
                tmp2 = matrix->start[i];
                delete_classic(matrix, i);
                delete_dyn(relations, tmp2);
                delete_dyn(smooth, tmp2);
                delete_classic(rel_weight, tmp2);
                tmp = 0;
                k = 0;

                while (k < matrix->len)
                {
                    if (matrix->start[k] == row_delimiter)
                    {
                        tmp++;
                        k++;
                    }

                    else if (matrix->start[k] == tmp2)
                    {
                        delete_classic(matrix, k);
                        weights.start[tmp] -= 1;

                        tmp3 = get_index(sorted, tmp);

                        bubble_sort_down(&sorted, weights, tmp3); // update the sorted list of line weights
                    }
                    else
                    {
                        if (matrix->start[k] > tmp2) matrix->start[k]--;
                        k++;
                    }
                }
            }
        }
    }

    unsigned long line_index, a, b, A, B, TMP, M;
    dyn_array_classic remainer;
    init_classic(&remainer);

    while (weights.start[sorted.start[0]] <= merge_bound)
    {
        i = 0;
        tmp2 = 0;

        while (tmp2 != sorted.start[0]){
            if (matrix->start[i] == row_delimiter) tmp2++;
            i++;
        }

        if (weights.start[sorted.start[0]] > 0)
        {
            pivot = matrix->start[i];

            for (unsigned long j = 1 ; j < weights.start[tmp2] ; j++)
            {
                tmp = matrix->start[i+j];
                if (rel_weight->start[tmp] < rel_weight->start[pivot]) pivot = tmp;
            }

            remainer.len = 0;
            for (unsigned long j = 0 ; j < weights.start[tmp2] ; j++)
            {
                if (matrix->start[i+j] != pivot)
                {
                    tmp = matrix->start[i+j];
                    append_classic(&remainer, tmp);
                    mpz_mul(relations->start[tmp], relations->start[tmp], relations->start[pivot]);
                    mpz_mul(smooth->start[tmp], smooth->start[tmp], smooth->start[pivot]);
                    mpz_mod(smooth->start[tmp], smooth->start[tmp], n);
                }
            }

            delete_dyn(relations, pivot);
            delete_dyn(smooth, pivot);
            line_index = 0;
            j = 0;

            while (line_index < weights.len)
            {
                a = 0;
                b = weights.start[line_index]-1;
                tmp = (a+b)>>1;
                while (a <= b && b < weights.start[line_index] && matrix->start[j+tmp] != pivot)
                {
                    if (matrix->start[j+tmp] < pivot) a = tmp+1;
                    else b = tmp-1;
                    tmp = (a+b)>>1;
                }

                if (a > b || b >= weights.start[line_index] || matrix->start[j+tmp] != pivot)
                {
                    j += a;
                    while (matrix->start[j] != row_delimiter)
                    {
                        if (matrix->start[j] > pivot) matrix->start[j]--;
                        j++;
                    }
                }

                else
                {
                    delete_classic(matrix, j+tmp);
                    weights.start[line_index]--;

                    for (size_t m = 0 ; m < remainer.len ; m++)
                    {
                        M = remainer.start[m];
                        A = 0;
                        B = weights.start[line_index]-1;
                        TMP = (A+B)>>1;
                        while (A <= B && B < weights.start[line_index] && matrix->start[j+TMP] != M)
                        {
                            if (matrix->start[j+TMP] < M) A = TMP+1;
                            else B = TMP-1;
                            TMP = (A+B)>>1;
                        }

                        if (A <= B && B < weights.start[line_index] && matrix->start[j+TMP] == M)
                        {
                            delete_classic(matrix, j+TMP);
                            weights.start[line_index]--;
                            rel_weight->start[M]--;
                        }
                        else
                        {
                            insert_classic(matrix, M, j+A);
                            weights.start[line_index]++;
                            rel_weight->start[M]++;
                        }
                    }

                    for (unsigned long z = 0 ; z < weights.start[line_index] ; z++)
                    {
                        if (matrix->start[j+z] > pivot) matrix->start[j+z]--;
                    }

                    tmp3 = get_index(sorted, line_index);

                    bubble_sort_down(&sorted, weights, tmp3); // update the sorted list of line weights
                    tmp3 = 0;
                    bubble_sort_up(&sorted, weights, tmp3); // update the sorted list of line weights

                    while (j < matrix->len && matrix->start[j] != row_delimiter) j++;
                }
                j++;
                line_index++;
            }
            delete_classic(rel_weight, pivot);
        }
        i = 0;
        tmp3 = 0;
        while (tmp3 != tmp2)
        {
            if (matrix->start[i] == row_delimiter) tmp3++;
            i++;
        }

        delete_classic(matrix, i);
        delete_classic(&weights, tmp2);

        unsigned long z = 0;
        while (z < sorted.len)
        {
            if (sorted.start[z] == tmp2) delete_classic(&sorted, z);
            else
            {
                if (sorted.start[z] > tmp2) sorted.start[z]--;
                z++;
            }
        }
    }

    free(sorted.start);
    sorted.start = NULL;

    free(weights.start);
    weights.start = NULL;
}