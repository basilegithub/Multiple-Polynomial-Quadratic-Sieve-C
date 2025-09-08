#include <gmp.h>
#include <stdlib.h>

#include "structures.h"

void reduce_relations(dyn_array* relations, dyn_array* smooth, dyn_array_classic* primes, mpz_t n)
{
    unsigned long count;
    unsigned long i = 0;
    unsigned long index = 0,index2 = 0;
    unsigned char tmp;
    unsigned char flag = 1;
    mpz_t tmp2;
    mpz_init(tmp2);
    while (flag)
    {
        flag = 0;
        for (unsigned long i = 0 ; i < primes->len ; i++)
        {
            count = 0;
            for (unsigned long j = 0 ; j < relations->len && count < 3 ; j++)
            {
                tmp = 0;
                mpz_set_ui(tmp2,*(primes->start+i));
                while (mpz_divisible_p(*(relations->start+j),tmp2))
                {
                    tmp ^= 1;
                    mpz_mul_ui(tmp2,tmp2,*(primes->start+i));
                }
                if (tmp)
                {
                    count++;
                    if (count == 1) index = j;
                    else index2 = j;
                }
            }
            if (count == 1)
            {
                delete_dyn_unsorted(relations,index);
                delete_dyn_unsorted(smooth,index);
                flag = 1;
            }
            else if (count == 2)
            {
                mpz_mul(*(relations->start+index),*(relations->start+index),*(relations->start+index2));
                mpz_mul(*(smooth->start+index),*(smooth->start+index),*(smooth->start+index2));
                mpz_mod(*(smooth->start+index),*(smooth->start+index),n);
                delete_dyn_unsorted(relations,index2);
                delete_dyn_unsorted(smooth,index2);
                flag = 1;
            }
        }
    }
    mpz_clear(tmp2);
}

void reduce_matrix(dyn_array_classic* matrix, dyn_array* relations, dyn_array* smooth, unsigned long limit, mpz_t n, dyn_array_classic* rel_weight, unsigned long merge_bound)
{
    dyn_array_classic weights;
    init_classic(&weights);
    dyn_array_classic sorted;
    init_classic(&sorted);
    unsigned long line_len = 0;
    unsigned long place = 0;
    unsigned long tmp, tmp2, tmp3 = 0;
    signed long lol;
    for (unsigned long i = 0 ; i < matrix->len ; i++)
    {
        if (*(matrix->start+i) == limit)
        {
            append_classic(&weights,line_len);
            if (weights.len == 1) append_classic(&sorted,place);
            else
            {
                signed long a = 0;
                signed long b = place-1;
                lol = (a+b)/2;
                tmp2 = *(sorted.start+lol);
                while (a <= b && *(weights.start+tmp2) != line_len)
                {
                    if (*(weights.start+tmp2) < line_len) a = lol+1;
                    else b = lol-1;
                    lol = (a+b)/2;
                    tmp2 = *(sorted.start+lol);
                }
                if (*(weights.start+tmp2) == line_len) insert_classic(&sorted,place,lol);
                else insert_classic(&sorted,place,a);
            }
            line_len = 0;
            place++;
        }
        else line_len++;
    }
    // weights are sorted
    unsigned long i,j,tmp4,k;
    unsigned long nb_lines;
    int flag = 1;
    unsigned long pivot,remain;
    while(flag)
    {
        place = 0;
        nb_lines = 0;
        flag = 0;
        i = 0;
        line_len = 0;
        while (i < matrix->len)
        {
            if (*(matrix->start+i) == limit)
            {
                if (line_len == 0)
                {
                    delete_classic(matrix,i);
                    flag = 1;
                    delete_classic(&weights,place);
                    k = 0;
                    while (k < sorted.len)
                    {
                        if (*(sorted.start+k) == place) delete_classic(&sorted,k);
                        else
                        {
                            if (*(sorted.start+k) > place) *(sorted.start+k) -= 1;
                            k++;
                        }
                    }
                }
                else if (line_len == 1)
                {
                    tmp = *(matrix->start+i-1);
                    delete_classic(matrix,i);
                    delete_classic(matrix,i-1);
                    delete_dyn(relations,tmp);
                    delete_dyn(smooth,tmp);
                    delete_classic(rel_weight,tmp);
                    delete_classic(&weights,place);
                    tmp2 = 0;
                    j = 0;
                    while (j < matrix->len)
                    {
                        if (*(matrix->start+j) == limit)
                        {
                            tmp2++;
                            j++;
                        }
                        else if (*(matrix->start+j) == tmp)
                        {
                            delete_classic(matrix,j);
                            *(weights.start+tmp2) -= 1;
                            for (unsigned long z = 0 ; z < sorted.len ; z++)
                            {
                                if (*(sorted.start+z) == tmp2)
                                {
                                    tmp3 = z;
                                    break;
                                }
                            }
                            while (tmp3 > 0 && *(weights.start+*(sorted.start+tmp3)) < *(weights.start+*(sorted.start+tmp3-1)))
                            {
                                tmp4 = *(sorted.start+tmp3-1);
                                *(sorted.start+tmp3-1) = *(sorted.start+tmp3);
                                *(sorted.start+tmp3) = tmp4;
                                tmp3--;
                            }
                        }
                        else
                        {
                            if (*(matrix->start+j) > tmp) *(matrix->start+j) -= 1;
                            j++;
                        }
                    }
                    k = 0;
                    while (k < sorted.len)
                    {
                        if (*(sorted.start+k) == place) delete_classic(&sorted,k);
                        else
                        {
                            if (*(sorted.start+k) > place) *(sorted.start+k) -= 1;
                            k++;
                        }
                    }
                    line_len = 0;
                    i = 0;
                    nb_lines = 0;
                    flag = 1;
                    place = 0;
                }
                else if (line_len == 2)
                {
                    pivot = *(matrix->start+i-2);
                    remain = *(matrix->start+i-1);
                    mpz_mul(*(relations->start+remain),*(relations->start+remain),*(relations->start+pivot));
                    mpz_mul(*(smooth->start+remain),*(smooth->start+remain),*(smooth->start+pivot));
                    mpz_mod(*(smooth->start+remain),*(smooth->start+remain),n);
                    delete_classic(matrix,i-2);
                    delete_classic(matrix,i-2);
                    delete_classic(matrix,i-2);
                    *(rel_weight->start+remain) -= 1;
                    delete_dyn(relations,pivot);
                    delete_dyn(smooth,pivot);
                    delete_classic(rel_weight,pivot);
                    delete_classic(&weights,place);

                    tmp2 = 0;
                    j = 0;
                    unsigned char has_pivot = 0;
                    while (j < matrix->len)
                    {
                        if (*(matrix->start+j) == limit)
                        {
                            if (has_pivot)
                            {
                                has_pivot = 0;
                                insert_classic(matrix,remain-1,j);
                                *(rel_weight->start+remain-1) += 1;
                                *(weights.start+tmp2) += 1;
                                for (unsigned long z = 0 ; z < sorted.len ; z++)
                                {
                                    if (*(sorted.start+z) == tmp2)
                                    {
                                        tmp3 = z;
                                        break;
                                    }
                                }
                                while (tmp3 < sorted.len-1 && *(weights.start+*(sorted.start+tmp3)) > *(weights.start+*(sorted.start+tmp3+1)))
                                {
                                    tmp4 = *(sorted.start+tmp3+1);
                                    *(sorted.start+tmp3+1) = *(sorted.start+tmp3);
                                    *(sorted.start+tmp3) = tmp4;
                                    tmp3++;
                                }
                                j++;
                            }
                            tmp2++;
                            j++;
                        }
                        else if (*(matrix->start+j) == pivot)
                        {
                            has_pivot = 1;
                            delete_classic(matrix,j);
                            *(weights.start+tmp2) -= 1;
                            for (unsigned long z = 0 ; z < sorted.len ; z++)
                            {
                                if (*(sorted.start+z) == tmp2)
                                {
                                    tmp3 = z;
                                    break;
                                }
                            }
                            while (tmp3 > 0 && *(weights.start+*(sorted.start+tmp3)) < *(weights.start+*(sorted.start+tmp3-1)))
                            {
                                tmp4 = *(sorted.start+tmp3-1);
                                *(sorted.start+tmp3-1) = *(sorted.start+tmp3);
                                *(sorted.start+tmp3) = tmp4;
                                tmp3--;
                            }
                        }
                        else if (*(matrix->start+j) == remain && has_pivot)
                        {
                            has_pivot = 0;
                            delete_classic(matrix,j);
                            *(rel_weight->start+remain-1) -= 1;
                            *(weights.start+tmp2) -= 1;
                            for (unsigned long z = 0 ; z < sorted.len ; z++)
                            {
                                if (*(sorted.start+z) == tmp2)
                                {
                                    tmp3 = z;
                                    break;
                                }
                            }
                            while (tmp3 > 0 && *(weights.start+*(sorted.start+tmp3)) < *(weights.start+*(sorted.start+tmp3-1)))
                            {
                                tmp4 = *(sorted.start+tmp3-1);
                                *(sorted.start+tmp3-1) = *(sorted.start+tmp3);
                                *(sorted.start+tmp3) = tmp4;
                                tmp3--;
                            }
                        }
                        else if (*(matrix->start+j) > remain && has_pivot)
                        {
                            has_pivot = 0;
                            insert_classic(matrix,remain-1,j);
                            *(rel_weight->start+remain-1) += 1;
                            *(weights.start+tmp2) += 1;
                            for (unsigned long z = 0 ; z < sorted.len ; z++)
                            {
                                if (*(sorted.start+z) == tmp2)
                                {
                                    tmp3 = z;
                                    break;
                                }
                            }
                            while (tmp3 < sorted.len-1 && *(weights.start+*(sorted.start+tmp3)) > *(weights.start+*(sorted.start+tmp3+1)))
                            {
                                tmp4 = *(sorted.start+tmp3+1);
                                *(sorted.start+tmp3+1) = *(sorted.start+tmp3);
                                *(sorted.start+tmp3) = tmp4;
                                tmp3++;
                            }
                            j++;
                        }
                        else
                        {
                            if (*(matrix->start+j) > pivot) *(matrix->start+j) -= 1;
                            j++;
                        }
                    }
                    k = 0;
                    while (k < sorted.len)
                    {
                        if (*(sorted.start+k) == place) delete_classic(&sorted,k);
                        else
                        {
                            if (*(sorted.start+k) > place) *(sorted.start+k) -= 1;
                            k++;
                        }
                    }


                    line_len = 0;
                    i = 0;
                    nb_lines = 0;
                    flag = 1;
                    place = 0;
                }
                else
                {
                    i++;
                    line_len = 0;
                    nb_lines++;
                    place++;
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
            flag = 1;
            place = 0;
            i = 0;
            while (place != *(sorted.start))
            {
                if (*(matrix->start+i) == limit) place++;
                i++;
            }
            if (*(matrix->start+i) == limit)
            {
                delete_classic(matrix,i);
                nb_lines--;
                delete_classic(&weights,place);
                for (unsigned long aa = 1 ; aa < sorted.len ; aa++)
                {
                    if (*(sorted.start+aa) > *(sorted.start)) *(sorted.start+aa) -= 1;
                }
                delete_classic_first(&sorted);
            }
            else
            {
                *(weights.start+place) -= 1;
                tmp2 = *(matrix->start+i);
                delete_classic(matrix,i);
                delete_dyn(relations,tmp2);
                delete_dyn(smooth,tmp2);
                delete_classic(rel_weight,tmp2);
                tmp = 0;
                k = 0;
                while (k < matrix->len)
                {
                    if (*(matrix->start+k) == limit)
                    {
                        tmp++;
                        k++;
                    }
                    else if (*(matrix->start+k) == tmp2)
                    {
                        delete_classic(matrix,k);
                        *(weights.start+tmp) -= 1;
                        for (unsigned long z = 0 ; z < sorted.len ; z++)
                        {
                            if (*(sorted.start+z) == tmp)
                            {
                                tmp3 = z;
                                break;
                            }
                        }
                        while (tmp3 > 0 && *(weights.start+*(sorted.start+tmp3)) < *(weights.start+*(sorted.start+tmp3-1)))
                        {
                            tmp4 = *(sorted.start+tmp3-1);
                            *(sorted.start+tmp3-1) = *(sorted.start+tmp3);
                            *(sorted.start+tmp3) = tmp4;
                            tmp3--;
                        }
                    }
                    else
                    {
                        if (*(matrix->start+k) > tmp2) *(matrix->start+k) -= 1;
                        k++;
                    }
                }
            }
        }
    }
    /*tmp = 0; tmp2 = 0; tmp3 = 0;
    while (tmp < matrix->len)
    {
        if (*(matrix->start+tmp) == limit)
        {
            if (tmp2 != *(weights.start+tmp3)) printf("t bz\n");
            tmp2 = 0;
            tmp3++;
        }
        else tmp2++;
        tmp++;
    }*/

    unsigned long line_index,a,b,A,B,TMP,M;
    while (*(weights.start+*(sorted.start)) <= merge_bound)
    {
        i = 0;
        tmp2 = 0;
        while (tmp2 != *(sorted.start)){
            if (*(matrix->start+i) == limit) tmp2++;
            i++;
        }
        if (*(weights.start+*(sorted.start)) > 0)
        {
            pivot = *(matrix->start+i);
            for (unsigned long j = 1 ; j < *(weights.start+tmp2) ; j++)
            {
                tmp = *(matrix->start+i+j);
                if (*(rel_weight->start+tmp) < *(rel_weight->start+pivot)) pivot = tmp;
            }
            dyn_array_classic remainer;
            init_classic(&remainer);
            for (unsigned long j = 0 ; j < *(weights.start+tmp2) ; j++)
            {
                if (*(matrix->start+i+j) != pivot)
                {
                    tmp = *(matrix->start+i+j);
                    append_classic(&remainer,tmp);
                    mpz_mul(*(relations->start+tmp),*(relations->start+tmp),*(relations->start+pivot));
                    mpz_mul(*(smooth->start+tmp),*(smooth->start+tmp),*(smooth->start+pivot));
                    mpz_mod(*(smooth->start+tmp),*(smooth->start+tmp),n);
                }
            }
            delete_dyn(relations,pivot);
            delete_dyn(smooth,pivot);
            line_index = 0;
            j = 0;
            while (line_index < weights.len)
            {
                a = 0;
                b = *(weights.start+line_index)-1;
                tmp = (a+b)>>1;
                while (a <= b && b < *(weights.start+line_index) && *(matrix->start+j+tmp) != pivot)
                {
                    if (*(matrix->start+j+tmp) < pivot) a = tmp+1;
                    else b = tmp-1;
                    tmp = (a+b)>>1;
                }
                if (a > b || b >= *(weights.start+line_index) || *(matrix->start+j+tmp) != pivot)
                {
                    j += a;
                    while (*(matrix->start+j) != limit)
                    {
                        if (*(matrix->start+j) > pivot) *(matrix->start+j) -= 1;
                        j++;
                    }
                }
                else
                {
                    delete_classic(matrix,j+tmp);
                    *(weights.start+line_index) -= 1;
                    for (unsigned long m = 0 ; m < remainer.len ; m++)
                    {
                        M = *(remainer.start+m);
                        A = 0;
                        B = *(weights.start+line_index)-1;
                        TMP = (A+B)>>1;
                        while (A <= B && B < *(weights.start+line_index) && *(matrix->start+j+TMP) != M)
                        {
                            if (*(matrix->start+j+TMP) < M) A = TMP+1;
                            else B = TMP-1;
                            TMP = (A+B)>>1;
                        }
                        if (A <= B && B < *(weights.start+line_index) && *(matrix->start+j+TMP) == M)
                        {
                            delete_classic(matrix,j+TMP);
                            *(weights.start+line_index) -= 1;
                            *(rel_weight->start+M) -= 1;
                        }
                        else
                        {
                            insert_classic(matrix,M,j+A);
                            *(weights.start+line_index) += 1;
                            *(rel_weight->start+M) += 1;
                        }
                    }
                    for (unsigned long z = 0 ; z < *(weights.start+line_index) ; z++)
                    {
                        if (*(matrix->start+j+z) > pivot) *(matrix->start+j+z) -= 1;
                    }
                    for (unsigned long z = 0 ; z < sorted.len ; z++)
                    {
                        if (*(sorted.start+z) == line_index)
                        {
                            tmp3 = z;
                            break;
                        }
                    }
                    while (tmp3 > 0 && *(weights.start+*(sorted.start+tmp3)) < *(weights.start+*(sorted.start+tmp3-1)))
                    {
                        tmp4 = *(sorted.start+tmp3-1);
                        *(sorted.start+tmp3-1) = *(sorted.start+tmp3);
                        *(sorted.start+tmp3) = tmp4;
                        tmp3--;
                    }
                    while (tmp3 < sorted.len-1 && *(weights.start+*(sorted.start+tmp3)) > *(weights.start+*(sorted.start+tmp3+1)))
                    {
                        tmp4 = *(sorted.start+tmp3+1);
                        *(sorted.start+tmp3+1) = *(sorted.start+tmp3);
                        *(sorted.start+tmp3) = tmp4;
                        tmp3++;
                    }
                    while (j < matrix->len && *(matrix->start+j) != limit) j++;
                }
                j++;
                line_index++;
            }
            delete_classic(rel_weight,pivot);
        }
        i = 0;
        tmp3 = 0;
        while (tmp3 != tmp2)
        {
            if (*(matrix->start+i) == limit) tmp3++;
            i++;
        }
        delete_classic(matrix,i);
        delete_classic(&weights,tmp2);
        unsigned long z = 0;
        while (z < sorted.len)
        {
            if (*(sorted.start+z) == tmp2) delete_classic(&sorted,z);
            else
            {
                if (*(sorted.start+z) > tmp2) *(sorted.start+z) -= 1;
                z++;
            }
        }
    }

    /*tmp = 0; tmp2 = 0; tmp3 = 0;
    while (tmp < matrix->len)
    {
        if (*(matrix->start+tmp) == limit)
        {
            if (tmp2 != *(weights.start+tmp3)) printf("t bz\n");
            tmp2 = 0;
            tmp3++;
        }
        else tmp2++;
        tmp++;
    }*/

    realloc(sorted.start,0);
    realloc(weights.start,0);
}