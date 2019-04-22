//
// Created by Смирнов Влад on 2019-02-24.
//

//#include "UCR_MDTW.hpp"

/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected � 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/

// todo how to launch: g++ UCR_MDTW.cpp -o mdtw
// todo how to launch: g++ UCR_DTW.cpp -o dtw
// todo how to launch: ./mdtw s_query.txt s_query.txt 5000 100
// todo how to launch: ./dtw s_query.txt s_query.txt 5000 100

#include "ucr_lib.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <iostream>
#include <string>
#include <map>
#include <climits>
#include <iomanip>
#include <cmath>

using namespace std;


/// If expected error happens, teminated the program.
void error(int id)
{
    if (id == 1)
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if (id == 2)
        printf("ERROR : File not Found!!!\n\n");
    else if (id == 3)
        printf("ERROR : Can't create Output File!!!\n\n");
    else if (id == 4)
    {
        printf("ERROR : Invalid Number of Arguments!!!\n");
        printf("Command Usage:  UCR_DTW.exe  data-file  query-file   m   R\n\n");
        printf("For example  :  UCR_DTW.exe  data.txt   query.txt   128  0.05\n");
    }
    exit(1);
}

double one_dim_dist(double x, double y, unsigned short distance_function)
{
    switch (distance_function)
    {
    case 1:
        return abs(x - y);
    case 2:
        return (x - y) * (x - y);
    case 3:
        return 0.0;

    default:
        cout << "unknowb distance function" << std::endl;
        return INF;
    }
}

double cosine_dist(double *x, double *y)
{
    double mul = 0.0;
    double d_x = 0.0;
    double d_y = 0.0 ;

    for(unsigned int i = 0; i < DIM; ++i)
    {
        mul += x[i] * y[i];
        d_x += x[i] * x[i];
        d_y += y[i] * y[i];
    }

    if (d_x == 0.0f || d_y == 0.0f)
        return 1;

    return 1 - abs(mul) / (sqrt(d_x) * sqrt(d_y));
}

// Need to use sqrt ??
double dist(double *x, double *y, unsigned short distance_function)
{
    double d = 0;
    if (distance_function == 3)
        return cosine_dist(x, y);

    for (unsigned int i = 0; i < DIM; i++)
    {
        d += one_dim_dist(x[i], y[i], distance_function);
    }

    return d;
}

double *dim_dist(double *x, double *y, double *d)
{

    for (unsigned int i = 0; i < DIM; i++)
    {
        d[i] = one_dim_dist(x[i], y[i], ESTIMATION_ORD);
    }

    return d;
}

double sum(double *x, int size)
{
    double dd = 0;

    for (unsigned int s = 0; s < size; s++)
    {
        dd += x[s];
    }

    return dd;
}

double *min(double *x, double *y)
{
    double xo = 0;
    double yo = 0;
    for (unsigned int i = 0; i < DIM; i++)
    {
        xo += x[i] * x[i];
        yo += y[i] * y[i];
    }

    return xo < yo ? x : y;
}

double *max(double *x, double *y)
{

    double xo = 0;
    double yo = 0;
    for (unsigned int i = 0; i < DIM; i++)
    {
        xo += x[i] * x[i];
        yo += y[i] * y[i];
    }

    return xo > yo ? x : y;
}


/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int comp(const void *a, const void *b)
{
    Index *x = (Index *) a;
    Index *y = (Index *) b;
    double v_x = 0, v_y = 0;

    for (unsigned int i = 0; i < DIM; i++)
    {
        v_x += pow(x->value[i], 2);
        v_y += pow(y->value[i], 2);
    }

    if ((v_y - v_x) > 0)
        return 1;
    else if ((v_y - v_x) < 0)
        return -1;
    else
        return 0;
}

/// Initial the queue at the begining step of envelop calculation
void init(deque *d, int capacity)
{
    d->capacity = capacity;
    d->size = 0;
    d->dq = (int *) malloc(sizeof(int) * d->capacity);
    d->f = 0;
    d->r = d->capacity - 1;
}

/// Init dynamic array
void init_array(double **&array, int size, int dim)
{
    array = (double **) calloc(size, sizeof(double *));

    if (array == NULL)
        error(1);

    for (unsigned int i = 0; i < size; i++)
    {
        array[i] = (double *) calloc(dim, sizeof(double));
    }
}

/// Free dynamic array
void destroy_array(double **&array, int size)
{
    for (unsigned int i = 0; i < size; i++)
    {
        free(array[i]);
    }

    free(array);
}

/// Destroy the queue
void destroy(deque *d)
{
    free(d->dq);
}

/// Insert to the queue at the back
void push_back(struct deque *d, int v)
{
    d->dq[d->r] = v;
    d->r--;
    if (d->r < 0)
        d->r = d->capacity - 1;
    d->size++;
}

/// Delete the current (front) element from queue
void pop_front(struct deque *d)
{
    d->f--;
    if (d->f < 0)
        d->f = d->capacity - 1;
    d->size--;
}

/// Delete the last element from queue
void pop_back(struct deque *d)
{
    d->r = (d->r + 1) % d->capacity;
    d->size--;
}

/// Get the value at the current position of the circular queue
int front(struct deque *d)
{
    int aux = d->f - 1;

    if (aux < 0)
        aux = d->capacity - 1;
    return d->dq[aux];
}

/// Get the value at the last position of the circular queueint back(struct deque *d)
int back(struct deque *d)
{
    int aux = (d->r + 1) % d->capacity;
    return d->dq[aux];
}

/// Check whether or not the queue is empty
int empty(struct deque *d)
{
    return d->size == 0;
}

/// Finding the envelop of min and max value for LB_Keogh
/// Implementation idea is intoruduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
void lower_upper_lemire(double **t, int len, int r, double **l)
{
    struct deque *dl;

    dl = (deque *) malloc(sizeof(deque) * 2 * DIM);

    for (unsigned int d = 0; d < 2 * DIM; d++)
    {

        init(&dl[d], 2 * r + 2);
        push_back(&dl[d], 0);
    }

    for (unsigned int i = 1; i < len; i++)
    {
        for (unsigned int d = 0; d < DIM; d++)
        {
            if (i > r)
            {
                l[i - r - 1][d * 2] = t[front(&dl[d * 2])][d];
                l[i - r - 1][d * 2 + 1] = t[front(&dl[d * 2 + 1])][d];

            }
            if (t[i][d] > t[(i - 1)][d])
            {
                pop_back(&dl[d * 2]);

                while (!empty(&dl[d * 2]) && t[i][d] > t[back(&dl[d * 2])][d])
                {
                    pop_back(&dl[d * 2]);
                }

            } else
            {
                pop_back(&dl[d * 2 + 1]);
                while (!empty(&dl[d * 2 + 1]) && t[i][d] < t[back(&dl[d * 2 + 1])][d])
                    pop_back(&dl[d * 2 + 1]);

            }

            push_back(&dl[d * 2], i);
            push_back(&dl[d * 2 + 1], i);
            if (i == 2 * r + 1 + front(&dl[d * 2]))
                pop_front(&dl[d * 2]);
            else if (i == 2 * r + 1 + front(&dl[d * 2 + 1]))
                pop_front(&dl[d * 2 + 1]);

        }
    }

    for (unsigned int i = len; i < len + r + 1; i++)
    {
        for (unsigned int s = 0; s < DIM; s++)
        {
            l[i - r - 1][s * 2] = t[front(&dl[s * 2])][s];
            l[i - r - 1][s * 2 + 1] = t[front(&dl[s * 2 + 1])][s];
            if (i - front(&dl[s * 2]) >= 2 * r + 1)
                pop_front(&dl[s * 2]);
            if (i - front(&dl[s * 2 + 1]) >= 2 * r + 1)
                pop_front(&dl[s * 2 + 1]);
        }

    }

    for (unsigned int s = 0; s < DIM; s++)
    {
        destroy(&dl[s * 2]);
        destroy(&dl[s * 2 + 1]);
    }

    free(dl);

}

/// Calculate quick lower bound
/// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
/// However, because of z-normalization the top and bottom cannot give siginifant benefits.
/// And using the first and last points can be computed in constant time.
/// The prunning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
double lb_kim_hierarchy(double **t, double **q, int j, int len, double *mean, double *std,
                        unsigned short dist_function, double bsf)
{
    /// 1 point at front and back
    double lb;
    double x0[DIM];
    double y0[DIM];
    double x1[DIM];
    double y1[DIM];
    double x2[DIM];
    double y2[DIM];

    double ddd[16][DIM];

    double *d;

    for (unsigned int s = 0; s < DIM; s++)
    {
        x0[s] = (t[j][s] - mean[s]) / std[s];
        y0[s] = (t[(len - 1 + j)][s] - mean[s]) / std[s];
    }
    lb = dist(x0, q[0], dist_function) + dist(y0, q[len - 1], dist_function);

    if (lb >= bsf) return lb;

    /// 2 points at front

    for (unsigned int s = 0; s < DIM; s++)
    {
        x1[s] = (t[(j + 1)][s] - mean[s]) / std[s];
    }

    d = min(dim_dist(x1, q[0], ddd[0]), dim_dist(x0, q[1], ddd[1]));
    d = min(d, dim_dist(x1, q[1], ddd[2]));
    lb += sum(d, DIM);

    if (lb >= bsf)
        return lb;

    /// 2 points at back

    for (unsigned int s = 0; s < DIM; s++)
    {
        y1[s] = (t[(len - 2 + j)][s] - mean[s]) / std[s];
    }

    d = min(dim_dist(y1, q[len - 1], ddd[3]), dim_dist(y0, q[len - 2], ddd[4])); // todo check
    d = min(d, dim_dist(y1, q[len - 2], ddd[5]));
    lb += sum(d, DIM);

    if (lb >= bsf)
        return lb;

    /// 3 points at front

    for (unsigned int s = 0; s < DIM; s++)
    {
        x2[s] = (t[(j + 2)][s] - mean[s]) / std[s];
    }

    d = min(dim_dist(x0, q[2], ddd[6]), dim_dist(x1, q[2], ddd[7]));
    d = min(d, dim_dist(x2, q[2], ddd[8]));
    d = min(d, dim_dist(x2, q[1], ddd[9]));
    d = min(d, dim_dist(x2, q[0], ddd[10]));
    lb += sum(d, DIM);

    if (lb >= bsf)
        return lb;

    /// 3 points at back

    for (unsigned int s = 0; s < DIM; s++)
    {
        y2[s] = (t[(len - 3 + j)][s] - mean[s]) / std[s];
    }
//    double y2 = (t[(len-3+j)] - mean) / std;
    d = min(dim_dist(y0, q[len - 3], ddd[11]), dim_dist(y1, q[len - 3], ddd[12]));
    d = min(d, dim_dist(y2, q[len - 3], ddd[13]));
    d = min(d, dim_dist(y2, q[len - 2], ddd[14]));
    d = min(d, dim_dist(y2, q[len - 1], ddd[15]));
    lb += sum(d, DIM);

    return lb;
}

/// LB_Keogh 1: Create Envelop for the query
/// Note that because the query is known, envelop can be created once at the begenining.
///
/// Variable Explanation,
/// order : sorted indices for the query.
/// uo, lo: upper and lower envelops for the query, which already sorted.
/// t     : a circular array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
double lb_keogh_cumulative(int *order, double **t, double **bounds, double **cb, int j,
                           int len, double *mean, double *std, unsigned short dist_function,
                           double best_so_far)
{
    double lb = 0;
    double x, d;

    for (unsigned int i = 0; i < len && lb < best_so_far; i++)
    {
        for (unsigned int s = 0; s < DIM; s++)
        {
            x = (t[(order[i] + j)][s] - mean[s]) / std[s];
            d = 0;

            if (x > bounds[i][s * 2])
            {
                d = one_dim_dist(x, bounds[i][s * 2], dist_function);
            } else if (x < bounds[i][s * 2 + 1])
            {
                d = one_dim_dist(x, bounds[i][2 * s + 1], dist_function);
            }

            lb += d;
            cb[order[i]][s] = d;
        }

    }
    return lb;
}

/// LB_Keogh 2: Create Envelop for the data
/// Note that the envelops have been created (in main function) when each data point has been read.
///
/// Variable Explanation,
/// tz: Z-normalized data
/// qo: sorted query
/// cb: (output) current bound at each position. Used later for early abandoning in DTW.
/// l,u: lower and upper envelop of the current data
double lb_keogh_data_cumulative(int *order, double **tz, double **qo, double **cb, double **l, int len, double *mean,
                                double *std, unsigned short dist_function, double best_so_far)
{
    double lb = 0;
    double ll[DIM * 2], d;

    for (unsigned int i = 0; i < len && lb < best_so_far; i++)
    {
        for (unsigned int s = 0; s < DIM; s++)
        {
            ll[s * 2] = (l[order[i]][s * 2] - mean[s]) / std[s];
            ll[s * 2 + 1] = (l[order[i]][s * 2 + 1] - mean[s]) / std[s];

            d = 0;
            if (qo[i][s] > ll[s * 2])
                d = one_dim_dist(qo[i][s], ll[s * 2], dist_function);
            else if (qo[i][s] < ll[s * 2 + 1])
            {
                d = one_dim_dist(qo[i][s], ll[s * 2 + 1], dist_function);
            }

            lb += d;
            cb[order[i]][s] = d;
        }
    }
    return lb;
}


// ################################################################
// ###########################  DEBUG   ###########################
// ################################################################

/// Print function for debugging
void print(double **x, int len)
{
    for (unsigned int i = 0; i < len; i++)
    {
        for (unsigned int j = 0; j < DIM; j++)
        {
            cout << " " << x[i][j];
        }
        cout << endl;
    }
}

// ################################################################
// ################################################################
// ################################################################

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double **A, double **B, double **cb, int m, int r,
           unsigned short dist_function, double bsf)
{

    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i, j, k;
    double x, y, z, min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (double *) malloc(sizeof(double) * (2 * r + 1));
    for (k = 0; k < 2 * r + 1; k++) cost[k] = INF;

    cost_prev = (double *) malloc(sizeof(double) * (2 * r + 1));
    for (k = 0; k < 2 * r + 1; k++) cost_prev[k] = INF;

    for (i = 0; i < m; i++)
    {
        k = max(0, r - i);
        min_cost = INF;

        for (j = max(0, i - r); j <= min(m - 1, i + r); j++, k++)
        {

            /// Initialize all row and column
            if ((i == 0) && (j == 0))
            {

                cost[k] = dist(A[0], B[0], dist_function); // todo metrics
                min_cost = cost[k];
                continue;
            }

            if ((j - 1 < 0) || (k - 1 < 0)) y = INF;
            else y = cost[k - 1];
            if ((i - 1 < 0) || (k + 1 > 2 * r)) x = INF;
            else x = cost_prev[k + 1];
            if ((i - 1 < 0) || (j - 1 < 0)) z = INF;
            else z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = min(min(x, y), z) + dist(A[i], B[j], dist_function); // todo metrics

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {
                min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if ((OPTIMIZATION == 1) && i + r < m - 1 && min_cost + sum(cb[i + r + 1], DIM) >= bsf)
        {
            free(cost);
            free(cost_prev);
            return min_cost + sum(cb[i + r + 1], DIM);
        }

        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];

    free(cost);
    free(cost_prev);
    return final_dtw;
}


