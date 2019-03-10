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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cstring>
#include <iostream>

//#define min(x,y) ((x)<(y)?(x):(y))
//#define max(x,y) ((x)>(y)?(x):(y))
//#define dist(x,y) ((x-y)*(x-y))

#define INF 1e20       //Pseudo Infitinte number for this code

#define DIM 2

using namespace std;

/// Data structure for sorting the query
typedef struct Index // todo MDTW
{
    double value[DIM];
    int index;
} Index;

/// Data structure (circular array) for finding minimum and maximum for LB_Keogh envolop
struct deque
{
    int *dq;
    int size, capacity;
    int f, r;
};

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

double dist(double *x, double *y)
{
    double d = 0;
    for (int i = 0; i < DIM; i++)
    {
        d += (x[i] - y[i]) * (x[i] - y[i]);
    }

    return d;
}

double *dim_dist(double *x, double *y, double *d)
{

    for (int i = 0; i < DIM; i++)
    {
        d[i] = (x[i] - y[i]) * (x[i] - y[i]);
    }

    return d;
}

double dist(double *x, int size)
{
    double d = 0;
    for (int i = 0; i < size; i++)
    {
        d += x[i] * x[i];
    }

    return d;
}

double dist(double x, double y)
{

    return (x - y) * (x - y);
}

double sum(double *x, int size)
{
    double dd = 0;

    for (int s = 0; s < size; s++)
    {
        dd += x[s];
    }

    return dd;
}

double *min(double *x, double *y)
{
    double xo = 0;
    double yo = 0;
    for (int i = 0; i < DIM; i++)
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
    for (int i = 0; i < DIM; i++)
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

    for (int i = 0; i < DIM; i++)
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

    for (int i = 0; i < size; i++)
    {
        array[i] = (double *) calloc(dim, sizeof(double));
    }
}

/// Free dynamic array
void destroy_array(double **&array, int size)
{
    for (int i = 0; i < size; i++)
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

    for (int d = 0; d < 2 * DIM; d++)
    {

        init(&dl[d], 2 * r + 2);
        push_back(&dl[d], 0);
    }

    for (int i = 1; i < len; i++)
    {
        for (int d = 0; d < DIM; d++)
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

    for (int i = len; i < len + r + 1; i++)
    {
        for (int s = 0; s < DIM; s++)
        {
            l[i - r - 1][s * 2] = t[front(&dl[s * 2])][s];
            l[i - r - 1][s * 2 + 1] = t[front(&dl[s * 2 + 1])][s];
            if (i - front(&dl[s * 2]) >= 2 * r + 1)
                pop_front(&dl[s * 2]);
            if (i - front(&dl[s * 2 + 1]) >= 2 * r + 1)
                pop_front(&dl[s * 2 + 1]);
        }

    }

    for (int s = 0; s < DIM; s++)
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
double lb_kim_hierarchy(double **t, double **q, int j, int len, double *mean, double *std, double bsf = INF)
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

    for (int s = 0; s < DIM; s++)
    {

        x0[s] = (t[j][s] - mean[s]) / std[s];
        y0[s] = (t[(len - 1 + j)][s] - mean[s]) / std[s];
    }

    lb = dist(x0, q[0]) + dist(y0, q[len - 1]);

    if (lb >= bsf) return lb;

    /// 2 points at front

    for (int s = 0; s < DIM; s++)
    {
        x1[s] = (t[(j + 1)][s] - mean[s]) / std[s];
    }


    d = min(dim_dist(x1, q[0], ddd[0]), dim_dist(x0, q[1], ddd[1]));
    d = min(d, dim_dist(x1, q[1], ddd[2]));
    lb += sum(d, DIM);

    if (lb >= bsf) return lb;

    /// 2 points at back

    for (int s = 0; s < DIM; s++)
    {
        y1[s] = (t[(len - 2 + j)][s] - mean[s]) / std[s];
    }

    d = min(dim_dist(y1, q[len - 1], ddd[3]), dim_dist(y0, q[len - 2], ddd[4])); // todo check
    d = min(d, dim_dist(y1, q[len - 2], ddd[5]));
    lb += sum(d, DIM);

    if (lb >= bsf) return lb;

    /// 3 points at front

    for (int s = 0; s < DIM; s++)
    {
        x2[s] = (t[(j + 2)][s] - mean[s]) / std[s];
    }

    d = min(dim_dist(x0, q[2], ddd[6]), dim_dist(x1, q[2], ddd[7]));
    d = min(d, dim_dist(x2, q[2], ddd[8]));
    d = min(d, dim_dist(x2, q[1], ddd[9]));
    d = min(d, dim_dist(x2, q[0], ddd[10]));
    lb += sum(d, DIM);

    if (lb >= bsf) return lb;

    /// 3 points at back

    for (int s = 0; s < DIM; s++)
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
double lb_keogh_cumulative(int *order, double **t, double **bounds, double **cb, int j, int len, double *mean, double *std,
                           double best_so_far = INF)
//double lb_keogh_cumulative(int* order, double *t, double *uo, double *lo, double *cb, int j, int len, double mean, double std, double best_so_far = INF)
{
    double lb = 0;
    double x, d;

    cout << "x = HUI"  << endl;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        for (int s = 0; s < DIM; s++)
        {
            x = (t[(order[i] + j)][s] - mean[s]) / std[s];
            d = 0;

            if (x > bounds[i][s * 2])
            {
                d = dist(x, bounds[i][s * 2]);
            } else if (x < bounds[i][s * 2 + 1])
            {
                d = dist(x, bounds[i][2 * s + 1]);
            }

//            if (s == 0)
//            {
////                cout << "d = " << d << endl;
////                cout << "s = " << s << endl;
//                cout << "x = " << x << endl;
//                cout << "x = " << x << endl;
////                cout << "u_bounds = " << bounds[i][0] << endl;
////                cout << "l_bounds = " << bounds[i][1] << endl;
//            }

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
                                double *std, double best_so_far = INF)
{
    double lb = 0;
    double ll[DIM * 2], d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        for (int s = 0; s < DIM; s++)
        {
            ll[s * 2] = (l[order[i]][s * 2] - mean[s]) / std[s];
            ll[s * 2 + 1] = (l[order[i]][s * 2 + 1] - mean[s]) / std[s];

            d = 0;
            if (qo[i][s] > ll[s * 2])
                d = dist(qo[i][s], ll[s * 2]);
            else if (qo[i][s] < ll[s * 2 + 1])
            {
                d = dist(qo[i][s], ll[s * 2 + 1]);
            }

            lb += d;
            cb[order[i]][s] = d;
        }
    }
    return lb;
}

/// Print function for debugging
void printArray(double *x, int len)
{
    for (int i = 0; i < len; i++)
        printf(" %6.2lf", x[i]);
    printf("\n");
}

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double **A, double **B, double **cb, int m, int r, double bsf = INF)
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

                cost[k] = dist(A[0], B[0]); // todo metrics
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
            cost[k] = min(min(x, y), z) + dist(A[i], B[j]); // todo metrics

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {
                min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (i + r < m - 1 && min_cost + sum(cb[i + r + 1], DIM) >= bsf)
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

/// Main Function
int main(int argc, char *argv[])
{
    FILE *fp;            /// data file pointer
    FILE *qp;            /// query file pointer
    double bsf;          /// best-so-far
    double **t, **q;       /// data array and query array
    int *order;          /// new order of the query
    double **l, **qo, **lo, **tz, **cb, **cb1, **cb2, **l_d;


    double d;
    long long i, j;
    double *ex, *ex2, *std, *m_ex, *m_ex2;
    int m = -1, r = -1;
    long long loc = 0;
    double t1, t2;
    int kim = 0, keogh = 0, keogh2 = 0;
    double dist = 0, lb_kim = 0, lb_k = 0, lb_k2 = 0;
    double **buffer, **l_buff;
    Index *Q_tmp;

    /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
    int EPOCH = 100000;

    /// If not enough input, display an error.
    if (argc <= 3)
        error(4);

    /// read size of the query
    if (argc > 3)
        m = atol(argv[3]);

    /// read warping windows
    if (argc > 4)
    {
        double R = atof(argv[4]);
        if (R <= 1)
            r = floor(R * m);
        else
            r = floor(R);
    }

    fp = fopen(argv[1], "r");
    if (fp == NULL)
        error(2);

    qp = fopen(argv[2], "r");
    if (qp == NULL)
        error(2);

    /// start the clock
    t1 = clock();


    /// malloc everything here

    init_array(q, m, DIM);

    init_array(qo, m, DIM);
    init_array(lo, m, DIM);

    order = (int *) malloc(sizeof(int) * m);
    if (order == NULL)
        error(1);

    Q_tmp = (Index *) malloc(sizeof(Index) * m);
    if (Q_tmp == NULL)
        error(1);

    init_array(l, m, 2 * DIM);
    init_array(cb, m, DIM);
    init_array(cb1, m, DIM);
    init_array(cb2, m, DIM);
    init_array(l_d, m, DIM);
    init_array(t, m * 2, DIM);
    init_array(tz, m, DIM);
    init_array(buffer, EPOCH, DIM);
    init_array(l_buff, EPOCH, 2 * DIM);


    /// Read query file
    bsf = INF;
    i = 0;
    j = 0;

    ex = (double *) calloc(DIM, sizeof(double));
    ex2 = (double *) calloc(DIM, sizeof(double));
    m_ex = (double *) calloc(DIM, sizeof(double));
    m_ex2 = (double *) calloc(DIM, sizeof(double));
    std = (double *) calloc(DIM, sizeof(double));

    memset(ex, 0, DIM);
    memset(ex2, 0, DIM);
    memset(std, 0, DIM);


    while (fscanf(qp, "%lf", &d) != EOF && i < m * DIM)
    {
        ex[i % DIM] += d;
        ex2[i % DIM] += d * d;
        q[i / DIM][i % DIM] = d;
        i++;
    }
    fclose(qp);

    /// Do z-normalize the query, keep in same array, q
    for (int s = 0; s < DIM; s++)
    {

        ex[s] = ex[s] / m;
        ex2[s] = ex2[s] / m;
        std[s] = sqrt(ex2[s] - ex[s] * ex[s]);

        for (i = 0; i < m; i++)
        {

            q[i][s] = (q[i][s] - ex[s]) / std[s];
        }
    }

    /// Create envelop of the query: lower envelop, l, and upper envelop, u
    lower_upper_lemire(q, m, r, l);

    /// Sort the query one time by abs(z-norm(q[i]))
    for (i = 0; i < m; i++)
    {
        for (int s = 0; s < DIM; s++)
        {
            Q_tmp[i].value[s] = q[i][s];
        }
        Q_tmp[i].index = i;
    }
    qsort(Q_tmp, m, sizeof(Index), comp);

    /// also create another arrays for keeping sorted envelop
    for (i = 0; i < m; i++)
    {

        int o = Q_tmp[i].index;
        order[i] = o;

        for (int s = 0; s < DIM; s++)
        {
            qo[i][s] = q[o][s];
            lo[i][s] = l[o][s];
            lo[i][s + 1] = l[o][s + 1];
        }
    }
    free(Q_tmp);

    /// Initial the cummulative lower bound
//    for( i=0; i<m; i++) // todo calloc???
//    {   cb[i]=0;
//        cb1[i]=0;
//        cb2[i]=0;
//    }

    i = 0;          /// current index of the data in current chunk of size EPOCH
    j = 0;          /// the starting index of the data in the circular array, t

    memset(ex, 0, DIM);
    memset(ex2, 0, DIM);
    memset(std, 0, DIM);

    bool done = false;
    int it = 0, ep = 0, k = 0;
    long long I;    /// the starting index of the data in current chunk of size EPOCH

    int b = 0;

    while (!done)
    {
        /// Read first m-1 points
        ep = 0;
        if (it == 0)
        {
            for (k = 0; k < (m - 1) * DIM; k++)
            {
                if (fscanf(fp, "%lf", &d) != EOF)
                {
                    buffer[k / DIM][k % DIM] = d;
                }
            }
        } else
        {
            for (k = 0; k < (m - 1) * DIM; k++)
            {
                buffer[k / DIM][k % DIM] = buffer[EPOCH - m + 1 + k / DIM][k % DIM];
            }
        }

        /// Read buffer of size EPOCH or when all data has been read.
        ep = m - 1;

        b = 0;

        while (ep < EPOCH)
        {
            for (int s = 0; s < DIM; s++)
            {
                if (fscanf(fp, "%lf", &d) == EOF)
                {
                    b = 1;
                    break;
                }
                buffer[ep][s] = d;
            }
            if (b == 1) break;
            ep++;
        }

        /// Data are read in chunk of size EPOCH.
        /// When there is nothing to read, the loop is end.
        if (ep <= (m - 1))
        {

            done = true;
        } else
        {

            lower_upper_lemire(buffer, ep, r, l_buff);


            /// Just for printing a dot for approximate a million point. Not much accurate.
            if (it % (1000000 / (EPOCH - m + 1)) == 0)
                fprintf(stderr, ".");

            /// Do main task here..
            for (int s = 0; s < DIM; s++)
            {
                ex[s] = 0;
                ex2[s] = 0;
                std[s] = 0;
            }

            for (i = 0; i < ep; i++)
            {

                for (int s = 0; s < DIM; s++)
                {
                    /// A bunch of data has been read and pick one of them at a time to use
                    d = buffer[i][s];

                    /// Calcualte sum and sum square
                    ex[s] += d;
                    ex2[s] += d * d;

                    /// t is a circular array for keeping current data
                    t[i % m][s] = d;


                    /// Double the size for avoiding using modulo "%" operator
                    t[(i % m) + m][s] = d;

                }

                /// Start the task when there are more than m-1 points in the current chunk
                if (i >= m - 1)
                {
                    for (int s = 0; s < DIM; s++)
                    {
                        m_ex[s] = ex[s] / m;
                        m_ex2[s] = ex2[s] / m;
                        std[s] = sqrt(m_ex2[s] - m_ex[s] * m_ex[s]);
                    }

                    /// compute the start location of the data in the current circular array, t
                    j = (i + 1) % m;
                    /// the start location of the data in the current chunk
                    I = i - (m - 1);

                    /// Use a constant lower bound to prune the obvious subsequence
                    lb_kim = lb_kim_hierarchy(t, q, j, m, m_ex, std, bsf);

                    if (lb_kim < bsf)
                    {
                        /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                        /// uo, lo are envelop of the query.
                            lb_k = lb_keogh_cumulative(order, t, lo, cb1, j, m, m_ex, std, bsf);
                        if (lb_k < bsf)
                        {
                            /// Take another linear time to compute z_normalization of t.
                            /// Note that for better optimization, this can merge to the previous function.
                            for (k = 0; k < m; k++)
                            {
                                for (int s = 0; s < DIM; s++)
                                {
                                    tz[k][s] = (t[(k + j)][s] - m_ex[s]) / std[s];
                                }
                            }

                            /// Use another lb_keogh to prune
                            /// qo is the sorted query. tz is unsorted z_normalized data.
                            /// l_buff, u_buff are big envelop for all data in this chunk
                            lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff + I, m, m_ex, std, bsf);

                            if (lb_k2 < bsf)
                            {
                                /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                                /// Note that cb and cb2 will be cumulative summed here.
                                if (lb_k > lb_k2)
                                {
                                    for (int s = 0; s < DIM; s++)
                                    {
                                        cb[m - 1][s] = cb1[m - 1][s];
                                    }
                                    for (k = m - 2; k >= 0; k--)
                                    {
                                        for (int s = 0; s < DIM; s++)
                                        {
                                            cb[k][s] = cb[k + 1][s] + cb1[k][s];
                                        }
                                    }
                                } else
                                {
                                    for (int s = 0; s < DIM; s++)
                                    {
                                        cb[m - 1][s] = cb2[m - 1][s];
                                    }
                                    for (k = m - 2; k >= 0; k--)
                                    {
                                        for (int s = 0; s < DIM; s++)
                                        {
                                            cb[k][s] = cb[k + 1][s] + cb2[k][s];
                                        }
                                    }
                                }

                                /// Compute DTW and early abandoning if possible
                                dist = dtw(tz, q, cb, m, r, bsf); // todo check performance

                                if (dist < bsf)
                                {   /// Update bsf
                                    /// loc is the real starting location of the nearest neighbor in the file
                                    bsf = dist;
                                    loc = (it) * (EPOCH - m + 1) + i - m + 1;
                                }
                            } else
                                keogh2++;
                        } else
                            keogh++;
                    } else
                        kim++;

                    /// Reduce obsolute points from sum and sum square
                    for (int s = 0; s < DIM; s++)
                    {
                        ex[s] -= t[j][s];
                        ex2[s] -= t[j][s] * t[j][s];
                    }
                }
            }

            /// If the size of last chunk is less then EPOCH, then no more data and terminate.
            if (ep < EPOCH)
                done = true;
            else
                it++;
        }
    }

    i = (it) * (EPOCH - m + 1) + ep;

//    free(q);
    destroy_array(q, m);
    destroy_array(qo, m);
    destroy_array(lo, m);
    free(order);
    destroy_array(l, m);
    destroy_array(cb, m);
    destroy_array(cb1, m);
    destroy_array(cb2, m);
    destroy_array(l_d, m);
    destroy_array(t, m * 2);
    destroy_array(tz, m);
    destroy_array(buffer, m);
    destroy_array(l_buff, m);

    fclose(fp);

    free(ex);
    free(m_ex);
    free(ex2);
    free(m_ex2);
    free(std);

    t2 = clock();
    printf("\n");

    /// Note that loc and i are long long.
    cout << "Location : " << loc << endl;
    cout << "Distance : " << sqrt(bsf) << endl;
    cout << "Data Scanned : " << i << endl;
    cout << "Total Execution Time : " << (t2 - t1) / CLOCKS_PER_SEC << " sec" << endl;

    /// printf is just easier for formating ;)
    printf("\n");
    printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i) * 100);
    printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i) * 100);
    printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i) * 100);
    printf("DTW Calculation     : %6.2f%%\n", 100 - (((double) kim + keogh + keogh2) / i * 100));
    return 0;
}
