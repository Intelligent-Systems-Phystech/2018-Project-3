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


/// Main Function
int main(int argc, char *argv[])
{
    FILE *fp;            /// data file pointer
    FILE *qp;            /// query file pointer
    double best_in_choosed;          /// best-so-far
    double **query, **search_query;       /// data array and query array
    int *order;          /// new order of the query
    double **l, **qo, **bounds, **tz, **cb, **cb1, **cb2, **l_d;

    double d;
    long long i, j;
    double *ex, *ex2, *std, *m_ex, *m_ex2;

    int query_len = -1;
    int warp_window = -1;
    unsigned int closest_series_num = 1;

    // corresponds to choosing internal distance function
    // 1: L_1 metric
    // 2: L_2 metric
    // 3: cosine similarity
    unsigned int distance_function = 2;

    long long current_location = 0;
    double t1, t2;
    int kim = 0, keogh = 0, keogh2 = 0;
    double dist = 0, lb_kim = 0, lb_k = 0, lb_k2 = 0;
    double **buffer, **l_buff;
    Index *Q_tmp;

    auto best_locations = map<double, long long, std::greater<double>> ();
    auto best_location_last_added = std::pair<long long, double> (LLONG_MIN, -1.0);

    /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for
    /// reducing the floating point error.
    int EPOCH = 100000;

    /// If not enough input, display an error.
    if (argc <= 4)
        error(4);

    /// read size of the query
    if (argc > 3)
        query_len = static_cast<int>(atol(argv[3]));

    /// read warping windows
    if (argc > 4)
    {
        double R = atof(argv[4]);
        if (R <= 1)
            warp_window = static_cast<int>(floor(R * query_len));
        else
            warp_window = static_cast<int>(floor(R));
    }

    if (argc > 5)
        closest_series_num = static_cast<unsigned int>(atol(argv[5]));
    
    if (argc > 6)    
        distance_function = static_cast<unsigned int>(atol(argv[6]));

    fp = fopen(argv[1], "r");
    if (fp == NULL)
        error(2);

    qp = fopen(argv[2], "r");
    if (qp == NULL)
        error(2);

    /// start the clock
    t1 = clock();

    /// malloc everything here

    init_array(search_query, query_len, DIM);
    init_array(qo, query_len, DIM);
    init_array(bounds, query_len, 2 * DIM);

    order = (int *) malloc(sizeof(int) * query_len);
    if (order == NULL)
        error(1);

    Q_tmp = (Index *) malloc(sizeof(Index) * query_len);
    if (Q_tmp == NULL)
        error(1);

    init_array(l, query_len, 2 * DIM);
    init_array(cb, query_len, DIM);
    init_array(cb1, query_len, DIM);
    init_array(cb2, query_len, DIM);
    init_array(l_d, query_len, DIM);
    init_array(query, query_len * 2, DIM);
    init_array(tz, query_len, DIM);
    init_array(buffer, EPOCH, DIM);
    init_array(l_buff, EPOCH, 2 * DIM);

    /// Read query file
    best_in_choosed = INF;
    i = 0;

    ex = (double *) calloc(DIM, sizeof(double));
    ex2 = (double *) calloc(DIM, sizeof(double));
    m_ex = (double *) calloc(DIM, sizeof(double));
    m_ex2 = (double *) calloc(DIM, sizeof(double));
    std = (double *) calloc(DIM, sizeof(double));

    memset(ex, 0, DIM);
    memset(ex2, 0, DIM);
    memset(std, 0, DIM);

    while (fscanf(qp, "%lf", &d) != EOF && i < query_len * DIM)
    {
        ex[i % DIM] += d;
        ex2[i % DIM] += d * d;
        search_query[i / DIM][i % DIM] = d;
        i++;
    }
    fclose(qp);

    /// Do z-normalize the query, keep in same array, q
    for (unsigned int s = 0; s < DIM; s++)
    {
        ex[s] = ex[s] / query_len;
        ex2[s] = ex2[s] / query_len;
        std[s] = sqrt(ex2[s] - ex[s] * ex[s]);

        for (i = 0; i < query_len; i++)
        {
            search_query[i][s] = (search_query[i][s] - ex[s]) / std[s];
        }
    }

    /// Create envelop of the query: lower envelop, l, and upper envelop, u
    lower_upper_lemire(search_query, query_len, warp_window, l);

    /// Sort the query one time by abs(z-norm(q[i]))
    for (i = 0; i < query_len; i++)
    {
        for (unsigned int s = 0; s < DIM; s++)
        {
            Q_tmp[i].value[s] = search_query[i][s];
        }
        Q_tmp[i].index = i;
    }
    qsort(Q_tmp, query_len, sizeof(Index), comp);

    /// also create another arrays for keeping sorted envelop
    for (i = 0; i < query_len; i++)
    {
        int o = Q_tmp[i].index;
        order[i] = o;

        for (unsigned int s = 0; s < DIM; s++)
        {
            qo[i][s] = search_query[o][s];
            bounds[i][2 * s] = l[o][2 * s];
            bounds[i][2 * s + 1] = l[o][2 * s + 1];
        }
    }
    free(Q_tmp);

    /// current index of the data in current chunk of size EPOCH
    /// the starting index of the data in the circular array, t

    memset(ex, 0, DIM);
    memset(ex2, 0, DIM);
    memset(std, 0, DIM);

    bool done = false;
    int it = 0, ep = 0, k = 0;
    int b = 0;
    long long I;    /// the starting index of the data in current chunk of size EPOCH

    while (!done)
    {
        /// Read first m-1 points
        if (it == 0)
        {
            for (k = 0; k < (query_len - 1) * DIM; k++)
            {
                if (fscanf(fp, "%lf", &d) != EOF)
                {
                    buffer[k / DIM][k % DIM] = d;
                }
            }
        } else
        {
            for (k = 0; k < (query_len - 1) * DIM; k++)
            {
                buffer[k / DIM][k % DIM] = buffer[EPOCH - query_len + 1 + k / DIM][k % DIM];
            }
        }

        /// Read buffer of size EPOCH or when all data has been read.
        ep = query_len - 1;

        b = 0;

        while (ep < EPOCH)
        {
            for (unsigned int s = 0; s < DIM; s++)
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
        if (ep <= (query_len - 1))
        {
            done = true;
        } else
        {
            lower_upper_lemire(buffer, ep, warp_window, l_buff);

            /// Just for printing a dot for approximate a million point. Not much accurate.
            if (it % (1000000 / (EPOCH - query_len + 1)) == 0)
                fprintf(stderr, ".");

            /// Do main task here..
            for (unsigned int s = 0; s < DIM; s++)
            {
                ex[s] = 0;
                ex2[s] = 0;
                std[s] = 0;
            }

            for (i = 0; i < ep; i++)
            {

                for (unsigned int s = 0; s < DIM; s++)
                {
                    /// A bunch of data has been read and pick one of them at a time to use
                    d = buffer[i][s];

                    /// Calcualte sum and sum square
                    ex[s] += d;
                    ex2[s] += d * d;

                    /// t is a circular array for keeping current data
                    query[i % query_len][s] = d;


                    /// Double the size for avoiding using modulo "%" operator
                    query[(i % query_len) + query_len][s] = d;

                }

                /// Start the task when there are more than m-1 points in the current chunk
                if (i >= query_len - 1)
                {
                    for (unsigned int s = 0; s < DIM; s++)
                    {
                        m_ex[s] = ex[s] / query_len;
                        m_ex2[s] = ex2[s] / query_len;
                        std[s] = sqrt(m_ex2[s] - m_ex[s] * m_ex[s]);
                    }

                    /// compute the start location of the data in the current circular array, t
                    j = (i + 1) % query_len;
                    /// the start location of the data in the current chunk
                    I = i - (query_len - 1);

                    if (OPTIMIZATION == 0)
                    {
                        for (k = 0; k < query_len; k++)
                        {
                            for (unsigned int s = 0; s < DIM; s++)
                            {
                                tz[k][s] = (query[(k + j)][s] - m_ex[s]) / std[s];
                            }
                        }

                        dist = dtw(tz, search_query, cb, query_len, warp_window, distance_function,
                            best_in_choosed);

                        if (best_locations.begin()->first > dist or best_locations.size() < closest_series_num)
                        {
                            current_location = (it) * (EPOCH - query_len + 1) + i - query_len + 1;

                            if (best_location_last_added.first <= current_location - query_len)
                            {
                                if (best_locations.size() == closest_series_num)
                                    best_locations.erase(best_locations.begin()->first);
                                best_locations[dist] = current_location;
                                best_location_last_added.first = current_location;
                                best_location_last_added.second = dist;
                            }
                            else if (dist < best_location_last_added.second )
                            {
                                best_locations.erase(best_location_last_added.second);
                                best_locations[dist] = current_location;
                                best_location_last_added.first = current_location;
                                best_location_last_added.second = dist;
                            }

                            best_in_choosed = best_locations.begin()->first;
                        }
                    } else
                    {
                        lb_kim = lb_kim_hierarchy(query, search_query, j, query_len, m_ex, std,
                            distance_function, best_in_choosed);

                        if (lb_kim < best_in_choosed)
                        {
                            /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                            /// uo, lo are envelop of the query.
                            lb_k = lb_keogh_cumulative(order, query, bounds, cb1, j, query_len, m_ex, std,
                                distance_function, best_in_choosed);
                            if (lb_k < best_in_choosed)
                            {
                                /// Take another linear time to compute z_normalization of t.
                                /// Note that for better optimization, this can merge to the previous function.
                                for (k = 0; k < query_len; k++)
                                {
                                    for (unsigned int s = 0; s < DIM; s++)
                                    {
                                        tz[k][s] = (query[(k + j)][s] - m_ex[s]) / std[s];
                                    }
                                }

                                /// Use another lb_keogh to prune
                                /// qo is the sorted query. tz is unsorted z_normalized data.
                                /// l_buff, u_buff are big envelop for all data in this chunk
                                lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff + I, query_len, m_ex, std,
                                    distance_function, best_in_choosed);

                                if (lb_k2 < best_in_choosed)
                                {
                                    /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                                    /// Note that cb and cb2 will be cumulative summed here.
                                    if (lb_k > lb_k2)
                                    {
                                        for (unsigned int s = 0; s < DIM; s++)
                                        {
                                            cb[query_len - 1][s] = cb1[query_len - 1][s];
                                        }
                                        for (k = query_len - 2; k >= 0; k--)
                                        {
                                            for (unsigned int s = 0; s < DIM; s++)
                                            {
                                                cb[k][s] = cb[k + 1][s] + cb1[k][s];
                                            }
                                        }
                                    } else
                                    {
                                        for (unsigned int s = 0; s < DIM; s++)
                                        {
                                            cb[query_len - 1][s] = cb2[query_len - 1][s];
                                        }
                                        for (k = query_len - 2; k >= 0; k--)
                                        {
                                            for (unsigned int s = 0; s < DIM; s++)
                                            {
                                                cb[k][s] = cb[k + 1][s] + cb2[k][s];
                                            }
                                        }
                                    }

                                    /// Compute DTW and early abandoning if possible
                                    // todo check performance
                                    dist = dtw(tz, search_query, cb, query_len, warp_window, distance_function,
                                        best_in_choosed);

                                    if (best_locations.begin()->first > dist or best_locations.size() < closest_series_num)
                                    {
                                        current_location = (it) * (EPOCH - query_len + 1) + i - query_len + 1;

                                        if (best_location_last_added.first <= current_location - query_len)
                                        {
                                            if (best_locations.size() == closest_series_num)
                                                best_locations.erase(best_locations.begin()->first);
                                            best_locations[dist] = current_location;
                                            best_location_last_added.first = current_location;
                                            best_location_last_added.second = dist;
                                        }
                                        else if (dist < best_location_last_added.second)
                                        {
                                            best_locations.erase(best_location_last_added.second);
                                            best_locations[dist] = current_location;
                                            best_location_last_added.first = current_location;
                                            best_location_last_added.second = dist;
                                        }

                                        best_in_choosed = best_locations.begin()->first;
                                    }
                                } else
                                    keogh2++;
                            } else
                                keogh++;
                        } else
                            kim++;
                    }
                            /// Use a constant lower bound to prune the obvious subsequence


                    /// Reduce obsolute points from sum and sum square
                    for (unsigned int s = 0; s < DIM; s++)
                    {
                        ex[s] -= query[j][s];
                        ex2[s] -= query[j][s] * query[j][s];
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

    i = (it) * (EPOCH - query_len + 1) + ep;

//    free(q);
    destroy_array(search_query, query_len);
    destroy_array(qo, query_len);
    destroy_array(bounds, query_len);
    free(order);
    destroy_array(l, query_len);
    destroy_array(cb, query_len);
    destroy_array(cb1, query_len);
    destroy_array(cb2, query_len);
    destroy_array(l_d, query_len);
    destroy_array(query, query_len * 2);
    destroy_array(tz, query_len);
    destroy_array(buffer, query_len);
    destroy_array(l_buff, query_len);

    fclose(fp);

    free(ex);
    free(m_ex);
    free(ex2);
    free(m_ex2);
    free(std);

    t2 = clock();
    printf("\n");

    std::cout << std::fixed;
    std::cout << std::setprecision(4);

    /// Note that loc and i are long long.
    cout << "Location : " << best_locations.rbegin()->second << endl;
    cout << "Distance : " << best_locations.rbegin()->first << endl;
    cout << "Data Scanned : " << i << endl;
    cout << "Total Execution Time : " << (t2 - t1) / CLOCKS_PER_SEC << " sec" << endl;

    /// printf is just easier for formating ;)
    printf("\n");
    printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i) * 100);
    printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i) * 100);
    printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i) * 100);
    printf("DTW Calculation     : %6.2f%%\n", 100 - (((double) kim + keogh + keogh2) / i * 100));

    cout << endl
         << std::left << std::setw(8) << "locate"
         << std::right << std::setw(8) << "dist" << endl
         << std::string(16, '-') << endl;

    for (auto it = best_locations.rbegin(); it != best_locations.rend(); ++it)
        cout << std::left << std::setw(8) << it->second
             << std::right << std::setw(8) << it->first << endl;

    return 0;
}