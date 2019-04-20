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

int main(int argc, char *argv[])
{
    FILE *fp;            /// data file pointer
    FILE *qp;            /// query file pointer
    double best_in_choosed;          /// best-so-far
    double **query, **search_query;       /// data array and query array
    double d;
    long long i, j;
    double *ex, *ex2, *std, *m_ex, *m_ex2;

    int query_len = -1;
    int warp_window = -1;

    // corresponds to choosing internal distance function
    // 1: L_1 metric
    // 2: L_2 metric
    // 3: cosine similarity
    unsigned int distance_function = 2;

    long long current_location = 0;
}