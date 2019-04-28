//
// Created by Смирнов Влад on 2019-04-28.
//

#include "DTW_matrix.h"

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define INF 1.79769e+308

using namespace std;

void getDtwMatrix(char* pathName);
double dtw(double** fRow, double** sRow);
double dist(double* a, double* b);


int dim = 3;
int distType = 0;
int num_size;
int seq_size = 0;
double** matrix;
double** dtwMatrix;
int halfWindow;



void getDtwMatrix(char* pathName) {

    /// vars initialization

    FILE *fp, *fw;

    fp = fopen(pathName,"r");
    if( fp == NULL ) {
        printf("ERROR : Cannot open file\n\n");
        exit(1);
    }


    int j = 0, k = 0, s = 0; /// [j][k][s]

    double d = 0.0;

    double** row;

    vector<double**> data;

    /// upload data

    cout << "start loading" << endl;

    int count = 0;

    while(fscanf(fp, "%lf", &d) != EOF) {

        count++;

        if (s == dim) {
            s = 0;
            k++;
        }

        if (k >= seq_size) {
            k = 0;
            s = 0;
        }

        if (k == 0 && s == 0) {
            row = (double**)malloc(sizeof(double*) * seq_size);
            for (int i = 0; i < seq_size; i++) {
                row[i] = (double*)malloc(sizeof(double) * dim);
            }
            data.push_back(row);
            j++;

        }

        row[k][s++] = d;
    }

    if (k < seq_size - 1 && s < dim - 1) {
        cout << "Last element is not full, pop it back! "<< endl;
        data.pop_back();
    }

    num_size = data.size();

    fclose(fp);

    cout << "end loading, line num = " << data.size() << endl;


    dtwMatrix = (double **) malloc(sizeof(double*) * num_size);
    for (int i = 0; i < num_size; i++) {
        dtwMatrix[i] = (double *) malloc(sizeof(double) * num_size);
    }

    matrix = (double **) malloc(sizeof(double*) * seq_size);
    for (int i = 0; i < seq_size; i++) {
        matrix[i] = (double *) malloc(sizeof(double) * seq_size);
    }

    /// calculate DTW Matrix

    cout << "start dtw calculation" << endl;

    for (int i = 0; i < num_size; i++) {

        dtwMatrix[i][i] = 0.0;


        for (int r = i + 1; r < num_size; r++) {

            dtwMatrix[i][r] = dtw(data.at(i), data.at(r));
            dtwMatrix[r][i] = dtwMatrix[i][r];
        }
    }

    cout << "end dtw calculation" << endl;

    /// print dtw Matrix

    fw = fopen("matrix_result.txt","w");

    for (int i = 0; i < num_size; i++) {

        for (int r = 0; r < num_size; r++) {

            fprintf(fw, "%10.2f\t", dtwMatrix[i][r]);
        }

        fprintf(fw, "\n\n");
    }

    fclose(fw);

    /// free resources

    for (int i = 0; i < num_size; i++) {
        row = data.at(i);
        for (s = 0; s < dim; s ++) {
            free(row[s]);
        }
        free(row);
    }

    for (int i = 0; i < num_size; i++) {
        free(dtwMatrix[i]);
    }

    free(dtwMatrix);

    for (int i = 0; i < seq_size; i++) {
        free(matrix[i]);
    }

    free(matrix);

    data.clear();
}

double dtw(double** fRow, double** sRow) {

    double d = 0.0;
    double m;
    int minJ, maxJ;

    int i, j;

    for (i = 0; i < seq_size; i++) {
        minJ = max(i - halfWindow, 0);
        maxJ = min(i + halfWindow, seq_size) - 1;
        for (j = minJ; j <= maxJ; j++) {

            if (i == 0 && j == 0) {
                m = 0.0;
            } else if (i == 0) {
                m = matrix[0][j - 1];
            } else if (j == 0) {
                m = matrix[i - 1][0];
            } else if (j == minJ) {
                m = min(matrix[i - 1][j - 1], matrix[i - 1][j]);
            } else if (j == maxJ) {
                m = min(matrix[i - 1][j - 1], matrix[i][j - 1]);
            } else {
                m = min(min(matrix[i - 1][j - 1], matrix[i - 1][j]), matrix[i][j - 1]);
            }

            matrix[i][j] = m + dist(fRow[i], sRow[j]);

        }
    }

    i = seq_size - 1;
    j = seq_size - 1;

    double minArr[3];

    while (true) {

        minJ = max(i - halfWindow, 0);
        maxJ = min(i + halfWindow, seq_size) - 1;

        d += matrix[i][j];

        if (i == 0 && j == 0) {
            break;
        } else if (i == 0) {
            minArr[0] = matrix[i][j - 1];
            minArr[1] = INF;
            minArr[2] = INF;
        } else if (j == 0) {
            minArr[0] = INF;
            minArr[1] = INF;
            minArr[2] = matrix[i - 1][j];
        } else if (j == minJ) {
            minArr[0] = INF;
            minArr[1] = matrix[i - 1][j - 1];
            minArr[2] = matrix[i - 1][j];
        } else if (j == maxJ) {
            minArr[0] = matrix[i][j - 1];
            minArr[1] = matrix[i - 1][j - 1];
            minArr[2] = INF;
        } else {
            minArr[0] = matrix[i][j - 1];
            minArr[1] = matrix[i - 1][j - 1];
            minArr[2] = matrix[i - 1][j];
        }


        if (minArr[0] > minArr[1])
            if (minArr[1] > minArr[2]) {
                i--;
            } else {
                i--;
                j--;
            }
        else if (minArr[0] > minArr[2])
            i--;
        else
            j--;

    }

    return d;
}

double dist(double* a, double* b) {

    double d = 0.0;

    switch (distType) {
        case 0: /// L2

            for (int i = 0; i < dim; i++) {
                d += (a[i] - b[i]) * (a[i] - b[i]);
            }

            break;

        case 1: /// L1

            for (int i = 0; i < dim; i++) {
                d += abs(a[i] - b[i]);
            }

            break;
    }

    return d;

}

/// argv : [1] path [2] dimensions [3] dist type {0 - L2, 1 - L1} [4] sequence length [5] window length {positive, even, less or equal to "sequence length"}
int main(  int argc , char *argv[] ) {

    /// If not enough input, display an error.
    if (argc < 5) {
        printf("ERROR : Not enough input! \nNecessary: [1] path [2] dimensions [3] dist type {0 - L2, 1 - L1} [4] sequence length [5] window length {positive, even, less or equal to \"sequence length\"}\n\n");
        exit(1);
    }

    cout << "d " << argv[0] << endl;
    cout << "d " << argv[1] << endl;
    cout << "d " << argv[2] << endl;
    cout << "d " << argv[3] << endl;
    cout << "d " << argv[4] << endl;
    cout << "d " << argv[5] << endl;


    dim = atol(argv[2]);
    distType = atol(argv[3]);
    seq_size = atol(argv[4]);
    if (argc > 5) {
        halfWindow = atol(argv[5]) / 2;
    } else {
        halfWindow = (seq_size + 1) / 2;
    }

    getDtwMatrix(argv[1]);

    return 0;
}
