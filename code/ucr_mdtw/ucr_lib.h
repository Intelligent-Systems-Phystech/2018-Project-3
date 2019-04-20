

#define INF 1e20       //Pseudo Infitinte number for this code
#define DIM 3
#define OPTIMIZATION 0
#define ESTIMATION_ORD 2

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

void error(int id);

double one_dim_dist(double x, double y, unsigned short distance_function);
double cosine_similarity(double *x, double *y);
double dist(double *x, double *y, unsigned short distance_function);
double *dim_dist(double *x, double *y, double *d);

double sum(double *x, int size);
double *min(double *x, double *y);
double *max(double *x, double *y);
void init_array(double **&array, int size, int dim);

int comp(const void *a, const void *b);
void init(deque *d, int capacity);
void destroy_array(double **&array, int size);
void destroy(deque *d);
void push_back(struct deque *d, int v);
void pop_front(struct deque *d);
void pop_back(struct deque *d);
int front(struct deque *d);
int back(struct deque *d);
int empty(struct deque *d);

// Bounds
void lower_upper_lemire(double **t, int len, int r, double **l);
double lb_kim_hierarchy(double **t, double **q, int j, int len, double *mean, double *std,
                        unsigned short dist_function, double bsf=INF);

double lb_keogh_cumulative(int *order, double **t, double **bounds, double **cb, int j,
    int len, double *mean, double *std, unsigned short dist_function, double best_so_far=INF);

double lb_keogh_data_cumulative(int *order, double **tz, double **qo, double **cb, double **l, int len, double *mean,
                                double *std, unsigned short dist_function, double best_so_far = INF);


double dtw(double **A, double **B, double **cb, int m, int r, unsigned short dist_function, double bsf=INF);


// ################################################################
// ###########################  DEBUG   ###########################
// ################################################################
void print(double **x, int len);
// ################################################################
// ################################################################
// ################################################################
