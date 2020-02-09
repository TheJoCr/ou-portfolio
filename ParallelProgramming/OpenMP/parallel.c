#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <math.h>

double g(double x, double y) {
    if(x == 0.0 || y == 0.0) {
        return 0;
    }
    if(x == 1.0) {
        return y;
    }
    if(y == 1.0) {
        return x;
    }
    printf("This should not happen x = %f, y = %f!\n", x, y);
    return 0;
}

double f(double x, double y) {
    return 0.0;
}

double** alloc2d(int rows, int cols) {
    double** arr2d = (double **) malloc(sizeof(double*) * rows );
    double* data = (double *) malloc(sizeof(double) * rows * cols );
    for(int i = 0; i < rows ; i++)
        arr2d[i] = &data[cols * i];
    return arr2d;
}

void free2d( double** arr2d ) {
    free( *arr2d );
    free( arr2d );
}

void print2d( double** arr2d, int rows, int cols ) {
    for(int i = 0; i < rows; i ++) {
        printf("[");
        for(int j = 0; j < cols; j++) {
            printf(" %.4f,", arr2d[i][j] );
        }
        printf("]\n");
    }
}

int main(int argc, char* args[]) {
    clock_t start = clock();

    double a = 0, b = 1;
    double c = 0, d = 1;
    int i, j, k;
    int num_t = atoi(args[2]);
    int m = atoi(args[1]); // choose a value for m
    int n = m ; // we choose the same value for n
    double h = (b - a)/n; // interval size on x-axis
    double l = (d - c)/m; // interval size on y-axis

    // init xs init ys us
    double x[n + 1];
    double y[m + 1];

    // Tricky initialization of u:
    double** u = alloc2d(n+1, m+1);
#pragma omp parallel for num_threads(num_t)
    for (i = 0; i <= n; i++) x[i] = a + i*h;// interval points on x-axis
#pragma omp parallel for num_threads(num_t)
    for (j = 0; j <= m; j++) y[j] = c + j*l; // interval points on y axis

    // we want to find u(xi, yj), for i=1, 2, n - 1, and j = 1, 2, ..., m - 1
    // In the following we denote u(xi, yj) as u[i][j].
    // The values of u[i][j] on the boundary of the rectangle are given as follows:
#pragma omp parallel for num_threads(num_t)
    for (i = 1; i <= n; i++) u[i][0] = g(x[i], y[0]);
#pragma omp parallel for num_threads(num_t)
    for (i = 1; i <= n; i++) u[i][m] = g(x[i], y[m]);
#pragma omp parallel for num_threads(num_t)
    for (j = 0; j <= m; j++) u[0][j] = g(x[0], y[j]);
#pragma omp parallel for num_threads(num_t)
    for (j = 0; j <= m; j++) u[n][j] = g(x[n], y[j]);

    // We have to determine the value of u[i][j] inside the rectangle
    // If m = n, this gives a linear system of (n - 1)2 equations: AW = b
    // The matrix A is a block tridiagonal matrix with matrix B along the diagonal and identity matrices on
    //  the super and sub diagonals.
    // The matrix B itself is a tridiagonal matrix with the value 4 along the main diagonal and the values -1
    //  on the super and sub diagonals
    // We use here Gauss-Seidel  iterative algorithm.

    double lambda = (h*h)/(l*l);
    double div = 1.0/(2 * (lambda + 1));
    double minhsq = -h*h;

    // initial guess
    for (i = 1; i < n; i++) {
        for (j = 1; j < m; j++) {
            u[i][j] = 0.0;
        }
    }

    // Start iteration
    double modulus = 1.0;
    int ittr_count = 0;
    double* intermediate = malloc( (n+1) * sizeof(double) );
    double* new_intermediate = malloc( (n+1) * sizeof(double) );
#pragma omp parallel for num_threads(num_t)
    for( i = 0; i < n+1; i ++) {
        intermediate[i] = 0;
        new_intermediate[i] = 0;
    }

    // Let's parallelize the beef of the work:
#pragma omp parallel num_threads(num_t) private(i, j, k)
    {
    while (modulus > 0.0000001) {
#pragma omp single
        ittr_count++;
        modulus = 0;
        double priv_modulus = 0;
        for (j = m - 1; j > 0; j--) {
            // The trick here is to recognize that only one term - here labelled 'down' - is actually needed.
            // We compute everything else, then it just becomes a prefix sum problem.
#pragma omp for schedule(static)
            for (k = n; k >= 1; k-- ) {
                double new_val = div *
                                 (minhsq * f(x[k], y[j])
                                  + lambda * u[k][j + 1] // RIGHT
                                  + u[k - 1][j] // UP
                                  + lambda * u[k][j - 1] // LEFT
                                 );
                intermediate[k] = new_val;
            }
#pragma omp single
            intermediate[n] = u[n][j]; // The first down term
            // Find the prefix sums of intermediate vals - this adds in the below term as we go.
            // Outer loop is the same at all sizes
            for (i = 1; i < n + 1; i *= 2) {
                double div_pow = pow(div, i);
            // But the inner loop can be split up.
#pragma omp for schedule(static)
                for (int r = n; r >= i; r--) {
                    new_intermediate[n - r] = intermediate[n - r] + div_pow * intermediate[n - r + i];
                }

#pragma omp for schedule(static)
                for (int r = 0; r < n; r++) {
                    intermediate[r] = new_intermediate[r];
                }
            }
#pragma omp for schedule(static)
            for (i = 1; i < n; i++) {
                priv_modulus += (intermediate[i] - u[i][j]) * (intermediate[i] - u[i][j]);
                u[i][j] = intermediate[i];
            }
        }
#pragma omp atomic
        modulus += priv_modulus;
#pragma omp barrier
    }
}
    clock_t end = clock();

//    printf("Final (after %d ittrs):\n", ittr_count);
//    print2d(u, n+1, m+1);

    printf("Took %f ms to solve with n = %d, p = %d\n", ((double) (end-start)) / (CLOCKS_PER_SEC / 1000.), n, num_t );
    free2d(u);
}
