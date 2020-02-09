#include <stdio.h>
#include <time.h>
#include <stdlib.h>

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

    int m = atoi(args[1]); // choose a value for m
    int n = m ; // we choose the same value for n
    double h = (b - a)/n; // interval size on x-axis
    double l = (d - c)/m; // interval size on y-axis

    // init xs init ys us
    double x[n + 1];
    double y[m + 1];

    // Tricky initialization of u:
    double** u = alloc2d(n+1, m+1);

    for (i = 0; i <= n; i++) x[i] = a + i*h;// interval points on x-axis
    for (j = 0; j <= m; j++) y[j] = c + j*l; // interval points on y axis

    // we want to find u(xi, yj), for i=1, 2, n - 1, and j = 1, 2, ..., m - 1
    // In the following we denote u(xi, yj) as u[i][j].
    // The values of u[i][j] on the boundary of the rectangle are given as follows:
    for (i = 1; i <= n; i++) u[i][0] = g(x[i], y[0]);
    for (i = 1; i <= n; i++) u[i][m] = g(x[i], y[m]);
    for (j = 0; j <= m; j++) u[0][j] = g(x[0], y[j]);
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
    while (modulus > 0.0000001) {
        ittr_count ++;
        modulus = 0;
        // Because we store values of G in the matrix, there is no need to make this complicated.
        for (j = m - 1; j > 0; j --) {
            for (k = n - 1; k >= 1; k --) {
                double new_val = div *
                          ( minhsq * f(x[k], y[j])
                            + u[k - 1][j] // UP
                            + lambda * u[k][j + 1] // RIGHT
                            + u[k + 1][j] // DOWN
                            + lambda * u[k][j - 1] // LEFT
                          );
                modulus += (new_val - u[k][j]) * (new_val - u[k][j]);
                u[k][j] = new_val;
            }
        }

    }
    clock_t end = clock();

//    printf("Final (after %d ittrs):\n", ittr_count);
//    print2d(u, n+1, m+1);

    printf("Took %f ms to calc (serial) with n = %d\n", ((double) (end-start)) / (CLOCKS_PER_SEC / 1000.), n );
    free2d(u);
}