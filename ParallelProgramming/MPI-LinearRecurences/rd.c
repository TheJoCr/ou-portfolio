#define MASTER_RANK 0
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>


typedef struct {
    double a;
    double b;
    double c;
    double d;
} Mat22;

typedef struct {
    double a;
    double b;
} Vec2;

const Vec2 Vec2_Identity = (Vec2) {0,0};
const Mat22 Mat22_Identity = (Mat22) {1, 0, 0, 1};

void solveLinearRecurrenceParallel(
        double Z0,
        const double *factors,
        const double *summands,
        double *output,
        int localSize
);

void solveLinearRecurrenceParallel_Mat(
        Vec2 Z0,
        const Mat22 *factors,
        Vec2 *output,
        int localSize
);

void solveLinearRecurrenceParallel_Reverse(
        double ZN,
        const double *factors,
        const double *summands,
        double *output,
        int localSize
);

void multMat2(
        const  Mat22 *a,
        const Mat22 *b,
        Mat22 *output
);

void multMatVec2(
        const Mat22 *a,
        const Vec2 *b,
        Vec2 *output
);

void addVec2 (
        const Vec2 *a,
        const Vec2 *b,
        Vec2 *output
);

void createMatType(
        MPI_Datatype *mpi_mat_type
);

void createVecType(
        MPI_Datatype *mpi_vec_type
);
// These are defined globablly for ease of use.
MPI_Datatype mpi_mat22_type, mpi_vec2_type;

char* getHeader(
        char* variable,
        int value
);

void printArr(
        char* header,
        char* varName,
        int localNumIntervals,
        double data[]
);

/**
 * Computes Ax + d
 */
void multMat2AndAddVec2(
        const Mat22 *a,
        const Vec2 *x,
        const Vec2 *d,
        Vec2 *output
);

int main (int argc, char* argv[])
{ /* main */
    int        number_of_processes;  /* Number of processes in this run    */
    int        my_rank;              /* Rank of this process               */
    int        mpi_error_code;       /* Error code returned by MPI call    */

    if( argc != 2 ) {
        fprintf( stderr, "Must supply two args! Usage: Homework3 numIntervals" );
        exit(1);
    }

    double startTime, endTime;
    if( my_rank == MASTER_RANK ) {
        startTime = MPI_Wtime();
    }

    int globalNumIntervals = (int) strtol( argv[1], NULL, 10 );

    mpi_error_code =
            MPI_Init(&argc, &argv);

    mpi_error_code =
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    mpi_error_code =
            MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    int iAmLastProc = my_rank + 1 == number_of_processes;
    int iAmFirstProc = my_rank == 0;

    if( my_rank == MASTER_RANK ) {
        //Start up logic
    }

    // Initialize global variables.
    createMatType( &mpi_mat22_type );
    createVecType( &mpi_vec2_type );

    // Grid size in x direction
    double h = 2.0 / (globalNumIntervals - 1);
    // Grid size in t direction
    double k = .01;

    int localNumIntervals = globalNumIntervals / number_of_processes;

    const double offDiagValue = -(k/(h*h));
    const double onDiagValue = 1 + 2 * (k/(h*h));

    double* as = malloc( localNumIntervals * sizeof(double) );
    double* bs = malloc( localNumIntervals * sizeof(double) );
    double* cs = malloc( localNumIntervals * sizeof(double) );
    double* ds = malloc( localNumIntervals * sizeof(double) );

    for( int i = 0; i < localNumIntervals; i ++) {
        as[i] = ( iAmFirstProc && i == 0 ) ? 0 : offDiagValue;
        bs[i] = onDiagValue;
        cs[i] = ( iAmLastProc && i + 1 == localNumIntervals ) ?
                              0 : offDiagValue;

        double x = (my_rank * localNumIntervals + i) * h;
        ds[i] = x * (2-x);
    }

    printArr( "Original Xs", "x", localNumIntervals, ds );

    // ====================== Solve for Ws ======================== //
    // First, we have to solve for our weights w, defined by the linear recurrence:
    // ( yi   ) = (   0      1   ) * ( yi-1 )
    // ( yi+1 )   (-ai/ci  bi/ci )   ( yi   )
    // Note that this does not depend on ds, so we can calculate this once.

    //Populate factors:
    Mat22 *mat_vFactors = malloc(localNumIntervals * sizeof(Mat22));

    int end = iAmLastProc ? localNumIntervals - 1 : localNumIntervals;
    for (int i = 0; i < end; i++) {
        mat_vFactors[i] = (Mat22) {               0,             1,
                                    - as[i] / cs[i], bs[i] / cs[i]};
    }

    Vec2 startVec = iAmFirstProc ? (Vec2) {1, bs[0] / cs[0]} : Vec2_Identity;

    Vec2 *vs = malloc(localNumIntervals * sizeof(Vec2));
    double *ws = malloc(localNumIntervals * sizeof(double));

    solveLinearRecurrenceParallel_Mat(
            startVec,
            mat_vFactors,
            vs,
            localNumIntervals
    );

    //Divide through to get our ws, where wi = yi/yi+1 (stored in vector vi)
    for (int i = 0; i < end; i++) {
        ws[i] = vs[i].a / vs[i].b;
    }

    if (iAmLastProc)
        ws[localNumIntervals - 1] = 0; //Not technically defined, but set it to 0 for ease.

    // These are no longer needed - they were only used in calculating ws.
    free(vs); free(mat_vFactors);

    // Allocate a set of calculation vectors used in intermediate steps.
    double *gFactors = malloc(localNumIntervals * sizeof(double));
    double *gSums = malloc(localNumIntervals * sizeof(double));
    double *gs = malloc(localNumIntervals * sizeof(double));

    double *xFactors = malloc(localNumIntervals * sizeof(double));
    double *xSums = malloc(localNumIntervals * sizeof(double));
    double *xs = malloc(localNumIntervals * sizeof(double));

    // Start the main loop:
    for( int r = 0; r < 10; r ++) {
        // There are three phases to solving a tri-diagonal system via
        // turning it into a set of linear recurrences. We've already
        // solved for Ws in an initial processing step, so we only have
        // to do the last 2 here:

        // ====================== Solve for Gs ======================== //
        // These satisfy the recurrence:
        // gi = (di - ai * g(i-1)) / (bi - ai * w(i-1))
        // g1 = d1 / b1
        for (int i = 1; i < localNumIntervals; i++) {
            double dividend = (bs[i] - as[i] * ws[i - 1]);

            gFactors[i] = -1 * as[i] / dividend;
            gSums[i] = ds[i] / dividend;
        }

        double prevW, g0;
        // Finding gi requires knowledge of the previous w -
        // this might reside on the previous processor.
        // Grab that here:
        if (!iAmLastProc) {
            MPI_Send(
                    &ws[localNumIntervals - 1], 1, MPI_DOUBLE, my_rank + 1,
                    0, MPI_COMM_WORLD
            );
        }
        if (!iAmFirstProc) {
            MPI_Recv(
                    &prevW, 1, MPI_DOUBLE, my_rank - 1,
                    0, MPI_COMM_WORLD, NULL
            );

            double base = (bs[0] - as[0] * prevW);

            gFactors[0] = -1 * as[0] / base;
            gSums[0] = ds[0] / base;
        } else {
            gFactors[0] = 1;
            gSums[0] = ds[0] / bs[0];
        }

        solveLinearRecurrenceParallel(
                g0,
                gFactors,
                gSums,
                gs,
                localNumIntervals
        );

        // ====================== Solve for Xs ======================== //

        // Finally, solve for your xs. These have a simple recurrence as follows:
        for (int i = 0; i < localNumIntervals; i++) {
            xFactors[i] = -ws[i];
            xSums[i] = gs[i];
        }
        // We have to go in reverse here.
        double xn = gs[localNumIntervals - 1];

        solveLinearRecurrenceParallel_Reverse(
                xn,
                xFactors,
                xSums,
                xs,
                localNumIntervals
        );

        printArr(getHeader("I", r + 1), "x", localNumIntervals, xs);

        // That's our result! Copy this into ds for the next iteration.
        memcpy(ds, xs, localNumIntervals * sizeof(double));
    }

    free(gFactors); free(gSums);
    free(xFactors); free(xSums);

    free(ws); free(gs); free(xs);
    free(as); free(bs); free(cs); free(ds);

    if( my_rank == MASTER_RANK ) {
        endTime = MPI_Wtime();
        printf("RD, %d, %d, %f\n",
               number_of_processes, globalNumIntervals, endTime - startTime );
    }

    mpi_error_code =
            MPI_Finalize();
}

void solveLinearRecurrenceParallel(
        const double Z0,
        const double *factors,
        const double *summands,
        double *output,
        int localSize
) {
    int number_of_processes, my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    double * leftOutput = malloc( localSize * sizeof(double) );
    double * leftM = malloc( localSize * sizeof(double) );

    double * thisOutput = malloc( localSize * sizeof(double) );
    double * thisM = malloc( localSize * sizeof(double) );

    double * nextOutput = malloc( localSize * sizeof(double) );
    double * nextM = malloc( localSize * sizeof(double) );

    for( int i = 0; i < localSize; i ++ ) {
        thisOutput[i] = summands[i];
        thisM[i] = factors[i];
    }

    if( my_rank == 0 ) {
        thisOutput[0] = Z0;
        thisM[0] = 1;
    }

    // Phase 1: Here, we send parts of our local data
    // In this phase, we only communicate with our immediate
    // neighbor.
    for( int i = 1; i <= localSize; i *= 2 ) {
        // Do the local updates
        for( int j = i; j < localSize; j ++ ) {
            nextOutput[j] = thisOutput[j] + thisM[j] * thisOutput[j - i];
            nextM[j] = thisM[j] * thisM[j-1];
        }

        //Pass the tail of our array to the right 1.
        if( my_rank < number_of_processes - 1 ) {
            MPI_Send(
                    &thisOutput[localSize - i], i, MPI_DOUBLE, my_rank + 1,
                    0, MPI_COMM_WORLD
            );
            MPI_Send(
                    &thisM[localSize - i], i, MPI_DOUBLE, my_rank + 1,
                    1, MPI_COMM_WORLD
            );
        }

        //Receive the tail of the array to the left and use it to update
        if( my_rank > 0 ) {
            MPI_Recv(
                    leftOutput, i, MPI_DOUBLE, my_rank - 1,
                    0, MPI_COMM_WORLD, NULL
            );
            MPI_Recv(
                    leftM, i, MPI_DOUBLE, my_rank - 1,
                    1, MPI_COMM_WORLD, NULL
            );

            for( int j = 0; j < i; j ++ ) {
                nextOutput[j] = thisOutput[j] + thisM[j] * leftOutput[j];
                nextM[j] = thisM[j] * leftM[j];
            }
        }
        else {
            for( int j = 0; j < i; j ++ ) {
                nextOutput[j] = thisOutput[j];
                nextM[j] = thisM[j];
            }
        }

        // Swap the next output and M with this one.
        double* temp = nextOutput;
        nextOutput = thisOutput;
        thisOutput = temp;

        temp = nextM;
        nextM = thisM;
        thisM = temp;
    }
    //These are only used in the first phase - after that they are no longer necessary.
    free( nextOutput );
    free( nextM );

    //Phase II. We now send full copies of our data to the processors on the right.
    for( int i = 2; i < number_of_processes; i *= 2) {
        if( my_rank + i < number_of_processes ) {
            // I need to send my whole array i spots to the right.
            MPI_Send( thisOutput, localSize, MPI_DOUBLE, my_rank + i,
                      0, MPI_COMM_WORLD );
            MPI_Send( thisM, localSize, MPI_DOUBLE, my_rank + i,
                      1, MPI_COMM_WORLD );
        }

        if( my_rank >= i ) {
            MPI_Recv( leftOutput, localSize, MPI_DOUBLE, my_rank - i,
                      0, MPI_COMM_WORLD, NULL );
            MPI_Recv( leftM, localSize, MPI_DOUBLE, my_rank - i,
                      1, MPI_COMM_WORLD, NULL );

            for( int j = 0; j < localSize; j ++ ) {
                thisOutput[j] = thisOutput[j] + thisM[j] * leftOutput[j];
                thisM[j] = thisM[j] * leftM[j];
            }
        }
    }

    for( int i = 0; i < localSize; i ++)
        output[i] = thisOutput[i];

    free( leftOutput );
    free( leftM );

    free( thisOutput );
    free( thisM );
}

void solveLinearRecurrenceParallel_Reverse(
        const double ZN,
        const double *factors,
        const double *summands,
        double *output,
        int localSize
) {
    int number_of_processes, my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    double * rightOutput = malloc( localSize * sizeof(double) );
    double * rightM = malloc( localSize * sizeof(double) );

    double * thisOutput = malloc( localSize * sizeof(double) );
    double * thisM = malloc( localSize * sizeof(double) );

    double * nextOutput = malloc( localSize * sizeof(double) );
    double * nextM = malloc( localSize * sizeof(double) );

    for( int i = 0; i < localSize; i ++ ) {
        thisOutput[i] = summands[i];
        thisM[i] = factors[i];
    }

    if( my_rank == number_of_processes - 1 ) {
        thisOutput[ localSize - 1 ] = ZN;
        thisM[ localSize - 1 ] = 1;
    }

    // Phase 1: Here, we send parts of our local data
    // In this phase, we only communicate with our immediate
    // neighbor.
    for( int i = 1; i <= localSize; i *= 2 ) {

        // Do the local updates
        for( int j = 0; j < localSize - i; j ++ ) {
            nextOutput[j] = thisOutput[j] + thisM[j] * thisOutput[j + i];
            nextM[j] = thisM[j] * thisM[j + i];
        }

        //Pass the head of our array (size i) to the left 1.
        if( my_rank > 0 ) {
            MPI_Send(
                    thisOutput, i, MPI_DOUBLE, my_rank - 1,
                    0, MPI_COMM_WORLD
            );
            MPI_Send(
                    thisM, i, MPI_DOUBLE, my_rank - 1,
                    1, MPI_COMM_WORLD
            );
        }
        //Receive the head of the array to the right and use it to update
        if( my_rank < number_of_processes - 1 ) {
            MPI_Recv(
                    rightOutput, i, MPI_DOUBLE, my_rank + 1,
                    0, MPI_COMM_WORLD, NULL
            );
            MPI_Recv(
                    rightM, i, MPI_DOUBLE, my_rank + 1,
                    1, MPI_COMM_WORLD, NULL
            );

            for( int j = localSize - 1; j >= localSize - i; j -- ) {
                int rightIndex = j - localSize + i; // i - 1 to 0
                nextOutput[j] = thisOutput[j] + thisM[j] * rightOutput[rightIndex];
                nextM[j] = thisM[j] * rightM[rightIndex];
            }
        }
        else {
            for( int j = localSize - 1; j >= localSize - i; j -- ) {
                nextOutput[j] = thisOutput[j];
                nextM[j] = thisM[j];
            }
        }

        // Swap the next output and M with this one.
        double* temp = nextOutput;
        nextOutput = thisOutput;
        thisOutput = temp;

        temp = nextM;
        nextM = thisM;
        thisM = temp;
    }

    //These are only used in the first phase - after that they are no longer necessary.
    free( nextOutput );
    free( nextM );

    //Phase II. We now send full copies of our data to the processors on the right.
    for( int i = 2; i < number_of_processes; i *= 2) {
        if( my_rank - i >= 0 ) {
            // I need to send my whole array i spots to the left.
            MPI_Send( thisOutput, localSize, MPI_DOUBLE, my_rank - i,
                      0, MPI_COMM_WORLD );
            MPI_Send( thisM, localSize, MPI_DOUBLE, my_rank - i,
                      1, MPI_COMM_WORLD );
        }

        if( my_rank + i < number_of_processes ) {
            MPI_Recv( rightOutput, localSize, MPI_DOUBLE, my_rank + i,
                      0, MPI_COMM_WORLD, NULL );
            MPI_Recv( rightM, localSize, MPI_DOUBLE, my_rank + i,
                      1, MPI_COMM_WORLD, NULL );

            for( int j = 0; j < localSize; j ++ ) {
                thisOutput[j] = thisOutput[j] + thisM[j] * rightOutput[j];
                thisM[j] = thisM[j] * rightM[j];
            }
        }
    }

    for( int i = 0; i < localSize; i ++)
        output[i] = thisOutput[i];

    free( rightOutput );
    free( rightM );

    free( thisOutput );
    free( thisM );
}

void solveLinearRecurrenceParallel_Mat(
        Vec2 Z0,
        const Mat22 *factors,
        Vec2 *output,
        int localSize
)
{
    int number_of_processes, my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    Vec2 * leftOutput = malloc( localSize * sizeof(Vec2) );
    Mat22 * leftM = malloc( localSize * sizeof(Mat22) );

    Vec2 * thisOutput = malloc( localSize * sizeof(Vec2) );
    Mat22 * thisM = malloc( localSize * sizeof(Mat22) );

    Vec2 * nextOutput = malloc( localSize * sizeof(Vec2) );
    Mat22 * nextM = malloc( localSize * sizeof(Mat22) );

    for( int i = 0; i < localSize; i ++ ) {
        thisOutput[i] = Vec2_Identity; // We only ever add zero, so no need to general case this.
        thisM[i] = factors[i];
    }

    if( my_rank == 0 ) {
        thisOutput[0] = Z0;
        thisM[0] = Mat22_Identity;
    }

    for( int i = 1; i <= localSize; i *= 2 ) {
        // Do the local updates
        for (int j = i; j < localSize; j++) {
            multMat2AndAddVec2( &thisM[j], &thisOutput[j - i], &thisOutput[j], &nextOutput[j] );
            multMat2( &thisM[j], &thisM[j - 1], &nextM[j] );
        }

        //Pass the tail of our array to the right 1.
        if (my_rank < number_of_processes - 1) {
            MPI_Send(
                    &thisOutput[localSize - i], i, mpi_vec2_type, my_rank + 1,
                    0, MPI_COMM_WORLD
            );
            MPI_Send(
                    &thisM[localSize - i], i, mpi_mat22_type, my_rank + 1,
                    1, MPI_COMM_WORLD
            );
        }

        //Receive the tail of the array to the left and use it to update
        if( my_rank > 0 ) {
            MPI_Recv(
                    leftOutput, i, mpi_vec2_type, my_rank - 1,
                    0, MPI_COMM_WORLD, NULL
            );
            MPI_Recv(
                    leftM, i, mpi_mat22_type, my_rank - 1,
                    1, MPI_COMM_WORLD, NULL
            );

            for( int j = 0; j < i; j ++ ) {
                multMat2AndAddVec2( &thisM[j], &leftOutput[j], &thisOutput[j], &nextOutput[j] );
                multMat2( &thisM[j], &leftM[j], &nextM[j] );
            }
        }
        else {
            for( int j = 0; j < i; j ++ ) {
                nextOutput[j] = thisOutput[j];
                nextM[j] = thisM[j];
            }
        }

        // Swap the next output and M with this one.
        Vec2* temp = nextOutput;
        nextOutput = thisOutput;
        thisOutput = temp;

        Mat22 *tempM = nextM;
        nextM = thisM;
        thisM = tempM;
    }

    free( nextOutput );
    free( nextM );

    //Phase II. We now send full copies of our data to the processors on the right.
    for( int i = 2; i < number_of_processes; i *= 2) {
        if( my_rank + i < number_of_processes ) {
            // I need to send my whole array i spots to the right.
            MPI_Send( thisOutput, localSize, mpi_vec2_type, my_rank + i,
                      0, MPI_COMM_WORLD );
            MPI_Send( thisM, localSize, mpi_mat22_type, my_rank + i,
                      1, MPI_COMM_WORLD );
        }

        if( my_rank >= i ) {
            MPI_Recv( leftOutput, localSize, mpi_vec2_type, my_rank - i,
                      0, MPI_COMM_WORLD, NULL );
            MPI_Recv( leftM, localSize, mpi_mat22_type, my_rank - i,
                      1, MPI_COMM_WORLD, NULL );

            for( int j = 0; j < localSize; j ++ ) {
                multMat2AndAddVec2( &thisM[j], &leftOutput[j], &thisOutput[j], &thisOutput[j] );
                multMat2( &thisM[j], &leftM[j], &thisM[j] );
            }
        }
    }

    // Copy into the output array.
    for( int i = 0; i < localSize; i ++)
        output[i] = thisOutput[i];

    // And free the local arrays that were used
    free( leftOutput );
    free( leftM );

    free( thisOutput );
    free( thisM );
}

/**
 * Computes Ax + d
 */
void multMat2AndAddVec2(
        const Mat22 *a,
        const Vec2 *x,
        const Vec2 *d,
        Vec2 *output
){
    Vec2 toAdd;
    multMatVec2( a, x, &toAdd );
    addVec2( d, &toAdd, output );
}

void multMat2(
        const Mat22 *a,
        const Mat22 *b,
        Mat22 *output
)
{
    //Use a temp variable to let the caller overwrite a or b with the product
    Mat22 temp;
    temp.a = a->a * b->a + a->b * b->c;
    temp.b = a->a * b->b + a->b * b->d;
    temp.c = a->c * b->a + a->d * b->c;
    temp.d = a->c * b->b + a->d * b->d;
    *output = temp;
}

void multMatVec2(
        const Mat22 *a,
        const Vec2 *b,
        Vec2 *output
)
{
    output->a = a->a * b->a + a->b * b->b;
    output->b = a->c * b->a + a->d * b->b;
}

void addVec2 (
        const Vec2 *a,
        const Vec2 *b,
        Vec2 *output
)
{
    output->a = a->a + b->a;
    output->b = a->b + b->b;
}

void createMatType(
        MPI_Datatype *mpi_mat_type
) {
    int blocklengths[] = {1,1,1,1};
    MPI_Datatype types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint     offsets[4];

    offsets[0] = offsetof(Mat22, a);
    offsets[1] = offsetof(Mat22, b);
    offsets[2] = offsetof(Mat22, c);
    offsets[3] = offsetof(Mat22, d);

    MPI_Type_create_struct(4, blocklengths, offsets, types, mpi_mat_type);
    MPI_Type_commit(mpi_mat_type);
}

void createVecType(
        MPI_Datatype *mpi_vec_type
) {
    int blocklengths[] = {1,1};
    MPI_Datatype types[] = {MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint     offsets[2];

    offsets[0] = offsetof(Vec2, a);
    offsets[1] = offsetof(Vec2, b);

    MPI_Type_create_struct(2, blocklengths, offsets, types, mpi_vec_type);
    MPI_Type_commit(mpi_vec_type);
}

void printArr(
        char* header,
        char* varName,
        int localNumIntervals,
        double data[]
) {
    int number_of_processes, my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    double *gatheredData = NULL;

    if (my_rank == MASTER_RANK) {
        gatheredData = malloc(localNumIntervals * number_of_processes * sizeof(double));
    }

    MPI_Gather(
            data, localNumIntervals, MPI_DOUBLE,
            gatheredData, localNumIntervals, MPI_DOUBLE,
            MASTER_RANK, MPI_COMM_WORLD
    );

    if (my_rank == MASTER_RANK) {
        //printf( "%s\n", header );
        for (int i = 0; i < number_of_processes * localNumIntervals; i++) {
            printf( "%f, ",
                    gatheredData[i] );
        }
        printf("\n");
        free( gatheredData );
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

char* getHeader(
        char* variable,
        int value
)
{
    char* buffer = malloc( 20 * sizeof(char) );
    sprintf( buffer, "=== %s=%d ===", variable, value);
    return buffer;
}