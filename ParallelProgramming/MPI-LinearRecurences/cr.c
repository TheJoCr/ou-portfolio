#define MASTER_RANK 0
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double a;
    double b;
    double c;
} ABCTuple;

void doEquationUpdate(
        ABCTuple *prev,
        ABCTuple *this,
        ABCTuple *next,
        ABCTuple *result
);

double getD(
        ABCTuple prev, double prevD,
        ABCTuple this, double thisD,
        ABCTuple next, double nextD
);

double solve (
    ABCTuple this,
    double thisD,
    double prevResult,
    double nextResult
);

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

void phase2LocalCR(
        int localNumIntervals,
        ABCTuple localEquations[],
        double localDs[],
        MPI_Datatype mpi_abc_type
);

void phase2TreeCR(
        int localNumIntervals,
        ABCTuple localEquations[],
        double localDs[],
        MPI_Datatype mpi_abc_type
);

void phase3TreeBackSub(
        int localNumIntervals,
        const ABCTuple localEquations[],
        const double localDs[],
        double next_u[]
);

void phase3LocalBackSub(
        int localNumIntervals,
        const ABCTuple localEquations[],
        const double localDs[],
        double next_u[]
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

    int globalNumIntervals = (int) strtol( argv[1], NULL, 10 );

    mpi_error_code =
            MPI_Init(&argc, &argv);

    mpi_error_code =
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    mpi_error_code =
            MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    double startTime, endTime;
    if( my_rank == MASTER_RANK ) {
        startTime = MPI_Wtime();
    }

    //First, an MPI type for the ABC Tuple
    int blocklengths[] = {1,1,1};
    MPI_Datatype types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Datatype mpi_abc_type;
    MPI_Aint     offsets[3];

    offsets[0] = offsetof(ABCTuple, a);
    offsets[1] = offsetof(ABCTuple, b);
    offsets[2] = offsetof(ABCTuple, c);

    MPI_Type_create_struct(3, blocklengths, offsets, types, &mpi_abc_type);
    MPI_Type_commit(&mpi_abc_type);

    // Grid size in x direction
    double h = 2.0 / (globalNumIntervals - 2);
    // Grid size in t direction
    double k = .01;

    int localNumIntervals = globalNumIntervals / number_of_processes;

    // First, every processor is in charge of generating their chunk of the
    // tridiagonal matrix and their chunk of u.

    // They are stored in rows of 3 parameters to one result - arrays of EqnVals
    ABCTuple *localEquations
            = malloc(localNumIntervals * sizeof( ABCTuple ));

    ABCTuple *newLocalEquations
            = malloc(localNumIntervals * sizeof( ABCTuple ));

    double *previous_u
            = malloc( localNumIntervals * sizeof( double ) );

    const double offDiagValue = -(k/(h*h));
    const double onDiagValue = 1 + 2 * (k/(h*h));

    for( int i = 0; i < localNumIntervals; i ++ ) {
        double x = ( my_rank * localNumIntervals + i ) * h;

        previous_u[i] = x * (2-x);

        localEquations[i].a = offDiagValue;
        localEquations[i].b = onDiagValue;
        localEquations[i].c = offDiagValue;
    }

    if( my_rank == 0 ) {
        localEquations[0].a = 0;
    }
    if( my_rank == number_of_processes - 1 ) {
        localEquations[localNumIntervals - 2].c = 0;
    }

    if( my_rank == number_of_processes - 1) {
        localEquations[localNumIntervals - 1].a = 0;
        localEquations[localNumIntervals - 1].b = 1;
        localEquations[localNumIntervals - 1].c = 0;

        previous_u[localNumIntervals - 1] = 0;
    }

    printArr( "Original Xs", "x", localNumIntervals, previous_u);

    double* localDs = malloc( localNumIntervals * sizeof(double) );
    double* next_u = malloc( localNumIntervals * sizeof(double) );
    memcpy( localDs, previous_u, localNumIntervals * sizeof(double) );
    for(int i = 0; i < localNumIntervals; i ++)
        next_u[i] = 0;

    for( int i = 0; i < 10; i ++ ) {

        for(int d = 0; d < localNumIntervals; d ++)
            newLocalEquations[d] = localEquations[d];

        // Part A: Locally Calculate d values using a similar procedure to I.A
        phase2LocalCR(localNumIntervals, newLocalEquations, localDs, mpi_abc_type);

        // Part B: Tree Pattern - do the computes across the whole processor tree now, similar to I.B
        phase2TreeCR(localNumIntervals, newLocalEquations, localDs, mpi_abc_type);

        // === PHASE III: BACKSUB ===
        // Part A: Back Sub across processors
        phase3TreeBackSub(localNumIntervals, newLocalEquations, localDs, next_u);

        // Part B: Back Sub locally within processors
        phase3LocalBackSub(localNumIntervals, newLocalEquations, localDs, next_u);

        // Finally, recombine the results.
        printArr(getHeader("I", i + 1), "x", localNumIntervals, next_u);

        for(int j = 0; j < localNumIntervals; j ++)
            localDs[j] = next_u[j];
    }

    free( localEquations );
    free( previous_u );
    free( localDs );
    free( next_u );


    if( my_rank == MASTER_RANK ) {
        endTime = MPI_Wtime();
        double totalTime = endTime - startTime;
        printf("CR, %d, %d, %f\n",
               number_of_processes, globalNumIntervals - 1, totalTime );
    }

    mpi_error_code =
            MPI_Finalize();
}

void doEquationUpdate(
        ABCTuple *prev,
        ABCTuple *this,
        ABCTuple *next,
        ABCTuple *result
) {
    double e = - this->a / prev->b;
    double f = - this->c / next->b;

    double a_prime = prev->a * e;
    double b_prime = this->b + prev->c * e + next->a * f;
    double c_prime = next->c * f;

    result->a = a_prime;
    result->b = b_prime;
    result->c = c_prime;
}

double getD(
        ABCTuple prev, double prevD,
        ABCTuple this, double thisD,
        ABCTuple next, double nextD
)
{
    double e = - this.a / prev.b;
    double f = - this.c / next.b;

    return thisD + prevD * e + nextD * f;
}


double solve (
        ABCTuple this,
        double thisD,
        double prevResult,
        double nextResult
)
{
    return (thisD - this.a * prevResult - this.c * nextResult ) / this.b;
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
        for (int i = 0; i < number_of_processes * localNumIntervals - 1; i++) {
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
    char* buffer = malloc( 15 * sizeof(char) );
    sprintf( buffer, "=== %s=%d ===", variable, value);
    return buffer;
}

void phase2LocalCR(
        int localNumIntervals,
        ABCTuple localEquations[],
        double localDs[],
        MPI_Datatype mpi_abc_type
) {
    MPI_Status s;
    int number_of_processes, my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    for( int j = 1; j < localNumIntervals; j *= 2 ) {
        int toMove = j - 1;
        ABCTuple abcTupleFromNextProc;
        double fromNextProc;

        // Pass it around.
        MPI_Sendrecv(
                &localDs[toMove], 1, MPI_DOUBLE, (number_of_processes + my_rank - 1 ) % number_of_processes, 0,
                &fromNextProc, 1, MPI_DOUBLE, (my_rank + 1) % number_of_processes, 0,
                MPI_COMM_WORLD, &s
        );

        // Make sure we also get the needed set of ABC coeffs
        MPI_Sendrecv(
                &localEquations[toMove], 1, mpi_abc_type, (number_of_processes + my_rank - 1 ) % number_of_processes, 0,
                &abcTupleFromNextProc, 1, mpi_abc_type, (my_rank + 1) % number_of_processes, 0,
                MPI_COMM_WORLD, &s
        );

        // Now update the values at every j things. Skip the last for now...
        int stride = j * 2;
        for( int i = stride - 1; i < localNumIntervals - j; i += stride) {
            // i is the center
            // i + j is the next one
            // i - j is the next one
            localDs[i] = getD(
                    localEquations[i - j], localDs[i-j],
                    localEquations[i], localDs[i],
                    localEquations[i + j], localDs[i + j]
            );

            doEquationUpdate(
                    &localEquations[i - j],
                    &localEquations[i],
                    &localEquations[i + j],
                    &localEquations[i]
            );
        }

        if( my_rank + 1 < number_of_processes) {
            //Here
            // numIntervals - j is the center
            // numIntervals - 2*j is the previous
            // And the one from the next processor next
            // Still we overwrite the previous value
            int i = localNumIntervals - 1;

            localDs[i] = getD(
                    localEquations[i - j], localDs[i - j],
                    localEquations[i], localDs[i],
                    abcTupleFromNextProc, fromNextProc
            );

            doEquationUpdate(
                    &localEquations[i - j],
                    &localEquations[i],
                    &abcTupleFromNextProc,
                    &localEquations[i]
            );
        }
    }
}

void phase2TreeCR(
        int localNumIntervals,
        ABCTuple localEquations[],
        double localDs[],
        MPI_Datatype mpi_abc_type
)
{
    MPI_Status s;
    int number_of_processes, my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    for( int j = 1; j < number_of_processes/2; j *= 2 )  {
        int stride = j * 2;
        int localIndex = localNumIntervals - 1;
        ABCTuple *localTuple = &localEquations[localIndex];
        double *localD = &localDs[localIndex];

        for( int i = j - 1; i < number_of_processes - stride; i += stride) {
            if( my_rank == i ) {
                //Send up
                MPI_Send(localTuple, 1, mpi_abc_type, my_rank + j, j, MPI_COMM_WORLD );
                MPI_Send(localD, 1, MPI_DOUBLE, my_rank + j, j, MPI_COMM_WORLD );
            }
        }
        for( int i = j - 1 + stride; i < number_of_processes; i += stride) {
            if( my_rank == i ) {
                //Send down
                MPI_Send(localTuple, 1, mpi_abc_type, my_rank - j, j, MPI_COMM_WORLD );
                MPI_Send(localD, 1, MPI_DOUBLE, my_rank - j, j, MPI_COMM_WORLD );
            }
        }

        if( my_rank % stride == stride - 1 &&
            my_rank != number_of_processes - 1) {

            ABCTuple below;
            ABCTuple above;

            double belowD;
            double aboveD;

            // Then I am a receiver. I get the messages
            // Receive from below
            MPI_Recv(&below, 1, mpi_abc_type, my_rank - j, j, MPI_COMM_WORLD, &s);
            MPI_Recv(&belowD, 1, MPI_DOUBLE, my_rank - j, j, MPI_COMM_WORLD, &s);

            // Receive from above
            MPI_Recv(&above, 1, mpi_abc_type, my_rank + j, j, MPI_COMM_WORLD, &s);
            MPI_Recv(&aboveD, 1, MPI_DOUBLE, my_rank + j, j, MPI_COMM_WORLD, &s);

            //And process
            localDs[localIndex] = getD(
                    below, belowD,
                    *localTuple, *localD,
                    above, aboveD
            );

            doEquationUpdate(
                    &below,
                    localTuple,
                    &above,
                    localTuple
            );
        }
    }
}

void phase3TreeBackSub(
        int localNumIntervals,
        const ABCTuple localEquations[],
        const double localDs[],
        double next_u[]
)
{
    MPI_Status s;
    int number_of_processes, my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    // The 'middle' answer is d/b
    if( my_rank == number_of_processes / 2 - 1 ) {
        //Then the last one = I have is the seminal value.
        // update my local values accordingly.
        next_u[localNumIntervals - 1] =
                localDs[localNumIntervals - 1] / localEquations[localNumIntervals - 1].b;
    }
    else if( number_of_processes == 1 ) {
        int i = localNumIntervals/2 - 1;
        next_u[i] = localDs[i] / localEquations[i].b;
    }



    // We do the inverse of the earlier tree structure, where senders receive and vice versa.
    for( int j = number_of_processes / 2; j > 0; j /= 2 )  {
        int stride = j * 2;
        // We always work with the end
        int localIndex = localNumIntervals - 1;
        ABCTuple localTuple = localEquations[localIndex];
        double localD = localDs[localIndex];

        for( int i = j - 1; i < number_of_processes; i += stride) {
            if( my_rank == i ) {
                double nextX = 0;
                double prevX = 0;

                if( i > j - 1) {
                    // Then we can receive from below
                    MPI_Recv(&prevX, 1, MPI_DOUBLE, my_rank - j, j, MPI_COMM_WORLD, &s);
                }
                if( i < number_of_processes - stride ) {
                    // Then we can receive from above
                    MPI_Recv(&nextX, 1, MPI_DOUBLE, my_rank + j, j, MPI_COMM_WORLD, &s);
                }

                // if we don't receive anything, we just let them chill out at 0.

                // Now find your local solution
                next_u[localIndex] = solve(
                        localTuple,
                        localD,
                        prevX,
                        nextX
                );
            }
        }
        if( my_rank % stride == stride - 1 &&
            my_rank != number_of_processes - 1) {
            // Then I am a sender. I get the messages
            // Send from below
            MPI_Send(&next_u[localIndex], 1, MPI_DOUBLE, my_rank + j, j, MPI_COMM_WORLD );
            // Send from above
            MPI_Send(&next_u[localIndex], 1, MPI_DOUBLE, my_rank - j, j, MPI_COMM_WORLD );
        }
    }
}


void phase3LocalBackSub(
        int localNumIntervals,
        const ABCTuple localEquations[],
        const double localDs[],
        double next_u[]
)
{
    MPI_Status s;
    int number_of_processes, my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    for( int j = localNumIntervals/2; j > 0; j /= 2 ) {
        // Now update the values at every j steps. Skip the first for now...
        int stride = j * 2;
        // Start and end indexes depend on where we are in the processor chain.
        int start = my_rank == 0 ?  stride + j - 1 : j - 1;
        int end = my_rank == number_of_processes - 1 ? localNumIntervals - 1 - j: localNumIntervals - 1;

        for( int i = start; i < end; i += stride) {
            // i is the center
            // i + j is the next one
            // i - j is the next one
            next_u[i] = solve(
                    localEquations[i],
                    localDs[i],
                    next_u[i-j],
                    next_u[i+j]
            );
        }

        //Special updates for the first and last proc, since they have 0s at the beginning and end.
        if( my_rank == 0 )
            next_u[j - 1] = solve(
                    localEquations[ j - 1 ],
                    localDs[ j - 1 ],
                    0,
                    next_u[ 2*j - 1 ]
            );

        if( my_rank == number_of_processes - 1)
            next_u[localNumIntervals - j - 1] = solve(
                    localEquations[localNumIntervals - j - 1],
                    localDs[localNumIntervals - j - 1],
                    next_u[localNumIntervals - 2*j - 1],
                    0
            );

        // Sometimes, we need information from the previous processor. Send that here:
        int toMove = localNumIntervals - 1;
        double fromPrevProc = 0;

        if( my_rank < number_of_processes - 1 ) {
            // I need to pass info forward
            MPI_Send( &next_u[toMove], 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD );
        }

        if( my_rank > 0 ) {
            // I need information from the previous
            MPI_Recv( &fromPrevProc, 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, &s );
        }

        //Here:
        // j - 1 is the center
        // the previous x came from the previous processor
        // 2*j - 1 is the next point
        int i = j - 1;
        next_u[i] = solve(
                localEquations[i],
                localDs[i],
                fromPrevProc,
                next_u[i + j]
        );
    }
}