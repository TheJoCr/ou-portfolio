#define MASTER_RANK 0
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

void randomMat(
        int numRows, int numCols, double ** mat );

void multMat(
        int numRows1, int numCols1, double **mat1,
        int numRows2, int numCols2, double **mat2,
        double **product );

void printMat(
        int numRows, int numCols, double **mat );

void multMatCumm(
        int numRows1, int numCols1, double **mat1,
        int numRows2, int numCols2, double **mat2,
        double **product);

void malloc2ddouble(double ***array, int n, int m);

int free2ddouble(double ***array);

int main (int argc, char* argv[])
{ /* main */
    int        number_of_processes;  /* Number of processes in this run    */
    int        my_rank;              /* Rank of this process               */
    int        mpi_error_code;       /* Error code returned by MPI call    */

    if( argc != 2 ) {
        fprintf( stderr, "Must supply two args! Usage: Homework2 matixSize" );
        exit(1);
    }

    int size = (int) strtol( argv[1], NULL, 10 );

    mpi_error_code =
            MPI_Init(&argc, &argv);

    mpi_error_code =
            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    mpi_error_code =
            MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    // Allocate memory that will only be populated on the master
    double **mat1 = NULL;
    double **mat2 = NULL;
    double **expectedResult = NULL;
    double **actualResult = NULL;

    double start, end;

    if( my_rank == MASTER_RANK ) {
        printf("Generating 2 random matrices of size %d\n", size);
        printf("First:\n");

        malloc2ddouble(&mat1, size, size);
        randomMat(size, size, mat1);
        //printMat( size, size, mat1);

        printf("Second:\n");
        malloc2ddouble(&mat2, size, size);
        randomMat(size, size, mat2);
        //printMat( size, size, mat2);

        // A square matrix times a square matrix produces a
        // same size matrix
        //if( size < 1025) {
        malloc2ddouble( &expectedResult, size, size);
            multMat(size, size, mat1,
                    size, size, mat2,
                    expectedResult);
        //}

        malloc2ddouble( &actualResult, size, size );

        printf("Finished generated matrices and calculating expected product.\n");
        //printMat( size, size, expectedResult );
        start = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Now, do the multiplication in parallel
#ifdef FOX
    // create a cartesian grid
    MPI_Comm comm_grid;
    MPI_Comm comm_myRow;
    MPI_Comm comm_myCol;
    int gridSize = (int) sqrt(number_of_processes);
    int gridSizes[] = {gridSize, gridSize};
    int wrap[] = {true,  true};
    int rowDim[] = {false, true};
    int colDim[] = {true, false};

    MPI_Cart_create(MPI_COMM_WORLD, 2, gridSizes, wrap, true, &comm_grid);
    MPI_Cart_sub(comm_grid, rowDim, &comm_myRow);
    MPI_Cart_sub(comm_grid, colDim, &comm_myCol);

    int blockSize = size / gridSize;

    MPI_Datatype chunkType, chunkTypeMaster, chunkTypeLocal;
    int globalSize[2] = {size, size};
    int subsizes[2] = {blockSize, blockSize};
    int starts[2] = {0,0};

    MPI_Type_create_subarray(2, globalSize, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &chunkType);
    MPI_Type_create_resized(chunkType, 0, blockSize * sizeof(double), &chunkTypeMaster);
    MPI_Type_commit(&chunkTypeMaster);

    MPI_Type_contiguous(blockSize * blockSize, MPI_DOUBLE, &chunkTypeLocal);
    MPI_Type_commit(&chunkTypeLocal);

    //With that created, we can scatter easily.

    int * scatterCounts = malloc( number_of_processes * sizeof(int) );
    int * scatterOffsets = malloc( number_of_processes * sizeof(int) );

    int offset = 0;
    for( int i = 0; i < gridSize; i ++) {
        for( int j = 0; j < gridSize; j ++) {
            int p = i * gridSize + j;
            scatterCounts[p] = 1;
            scatterOffsets[p] = offset;
            offset ++;
        }
        offset += (blockSize - 1) * gridSize;
    }

    double **myChunkOfMat1, **currentChunkOfMat1, **currentChunkOfMat2, **currentResult;
    malloc2ddouble( &myChunkOfMat1, blockSize, blockSize);
    malloc2ddouble( &currentChunkOfMat1, blockSize, blockSize);
    malloc2ddouble( &currentChunkOfMat2, blockSize, blockSize);
    malloc2ddouble( &currentResult, blockSize, blockSize);

    for(int i = 0; i < blockSize; i ++) {
        for(int j = 0; j < blockSize; j++) {
            currentResult[i][j] = 0;
        }
    }
    //First, send everyone their chunks of the matrix.
    void* mat1Ref = mat1 == NULL ? NULL : &(mat1[0][0]);
    void* mat2Ref = mat2 == NULL ? NULL : &(mat2[0][0]);
    MPI_Scatterv(mat1Ref, scatterCounts, scatterOffsets, chunkTypeMaster, &(myChunkOfMat1[0][0]), 1, chunkTypeLocal, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Scatterv(mat2Ref, scatterCounts, scatterOffsets, chunkTypeMaster, &(currentChunkOfMat2[0][0]), 1, chunkTypeLocal, MASTER_RANK, MPI_COMM_WORLD);

//    for( int p = 0; p < number_of_processes; p ++) {
//        if( my_rank == p) {
//            printf("Local mat1 (Rank %d):\n", p);
//            printMat(blockSize, blockSize, myChunkOfMat1);
//        }
//        MPI_Barrier( MPI_COMM_WORLD);
//    }

    // Get locations for sending and receiving chunks of B
    int myRow, myCol, prevInCol, nextInCol;
    MPI_Cart_shift(comm_myCol, 0, 1, &prevInCol, &nextInCol );
    MPI_Comm_rank(comm_myRow, &myCol);
    MPI_Comm_rank(comm_myCol, &myRow);

    MPI_Barrier(MPI_COMM_WORLD);

    for( int i = 0; i < gridSize; i ++) {
        // First - Broadcast A across the row
        int root = (myRow - i + gridSize) % gridSize;
        if( myCol == root ) {
            // Copy from myChunkOfMat1 into the current - that's what we're all using this time
            memcpy(&(currentChunkOfMat1[0][0]), &(myChunkOfMat1[0][0]), sizeof(double) * blockSize * blockSize);

        }
        MPI_Bcast( &(currentChunkOfMat1[0][0]), 1, chunkTypeLocal, root, comm_myRow);

        // Second - Do the multiplication by local chunk of B and add it to the running result
        multMatCumm(
                blockSize, blockSize, currentChunkOfMat1,
                blockSize, blockSize, currentChunkOfMat2,
                currentResult );

        if(my_rank == 0) {
            printf("Multiplied:\n");
            printMat( blockSize, blockSize, currentChunkOfMat1);
            printf("By:\n");
            printMat( blockSize, blockSize, currentChunkOfMat2);
        }

        // Third - Send/Recv_replace Chunk of B up.
        MPI_Status status;
        MPI_Sendrecv_replace(&(currentChunkOfMat2[0][0]), 1, chunkTypeLocal, nextInCol, 0, prevInCol, 0, comm_myCol, &status);
    }

    if( my_rank == 0 ) {
        printf("Gathering everything\n");
        printMat( blockSize, blockSize, currentResult);
    }

    //Finally, gather everything in the master
    void *actualRef = actualRef == NULL ? NULL : &(actualResult[0][0]);
    MPI_Gatherv(&(currentResult[0][0]), 1, chunkTypeLocal, actualRef, scatterCounts, scatterOffsets, chunkTypeMaster, MASTER_RANK, MPI_COMM_WORLD);
#else
    // And local data for doing manipulations on each processor
    int localNumCols = size/number_of_processes;
    double localCols[size][localNumCols];
    double localResult[size][localNumCols];

    int localNumRows = size/number_of_processes;
    double localRows[localNumRows][size];
    double newLocalRows[localNumRows][size];

    MPI_Datatype col1, col2, coltypeMaster, coltypeLocal, rowtype;

    // This is a type for sending rows that will allow us to
    // scatter and gather them
    MPI_Type_contiguous(size, MPI_DOUBLE, &rowtype);
    MPI_Type_commit(&rowtype);

    // And a type for sending columns. Note the size trick
    MPI_Type_vector(size, 1, size, MPI_DOUBLE, &col1);
    MPI_Type_create_resized(col1, 0, 1*sizeof(double), &coltypeMaster);
    MPI_Type_commit(&coltypeMaster);

    MPI_Type_vector(size, 1, localNumCols, MPI_DOUBLE, &col2);
    MPI_Type_create_resized(col2, 0, 1*sizeof(double), &coltypeLocal);
    MPI_Type_commit(&coltypeLocal);

    // First, distribute the rows of A using the type defined above
    // Note that only the sender uses the row type - receivers just
    // get the data as if it were a vector.
    MPI_Scatter(&(mat1[0][0]), localNumRows, rowtype,
                &(localRows[0][0]), localNumRows, rowtype,
                MASTER_RANK, MPI_COMM_WORLD );

    // Then the cols of B using the same trick.
    MPI_Scatter(&(mat2[0][0]), localNumCols, coltypeMaster,
                &(localCols[0][0]), localNumCols, coltypeLocal,
                MASTER_RANK, MPI_COMM_WORLD );

    for(int i = 0; i < number_of_processes; i ++) {
        // Do the local multiplication into the correct location in the local product
        // Everyone starts with row my_rank and column my_rank
        // After every increment they have row (my_rank + i) and column my_rank
        int rowOffset = (my_rank - i + number_of_processes) % number_of_processes ;
        rowOffset *= localNumRows;
        multMat(
                localNumRows, size, localRows,
                size, localNumCols, localCols,
                &localResult[rowOffset]
        );

        //And send/receive the message to the next processor
        MPI_Status status;

        int nextProc = (my_rank + 1) % number_of_processes;
        int prevProc = (my_rank - 1) % number_of_processes;
        MPI_Sendrecv_replace(
                &(localRows[0][0]), localNumRows, rowtype,
                nextProc, 0, prevProc, 0,
                MPI_COMM_WORLD, &status
        );
    }

    //Concatenate the resulting matrices
    MPI_Gather(
            &(localResult[0][0]), localNumCols, coltypeLocal,
            &(actualResult[0][0]), localNumCols, coltypeMaster,
            MASTER_RANK, MPI_COMM_WORLD
    );
#endif
    if( my_rank == MASTER_RANK ) {
        end = MPI_Wtime();
        //printf( "Distributed result:\n" );
        //printMat(size, size, actualResult);
        printf( "Finished computing. took %.4f seconds\n", end - start );
        //verify the results.
        if( size < 1025 ) {
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    double delta = fabs(expectedResult[i][j] - actualResult[i][j]);
                    if (delta > .01) {
                        printf("Does not match expected!\n");
                        exit(1);
                    }
                }
            }
        }
        printf("Results Verified.\n");
    }

    mpi_error_code =
            MPI_Finalize();
}


void randomMat(
        int numRows, int numCols, double ** mat ) {
    for(int r = 0; r < numRows; r++) {
        for(int c= 0; c < numCols; c ++) {
            //uniform between -1 and 1
            mat[r][c] = 2 * ((float)rand() / (float)RAND_MAX) - 1;
            //mat[r][c] = rand() % 4;
        }
    }
}

void multMat(
        int numRows1, int numCols1, double **mat1,
        int numRows2, int numCols2, double **mat2,
        double ** product )
{
    for(int r = 0; r < numRows1; r++) {
        for(int c = 0; c < numCols2; c ++) {
            product[r][c] = 0;
            for( int i = 0; i < numRows2; i++){
                product[r][c] += mat1[r][i] * mat2[i][c];
            }
        }
    }
}

void multMatCumm(
        int numRows1, int numCols1, double ** mat1,
        int numRows2, int numCols2, double ** mat2,
        double ** product )
{
    for(int r = 0; r < numRows1; r++) {
        for(int c = 0; c < numCols2; c ++) {
            for( int i = 0; i < numRows2; i++){
                product[r][c] += mat1[r][i] * mat2[i][c];
            }
        }
    }
}

void printMat(
        int numRows, int numCols, double ** mat )
{
    for(int r = 0; r < numRows; r++) {
        printf("[");
        for(int c= 0; c < numCols; c ++) {
            printf( "%5.2f", mat[r][c]);
            if( c + 1 != numCols ) {
                printf(" ");
            }
        }
        printf("]\n");
    }
}

void malloc2ddouble(double ***array, int n, int m) {
    // allocate contiguous items
    double *p = (double *)malloc(n*m*sizeof(double));

    // allocate row pointers into the memory
    (*array) = (double **)malloc(n*sizeof(double*));

    // and point them
    for (int i=0; i<n; i++) {
        (*array)[i] = &(p[i * m]);
    }
}

int free2ddouble(double ***array) {
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}
