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
