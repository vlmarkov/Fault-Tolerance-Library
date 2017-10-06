#include <stdio.h>

#include <signal.h>

#include <mpi.h>
#include <mpi-ext.h>


void ulfm_error_handler(MPI_Comm* comm, int* err, ...)
{
    int rank, size, len;
    char errstr[MPI_MAX_ERROR_STRING];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Error_string(*err, errstr, &len);
    fprintf(stderr, "Rank %04d/%04d: Notified of error %s\n", rank, size, errstr);
}

void iterative_test(MPI_Comm comm)
{
    int commsize;
    int rank;

    MPI_Comm comm_spare;

    MPI_Comm_size(comm, &commsize);
    MPI_Comm_rank(comm, &rank);

    printf("I'am rank %04d/%04d, starting...\n", rank, commsize);

    int i;
    for (i = 0; i < 10; i++)
    {
        if ((i == 1) && (rank == 1))
        {
            printf("I'am rank %04d: committing suicide\n", rank);
            raise(SIGKILL);
        }

        if ((i == 2) && (rank == 2))
        {
            printf("I'am rank %04d: committing suicide\n", rank);
            raise(SIGKILL);            
        }

        int rc = MPI_Barrier(comm);

        if (MPI_ERR_PROC_FAILED == rc)
        {
            MPIX_Comm_revoke(comm);
            /*
             * About to leave: let's be sure that everybody
             * received the same information
             */
            int allsucceeded = (rc == MPI_SUCCESS);
            MPIX_Comm_agree(comm, &allsucceeded);
            if (!allsucceeded)
            {
                /*
                 * We plan to join the shrink, thus the communicator
                 * should be marked as revoked
                 */
                //MPIX_Comm_revoke(comm);
                MPIX_Comm_shrink(comm, &comm_spare);
                MPI_Comm_free(&comm); // Release the revoked communicator
                comm = comm_spare;
            }
            
            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &commsize);
        }
    }

    printf("I'am rank %04d/%04d, finishing...\n", rank, commsize);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    /*
     * Create a new error handler for MPI_COMM_WORLD
     * This overrides the default MPI_ERRORS_ARE_FATAL so that ranks in this
     * communicator will not automatically abort if a failure occurs.
     */
    MPI_Errhandler ulfm_eh;
    MPI_Comm_create_errhandler(ulfm_error_handler, &ulfm_eh);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, ulfm_eh);

    iterative_test(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
