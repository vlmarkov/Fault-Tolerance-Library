#include <math.h>
#include <stdio.h>

#include <mpi.h>

#ifdef ULFM_SUPPORT
#include <mpi-ext.h>
#endif /* ULFM_SUPPORT */

const double eps = 1E-6;
const int n0 = 100;

#ifdef ULFM_SUPPORT
void ulfm_error_handler(MPI_Comm* comm, int* err, ...)
{
    int rank, size, len;
    char errstr[MPI_MAX_ERROR_STRING];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Error_string(*err, errstr, &len);
    fprintf(stderr, "Rank %d / %d: Notified of error %s\n", rank, size, errstr);
}
#endif /* ULFM_SUPPORT */

double func(double x) { return exp(-x * x); }

#ifdef ULFM_SUPPORT
void iterative_solver(MPI_Comm comm)
{
    const double a = -4.0;
    const double b = 4.0;

    int commsize, rank;    
    int n = n0, k, k_iter, n_iter;
    double sq[2] = { 0.0 }, delta = 1;
    
    MPI_Comm comm_spare;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &commsize);

    int test_iter = 0;

    for (k = 0; delta > eps; n *= 2, k ^= 1)
    {
        iteration_start:

        k_iter = k; // Rember previos value
        n_iter = n; // Rember previos value
        
        int points_per_proc = n / commsize;
        int lb = rank * points_per_proc;
        int ub = (rank == commsize - 1) ? (n - 1) : (lb + points_per_proc - 1);

        double h = (b - a) / n;
        double s = 0.0;

        for (int i = lb; i<= ub; i++)
        {
            s += func(a + h * (i + 0.5)); 
        }

        if (rank == 0 && test_iter == 0)
        {
            raise(SIGKILL);
        }
        test_iter++;

        int rc = MPI_Allreduce(&s, &sq[k], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (MPI_ERR_PROC_FAILED == rc)
        {
            MPIX_Comm_revoke(comm);

            int allsucceeded = (rc == MPI_SUCCESS);
            MPIX_Comm_agree(comm, &allsucceeded);
            if (!allsucceeded)
            {
                /*
                 * We plan to join the shrink, thus the communicator
                 * should be marked as revoked
                 */
                MPIX_Comm_shrink(comm, &comm_spare);
                MPI_Comm_free(&comm); // Release the revoked communicator
                comm = comm_spare;
            }
            
            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &commsize);

            k = k_iter; // Set previos value
            n = n_iter; // Set previos value
            goto iteration_start;
        }

        sq[k] *= h;

        if (n > n0)
        {
            delta = fabs(sq[k] -sq[k ^ 1]) / 3.0;
        }
    }

    if (rank == 0)
    {
        printf("Result Pi: %.12f; Runge rule: EPS %e, n %d\n", sq[k] * sq[k], eps, n / 2);
    }
}
#else
void iterative_solver(MPI_Comm comm)
{
    int commsize;
    int rank;

    MPI_Comm_size(comm, &commsize);
    MPI_Comm_rank(comm, &rank);

    const double a = -4.0;
    const double b = 4.0;

    int n = n0, k;
    double sq[2], delta = 1;

    for (k = 0; delta > eps; n *= 2, k ^= 1)
    {
        int points_per_proc = n / commsize;
        int lb = rank * points_per_proc;
        int ub = (rank == commsize - 1) ? (n - 1) : (lb + points_per_proc - 1);

        double h = (b - a) / n;
        double s = 0.0;

        for (int i = lb; i<= ub; i++)
        {
            s += func(a + h * (i + 0.5)); 
        }

        MPI_Allreduce(&s, &sq[k], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        sq[k] *= h;

        if (n > n0)
        {
            delta = fabs(sq[k] -sq[k ^ 1]) / 3.0;
        }
    }

    if (rank == 0)
    {
        printf("Result Pi: %.12f; Runge rule: EPS %e, n %d\n", sq[k] * sq[k], eps, n / 2);
    }
}
#endif /* ULFM_SUPPORT */

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

#ifdef ULFM_SUPPORT
    /*
     * Create a new error handler for MPI_COMM_WORLD
     * This overrides the default MPI_ERRORS_ARE_FATAL so that ranks in this
     * communicator will not automatically abort if a failure occurs.
     */
    MPI_Errhandler ulfm_eh;
    MPI_Comm_create_errhandler(ulfm_error_handler, &ulfm_eh);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, ulfm_eh);
#endif /* ULFM_SUPPORT */

    iterative_solver(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
