#include "mpi_ft.h"

int nprimes = 0;
int a = 1;
int b = 10000000;
int CHECKPOINT = 0;
int RESTORE = 0;

/* 
 * is_prime_number: Returns 1 if n is a prime number and 0 otherwise. 
 *                  This function uses trial division primality test.
 */
int is_prime_number(int n)
{
    int limit = sqrt(n) + 1;
    for (int i = 2; i <= limit; i++) {
        if (n % i == 0) {
            return 0;
        }
    }

    return (n > 1) ? 1 : 0;
}

int count_prime_numbers_par(int a, int b)
{
    //int nprimes = 0;
    int nprimes_global = 0;

    /* Count '2' as a prime number */
    int commsize = get_comm_size();
    int rank = get_comm_rank();

    if (a <= 2) {
        a = 2;

        if (rank == 0) {
            nprimes = 1;
        }
    }

    for (int i = a + rank; i <= b; i += commsize) {
        if (i % 2 > 0 && is_prime_number(i)) {
            nprimes++;
        }
        sleep(1);
    }

    MPI_Reduce(&nprimes, &nprimes_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    return nprimes_global;
}

double run_parallel()
{
    count_prime_numbers_par(a, b);

    return 0.0;
}

void signal_handler(int sig)
{
    printf("Received Signal: rank %d, signal %d\n", get_comm_rank(), sig);

    if (CHECKPOINT == 1) {
        make_snapshot(nprimes);
    }

    printf("Rank: %d, make checkpoint\n", get_comm_rank());

    MPI_Finalize();
    exit(1);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    if (argc == 1) {
        printf("usage: %s [checkpoint] [restore]\n", argv[0]);
        MPI_Finalize();
        return 0;
    }

    if (argc > 1) {
        if (strcmp(argv[1], "checkpoint") == 0) {
            CHECKPOINT = 1;
        } else if (strcmp(argv[1], "restore") == 0) {
            RESTORE = 1;
            nprimes = get_snapshot();
        }
    }

    signal(SIGINT,  signal_handler);
    signal(SIGTERM, signal_handler);

    run_parallel();

    MPI_Finalize();
    return 0;
}

