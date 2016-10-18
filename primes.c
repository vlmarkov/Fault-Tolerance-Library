/*
 * Prime numbers balanced
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#include <mpi.h>

const int a = 1;
const int b = 10000000;

enum {
    OPTION_INVALID    = -1,
    OPTION_CHECKPOINT = 0,
    OPTION_RESTORE    = 1
};

int PROCESS_OPTION = OPTION_INVALID;

int get_comm_rank();
int get_comm_size();

int is_prime_number(int n);

int count_prime_numbers(int a, int b);
int count_prime_numbers_par(int a, int b);

int count_prime_numbers_(int a, int b);
int count_prime_numbers_par_(int a, int b);

double run_serial();
double run_parallel();

MPI_File snapshot;

void create_file()
{
    char file_name[256] = { 0 };
    int myrank = get_comm_rank();

    sprintf(file_name,"%s_%d_snapshot.txt",__FILE__, myrank);

    MPI_File_open( MPI_COMM_WORLD,
                    file_name,
                    MPI_MODE_CREATE|MPI_MODE_WRONLY,
                    MPI_INFO_NULL,
                    &snapshot);
}

void make_snapshot(int value, int iter)
{
    MPI_Status status;
    char buf[256] = { 0 };
    sprintf(buf, "%d %d\n", value, iter);

    MPI_File_write(snapshot, buf, strlen(buf), MPI_CHAR, &status);
}

int get_comm_rank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int get_comm_size()
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

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

int count_prime_numbers(int a, int b)
{
    int nprimes = 0;

    /* Count '2' as a prime number */
    if (a <= 2) {
        nprimes = 1;
        a = 2;
    }

    /* Shift 'a' to odd number */
    if (a % 2 == 0) {
        a++;
    }

    /* Loop over odd numbers: a, a + 2, a + 4, ... , b */
    for (int i = a; i <= b; i++) {
        if (i % 2 > 0 && is_prime_number(i)) {
            nprimes++;
        }
    }

    return nprimes;
}

int count_prime_numbers_par(int a, int b)
{
    int nprimes        = 0;
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
    }

    MPI_Reduce(&nprimes, &nprimes_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return nprimes_global;
}

int count_prime_numbers_par_(int a, int b)
{
    int nprimes        = 0;
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

    create_file();

    for (int i = a + rank; i <= b; i += commsize) {
        if (i % 2 > 0 && is_prime_number(i)) {
            nprimes++;

            if ((nprimes % 50000) == 0) {
                make_snapshot(nprimes, i);
            }

        }
    }

    make_snapshot(nprimes, -1);
    MPI_File_close(&snapshot);

    MPI_Reduce(&nprimes, &nprimes_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return nprimes_global;
}


double run_serial()
{
    double t = MPI_Wtime();
    int n = count_prime_numbers(a, b);
    t = MPI_Wtime() - t;

    printf("Result (serial): %d\n", n);
    return t;
}


double run_parallel()
{
    int n;
    double t = MPI_Wtime();
    if (PROCESS_OPTION != OPTION_INVALID) {
        n = count_prime_numbers_par_(a, b);
    } else {
        n = count_prime_numbers_par(a, b);
    }
    t = MPI_Wtime() - t;

    printf("Process %d/%d execution time: %.6f\n", get_comm_rank(), get_comm_size(), t);

    if (get_comm_rank() == 0) {
        printf("Result (parallel): %d\n", n);
    }

    double tmax;
    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return tmax;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    if (argc == 1) {
        printf("usage: %s [checkpoint] [restore]\n", argv[0]);
        MPI_Finalize();
        return 0;
    }

    if (strcmp(argv[1], "checkpoint") == 0) {
        PROCESS_OPTION = OPTION_CHECKPOINT;
    } else if (strcmp(argv[1], "restore") == 0) {
        PROCESS_OPTION = OPTION_RESTORE;
    }

    // Start serial version
    double tserial = 0;
    if (get_comm_rank() == 0) {
        tserial = run_serial();
    }

    // Start parallel version
    double tparallel = run_parallel();

    if (get_comm_rank() == 0) {
        printf("Count prime numbers on [%d, %d]\n", a, b);
        printf("Execution time (serial): %.6f\n", tserial);
        printf("Execution time (parallel): %.6f\n", tparallel);
        printf("Speedup (processes %d): %.2f\n", get_comm_size(), tserial / tparallel);
    }

    MPI_Finalize();
    return 0;
}