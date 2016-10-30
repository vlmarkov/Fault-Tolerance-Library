/*
 * Prime numbers balanced
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#include <mpi.h>


/*****************************************************************************/
/* Global variables                                                          */
/*****************************************************************************/
const int a = 1;
const int b = 10000000;

enum {
    OPTION_INVALID    = -1,
    OPTION_CHECKPOINT = 0,
    OPTION_RESTORE    = 1
};

enum {
   STATE_1 = 1,
   STATE_2 = 2
};

int PROCESS_OPTION = OPTION_INVALID;
MPI_File snapshot;


/*****************************************************************************/
/* Control flow-macros                                                       */
/*****************************************************************************/
#define CHECKPOINT_INIT() create_file();
#define CHECKPOINT_FINALIZE() MPI_File_close(&snapshot);

#define CHECKPOINT_SET(name) name ## _checkpoint:
#define CHECKPOINT_GOTO(name) goto name ## _checkpoint


/*****************************************************************************/
/* Prototypes                                                                */
/*****************************************************************************/
int get_comm_rank();
int get_comm_size();

int is_prime_number(int n);

int count_prime_numbers(int a, int b);
int count_prime_numbers_par(int a, int b);

int count_prime_numbers_(int a, int b);
int count_prime_numbers_par_(int a, int b);

double run_serial();
double run_parallel();


void create_file()
{
    char file_name[256] = { 0 };

    sprintf(file_name,"%s_%d_snapshot.txt",__FILE__, get_comm_rank());

    MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  MPI_INFO_NULL,&snapshot);
}

void make_snapshot(char *data, int state)
{
    MPI_Status status;
    char buf[256] = { 0 };

    sprintf(buf, "%s:%d\n", data, state);
    MPI_File_write(snapshot, buf, strlen(buf), MPI_CHAR, &status);
}

void argument_parse(char *data, int *nprimes, int *i)
{
    // get the first token
    char *token = strtok(data, " ");
    *nprimes = atoi(token);

    // walk through other tokens
    token = strtok(NULL, " ");

    *i = atoi(token);
    *i += get_comm_size();
}

int get_snapshot(char *data)
{
    int state;
    int myrank          = get_comm_rank();
    char *buf           = NULL;
    size_t lenght       = 0;
    char file_name[256] = { 0 };

    sprintf(file_name,"%s_%d_snapshot.txt",__FILE__, myrank);

    FILE *snapshot = fopen(file_name, "r");
    if (snapshot) {
        while ((getline(&buf, &lenght, snapshot)) != -1) {
            // ...
        }

        // get the first token
        char *token = strtok(buf, ":");
        strcpy(data, token);

        // walk through other tokens
        token  = strtok(NULL, ":");
        state = atoi(token);

        if (buf) {
            free(buf);
        }
        fclose(snapshot);
        return state;
    }

    return -1;
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
    int i              = 0;
    int nprimes        = 0;
    int nprimes_global = 0;

    char checkpoint_data[256];

    /* Count '2' as a prime number */
    int commsize = get_comm_size();
    int rank = get_comm_rank();

    if (a <= 2) {
        a = 2;
        if (rank == 0) {
            nprimes = 1;
        }
    }

    i = a + rank;

    CHECKPOINT_INIT();


    if (PROCESS_OPTION == OPTION_RESTORE) {
        int state = get_snapshot(checkpoint_data);
        switch (state) 
        {
            case 1:
                argument_parse(checkpoint_data, &nprimes, &i);
                printf("state 1 rank = %d : nprimes %d i %d\n", rank, nprimes, i);
                CHECKPOINT_GOTO(STATE_1);
                break;

            case 2:
                argument_parse(checkpoint_data, &nprimes, &i);
                printf("state 2 rank = %d : nprimes %d i %d\n", rank, nprimes, i);
                CHECKPOINT_GOTO(STATE_2);

                break;
            default:
                fprintf(stderr, "can't read checkpoint\n");
                //exit(1);
        }

    }

    CHECKPOINT_SET(STATE_1);

    for ( ; i <= b; i += commsize) {
        if (i % 2 > 0 && is_prime_number(i)) {
            nprimes++;

            if ((nprimes % 50000) == 0) {
                sprintf(checkpoint_data, "%d %d", nprimes, i);
                make_snapshot(checkpoint_data, STATE_1);
                //MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
            }

        }
    }

    sprintf(checkpoint_data, "%d %d", nprimes, i);
    make_snapshot(checkpoint_data, STATE_2);

    CHECKPOINT_SET(STATE_2);
    CHECKPOINT_FINALIZE();

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
        printf("Result (parallel)        : %d\n", n);
    }

    double tmax;
    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return tmax;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    if (argc < 2) {
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
        printf("Execution time (serial)  : %.6f\n", tserial);
        printf("Execution time (parallel): %.6f\n", tparallel);
        printf("Speedup (processes %d)   : %.2f\n", get_comm_size(), tserial / tparallel);
    }

    MPI_Finalize();
    return 0;
}
