/*****************************************************************************/
/* C - check                                                                 */
/* P - point                                                                 */
/* L - library                                                               */
/*****************************************************************************/

#include "checkpoint_lib.h"

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>

#include <sys/stat.h>

#include <inttypes.h>

#include "zlib.h"

/*****************************************************************************/
/* Global variables                                                          */
/*****************************************************************************/
enum {
    CPL_CHECKPOINT_MODE = 1,
    CPL_RECOVERY_MODE   = 2
};


void   **cpl_checkpoint_table;

double cpl_start_time    = 0.0;
double cpl_start_time_local = 0.0;

double cpl_save_time       = 0.0;
double cpl_save_time_local = 0.0;

int cpl_snapshot_counter = 0;

int cpl_size    = 0;
int cpl_time    = 0;
int cpl_counter = 0;

int cpl_run_options = 0;


static double wtime_()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

static int get_comm_rank_()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

static int get_comm_size_()
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}


/*****************************************************************************/
/* Work with files                                                           */
/*****************************************************************************/
void CPL_FILE_OPEN(MPI_File *snapshot, int phase)
{
    MPI_Status status;

    static int counter  = 0;
    char file_name[256] = { 0 };
    char file_path[256] = { 0 };
    /* 
     * 1_1 => [PHASE_OF_CALCULATION]_[COUNTER]
     *
     * PHASE_OF_CALCULATION - global state, need for 'goto'
     * COUNTER              - each PHASE_OF_CALCULATION could reach many times
     */

    sprintf(file_path,"%s/%d", SNAPSHOT_DIR_NAME, get_comm_rank_());
    mkdir(file_path, 0777);

    sprintf(file_name,"%s/%d/%d_%d", SNAPSHOT_DIR_NAME, get_comm_rank_(), phase, counter++);

    MPI_File_open(MPI_COMM_WORLD,
                  file_name,
                  MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  MPI_INFO_NULL, snapshot);

    cpl_save_time_local = wtime_();
}

void CPL_FILE_CLOSE(MPI_File *snapshot)
{
    MPI_Status status;

    cpl_snapshot_counter += 1;

    cpl_save_time_local = wtime_() - cpl_save_time_local;
    cpl_save_time += cpl_save_time_local;

    double elapsed_time = wtime_() - cpl_start_time_local;

    // Write meta-data
    MPI_File_write(*snapshot, &elapsed_time, 1, MPI_DOUBLE, &status);
    MPI_File_write(*snapshot, &cpl_save_time, 1, MPI_DOUBLE, &status);
    MPI_File_write(*snapshot, &cpl_snapshot_counter, 1, MPI_INT, &status);
    MPI_File_write(*snapshot, &INTEGRITY_SNAPSHOT, strlen(INTEGRITY_SNAPSHOT), MPI_CHAR, &status);

    MPI_File_close(snapshot);
}

void CPL_SAVE_SNAPSHOT(MPI_File file, void *data, int n, MPI_Datatype type)
{
    MPI_Status status;
    MPI_File_write(file, data, n, type, &status);
}

FILE *CPL_OPEN_SNAPSHOT(char *file_name, char *mode)
{
    FILE *file = fopen(file_name, mode);
    if (!file) {
        fprintf(stderr, "Can't read %s\n", file_name);
        exit(1);
    }

    long offset = sizeof(cpl_start_time) \
                  + sizeof(cpl_save_time) \
                  + sizeof(cpl_snapshot_counter) \
                  + strlen(INTEGRITY_SNAPSHOT); 

    // Move file to the end of file
    fseek(file, -offset, SEEK_END);

    // copy the file into the buffer:
    fread(&cpl_start_time, sizeof(double), 1, file);
    fread(&cpl_save_time, sizeof(double), 1, file);
    fread(&cpl_snapshot_counter, sizeof(int), 1, file);

    // Move file ptr to the begin of file
    rewind(file);

    return file;
}

/*****************************************************************************/
/* Get last snapshot file                                                    */
/*****************************************************************************/
int CPL_GET_SNAPSHOT(char *last_checkpoint)
{
    int myrank = get_comm_rank_();

    char tmp_checkpoint[256] = { 0 };

    char * line = NULL;
    size_t len = 0;
    int i, j;

    FILE *file = fopen(INTEGRITY_SNAPSHOT_FILE, "r");
    if (!file) {
        fprintf(stderr, "Can't read from %s\n", INTEGRITY_SNAPSHOT_FILE);
        exit(1);
    }

    while ((getline(&line, &len, file)) != -1) {
        char tmp[10] = { 0 };
        sprintf(tmp, "%c", line[0]);
        if (atoi(tmp) == myrank) {
            strcpy(last_checkpoint, line);
        }
    }
    fclose(file);

    for (i = 0; i < strlen(last_checkpoint); i++) {
        if (last_checkpoint[i] == '=') {
            break;
        }
    }

    i += 1;

    for (j = 0; i < strlen(last_checkpoint); j++, i++) {
        tmp_checkpoint[j] = last_checkpoint[i];
    }

    strcpy(last_checkpoint, tmp_checkpoint);

    last_checkpoint[strlen(last_checkpoint) - 1] = '\0';

    printf("Rank %d, file %s, phase %d\n", myrank, last_checkpoint, last_checkpoint[0] - '0');

    return last_checkpoint[0] - '0';
}


/*****************************************************************************/
/* Get checkpoint index by checkpoint name                                   */
/*****************************************************************************/
int get_checkpoint_idx_by_name_(void **table, int size, void *name)
{
    int i;
    for (i = 0; i < size; i++) {
        if (table[i] == name)
            break;
    }
    return i;
}


/*****************************************************************************/
/* Init                                                                      */
/*****************************************************************************/
void **init_table_(int size)
{
    void **jump_table = (void**) malloc (sizeof(void*) * size);
    if (!jump_table) {
        fprintf(stderr, "[%s] Can't allocate memory\n", __FUNCTION__);
        exit(1);
    }

    mkdir(SNAPSHOT_DIR_NAME, 0777);

    return jump_table;
}

void CPL_INIT(int size, int argc, char *argv[])
{
    cpl_counter          = 0;
    cpl_size             = size;
    cpl_start_time_local = wtime_();
    cpl_checkpoint_table = init_table_(cpl_size);

    if (get_comm_rank_() == 0) {
        printf("\n[CPL_LIBRARY] note, use: %s [args] [recovery]\n", argv[0]);
    }

    while (argc-->0) {
        if (strcmp(argv[argc], "recovery") == 0) {
            if (get_comm_rank_() == 0) {
                printf("[CPL_LIBRARY] running options 'recovery'\n");
            }
            cpl_run_options = CPL_RECOVERY_MODE;
            break;
        }
    }
}

void CPL_FINALIZE()
{
    int i;

    // Zero process print statistic
    if (get_comm_rank_() == 0) {
        printf("\n");
        // Print up border
        for (i = 0; i < 80; i++) {
            printf("*");
        }
        printf("\n");
        cpl_start_time +=  wtime_() - cpl_start_time_local;
        printf("[CPL_LIBRARY] common elapsed time   = %f sec\n", cpl_start_time);

        double time_one_save = cpl_save_time / (double)cpl_snapshot_counter;
        cpl_save_time = time_one_save * (double)cpl_snapshot_counter;

        printf("[CPL_LIBRARY] elapsed time one save = %f sec\n", time_one_save);
        printf("[CPL_LIBRARY] elapsed time all save = %f sec\n", cpl_save_time);
        printf("[CPL_LIBRARY] amount of snapshots   = %d\n", cpl_snapshot_counter);

        // Print down border
        for (i = 0; i < 80; i++) {
            printf("*");
        }
        printf("\n");
    }
}

int IS_CPL_RECOVERY_MODE()
{
    if (cpl_run_options == CPL_RECOVERY_MODE)
        return 1;
    else
        return 0;
}

void *CPL_COMRESS_DATA_BLOCK(void *data, int size, MPI_Datatype type)
{
/*
    void *compressed_data = NULL;
    // TODO
    if (!compressed_data) {
        fprintf(stderr, "[%s] Can't allocate memory\n", __FUNCTION__);
        exit(1);
    }
    return data;
*/
}

int CPL_IS_DATA_DIFF(struct DeltaCP *buffer, void * data, int size, MPI_Datatype type, int block_idx)
{
    // TODO

    buffer->block  = block_idx;
    buffer->size   = size;
    buffer->type   = type;
    buffer->offset = 0;
    buffer->data   = data;

    return 1;
}

void CPL_SAVE_SNAPSHOT_DELTA(MPI_File file, struct DeltaCP data)
{
    MPI_Status status;
    char string[256] = { 0 };

    if (data.type == MPI_INT) {
        sprintf(string, "\nblock %d\nsize %d\ntype 1\n", data.block, data.size);
    } else if (data.type == MPI_DOUBLE) {
        sprintf(string, "\nblock %d\nsize %d\ntype 2\n", data.block, data.size);
    }

    MPI_File_write(file, string, strlen(string), MPI_CHAR, &status);
    MPI_File_write(file, data.data, data.size, data.type, &status);
}

/*
 * Compressing double-precision numbers by erasing the last 29 bits.
 * In that case we sacrifice number precision,
 * but helps to compress that data array.
 * 29 bits, because ieee 754 number representation (sign, exponent, mantissa)
 */
static double double_trunc(double x)
{
    int zerobits       = 29;
    uint64_t mask      = -(1LL << zerobits);
    uint64_t floatbits = (*((uint64_t*)(&x)));
    floatbits          &= mask;

    return *((double*)(&floatbits));
}


void CPL_SAVE_SNAPSHOT_DELTA_COMRESSED(MPI_File file,
                                       void *data,
                                       int size,
                                       MPI_Datatype type,
                                       int block_idx)
{
    struct DeltaCP buffer;

    uLong offset           = 12;
    uLong data_size        = 0;
    uLongf *compress_size  = NULL;
    void * compressed_data = NULL;

    if (type == MPI_INT) {
        data_size = size * sizeof(int);
        compressed_data = (int *) malloc((sizeof(int) * size) + offset);
        offset += data_size;
    } else if (type == MPI_DOUBLE) {
        data_size = size * sizeof(double);
        compressed_data = (double *) malloc((sizeof(double) * size) + offset);
        offset += data_size;
    }

    if (!compressed_data) {
        fprintf(stderr, "[%s] Can't allocate memory\n",__FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    compress_size = &offset;

    if (type == MPI_DOUBLE) {
        int i;
        double *trunc_ptr = (double *)data;

        for (i = 0; i < size; i++) {
            trunc_ptr[i] = double_trunc(trunc_ptr[i]);
        }
    }

    if (compress((Bytef*)compressed_data, compress_size, (Bytef*)data, data_size) != Z_OK) {
        fprintf(stderr, "[%s] Error in compress\n",__FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (CPL_IS_DATA_DIFF(&buffer, compressed_data, (int)*compress_size, MPI_CHAR, block_idx)) {
        CPL_SAVE_SNAPSHOT_DELTA(file, buffer);
    }

    free(compressed_data);
}