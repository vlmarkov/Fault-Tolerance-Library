#include "ulcp.h"

#include <mpi.h>
#include <mpi-ext.h> // ULFM support

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "zlib.h"

static int ulcp_is_data_packed(ulcp_t *buffer,
                               void *data,
                               int size,
                               MPI_Datatype type,
                               int block_idx)
{
    // TODO

    buffer->block  = block_idx;
    buffer->size   = size;
    buffer->type   = type;
    buffer->offset = 0;
    buffer->data   = data;

    return 1;
}

/*
 * Very WIP
 */
static void ulcp_get_double_delta(double *a, double *b, double *delta)
{
    int *part_of_double_a     = (int *)a;
    int *part_of_double_b     = (int *)b;
    int *part_of_double_delta = (int *)delta;

    *part_of_double_delta++ = *part_of_double_a++ ^ *part_of_double_b++;
    *part_of_double_delta   = *part_of_double_a   ^ *part_of_double_b;
}

void ulcp_snapshot_delta_save(MPI_File file, ulcp_t data)
{
    MPI_Status status;
    char string[256] = { 0 };

    if (data.type == MPI_INT)
    {
        sprintf(string, "\nblock %d\nsize %d\ntype 1\n", data.block, data.size);
    } 
    else if (data.type == MPI_DOUBLE)
    {
        sprintf(string, "\nblock %d\nsize %d\ntype 2\n", data.block, data.size);
    }

    MPI_File_write(file, string, strlen(string), MPI_CHAR, &status);
    MPI_File_write(file, data.data, data.size, data.type, &status);
}

void ulcp_snapshot_delta_save_compressed(MPI_File file,
                                         void *data,
                                         int size,
                                         MPI_Datatype type,
                                         int block_idx)
{
    ulcp_t buffer;

    uLong offset           = 12;
    uLong data_size        = 0;
    uLongf *compress_size  = NULL;
    void * compressed_data = NULL;

    if (type == MPI_INT)
    {
        data_size = size * sizeof(int);
        compressed_data = (int *)malloc((sizeof(int) * size) + offset);
        offset += data_size;
    } 
    else if (type == MPI_DOUBLE)
    {
        data_size = size * sizeof(double);
        compressed_data = (double *)malloc((sizeof(double) * size) + offset);
        offset += data_size;
    }

    if (!compressed_data)
    {
        fprintf(stderr, "[%s] Can't allocate memory\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    compress_size = &offset;

    int i = 0;
    double *new_snapshot = (double *)data;
    double *delta_snapshot = (double *)malloc(sizeof(double) * size);

    for (i = 0; i < size; i++)
    {
        ulcp_get_double_delta(&ulcp_base_snapshot[i], &new_snapshot[i], &delta_snapshot[i]);
        ulcp_base_snapshot[i] = new_snapshot[i];
    }

    if (compress((Bytef*)compressed_data, compress_size, (Bytef*)delta_snapshot, data_size) != Z_OK)
    {
        fprintf(stderr, "[%s] Error in compress\n",__FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (ulcp_is_data_packed(&buffer, compressed_data, (int)*compress_size, MPI_CHAR, block_idx))
    {
        ulcp_snapshot_delta_save(file, buffer);
    }

    free(delta_snapshot);
    free(compressed_data);

    printf("<%s><%d>\n", __FUNCTION__, __LINE__);
}

void ulcp_snapshot_save_compressed(MPI_File file,
                                   void *data,
                                   int size,
                                   MPI_Datatype type,
                                   int block_idx)
{
    ulcp_t buffer;

    uLong offset          = 12;
    uLong data_size       = 0;
    uLongf *compress_size = NULL;
    void *compressed_data = NULL;

    if (type == MPI_INT)
    {
        data_size = size * sizeof(int);
        compressed_data = (int *)malloc((sizeof(int) * size) + offset);
        offset += data_size;
    } 
    else if (type == MPI_DOUBLE)
    {
        data_size = size * sizeof(double);
        compressed_data = (double *)malloc((sizeof(double) * size) + offset);
        offset += data_size;
    }

    if (!compressed_data)
    {
        fprintf(stderr, "[%s] Can't allocate memory\n",__FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    compress_size = &offset;

    if (compress((Bytef*)compressed_data, compress_size, (Bytef*)data, data_size) != Z_OK)
    {
        fprintf(stderr, "[%s] Error in compress\n",__FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (ulcp_is_data_packed(&buffer, compressed_data, (int)*compress_size, MPI_CHAR, block_idx))
    {
        ulcp_snapshot_delta_save(file, buffer);
    }

    free(compressed_data);
}

void ulcp_snapshot_set_diff(MPI_Datatype type, int size)
{
    int i = 0;

    if (type == MPI_DOUBLE)
    {
        ulcp_base_snapshot = (double *)malloc(sizeof(double) * size);
    }

    for (i = 0; i < size; i++)
    {
        ulcp_base_snapshot[i] = 0.0;
    }
}
