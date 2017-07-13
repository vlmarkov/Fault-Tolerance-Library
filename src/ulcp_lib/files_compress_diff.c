#include "ulcp.h"

#include <mpi.h>
#ifdef ULFM_SUPPORT
#include <mpi-ext.h> // ULFM support
#endif /* ULFM_SUPPORT */

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
static void ulcp_get_delta_double(double *a, double *b, double *delta)
{
    int *part_of_double_a     = (int *)a;
    int *part_of_double_b     = (int *)b;
    int *part_of_double_delta = (int *)delta;

    *part_of_double_delta++ = *part_of_double_a++ ^ *part_of_double_b++;
    *part_of_double_delta   = *part_of_double_a   ^ *part_of_double_b;
}

static void ulcp_get_delta_int(int *a, int *b, int *delta)
{
    *delta = *a ^ *b;
}

/*
static void ulcp_get_delta_complex(void *a, void *b, void *delta, size_t size)
{
    int i;
    int *part_of_a     = (int *)a;
    int *part_of_b     = (int *)b;
    int *part_of_delta = (int *)delta;

    for (i = 0; i < size; i++)
    {        
        *part_of_delta++ = *part_of_a++ ^ *part_of_b++;
    }
}
*/

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
    void *data, int size, MPI_Datatype type, int block_idx)
{
    int i;
    ulcp_t buffer;

    uLong offset           = 12;
    uLong data_size        = 0;
    uLongf *compress_size  = NULL;
    void * compressed_data = NULL;

    if (type == MPI_INT)
    {
        data_size = size * sizeof(int);
        compressed_data = (void *)malloc((sizeof(int) * size) + offset);
        offset += data_size;
    } 
    else if (type == MPI_DOUBLE)
    {
        data_size = size * sizeof(double);
        compressed_data = (void *)malloc((sizeof(double) * size) + offset);
        offset += data_size;
    }

    if (!compressed_data)
    {
        fprintf(stderr, "[ULCP] [%s] Error in compress\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    compress_size = &offset;

    void *delta_snapshot = NULL;
    void *new_snapshot   = NULL;

    if (type == MPI_INT)
    {
        new_snapshot   = (int *)data;
        delta_snapshot = (int *)malloc(sizeof(int) * size);

        int *int_ulcp_base_snapshot = (int *)ulcp_base_snapshot;
        int *int_new_snapshot = (int *)new_snapshot;
        int *int_delta_snapshot = (int *)delta_snapshot;

        for (i = 0; i < size; i++)
        {

            ulcp_get_delta_int(&int_ulcp_base_snapshot[i],
                &int_new_snapshot[i],
                &int_delta_snapshot[i]);

            int_ulcp_base_snapshot[i] = int_new_snapshot[i];
        }
    }
    else if (type == MPI_DOUBLE)
    {
        new_snapshot   = (double *)data;
        delta_snapshot = (double *)malloc(sizeof(double) * size);
        
        double *double_ulcp_base_snapshot = (double *)ulcp_base_snapshot;
        double *double_new_snapshot = (double *)new_snapshot;
        double *double_delta_snapshot = (double *)delta_snapshot;

        for (i = 0; i < size; i++)
        {
            ulcp_get_delta_double(&double_ulcp_base_snapshot[i],
                &double_new_snapshot[i],
                &double_delta_snapshot[i]);

            double_ulcp_base_snapshot[i] = double_new_snapshot[i];
        }
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

void ulcp_snapshot_save_compressed_complex(MPI_File file,
                                   void *data,
                                   int size,
                                   int sizeof_type,
                                   int block_idx)
{
    ulcp_t buffer;

    uLong offset          = 12;
    uLong data_size       = 0;
    uLongf *compress_size = NULL;
    void *compressed_data = NULL;

    data_size = size * sizeof(int);
    compressed_data = (int *)malloc((sizeof_type * size) + offset);
    offset += data_size;

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

    if (ulcp_is_data_packed(&buffer,
                            compressed_data,
                            (int)*compress_size,
                            MPI_CHAR,
                            block_idx))
    {
        ulcp_snapshot_delta_save(file, buffer);
    }

    free(compressed_data);
}

void ulcp_snapshot_set_diff(ulcp_action_t * action)
{
    int i;
    int size_of_elem = 0;

    if (!action)
    {
        fprintf(stderr, "[ULCP] Bad pointer in %s\n", __FUNCTION__);
        exit(-1);
    }

    if (action->mode == ULCP_SET_MODE_SIMPLE)
    {
        if (action->mpi_type == MPI_DOUBLE)
        {
            ulcp_base_snapshot = (double *)malloc(sizeof(double) * action->size);
            size_of_elem = sizeof(double);
        }
    }
    else if (action->mode == ULCP_SET_MODE_COMPLEX)
    {
        ulcp_base_snapshot = (void *)malloc(action->user_type * action->size);
        size_of_elem = action->user_type;
    }

    if (!ulcp_base_snapshot)
    {
        fprintf(stderr, "[ULCP] Memory allocation error in %s\n", __FUNCTION__);
        exit(-1);
    }

    void * tmp_ptr = ulcp_base_snapshot;

    for (i = 0; i < action->size; i++)
    {
        bzero(tmp_ptr, size_of_elem);
        tmp_ptr += size_of_elem;
    }
}
