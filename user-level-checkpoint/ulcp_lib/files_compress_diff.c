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
#include "zfp.h"

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


static void ulcp_get_delta_complex(void *a, void *b, void *delta, int size)
{
    int i;
    int *part_of_a     = (int *)a;
    int *part_of_b     = (int *)b;
    int *part_of_delta = (int *)delta;

    size = size / sizeof(int);

    for (i = 0; i < size; i++)
    {        
        *part_of_delta++ = *part_of_a++ ^ *part_of_b++;
    }
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
        fprintf(stderr, "[ULCP] Memory allocation error in %s, by rank %d\n",
            __FUNCTION__, ulcp_get_comm_rank());
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
        fprintf(stderr, "[ULCP] Can't compress in %s, by rank %d\n",
            __FUNCTION__, ulcp_get_comm_rank());
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (ulcp_is_data_packed(&buffer, compressed_data, (int)*compress_size, MPI_CHAR, block_idx))
    {
        ulcp_snapshot_delta_save(file, buffer);
    }

    free(delta_snapshot);
    free(compressed_data);
}

void ulcp_snapshot_delta_save_compressed_complex(ulcp_action_save_t * action)
{
    int i;
    ulcp_t buffer;

    uLong offset          = 12;
    uLong data_size       = 0;
    uLongf *compress_size = NULL;
    void *compressed_data = NULL;
    void *delta_snapshot  = NULL;
    void *new_snapshot    = NULL;

    if (!action)
    {
        fprintf(stderr, "[ULCP] Bad pointer in %s, by rank %d\n",
            __FUNCTION__, ulcp_get_comm_rank());
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    data_size = action->size * action->size_of;

    compressed_data = (void *)malloc(data_size + offset);
    if (!compressed_data)
    {
        fprintf(stderr, "[ULCP] Memory allocation error in %s, by rank %d\n",
            __FUNCTION__, ulcp_get_comm_rank());
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    offset += data_size;
    compress_size = &offset;

    new_snapshot   = action->data;
    delta_snapshot = malloc(data_size);
    if (!delta_snapshot)
    {
        fprintf(stderr, "[ULCP] Memory allocation error in %s, by rank %d\n",
            __FUNCTION__, ulcp_get_comm_rank());
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    char *char_ulcp_base_snapshot = (char *)ulcp_base_snapshot;
    char *char_new_snapshot       = (char *)new_snapshot;
    char *char_delta_snapshot     = (char *)delta_snapshot;

    for (i = 0; i < action->size; i++)
    {

        ulcp_get_delta_complex(char_ulcp_base_snapshot,
                               char_new_snapshot,
                               char_delta_snapshot,
                               action->size_of);

        memcpy(char_ulcp_base_snapshot, char_new_snapshot, action->size_of);

        char_ulcp_base_snapshot += action->size_of;
        char_new_snapshot       += action->size_of;
        char_delta_snapshot     += action->size_of;
    }

    if (compress((Bytef*)compressed_data, compress_size,
                    (Bytef*)delta_snapshot, data_size) != Z_OK)
    {
        fprintf(stderr, "[ULCP] Can't compress in %s, by rank %d\n",
            __FUNCTION__, ulcp_get_comm_rank());
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    if (ulcp_is_data_packed(&buffer, compressed_data,
            (int)*compress_size, MPI_CHAR, action->delta_idx))
    {
        ulcp_snapshot_delta_save(action->file, buffer);
    }

    free(delta_snapshot);
    free(compressed_data);
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
        fprintf(stderr, "[ULCP] <%s> Can't allocate memory\n",__FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    compress_size = &offset;

    if (compress((Bytef*)compressed_data, compress_size, (Bytef*)data, data_size) != Z_OK)
    {
        fprintf(stderr, "[ULCP] <%s> Error in compress\n",__FUNCTION__);
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
    if (!action)
    {
        fprintf(stderr, "[ULCP] Bad pointer in %s\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    ulcp_base_snapshot = (void *)malloc(action->size_of * action->size);
    if (!ulcp_base_snapshot)
    {
        fprintf(stderr, "[ULCP] Memory allocation error in %s\n", __FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Memory zero out
    memset(ulcp_base_snapshot, 0, action->size);

    fprintf(stdout, "[ULPC] Succeffuly set up base checkpoint in %s(), by rank %d\n",
        __FUNCTION__, ulcp_get_comm_rank());
}


// Compress or decompress array
static int ulcp_zfp_compress(MPI_File  file,
                             double   *array,
                             int       size,
                             double    tolerance,
                             int       block_idx)
{
    int status = 0;    // return value: 0 = success
    zfp_type type;     // array scalar type
    zfp_field* field;  // array meta data
    zfp_stream* zfp;   // compressed stream
    void* buffer;      // storage for compressed stream
    size_t bufsize;    // byte size of compressed buffer
    bitstream* stream; // bit stream to write to or read from
    size_t zfpsize;    // byte size of compressed stream

    // Allocate meta data for the 1D array a[size]
    type  = zfp_type_double;
    field = zfp_field_1d(array, type, size);

    // Allocate meta data for a compressed stream
    zfp = zfp_stream_open(NULL);

    // Set compression mode and parameters via one of three functions
    //zfp_stream_set_rate(zfp, rate, type, 3, 0);
    //zfp_stream_set_precision(zfp, precision);
    zfp_stream_set_accuracy(zfp, tolerance);

    // Allocate buffer for compressed data
    bufsize = zfp_stream_maximum_size(zfp, field);
    buffer  = malloc(bufsize);

    // Associate bit stream with allocated buffer
    stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

    // Compress array and output compressed stream
    zfpsize = zfp_compress(zfp, field);
    if (!zfpsize)
    {
        fprintf(stderr, "[USLP] <%s> Compression failed\n", __FUNCTION__);
        status = 1;
    } 
    else
    {
        ulcp_t buffer_container;
        if (ulcp_is_data_packed(&buffer_container, buffer, zfpsize, MPI_CHAR, block_idx))
        {
            ulcp_snapshot_delta_save(file, buffer_container);
        }
    }

    // Clean up
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    free(buffer);

    return status;
}

void ulcp_snapshot_delta_save_zfp_compressed(MPI_File     file,
                                             void        *data,
                                             int          size,
                                             MPI_Datatype type,
                                             int          block_idx)
{
    int i;

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

    if (ulcp_zfp_compress(file, delta_snapshot, size, 1e-3, block_idx) != 0)
    {
        fprintf(stderr, "[ULCP] <%s> Error in compress\n",__FUNCTION__);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);        
    }

    free(delta_snapshot);
}
