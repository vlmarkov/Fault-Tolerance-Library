#ifndef _ULCP_FILES_COMPRESS_DIFF_H_
#define _ULCP_FILES_COMPRESS_DIFF_H_

#include <mpi.h>
#include <mpi-ext.h>

#include "ulcp_lib.h"

void ulcp_snapshot_delta_save(MPI_File file, struct DeltaCP data);

void ulcp_snapshot_delta_save_compressed(MPI_File file,
                                       void *data,
                                       int size,
                                       MPI_Datatype type,
                                       int block_idx);

void ulcp_snapshot_save_compressed(MPI_File file,
                                 void *data,
                                 int size,
                                 MPI_Datatype type,
                                 int block_idx);

void ulcp_snapshot_set_diff(MPI_Datatype type, int size);

#endif /* _ULCP_FILES_COMPRESS_DIFF_H_ */
