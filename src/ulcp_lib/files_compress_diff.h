#ifndef _ULCP_FILES_COMPRESS_DIFF_H_
#define _ULCP_FILES_COMPRESS_DIFF_H_

#include <mpi.h>
#ifdef ULFM_SUPPORT
#include <mpi-ext.h> // ULFM support
#endif /* ULFM_SUPPORT */


void ulcp_snapshot_delta_save(MPI_File file, ulcp_t data);

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

void ulcp_snapshot_save_compressed_complex(MPI_File file,
                                   void *data,
                                   int size,
                                   int sizeof_type,
                                   int block_idx);

void ulcp_snapshot_set_diff(ulcp_action_t * action);

#endif /* _ULCP_FILES_COMPRESS_DIFF_H_ */
