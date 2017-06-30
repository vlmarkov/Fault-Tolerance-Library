#ifndef _ULCP_FILES_ROUTINE_H_
#define _ULCP_FILES_ROUTINE_H_

#include <mpi.h>
#include <mpi-ext.h>

#include <stdio.h>

void ulcp_open_file(MPI_File* file, int snapshot_phase);
void ulcp_close_file(MPI_File* file);
void ulcp_snapshot_save(MPI_File file, void *data, int size, MPI_Datatype type);
int ulcp_get_snapshot(char *last_snapshot);
int ulcp_get_snapshot_idx_by_name(void **table, int size, void *name);
FILE* ulcp_open_snapshot(char *file_name, char *mode);

#endif /* _ULCP_FILES_ROUTINE_H_ */
