#ifndef _CHECKPOINT_LIB_H_
#define _CHECKPOINT_LIB_H_

#include <stdio.h>
#include <signal.h>
#include <sys/time.h>

#include <mpi.h>

#include <mpi.h>
#include <mpi-ext.h>

#define SNAPSHOT_DIR_NAME "snapshot"
#define INTEGRITY_SNAPSHOT "\n=end_of_file="
#define INTEGRITY_SNAPSHOT_FILE "integrity_file.txt"

/*****************************************************************************/
/* Global variables                                                          */
/*****************************************************************************/
extern void   **cpl_checkpoint_table;

extern int cpl_size;
extern int cpl_counter;

struct DeltaCP
{
    int            block;  // checkpoint block number
    int            offset;
    int            size;   // data size
    MPI_Datatype   type;   // data type
    void           *data;  // data itself
    struct DeltaCP *next;
};


/*****************************************************************************/
/* Initializing checkpoint library macros                                    */
/*****************************************************************************/
void CPL_INIT(int size, int argc, char *argv[]);
void CPL_FINALIZE();

/*****************************************************************************/
/* Declaration checkpoint macros                                             */
/*****************************************************************************/

/*
 * Description:
 * name - checkpoint_one, checkpoint_two, etc
 *      - phase_one, phase_two, etc
 *      - one, two, three, etc 
 */

#define CPL_DECLARATE_CHECKPOINT(name)                                        \
    cpl_checkpoint_table[cpl_counter++] = name;                               \


/*****************************************************************************/
/* Control flow-macros                                                       */
/*****************************************************************************/
#define CPL_GO_TO_CHECKPOINT(idx)                                             \
    goto *cpl_checkpoint_table[idx];                                          \


#define CPL_SET_CHECKPOINT(checkpoint_name)                                   \
    checkpoint_name :                                                         \

#define CPL_SAVE_STATE(checkpoint, user_save_callback)                        \
    user_save_callback(get_checkpoint_idx_by_name_(cpl_checkpoint_table,      \
                                                       cpl_size, checkpoint));\

int get_checkpoint_idx_by_name_(void **table, int size, void *name);


/*****************************************************************************/
/* Checkpoint-save functions                                                 */
/*****************************************************************************/
void CPL_FILE_OPEN(MPI_File *snapshot, int phase);
void CPL_FILE_CLOSE(MPI_File *snapshot);

void CPL_SAVE_SNAPSHOT(MPI_File file, void *data, int size, MPI_Datatype type);
int CPL_GET_SNAPSHOT(char *file_name);
FILE *CPL_OPEN_SNAPSHOT(char *file_name, char *mode);
int CPL_IS_DATA_PACK(struct DeltaCP *buf,
                     void * src,
                     int size,
                     MPI_Datatype type,
                     int delta_idx);

void CPL_SAVE_SNAPSHOT_DELTA(MPI_File file, struct DeltaCP data);

void CPL_SAVE_SNAPSHOT_DELTA_COMRESSED(MPI_File file,
                                       void *data,
                                       int size,
                                       MPI_Datatype type,
                                       int delta_idx);

void CPL_SAVE_SNAPSHOT_COMRESSED(MPI_File file,
                                 void *data,
                                 int size,
                                 MPI_Datatype type,
                                 int delta_idx);

/*****************************************************************************/
/* Run options functions                                                     */
/*****************************************************************************/
int IS_CPL_CHECKPOINT_MODE();
int IS_CPL_RECOVERY_MODE();

/*
 * Very WIP
 */
void get_delta_double(double *a, double *b, double *delta);
void CPL_SET_DIFF_SNAPSHOT(MPI_Datatype type, int size);

#endif /* _CHECKPOINT_LIB_H_ */
