#ifndef _CHECKPOINT_LIB_H_
#define _CHECKPOINT_LIB_H_

/*****************************************************************************/
/* C - check                                                                 */
/* P - point                                                                 */
/* L - library                                                               */
/*****************************************************************************/

#include <stdio.h>
#include <signal.h>
#include <sys/time.h>

#include <mpi.h>

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

/* 
 * Decription:
 * size - checkpoints numbers
 * time - in seconds for timer 
 */

#define CPL_INIT(size, time, argc, argv)                                      \
    cpl_init(size, time, argc, argv);                                         \


#define CPL_FINALIZE()                                                        \
    cpl_finalize();                                                           \


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


/*****************************************************************************/
/* Checkpoint-save macros                                                    */
/*****************************************************************************/
#define CPL_FILE_OPEN(file, phase)                                            \
    open_snapshot_file_(file, phase);                                         \


#define CPL_FILE_CLOSE(file)                                                  \
    close_snapshot_file_(file);                                               \


#define CPL_SAVE_SNAPSHOT(file, data, n, type)                                \
    write_to_snapshot_(file, data, n, type);                                  \


#define CPL_GET_SNAPSHOT(snapshot)                                            \
    get_last_snapshot_(snapshot);                                             \


#define CPL_SAVE_STATE(checkpoint, user_save_callback)                        \
    user_save_callback(get_checkpoint_idx_by_name_(cpl_checkpoint_table,      \
                                                       cpl_size, checkpoint));\


#define CPL_OPEN_SNAPSHOT(file_name, mode)                                    \
    cpl_open_file(file_name, mode);                                           \

/*****************************************************************************/
/* Timer                                                                     */
/*****************************************************************************/
#define CPL_TIMER_INIT()                                                      \
    timer_init_();                                                            \


#define CPL_TIMER_STOP()                                                      \
    timer_stop_();                                                            \


/*****************************************************************************/
/* Prototypes                                                                */
/*****************************************************************************/
void   **init_table_(int size);

void cpl_init(int size, double time, int argc, char *argv[]);
void cpl_finalize();


double wtime_();

void timer_init_();
void timer_stop_();

FILE *cpl_open_file(char *file_name, char *mode);
void write_to_snapshot_(MPI_File file, void *data, int n, MPI_Datatype type);

int get_last_snapshot_(char *last_checkpoint);
int get_checkpoint_idx_by_name_(void **table, int size, void *name);

void open_snapshot_file_(MPI_File *snapshot, int phase);
void close_snapshot_file_(MPI_File *snapshot);

/*****************************************************************************/
/* Run options functions                                                     */
/*****************************************************************************/
int IS_CPL_CHECKPOINT_MODE();
int IS_CPL_RECOVERY_MODE();


#define CPL_SAVE_SNAPSHOT_DELTA(file, buffer)                                 \
    cpl_save_snapshot_delta(file, buffer);                                    \

#define CPL_IS_DATA_DIFF(buffer, source, size, type, idx)                     \
    cpl_is_data_diff(buffer, source, size, type, idx);                        \

int cpl_is_data_diff(struct DeltaCP *buf, void * src, int size, MPI_Datatype type, int delta_idx);
void cpl_save_snapshot_delta(MPI_File file, struct DeltaCP data);

#endif /* _CHECKPOINT_LIB_H_ */
